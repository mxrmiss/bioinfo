#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: venn.R (v1.6 Layout Polish)
# Description: 
#   - [Visual] Threshold info moved to Bottom Caption (Single line).
#   - [Config] Toggles for Title and Caption (Legend).
#   - [Config] Custom Title support.
#   - [Visual] Pastel colors, Borderless, Clean style.
# ==============================================================================

# ==============================================================================
# Script: venn.R
# Purpose:
#   - 从多个差异表达（DEG）结果表中提取“显著基因集合”，绘制交集图：
#       * 当集合数 <= 4：绘制 Venn Diagram（ggvenn）
#       * 当集合数  > 4：绘制 UpSet Plot（ComplexHeatmap）
#   - 同时导出交集与各集合特异（unique）基因列表，便于后续功能富集、候选基因筛选等分析。
#
# Inputs (输入原材料):
#   1) DEG tables (TSV) placed under IN_DIR (default: "input/")
#      - 文件格式：制表符分隔 .tsv
#      - 关键表头（必须一模一样）：gene_id, log2fc, p_adjust
#      - 语义要求：
#          * gene_id   : 基因 ID（建议为 gene-level ID，且不同文件必须同一套 ID 体系）
#          * log2fc    : log2 fold change（可正可负）
#          * p_adjust  : 多重校正后的 P 值（FDR / padj；支持科学计数法）
#
#   2) File selection mode (两种输入文件选择方式):
#      - Manual mode: 在 TARGET_FILES 填入需要处理的文件名（相对 input/ 目录）
#          * 例如：TARGET_FILES <- c("foot_vs_mantle_DEG_all.tsv", "foot_vs_gill_DEG_all.tsv")
#      - Auto-scan mode: TARGET_FILES 为空时，自动扫描 input/ 下所有 .tsv
#          * 注意：会自动跳过文件名包含 tpm|count|sample|meta|template 的 .tsv
#
# Parameters (阈值与显示控制):
#   - FC_CUTOFF : |log2fc| 的阈值（默认 1.0）
#   - P_CUTOFF  : p_adjust 的阈值（默认 0.05）
#   - SHOW_PLOT_TITLE / CUSTOM_PLOT_TITLE : 控制标题显示与自定义标题
#   - SHOW_CAPTION : 是否在图底部显示阈值说明（caption）
#   - SET_COLORS : 集合填充色（最多支持 4 组 Venn，超过则用 UpSet）
#
# Outputs (输出产物):
#   输出目录：OUT_DIR（默认 output/venn/）
#   1) 图文件：
#      - 当集合数 <= 4：
#          * Venn_*.pdf
#          * Venn_*.png
#      - 当集合数  > 4：
#          * UpSet_Combined_Sets.pdf
#          * UpSet_Combined_Sets.png
#   2) 基因列表：
#      - Genes_Intersection_ALL.txt
#          * 所有集合共同交集（core genes）
#      - Genes_Unique_<SetName>.txt
#          * 每个集合的特异基因（相对于其它所有集合的 setdiff）
#
# Interpretation Notes (结果解读注意事项):
#   - 本脚本默认用 abs(log2fc) > FC_CUTOFF，因此“上调/下调”会混在同一个集合里做交集。
#     若要分别做“共同上调”或“共同下调”的交集：
#       * 建议输入文件层面就分开（例如使用 DEG_up.tsv / DEG_down.tsv），或在筛选逻辑中改条件。
#   - 交集分析的前提是各文件 gene_id 可比且一致：
#       * 不要混用 transcript_id 与 gene_id
#       * 不要混用不同物种的原始基因 ID（除非你已经做了同源映射/统一 ID）
#   - 当输入集合 > 4 时需要安装 ComplexHeatmap（否则会报错并停止）。
#
# Practical Tips (实用建议):
#   - 为了更稳定/更“生物学一致”的交集，建议输入使用统一标准产生的 DEG 结果表
#     （例如同一个 DE 流程、同一版本注释、同一过滤策略）。
#   - 文件命名会影响集合名称：脚本会自动从文件名中去掉 "_DEG_all" 与 "_vs_D0"；
#     若你希望图例更清晰，建议在文件名里体现 contrast 信息但尽量简短。
# ==============================================================================


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tools)
  
  if (!require("ggvenn", quietly = TRUE)) {
    stop("[ERROR] 'ggvenn' is missing. Install via: install.packages('ggvenn')")
  }
  library(ggvenn)
  
  has_complex_heatmap <- require("ComplexHeatmap", quietly = TRUE)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Input Control ---
TARGET_FILES <- c()

# --- Thresholds ---
FC_CUTOFF   <- 1.0    
P_CUTOFF    <- 0.05   

# --- Display Control ---
SHOW_PLOT_TITLE   <- TRUE       # Toggle Top Title
SHOW_CAPTION      <- TRUE       # Toggle Bottom Threshold Info
CUSTOM_PLOT_TITLE <- ""         # Override title string

# --- Aesthetics ---
FONT_FAMILY  <- "Arial"
EDGE_SIZE    <- 0.3             
# Fresh Pastel Palette
SET_COLORS   <- c("#8DD3C7", "#FB8072", "#BEBADA", "#80B1D3") 
PLOT_MARGIN  <- 20              

#Paths
IN_DIR       <- "input"
OUT_DIR      <- file.path("output", "venn")
TEMPLATE_DIR <- "templates"
TEMPLATE_FILE<- "volcano_input_template.tsv" 

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  tmpl_path <- file.path(TEMPLATE_DIR, TEMPLATE_FILE)
  if (!file.exists(tmpl_path)) {
    example_data <- data.frame(
      gene_id  = c("Gene_A", "Gene_B", "Gene_C"),
      log2fc   = c(2.5, -3.1, 0.2),
      p_adjust = c(0.0001, 1.2e-10, 0.8),
      stringsAsFactors = FALSE
    )
    write_tsv(example_data, tmpl_path)
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic Venn/UpSet Plotter v1.6 (Layout Polish)\n")
cat("================================================================\n")

init_template_system()

# --- Step 1: Identify Files ---
files_to_process <- c()

if (length(TARGET_FILES) > 0) {
  for (f in TARGET_FILES) {
    f_path <- file.path(IN_DIR, f)
    if (file.exists(f_path)) files_to_process <- c(files_to_process, f_path)
    else cat(sprintf("   [WARN] Specified file not found: %s\n", f))
  }
} else {
  cat("Mode: Auto-scan 'input/' for DEG tables...\n")
  all_files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
  files_to_process <- all_files[!grepl("tpm|count|sample|meta|template", basename(all_files), ignore.case = TRUE)]
}

if (length(files_to_process) < 2) stop("[ERROR] Need at least 2 DEG files.")

cat(sprintf(">> Found %d files. Extracting significant genes...\n", length(files_to_process)))

# --- Step 2: Load & Filter Data ---
gene_sets <- list()

for (f_path in files_to_process) {
  fname <- basename(f_path)
  set_name <- file_path_sans_ext(fname)
  set_name <- gsub("_DEG_all", "", set_name)
  set_name <- gsub("_vs_D0", "", set_name) 
  
  df <- tryCatch(read_tsv(f_path, show_col_types = FALSE), error = function(e) NULL)
  
  if (is.null(df)) next
  if (!all(c("gene_id", "log2fc", "p_adjust") %in% colnames(df))) {
    cat(sprintf("   [SKIP] Invalid headers in %s\n", fname))
    next
  }
  
  sig_genes <- df %>% 
    filter(p_adjust < P_CUTOFF & abs(log2fc) > FC_CUTOFF) %>%
    pull(gene_id)
  
  cat(sprintf("   + %s: %d sig genes\n", set_name, length(sig_genes)))
  
  if (length(sig_genes) > 0) {
    gene_sets[[set_name]] <- sig_genes
  }
}

n_sets <- length(gene_sets)
if (n_sets < 2) stop("[ERROR] Not enough valid gene sets (<2).")

# --- Step 3: Plotting Logic ---

# [FIX] Single Line Caption at Bottom
caption_text <- NULL
if (SHOW_CAPTION) {
  caption_text <- sprintf("Thresholds: P.adjust < %s & |log2FC| > %s", P_CUTOFF, FC_CUTOFF)
}

# Title Logic
final_title <- NULL
if (SHOW_PLOT_TITLE) {
  if (CUSTOM_PLOT_TITLE != "") {
    final_title <- CUSTOM_PLOT_TITLE
  } else {
    final_title <- "Venn Analysis"
  }
}

if (n_sets <= 4) {
  # >>> VENN DIAGRAM MODE <<<
  cat(">> Mode: Pastel Venn Diagram (<= 4 sets)\n")
  
  p <- ggvenn(
    gene_sets, 
    fill_color = SET_COLORS[1:n_sets],
    fill_alpha = 0.5,
    stroke_color = "white",
    stroke_size = 0.5,
    set_name_size = 6,
    text_size = 4
  ) +
  # [FIX] Title Top, Caption Bottom
  labs(title = final_title, caption = caption_text) +
  theme(
    text = element_text(family = FONT_FAMILY),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    # [FIX] Center align caption at bottom
    plot.caption = element_text(hjust = 0.5, size = 12, color = "grey30", margin = margin(t = 10)),
    plot.margin = margin(t = PLOT_MARGIN, r = PLOT_MARGIN, b = PLOT_MARGIN, l = PLOT_MARGIN, unit = "pt"),
    legend.position = "bottom" # Just in case
  )
  
  out_name <- paste0("Venn_", paste(names(gene_sets), collapse="_"))
  if (nchar(out_name) > 50) out_name <- "Venn_Combined_Sets"
  
  tryCatch({
    ggsave(file.path(OUT_DIR, paste0(out_name, ".pdf")), p, width = 8, height = 8, device = cairo_pdf)
    ggsave(file.path(OUT_DIR, paste0(out_name, ".png")), p, width = 8, height = 8, dpi = 300, bg = "white")
    cat("   [SUCCESS] Venn Diagram Saved.\n")
  }, error = function(e) cat(sprintf("   [ERROR] Save failed: %s\n", e$message)))
  
} else {
  # >>> UPSET PLOT MODE <<<
  cat(">> Mode: UpSet Plot (> 4 sets)\n")
  if (!has_complex_heatmap) stop("[ERROR] >4 sets requires 'ComplexHeatmap'.")
  
  m <- make_comb_mat(gene_sets)
  
  # For UpSet, we append caption to title because it handles layout differently
  upset_title <- final_title
  if (!is.null(caption_text)) {
    upset_title <- paste0(upset_title, "\n", caption_text)
  }
  
  out_name <- "UpSet_Combined_Sets"
  
  tryCatch({
    pdf(file.path(OUT_DIR, paste0(out_name, ".pdf")), width = 10, height = 6)
    draw(UpSet(m, column_title = upset_title))
    dev.off()
    
    png(file.path(OUT_DIR, paste0(out_name, ".png")), width = 3000, height = 1800, res = 300)
    draw(UpSet(m, column_title = upset_title))
    dev.off()
    cat("   [SUCCESS] UpSet Plot Saved.\n")
  }, error = function(e) { if (dev.cur() > 1) dev.off(); cat(sprintf("   [ERROR] Save failed: %s\n", e$message)) })
}

# --- Step 4: Export Intersections ---
cat(">> Exporting Intersection Lists...\n")
core_genes <- Reduce(intersect, gene_sets)
if (length(core_genes) > 0) {
  writeLines(core_genes, file.path(OUT_DIR, "Genes_Intersection_ALL.txt"))
}
for (set_name in names(gene_sets)) {
  others <- gene_sets[names(gene_sets) != set_name]
  unique_genes <- setdiff(gene_sets[[set_name]], unlist(others))
  if (length(unique_genes) > 0) {
    writeLines(unique_genes, file.path(OUT_DIR, paste0("Genes_Unique_", set_name, ".txt")))
  }
}

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
