#!/usr/bin/env Rscript
# ==============================================================================
# Script: upset.R
# Purpose:
#   - 从多个差异表达（DEG）结果表中提取“显著基因集合”，绘制 UpSet Plot（ComplexHeatmap）
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
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tools)
  library(grid)

  if (!require("ComplexHeatmap", quietly = TRUE)) {
    stop("[ERROR] 'ComplexHeatmap' is missing. Install it first.")
  }
  library(ComplexHeatmap)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# 指定要分析的文件名（只写文件名，不写路径）；留空则自动扫描 IN_DIR 下所有符合条件的 tsv
TARGET_FILES <- c()

# 输入目录
IN_DIR <- "input"

# 输出目录
OUT_DIR <- file.path("output", "upset")

# |log2FC| 阈值
FC_CUTOFF <- 1.0

# adjusted P value / FDR 阈值
P_CUTOFF <- 0.05

# 是否显示主标题
SHOW_PLOT_TITLE <- FALSE

# 自定义标题；留空则使用默认标题
CUSTOM_PLOT_TITLE <- ""

# 是否显示底部阈值说明
SHOW_CAPTION <- TRUE

# 字体名称；默认使用 sans，最稳
FONT_FAMILY <- "sans"

# 顶部交集柱颜色
UPSET_MAIN_BAR_COLOR <- "#F79D93"

# 右侧集合大小柱颜色
UPSET_SET_BAR_COLOR <- "#9FD5CB"

# 选中点与连接线颜色
UPSET_ACTIVE_DOT_COLOR <- "#F79D93"

# 背景点颜色
UPSET_BG_DOT_COLOR <- "#E3E6EB"

# 背景连接线颜色
UPSET_LINE_COLOR <- "#A9AFB8"

# 文字颜色
UPSET_TEXT_COLOR <- "#3F4650"

# 坐标轴文字颜色
UPSET_AXIS_COLOR <- "#6E7681"

# 柱边框颜色；NA 表示不画边框
UPSET_BORDER_COLOR <- NA

# PDF 宽度（英寸）
UPSET_PDF_WIDTH <- 10

# PDF 高度（英寸）
UPSET_PDF_HEIGHT <- 6

# PNG 宽度（像素）
UPSET_PNG_WIDTH <- 3000

# PNG 高度（像素）
UPSET_PNG_HEIGHT <- 1800

# PNG 分辨率
UPSET_PNG_RES <- 300

# 顶部柱图区高度
UPSET_TOP_HEIGHT <- unit(3.4, "cm")

# 右侧集合柱区域宽度
UPSET_RIGHT_WIDTH <- unit(2.8, "cm")

# 点大小
UPSET_DOT_SIZE <- unit(3.2, "mm")

# 连接线宽度
UPSET_LINE_WIDTH <- 1.2

# 标题字号
UPSET_TITLE_SIZE <- 15

# 底部说明字号
UPSET_CAPTION_SIZE <- 16

# 集合名/注释字号
UPSET_LABEL_SIZE <- 16

# 坐标轴字号
UPSET_AXIS_SIZE <- 14

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

ensure_outdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

safe_dev_off <- function() {
  while (dev.cur() > 1) dev.off()
}

build_upset_plot <- function(m, final_title = NULL) {
  set_ord  <- order(set_size(m), decreasing = TRUE)
  comb_ord <- order(comb_size(m), decreasing = TRUE)

  top_anno <- upset_top_annotation(
    m,
    height = UPSET_TOP_HEIGHT,
    gp = gpar(fill = UPSET_MAIN_BAR_COLOR, col = UPSET_BORDER_COLOR),
    annotation_name_gp = gpar(
      fontfamily = FONT_FAMILY,
      fontsize = UPSET_LABEL_SIZE,
      fontface = "bold",
      col = UPSET_TEXT_COLOR
    ),
    axis_param = list(
      gp = gpar(
        fontfamily = FONT_FAMILY,
        fontsize = UPSET_AXIS_SIZE,
        col = UPSET_AXIS_COLOR
      ),
      side = "left"
    )
  )

  right_anno <- upset_right_annotation(
    m,
    width = UPSET_RIGHT_WIDTH,
    gp = gpar(fill = UPSET_SET_BAR_COLOR, col = UPSET_BORDER_COLOR),
    annotation_name_gp = gpar(
      fontfamily = FONT_FAMILY,
      fontsize = UPSET_LABEL_SIZE,
      fontface = "bold",
      col = UPSET_TEXT_COLOR
    ),
    axis_param = list(
      gp = gpar(
        fontfamily = FONT_FAMILY,
        fontsize = UPSET_AXIS_SIZE,
        col = UPSET_AXIS_COLOR
      )
    )
  )

  UpSet(
    m,
    set_order = set_ord,
    comb_order = comb_ord,
    top_annotation = top_anno,
    right_annotation = right_anno,
    comb_col = UPSET_ACTIVE_DOT_COLOR,
    pt_size = UPSET_DOT_SIZE,
    lwd = UPSET_LINE_WIDTH,
    bg_col = UPSET_LINE_COLOR,
    bg_pt_col = UPSET_BG_DOT_COLOR,
    row_names_side = "left",
    row_names_gp = gpar(
      fontfamily = FONT_FAMILY,
      fontsize = UPSET_LABEL_SIZE,
      col = UPSET_TEXT_COLOR
    ),
    column_title = final_title,
    column_title_gp = gpar(
      fontfamily = FONT_FAMILY,
      fontsize = UPSET_TITLE_SIZE,
      fontface = "bold",
      col = UPSET_TEXT_COLOR
    )
  )
}

draw_caption <- function(caption_text) {
  if (!is.null(caption_text) && nzchar(caption_text)) {
    grid.text(
      label = caption_text,
      x = unit(0.5, "npc"),
      y = unit(0.02, "npc"),
      just = c("center", "bottom"),
      gp = gpar(
        fontfamily = FONT_FAMILY,
        fontsize = UPSET_CAPTION_SIZE,
        col = UPSET_AXIS_COLOR
      )
    )
  }
}

save_upset_plot <- function(ht, pdf_file, png_file, caption_text = NULL) {
  pdf(
    file = pdf_file,
    width = UPSET_PDF_WIDTH,
    height = UPSET_PDF_HEIGHT,
    family = FONT_FAMILY
  )
  draw(ht, padding = unit(c(8, 4, 10, 6), "mm"))
  draw_caption(caption_text)
  dev.off()

  png(
    filename = png_file,
    width = UPSET_PNG_WIDTH,
    height = UPSET_PNG_HEIGHT,
    res = UPSET_PNG_RES,
    bg = "white"
  )
  draw(ht, padding = unit(c(8, 4, 10, 6), "mm"))
  draw_caption(caption_text)
  dev.off()
}

# ==============================================================================
# 3. MAIN
# ==============================================================================

ensure_outdir(OUT_DIR)

cat("===============================================================\n")
cat(" Magic UpSet Plotter\n")
cat("===============================================================\n")

files_to_process <- c()

if (length(TARGET_FILES) > 0) {
  for (f in TARGET_FILES) {
    f_path <- file.path(IN_DIR, f)
    if (file.exists(f_path)) {
      files_to_process <- c(files_to_process, f_path)
    } else {
      cat(sprintf("[WARN] Specified file not found: %s\n", f))
    }
  }
} else {
  cat("Mode: Auto-scan 'input/' for DEG tables...\n")
  all_files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
  files_to_process <- all_files[!grepl("tpm|count|sample|meta", basename(all_files), ignore.case = TRUE)]
}

if (length(files_to_process) < 2) {
  stop("[ERROR] Need at least 2 DEG files.")
}

cat(sprintf("Found %d files. Extracting significant genes...\n", length(files_to_process)))

gene_sets <- list()

for (f_path in files_to_process) {
  fname <- basename(f_path)
  set_name <- file_path_sans_ext(fname)
  set_name <- gsub("_DEG_all", "", set_name)
  set_name <- gsub("_vs_D0", "", set_name)

  df <- tryCatch(
    read_tsv(f_path, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (is.null(df)) {
    cat(sprintf("[SKIP] Failed to read: %s\n", fname))
    next
  }

  if (!all(c("gene_id", "log2fc", "p_adjust") %in% colnames(df))) {
    cat(sprintf("[SKIP] Invalid headers in %s\n", fname))
    next
  }

  sig_genes <- df %>%
    filter(p_adjust < P_CUTOFF & abs(log2fc) > FC_CUTOFF) %>%
    pull(gene_id) %>%
    unique()

  cat(sprintf("  + %s: %d sig genes\n", set_name, length(sig_genes)))

  if (length(sig_genes) > 0) {
    gene_sets[[set_name]] <- sig_genes
  }
}

n_sets <- length(gene_sets)
if (n_sets < 2) {
  stop("[ERROR] Not enough valid gene sets (<2).")
}

caption_text <- NULL
if (SHOW_CAPTION) {
  caption_text <- sprintf("Thresholds: FDR < %s and |log2FC| > %s", P_CUTOFF, FC_CUTOFF)
}

final_title <- NULL
if (SHOW_PLOT_TITLE) {
  if (CUSTOM_PLOT_TITLE != "") {
    final_title <- CUSTOM_PLOT_TITLE
  } else {
    final_title <- "UpSet Analysis"
  }
}

cat("Mode: Royal UpSet Plot\n")

m <- make_comb_mat(gene_sets)
ht <- build_upset_plot(m, final_title = final_title)

out_name <- "UpSet_Combined_Sets"
pdf_file <- file.path(OUT_DIR, paste0(out_name, ".pdf"))
png_file <- file.path(OUT_DIR, paste0(out_name, ".png"))

tryCatch({
  save_upset_plot(ht, pdf_file, png_file, caption_text = caption_text)
  cat(sprintf("Saved PDF: %s\n", pdf_file))
  cat(sprintf("Saved PNG: %s\n", png_file))
}, error = function(e) {
  safe_dev_off()
  stop(sprintf("[ERROR] Save failed: %s", e$message))
})

cat("Exporting Intersection Lists...\n")

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

cat("---------------------------------------------------------------\n")
cat("All processing finished.\n")
