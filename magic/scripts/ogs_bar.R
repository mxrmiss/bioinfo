#!/usr/bin/env Rscript
# ==============================================================================
# MAGIC PLOTTING KIT - ogs_bar_evolution_full.R
# [功能] 绘制 Figure 1B 全套：Orthogroups 占比图 + 基因数堆叠图
#
# 【分类方案 · 5 类家族演化视角】
#   基于 OrthoFinder 的两个文件：
#     1. Orthogroups.GeneCount.tsv (聚类成功的家族)
#     2. Orthogroups_UnassignedGenes.tsv (未聚类的孤儿基因)
#
#   我们将数据划分为 5+2 种逻辑状态，用于两张不同的图：
#
#   [核心分类定义]
#     A. Single-copy orthologs: 所有物种都在，且全是 1 拷贝。
#     B. Multiple-copy orthologs: 所有物种都在，至少有一个 > 1。
#     C. Unique paralogs: 仅 1 个物种有 (Species-specific)。
#     D. Other orthologs: 物种数在 (1, 总数) 之间 (Patchy)。
#
#   [图 1: Orthogroups 占比图] (单位: 家族数)
#     - 包含: A, B, C, D
#     - 补充: Absent (对于 C 和 D，如果某物种没有，则记为 Absent)
#     * 不包含 Unclustered (因为它们不是家族)
#
#   [图 2: Gene 计数图] (单位: 基因数)
#     - 包含: A, B, C, D
#     - 补充: Unclustered genes (从 UnassignedGenes.tsv 读取)
#     * 不包含 Absent (因为基因数为 0)
#
# [输出]
#   1. output/.../ogs_bar_evo_OGs.pdf/.png   (占比图)
#   2. output/.../ogs_bar_evo_Genes.pdf/.png (基因数图)
# ==============================================================================

# --- A. 用户配置区 (User Configuration) ---
CONFIG <- list(
  # 输入文件
  INPUT_COUNTS_FILE     = "Orthogroups.GeneCount.tsv",
  INPUT_UNASSIGNED_FILE = "Orthogroups_UnassignedGenes.tsv",

  # 输出目录 & 文件前缀
  OUTPUT_DIR_NAME   = "ogs_bar_evo_full",
  OUTPUT_FILE_STEM  = "ogs_bar_evolution",

  # 标题控制
  SHOW_TITLE_OGS    = TRUE,
  TITLE_TEXT_OGS    = "Orthogroups composition per species",
  
  SHOW_TITLE_GENES  = TRUE,
  TITLE_TEXT_GENES  = "Gene counts per orthogroup category",

  # 物种名字是否缩写
  USE_SHORT_NAME    = TRUE,

  # 字体与字号
  FONT_FAMILY       = "sans",
  FONT_SIZE_AXIS    = 10,
  FONT_SIZE_TITLE   = 12,
  FONT_SIZE_LEGEND  = 9
)

# 统一配色方案（皇家马卡龙配色版；保留 CC 透明度）
COLOR_MAP <- c(
  "Single-copy orthologs"    = "#95C8F2CC", # Royal SkyBlue
  "Multiple-copy orthologs"  = "#F79D93CC", # Royal Coral
  "Unique paralogs"          = "#9FD5CBCC", # Royal Mint
  "Other orthologs"          = "#A99BEFCC", # Royal Lavender
  "Absent"                   = "#FBE3A8CC", # Royal Cream (原本灰色→改为温柔浅黄)
  "Unclustered genes"        = "#F6CD96CC"  # Royal Light Orange
)

# --- B. 环境加载 & 路径解析 ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(readr)
})

MAGIC_ROOT  <- getwd()
INPUT_PATH  <- file.path(MAGIC_ROOT, "input")
OUTPUT_PATH <- file.path(MAGIC_ROOT, "output", CONFIG$OUTPUT_DIR_NAME)

if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)

FILE_COUNTS     <- file.path(INPUT_PATH, CONFIG$INPUT_COUNTS_FILE)
FILE_UNASSIGNED <- file.path(INPUT_PATH, CONFIG$INPUT_UNASSIGNED_FILE)

if (!file.exists(FILE_COUNTS)) stop(paste0("错误: 找不到文件 ", FILE_COUNTS))
if (!file.exists(FILE_UNASSIGNED)) stop(paste0("错误: 找不到文件 ", FILE_UNASSIGNED))

message("✅ [magic] 环境检查通过，准备开始处理数据...")

# --- C. 工具函数 ---

clean_and_format_label <- function(raw_names, abbreviate = TRUE) {
  sapply(raw_names, function(x) {
    if (is.na(x) || x == "") return(x)
    clean_x <- gsub("'", "", x)
    if (abbreviate) {
      gsub("^([A-Za-z])[^_]+_([a-zA-Z0-9]+).*", "\\1.\\2", clean_x)
    } else {
      gsub("_", " ", clean_x)
    }
  })
}

# --- D. 核心处理逻辑 ---

process_data <- function(df_counts, df_unassigned) {
  
  # 1. 准备矩阵
  sp_cols <- setdiff(names(df_counts), c("Orthogroup", "Total"))
  n_sp    <- length(sp_cols)
  mat_counts <- as.matrix(df_counts[, sp_cols])
  total_ogs  <- nrow(df_counts)
  
  # 2. 对 OG 进行分类打标 (Vectorized)
  n_present <- rowSums(mat_counts > 0)
  max_val   <- apply(mat_counts, 1, max)
  
  # A. Single-copy: 全员存在且最大拷贝为1
  is_single <- (n_present == n_sp) & (max_val == 1)
  # B. Multi-copy: 全员存在且最大拷贝>1
  is_multi  <- (n_present == n_sp) & (max_val > 1)
  # C. Unique: 仅1个物种存在
  is_unique <- (n_present == 1)
  # D. Other: 介于中间
  is_other  <- (n_present > 1) & (n_present < n_sp)
  
  # 结果容器
  list_ogs   <- list() # 存图1数据
  list_genes <- list() # 存图2数据
  
  # 3. 遍历物种计算
  for (sp in sp_cols) {
    col_data <- df_counts[[sp]]
    
    # === 计算图 1 数据 (OGs 数量) ===
    # 注意：Category 是指该 OG 的全局属性。
    # 但如果该物种在这个 OG 里没基因 (count==0)，那对他来说就是 Absent。
    
    # 3.1 统计各状态下，该物种 Count > 0 的家族数
    n_og_single <- sum(is_single & (col_data > 0)) # 恒等于 sum(is_single)
    n_og_multi  <- sum(is_multi  & (col_data > 0)) # 恒等于 sum(is_multi)
    n_og_unique <- sum(is_unique & (col_data > 0)) # 只统计属于自己的 Unique
    n_og_other  <- sum(is_other  & (col_data > 0)) # 统计自己参与了的 Other
    
    # 3.2 剩下的就是 Absent
    n_og_absent <- total_ogs - (n_og_single + n_og_multi + n_og_unique + n_og_other)
    
    list_ogs[[sp]] <- data.frame(
      Species = sp,
      Category = c("Single-copy orthologs", "Multiple-copy orthologs", 
                   "Unique paralogs", "Other orthologs", "Absent"),
      CountOG = c(n_og_single, n_og_multi, n_og_unique, n_og_other, n_og_absent)
    )
    
    # === 计算图 2 数据 (Gene 数量) ===
    
    # 3.3 统计 GeneCount 表中的基因
    n_gene_single <- sum(col_data[is_single])
    n_gene_multi  <- sum(col_data[is_multi])
    n_gene_unique <- sum(col_data[is_unique])
    n_gene_other  <- sum(col_data[is_other])
    
    # 3.4 统计 UnassignedGenes 表中的孤儿基因
    #     需清洗 "|" 及后续内容
    if (sp %in% names(df_unassigned)) {
        raw_un <- df_unassigned[[sp]]
        clean_un <- gsub("\\|.*", "", ifelse(is.na(raw_un), "", raw_un))
        n_gene_unclustered <- sum(clean_un != "")
    } else {
        n_gene_unclustered <- 0
    }

    list_genes[[sp]] <- data.frame(
      Species = sp,
      Category = c("Single-copy orthologs", "Multiple-copy orthologs", 
                   "Unique paralogs", "Other orthologs", "Unclustered genes"),
      GeneCount = c(n_gene_single, n_gene_multi, n_gene_unique, 
                    n_gene_other, n_gene_unclustered)
    )
  }
  
  return(list(
    ogs   = bind_rows(list_ogs),
    genes = bind_rows(list_genes),
    total_ogs = total_ogs
  ))
}

# --- E. 数据读取与处理 ---

message("⏳ [magic] 正在读取数据...")
counts_raw     <- read_tsv(FILE_COUNTS, show_col_types = FALSE)
unassigned_raw <- read_tsv(FILE_UNASSIGNED, show_col_types = FALSE)

res <- process_data(counts_raw, unassigned_raw)
plot_data_ogs   <- res$ogs
plot_data_genes <- res$genes
total_ogs       <- res$total_ogs

# --- F. 格式化 (排序与标签) ---

# 1. 确定物种顺序 (按原表顺序倒序)
species_order_raw <- setdiff(names(counts_raw), c("Orthogroup", "Total"))
species_order_fmt <- clean_and_format_label(species_order_raw, abbreviate = CONFIG$USE_SHORT_NAME)

# 2. 应用格式化
plot_data_ogs$Species_Label   <- clean_and_format_label(plot_data_ogs$Species, abbreviate = CONFIG$USE_SHORT_NAME)
plot_data_genes$Species_Label <- clean_and_format_label(plot_data_genes$Species, abbreviate = CONFIG$USE_SHORT_NAME)

plot_data_ogs$Species_Label   <- factor(plot_data_ogs$Species_Label, levels = rev(species_order_fmt))
plot_data_genes$Species_Label <- factor(plot_data_genes$Species_Label, levels = rev(species_order_fmt))

# 3. 设置 Category 的堆叠顺序 (Levels)
# 图 1 顺序: Absent 最上/下? 通常 Absent 放顶层或底层。这里放顶层。
levels_ogs <- c("Absent", "Other orthologs", "Unique paralogs", 
                "Multiple-copy orthologs", "Single-copy orthologs")
plot_data_ogs$Category <- factor(plot_data_ogs$Category, levels = levels_ogs)

# 图 2 顺序: Unclustered 放顶层
levels_genes <- c("Unclustered genes", "Other orthologs", "Unique paralogs", 
                  "Multiple-copy orthologs", "Single-copy orthologs")
plot_data_genes$Category <- factor(plot_data_genes$Category, levels = levels_genes)

# 计算百分比 (图1专用)
plot_data_ogs <- plot_data_ogs %>%
  mutate(Percent = CountOG / total_ogs * 100)

# --- G. 绘图 ---

# === 绘制图 1: Orthogroups Composition (百分比) ===
message("🎨 [magic] 绘制图 1: Orthogroups 占比图...")

p1 <- ggplot(plot_data_ogs, aes(x = Percent, y = Species_Label, fill = Category)) +
  geom_col(width = 0.75, position = position_stack()) +
  scale_fill_manual(values = COLOR_MAP) +
  # 修复警告: 稍微放宽 limits 防止浮点数被切，并用 expand 贴边
  scale_x_continuous(limits = c(0, 100.1), expand = c(0, 0), breaks = seq(0, 100, 20)) +
  labs(
    x = "Proportion of orthogroups (%)",
    y = NULL,
    title = if (isTRUE(CONFIG$SHOW_TITLE_OGS)) CONFIG$TITLE_TEXT_OGS else NULL
  ) +
  theme_classic(base_family = CONFIG$FONT_FAMILY) +
  theme(
    axis.text.x  = element_text(size = CONFIG$FONT_SIZE_AXIS, color = "black"),
    axis.text.y  = element_text(size = CONFIG$FONT_SIZE_AXIS, face = "italic", color = "black"),
    axis.title.x = element_text(size = CONFIG$FONT_SIZE_AXIS + 1, color = "black"),
    plot.title   = element_text(size = CONFIG$FONT_SIZE_TITLE, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text     = element_text(size = CONFIG$FONT_SIZE_LEGEND),
    legend.title    = element_blank(), # 去掉 Category 标题更清爽
    plot.margin     = margin(t = 10, r = 20, b = 10, l = 10)
  )

# === 绘制图 2: Gene Counts (绝对数量) ===
message("🎨 [magic] 绘制图 2: Gene 数量堆叠图...")

p2 <- ggplot(plot_data_genes, aes(x = GeneCount, y = Species_Label, fill = Category)) +
  geom_col(width = 0.75, position = position_stack()) +
  scale_fill_manual(values = COLOR_MAP) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    x = "Number of genes",
    y = NULL,
    title = if (isTRUE(CONFIG$SHOW_TITLE_GENES)) CONFIG$TITLE_TEXT_GENES else NULL
  ) +
  theme_classic(base_family = CONFIG$FONT_FAMILY) +
  theme(
    axis.text.x  = element_text(size = CONFIG$FONT_SIZE_AXIS, color = "black"),
    axis.text.y  = element_text(size = CONFIG$FONT_SIZE_AXIS, face = "italic", color = "black"),
    axis.title.x = element_text(size = CONFIG$FONT_SIZE_AXIS + 1, color = "black"),
    plot.title   = element_text(size = CONFIG$FONT_SIZE_TITLE, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text     = element_text(size = CONFIG$FONT_SIZE_LEGEND),
    legend.title    = element_blank(),
    plot.margin     = margin(t = 10, r = 20, b = 10, l = 10)
  )

# --- H. 保存输出 ---

# 保存图 1
f1_pdf <- file.path(OUTPUT_PATH, paste0(CONFIG$OUTPUT_FILE_STEM, "_OGs.pdf"))
f1_png <- file.path(OUTPUT_PATH, paste0(CONFIG$OUTPUT_FILE_STEM, "_OGs.png"))
ggsave(f1_pdf, p1, width = 9, height = 6, device = cairo_pdf)
ggsave(f1_png, p1, width = 9, height = 6, dpi = 300, bg = "white")

# 保存图 2
f2_pdf <- file.path(OUTPUT_PATH, paste0(CONFIG$OUTPUT_FILE_STEM, "_Genes.pdf"))
f2_png <- file.path(OUTPUT_PATH, paste0(CONFIG$OUTPUT_FILE_STEM, "_Genes.png"))
ggsave(f2_pdf, p2, width = 9, height = 6, device = cairo_pdf)
ggsave(f2_png, p2, width = 9, height = 6, dpi = 300, bg = "white")

message("✅ [magic] 已生成图1 (OGs): ", f1_pdf)
message("✅ [magic] 已生成图2 (Genes): ", f2_pdf)
message("✨ [magic] 全部绘图任务完成！")
