#!/usr/bin/env Rscript
# ==============================================================================
# MAGIC PLOTTING KIT - ogs_bar.R
# [功能] Orthogroups 组成堆叠柱状图 (Figure 1B)
#
# 【方案 B · 全集视角 + 严格 SCO 定义】
#   分母 = 全部 orthogroups 数（矩阵总行数）
#
#   先在全矩阵上定义严格单拷贝家族（Strict Single-copy Orthogroups, SCO）：
#     - 对所有物种列来说，该行的计数都恰好为 1
#     - 即整个 OG 在所有物种中均为 1 拷贝
#
#   对每个物种、每个 Orthogroup 按如下规则划分 4 类（互斥且覆盖）：
#     1. Single-copy
#        - 该 Orthogroup 是严格 SCO（全物种均为 1 拷贝）
#        - 且该物种在此 OG 中 counts_sp > 0（必然为 1）
#     2. Species-specific
#        - 该物种 counts_sp > 0
#        - 且「其他所有物种 counts == 0」
#        - 且该 OG 不是 strict SCO
#        - 即真·物种特有家族
#     3. Multi-copy
#        - 该物种 counts_sp > 0
#        - 且该 OG 不是严格 SCO，也不是 species-specific
#        - 即「所有非 SCO、非特有家族」在该物种上的贡献
#     4. Absent
#        - 该物种 counts_sp == 0
#
#   图 1（OG 视角）：
#     - 对每个物种，用「四类 OG 数 / 全 OG 数」得到百分比，画堆叠条。
#
#   图 2（基因视角）：
#     - 对每个物种，用同一分类，对基因数进行汇总：
#         * Single-copy 基因数 = 严格 SCO 中该物种的基因数（= SCO 家族数）
#         * Species-specific 基因数 = species-specific OG 中该物种 counts_sp 之和
#         * Multi-copy 基因数 = multi-copy OG 中该物种 counts_sp 之和
#       （Absent 对该物种来说没有基因，不再单独画一块）
#
# [输出]
#   若不改 CONFIG：
#     output/ogs_bar/ogs_bar.pdf        （OG 比例图）
#     output/ogs_bar/ogs_bar.png
#     output/ogs_bar/ogs_bar_genes.pdf  （基因数图）
#     output/ogs_bar/ogs_bar_genes.png
#   也可以在 CONFIG 中自定义输出目录名与文件前缀。
# ==============================================================================

# --- A. 用户配置区 (User Configuration) ---
CONFIG <- list(
  # 输入文件（固定为 OrthoFinder 的基因家族计数表）
  INPUT_COUNTS_FILE = "Orthogroups.GeneCount.tsv",

  # 输出目录 & 文件前缀：
  # - 若留空 ""，则自动使用脚本名（去掉 .R），例如 ogs_bar
  # - 若填写非空字符串，则使用自定义值
  OUTPUT_DIR_NAME   = "",   # 例如 "cafe_tree"；留空则用脚本名
  OUTPUT_FILE_STEM  = "",   # 例如 "Fig1B_orthogroups_composition"；留空则用脚本名

  # 标题控制（OG 比例图）
  SHOW_TITLE_OGS    = TRUE,  # TRUE 显示标题，FALSE 不显示标题
  TITLE_TEXT_OGS    = "Orthogroups composition per species (strict SCO-based)",

  # 标题控制（基因数图）
  SHOW_TITLE_GENES  = TRUE,
  TITLE_TEXT_GENES  = "Gene counts per orthogroup category",

  # 物种名字是否缩写为 C.gigas 形式
  USE_SHORT_NAME    = TRUE,

  # 字体与字号
  FONT_FAMILY       = "sans",
  FONT_SIZE_AXIS    = 10,
  FONT_SIZE_TITLE   = 12,
  FONT_SIZE_LEGEND  = 9
)

# 顶刊风格配色（4 类：Single / Multi / Species-specific / Absent）
COLOR_ORTHO <- c(
  "Single-copy"      = "#4DBBD5CC",
  "Multi-copy"       = "#E64B35CC",
  "Species-specific" = "#00A087CC",
  "Absent"           = "#B09C85CC"
)

# --- B. 环境加载 & 路径、脚本名解析 ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)  # 包括 dplyr / tidyr / readr 等
  library(readr)
})

# 获取当前脚本名（从 --file 参数解析），默认 "ogs_bar"
get_script_basename <- function(default = "ogs_bar") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path   <- sub("^--file=", "", file_arg[1])
    base          <- basename(script_path)
    base_no_ext   <- sub("\\.R$", "", base, ignore.case = TRUE)
    if (nzchar(base_no_ext)) return(base_no_ext)
  }
  return(default)
}

SCRIPT_BASENAME <- get_script_basename("ogs_bar")

# 若 CONFIG 中未指定则使用脚本名作为默认
OUTPUT_DIR_NAME  <- CONFIG$OUTPUT_DIR_NAME
if (is.null(OUTPUT_DIR_NAME) || OUTPUT_DIR_NAME == "") {
  OUTPUT_DIR_NAME <- SCRIPT_BASENAME
}

OUTPUT_FILE_STEM <- CONFIG$OUTPUT_FILE_STEM
if (is.null(OUTPUT_FILE_STEM) || OUTPUT_FILE_STEM == "") {
  OUTPUT_FILE_STEM <- SCRIPT_BASENAME
}

MAGIC_ROOT  <- getwd()
INPUT_PATH  <- file.path(MAGIC_ROOT, "input")
OUTPUT_PATH <- file.path(MAGIC_ROOT, "output", OUTPUT_DIR_NAME)

if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)

FILE_COUNTS <- file.path(INPUT_PATH, CONFIG$INPUT_COUNTS_FILE)

if (!file.exists(FILE_COUNTS)) {
  stop(paste0("错误: 找不到 Orthogroups.GeneCount 表：", FILE_COUNTS))
}

message("✅ [magic][ogs_bar] 环境检查通过，开始读取 Orthogroups.GeneCount.tsv ...")

# --- C. 工具函数 ---

# 物种名格式化，和 cafe_tree.R 逻辑保持一致
clean_and_format_label <- function(raw_names, abbreviate = TRUE) {
  sapply(raw_names, function(x) {
    if (is.na(x) || x == "") return(x)
    clean_x <- gsub("'", "", x)
    if (abbreviate) {
      # 例如 Branchiostoma_floridae -> B.floridae
      gsub("^([A-Za-z])[^_]+_([a-zA-Z0-9]+).*", "\\1.\\2", clean_x)
    } else {
      gsub("_", " ", clean_x)
    }
  })
}

# 【全集视角 + 严格 SCO】按物种统计：
#   - 家族数（OG 数）四分类
#   - 基因数（三分类：Single / Multi / Species-specific）
classify_orthogroups_global_strict_sco <- function(df) {
  # 所有物种列（去掉 Orthogroup / Total）
  sp_cols  <- setdiff(names(df), c("Orthogroup", "Total"))
  mat_sp   <- as.matrix(df[, sp_cols, drop = FALSE])

  # 严格 Single-copy Orthogroups：
  # 每一行所有物种计数都恰好为 1
  strict_sco <- apply(mat_sp, 1, function(x) all(x == 1))

  ogs_list   <- list()
  genes_list <- list()
  total_ogs  <- nrow(df)  # 分母 = 全部 OG 数

  for (sp in sp_cols) {
    counts_sp <- df[[sp]]
    others    <- df[, setdiff(sp_cols, sp), drop = FALSE]

    # 当前 OG 在其他物种中是否全部为 0
    others_all_zero <- rowSums(others > 0) == 0

    # 基础状态
    is_absent  <- counts_sp == 0
    is_present <- counts_sp > 0

    # 四类分类（互斥）：
    # 1) Single-copy：严格 SCO 中该物种的贡献
    is_single <- is_present & strict_sco

    # 2) Species-specific：该物种有，其他全部 0，且不是 SCO
    is_spec   <- is_present & others_all_zero & !strict_sco

    # 3) Multi-copy：其余所有 present 且非 SCO、非特有
    is_multi  <- is_present & !is_single & !is_spec

    # 4) Absent：is_absent

    # ---- 家族数（OG 数） ----
    n_single  <- sum(is_single)
    n_multi   <- sum(is_multi)
    n_spec    <- sum(is_spec)
    n_absent  <- sum(is_absent)

    ogs_list[[sp]] <- data.frame(
      Species  = sp,
      Category = c("Single-copy", "Multi-copy", "Species-specific", "Absent"),
      CountOG  = c(n_single, n_multi, n_spec, n_absent)
    )

    # ---- 基因数（只统计有基因的三类）----
    gene_single <- sum(counts_sp[is_single])  # 对 strict SCO 来说等于 n_single
    gene_spec   <- sum(counts_sp[is_spec])
    gene_multi  <- sum(counts_sp[is_multi])

    genes_list[[sp]] <- data.frame(
      Species   = sp,
      Category  = c("Single-copy", "Multi-copy", "Species-specific"),
      GeneCount = c(gene_single, gene_multi, gene_spec)
    )
  }

  list(
    ogs   = bind_rows(ogs_list),
    genes = bind_rows(genes_list),
    total_ogs = total_ogs
  )
}

# --- D. 数据读取与整理 ---

counts_raw <- read_tsv(FILE_COUNTS, show_col_types = FALSE)

message("ℹ️  [magic][ogs_bar] 读取到矩阵维度：",
        nrow(counts_raw), " orthogroups × ",
        ncol(counts_raw) - 2, " species.")

class_res      <- classify_orthogroups_global_strict_sco(counts_raw)
bar_ogs_raw    <- class_res$ogs
bar_genes_raw  <- class_res$genes
total_ogs      <- class_res$total_ogs

# 物种名格式化
bar_ogs_raw$Species_raw   <- bar_ogs_raw$Species
bar_genes_raw$Species_raw <- bar_genes_raw$Species

# 保持物种顺序与原表列顺序一致（与树对应），倒序以便从上到下顺序一致
species_order_raw <- setdiff(names(counts_raw), c("Orthogroup", "Total"))
species_order_fmt <- clean_and_format_label(
  species_order_raw,
  abbreviate = CONFIG$USE_SHORT_NAME
)

bar_ogs_raw$Species <- clean_and_format_label(
  bar_ogs_raw$Species_raw,
  abbreviate = CONFIG$USE_SHORT_NAME
)
bar_genes_raw$Species <- clean_and_format_label(
  bar_genes_raw$Species_raw,
  abbreviate = CONFIG$USE_SHORT_NAME
)

bar_ogs_raw$Species <- factor(
  bar_ogs_raw$Species,
  levels = rev(species_order_fmt)
)
bar_genes_raw$Species <- factor(
  bar_genes_raw$Species,
  levels = rev(species_order_fmt)
)

# ---- OG 比例图数据 ----
bar_ogs <- bar_ogs_raw %>%
  mutate(Percent = CountOG / total_ogs * 100)

bar_ogs$Category <- factor(
  bar_ogs$Category,
  levels = c("Single-copy", "Multi-copy", "Species-specific", "Absent")
)

# ---- 基因数图数据 ----
bar_genes <- bar_genes_raw
bar_genes$Category <- factor(
  bar_genes$Category,
  levels = c("Single-copy", "Multi-copy", "Species-specific")
)

# --- E. 绘图 ---

## 1) Orthogroups 比例图（按家族数）
message("🎨 [magic][ogs_bar] 绘制 Orthogroups 组成堆叠图 (strict SCO-based, 全 OG 视角) ...")

plot_title_ogs <- if (isTRUE(CONFIG$SHOW_TITLE_OGS)) CONFIG$TITLE_TEXT_OGS else NULL

p_bar_ogs <- ggplot(bar_ogs, aes(x = Percent, y = Species, fill = Category)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = COLOR_ORTHO, name = "Orthogroups") +
  scale_x_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    x = "Proportion of orthogroups (%)",
    y = NULL,
    title = plot_title_ogs
  ) +
  theme_classic(base_family = CONFIG$FONT_FAMILY) +
  theme(
    axis.text.x  = element_text(
      size   = CONFIG$FONT_SIZE_AXIS,
      colour = "black"
    ),
    axis.text.y  = element_text(
      size   = CONFIG$FONT_SIZE_AXIS,
      face   = "italic",
      colour = "black"
    ),
    axis.title.x = element_text(
      size   = CONFIG$FONT_SIZE_AXIS + 1,
      colour = "black"
    ),
    plot.title   = if (isTRUE(CONFIG$SHOW_TITLE_OGS)) {
      element_text(
        size   = CONFIG$FONT_SIZE_TITLE,
        face   = "bold",
        hjust  = 0.5
      )
    } else {
      element_blank()
    },
    legend.title = element_text(
      size = CONFIG$FONT_SIZE_LEGEND,
      face = "bold"
    ),
    legend.text  = element_text(size = CONFIG$FONT_SIZE_LEGEND),
    legend.position   = "right",
    legend.box.margin = margin(t = 5, r = 5, b = 5, l = 5),
    plot.margin       = margin(t = 10, r = 20, b = 10, l = 10)
  )

## 2) 基因数构成图
message("🎨 [magic][ogs_bar] 绘制基因数构成柱状图 (Single/Multi/Species-specific) ...")

plot_title_genes <- if (isTRUE(CONFIG$SHOW_TITLE_GENES)) CONFIG$TITLE_TEXT_GENES else NULL

p_bar_genes <- ggplot(bar_genes, aes(x = GeneCount, y = Species, fill = Category)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = COLOR_ORTHO[c("Single-copy", "Multi-copy", "Species-specific")],
    name   = "Genes"
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  labs(
    x = "Number of genes",
    y = NULL,
    title = plot_title_genes
  ) +
  theme_classic(base_family = CONFIG$FONT_FAMILY) +
  theme(
    axis.text.x  = element_text(
      size   = CONFIG$FONT_SIZE_AXIS,
      colour = "black"
    ),
    axis.text.y  = element_text(
      size   = CONFIG$FONT_SIZE_AXIS,
      face   = "italic",
      colour = "black"
    ),
    axis.title.x = element_text(
      size   = CONFIG$FONT_SIZE_AXIS + 1,
      colour = "black"
    ),
    plot.title   = if (isTRUE(CONFIG$SHOW_TITLE_GENES)) {
      element_text(
        size   = CONFIG$FONT_SIZE_TITLE,
        face   = "bold",
        hjust  = 0.5
      )
    } else {
      element_blank()
    },
    legend.title = element_text(
      size = CONFIG$FONT_SIZE_LEGEND,
      face = "bold"
    ),
    legend.text  = element_text(size = CONFIG$FONT_SIZE_LEGEND),
    legend.position   = "right",
    legend.box.margin = margin(t = 5, r = 5, b = 5, l = 5),
    plot.margin       = margin(t = 10, r = 20, b = 10, l = 10)
  )

# --- F. 保存 ---

file_pdf_ogs   <- file.path(OUTPUT_PATH,
                            paste0(OUTPUT_FILE_STEM, ".pdf"))
file_png_ogs   <- file.path(OUTPUT_PATH,
                            paste0(OUTPUT_FILE_STEM, ".png"))
file_pdf_genes <- file.path(OUTPUT_PATH,
                            paste0(OUTPUT_FILE_STEM, "_genes.pdf"))
file_png_genes <- file.path(OUTPUT_PATH,
                            paste0(OUTPUT_FILE_STEM, "_genes.png"))

ggsave(file_pdf_ogs,   p_bar_ogs,   width = 8, height = 6, device = cairo_pdf)
ggsave(file_png_ogs,   p_bar_ogs,   width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file_pdf_genes, p_bar_genes, width = 8, height = 6, device = cairo_pdf)
ggsave(file_png_genes, p_bar_genes, width = 8, height = 6, dpi = 300, bg = "white")

message("✅ [magic][ogs_bar] 已生成：", file_pdf_ogs)
message("✅ [magic][ogs_bar] 已生成：", file_png_ogs)
message("✅ [magic][ogs_bar] 已生成：", file_pdf_genes)
message("✅ [magic][ogs_bar] 已生成：", file_png_genes)
message("✨ [magic][ogs_bar] Orthogroups 组成（家族数 + 基因数）双视角图绘制完成。")

