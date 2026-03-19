#!/usr/bin/env Rscript
# =============================================================================
# plot_19gene_heatmap.R
# -----------------------------------------------------------------------------
# 功能：
# 1. 读取 19genes.tsv（基因注释表）和 19genes_tpm.tsv（表达矩阵）
# 2. 通过 Gene_ID 匹配表达矩阵
# 3. 绘制“列按组织分组 + 行按 Module 分组，但不聚类”的 grouped heatmap
# 4. 保留重复样本，但不显示重复样本名称；每个组织仅显示一次组织名
# 5. 行名显示 Gene_symbol，放在右侧
# 6. 左侧不显示 Module 文字标签，而使用颜色条表示 Module 分组
# 7. 顶部不显示组织颜色条，也不显示组织分组图例
# 8. 同时输出高分辨率 PNG 和矢量 PDF，适合论文作图
#
# 说明：
# - 不使用命令行传参；所有参数都在下方“用户参数设置区”统一修改
# - 默认目标字体为 Arial；若系统不可用，则自动回退到更稳妥的 sans 字体链
# - 字体策略走“路线1”：强制使用 cairo / ragg 设备链，不走基础 pdf()
# - 运行时会输出简明屏幕信息，便于判断进度与报错位置
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# =============================================================================
#                           用户参数设置区（皇上改这里）
# =============================================================================

# ---------------------------- 1. 输入输出路径设置 ----------------------------
# 基因注释表（必须为制表符分隔）
# 需要至少包含以下列：
# Module, Gene_symbol, Gene_ID
GENE_INFO_TSV <- "19genes.tsv"

# 表达矩阵表（必须为制表符分隔）
# 第一列应为 Gene_ID；后面各列是表达量（如 TPM）
# 允许列名重复，例如 foot foot foot gill gill ...
EXPR_TSV <- "19genes_tpm.tsv"

# 输出目录
OUTDIR <- "heatmap_out"

# 输出文件名前缀
OUT_PREFIX <- "Sc_19genes_grouped_heatmap"


# ---------------------------- 2. 输出文件类型 ----------------------------
# 是否输出 PNG
OUTPUT_PNG <- TRUE

# 是否输出 PDF
OUTPUT_PDF <- TRUE


# ---------------------------- 3. 图片尺寸与清晰度 ----------------------------
# PNG 宽度（单位：英寸）
PNG_WIDTH <- 10

# PNG 高度（单位：英寸）
PNG_HEIGHT <- 9

# PNG 分辨率（单位：dpi）
# 一般论文图推荐 300 或 600
PNG_DPI <- 600

# PDF 宽度（单位：英寸）
PDF_WIDTH <- 10

# PDF 高度（单位：英寸）
PDF_HEIGHT <- 9


# ---------------------------- 4. 标题设置 ----------------------------
# 是否显示主标题
SHOW_MAIN_TITLE <- FALSE

# 主标题文本
MAIN_TITLE <- "Grouped heatmap of 19 candidate genes across six tissues"

# 主标题字号
MAIN_TITLE_FONT_SIZE <- 14

# 主标题是否加粗
MAIN_TITLE_BOLD <- TRUE


# ---------------------------- 5. 字体设置 ----------------------------
# 默认目标字体
# 说明：
# 1）脚本会优先尝试检测系统中是否有 Arial
# 2）若找不到，则自动回退到 Liberation Sans / DejaVu Sans / sans
# 3）本脚本只走 cairo / ragg 路线，不使用基础 pdf() 设备
FONT_FAMILY_TARGET <- "sans"

# 行标签（Gene_symbol）字号
ROW_NAME_FONT_SIZE <- 16

# 行标签是否加粗
ROW_NAME_BOLD <- FALSE

# 顶部组织名称字号（foot / gill / intest / lips / mantle / siphon）
COLUMN_GROUP_TITLE_FONT_SIZE <- 18

# 顶部组织名称是否加粗
COLUMN_GROUP_TITLE_BOLD <- FALSE

# 左侧 Module 颜色条标题字号
# 当前默认不显示标题，仅保留颜色条本体；此参数保留给将来扩展
MODULE_BAR_TITLE_FONT_SIZE <- 10

# 图例标题字号
LEGEND_TITLE_FONT_SIZE <- 16

# 图例标签字号
LEGEND_LABEL_FONT_SIZE <- 14


# ---------------------------- 6. 热图数值设置 ----------------------------
# 是否按行做 z-score
# TRUE  = 每个基因在所有样本中做标准化，更适合看组织间相对模式
# FALSE = 使用原始表达值作图
USE_ROW_ZSCORE <- TRUE

# z-score 截断下限
ZSCORE_MIN <- -2

# z-score 截断上限
ZSCORE_MAX <- 2

# 热图图例标题
HEATMAP_LEGEND_TITLE <- "z-score"


# ---------------------------- 7. 行列结构设置 ----------------------------
# 行名显示在右侧（按皇上要求）
ROW_NAMES_SIDE <- "right"

# 图例放左侧（按皇上要求）
HEATMAP_LEGEND_SIDE <- "left"
ANNOTATION_LEGEND_SIDE <- "left"

# 是否显示左侧 Module 颜色条
SHOW_MODULE_COLOR_BAR <- TRUE

# 左侧 Module 颜色条宽度（单位：mm）
MODULE_COLOR_BAR_WIDTH_MM <- 4

# 是否显示模块之间的分隔
SHOW_MODULE_SEPARATORS <- TRUE

# 模块之间的间距（单位：mm）
ROW_GAP_MM <- 1.5

# 不同组织组之间的间距（单位：mm）
COLUMN_GAP_MM <- 2.5

# 是否显示单元格边框
SHOW_CELL_BORDER <- FALSE

# 单元格边框颜色
CELL_BORDER_COLOR <- "white"

# 单元格边框线宽
CELL_BORDER_LWD <- 0.4

# 单个色块缩放比例（0~1之间）
# 数值越小，格子之间留白越明显，越像上传图片中的“独立圆角方块”
TILE_SHRINK <- 0.85

# 圆角半径（单位：mm）
# 数值越大，色块四角越圆润
ROUNDED_TILE_RADIUS_MM <- 1.8


# ---------------------------- 8. 行列分组顺序设置 ----------------------------
# 组织顺序（必须与表达矩阵列名中的组织名称一致）
TISSUE_ORDER <- c("foot", "gill", "intest", "lips", "mantle", "siphon")

# Module 顺序
MODULE_ORDER <- c(
  "Neural regulation / membrane excitability",
  "Contractile execution / calcium-associated regulation",
  "Extracellular interface / ECM-related interaction",
  "Glycosylation / surface modification / secretion-related",
  "Metabolic support / energetic support"
)

# 是否严格按照 19genes.tsv 当前顺序显示每个模块内部基因
# TRUE = 完全按文件中的顺序
# FALSE = 仍不聚类，但可按 Gene_symbol 字母排序
KEEP_GENE_ORDER_AS_INPUT <- TRUE


# ---------------------------- 9. 配色设置 ----------------------------
# 热图低值 / 中间 / 高值颜色
HEATMAP_COLORS <- c("#3B8BC2", "white", "#E76F51")

# 左侧 Module 颜色条颜色
# 注意：
# 这里不再显示文字标签，也不显示图例
# 因此颜色需要在图注中解释
MODULE_COLORS <- c(
  "Neural regulation / membrane excitability" = "#F79D93",
  "Contractile execution / calcium-associated regulation" = "#95C8F2",
  "Extracellular interface / ECM-related interaction" = "#F6CD96",
  "Glycosylation / surface modification / secretion-related" = "#9FD5CB",
  "Metabolic support / energetic support" = "#C3A3F5"
)


# ---------------------------- 10. 版式与留白 ----------------------------
# 热图主体四周留白（单位：mm）
# 顺序：上、右、下、左
PLOT_PADDING_MM <- c(4, 4, 4, 4)

# 热图对象内部名称（一般无需修改）
HEATMAP_NAME <- "expr"


suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")
stop_with <- function(...) stop(sprintf(...), call. = FALSE)
`%||%` <- function(a, b) if (!is.null(a)) a else b
make_bold_face <- function(flag) if (isTRUE(flag)) "bold" else "plain"

row_zscore <- function(mat) {
  z <- t(apply(mat, 1, function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

clamp_matrix <- function(mat, low, high) {
  mat[mat < low] <- low
  mat[mat > high] <- high
  mat
}

open_png <- function(path, w, h, res) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = w, height = h, units = "in", res = res)
  } else {
    png(path, width = w * res, height = h * res, res = res, type = "cairo")
  }
}

open_pdf <- function(path, w, h) {
  if (capabilities("cairo")) {
    cairo_pdf(path, width = w, height = h)
  } else {
    pdf(path, width = w, height = h, useDingbats = FALSE)
  }
}

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(GENE_INFO_TSV)) stop_with("找不到文件：%s", GENE_INFO_TSV)
if (!file.exists(EXPR_TSV)) stop_with("找不到文件：%s", EXPR_TSV)

msg("Read input")

gene_info <- tryCatch(
  read.delim(GENE_INFO_TSV, header = TRUE, sep = "\t", check.names = FALSE, quote = ""),
  error = function(e) stop_with("读取基因注释表失败：%s", e$message)
)

req_cols <- c("Module", "Gene_symbol", "Gene_ID")
miss_cols <- setdiff(req_cols, colnames(gene_info))
if (length(miss_cols) > 0) stop_with("19genes.tsv 缺少必要列：%s", paste(miss_cols, collapse = ", "))

gene_info$Module <- trimws(gene_info$Module)
gene_info$Gene_symbol <- trimws(gene_info$Gene_symbol)
gene_info$Gene_ID <- trimws(gene_info$Gene_ID)

if (anyDuplicated(gene_info$Gene_ID) > 0) {
  dup_ids <- unique(gene_info$Gene_ID[duplicated(gene_info$Gene_ID)])
  stop_with("19genes.tsv 中存在重复 Gene_ID：%s", paste(dup_ids, collapse = ", "))
}

if (!all(gene_info$Module %in% MODULE_ORDER)) {
  bad <- unique(gene_info$Module[!gene_info$Module %in% MODULE_ORDER])
  stop_with("19genes.tsv 中存在未纳入 MODULE_ORDER 的 Module 名称：%s", paste(bad, collapse = " | "))
}

gene_info$Module <- factor(gene_info$Module, levels = MODULE_ORDER)
gene_info$.__ord__ <- seq_len(nrow(gene_info))

if (isTRUE(KEEP_GENE_ORDER_AS_INPUT)) {
  gene_info <- gene_info[order(gene_info$Module, gene_info$.__ord__), , drop = FALSE]
} else {
  gene_info <- gene_info[order(gene_info$Module, gene_info$Gene_symbol), , drop = FALSE]
}
gene_info$.__ord__ <- NULL

expr_df <- tryCatch(
  read.delim(EXPR_TSV, header = TRUE, sep = "\t", check.names = FALSE, quote = ""),
  error = function(e) stop_with("读取表达矩阵失败：%s", e$message)
)

if (ncol(expr_df) < 2) stop_with("表达矩阵至少应包含 2 列：第一列 Gene_ID + 至少 1 列表达值")

colnames(expr_df)[1] <- "Gene_ID"
expr_df$Gene_ID <- trimws(expr_df$Gene_ID)

sample_tissues <- trimws(colnames(expr_df)[-1])
if (!all(sample_tissues %in% TISSUE_ORDER)) {
  bad_tissues <- unique(sample_tissues[!sample_tissues %in% TISSUE_ORDER])
  stop_with("表达矩阵中检测到未在 TISSUE_ORDER 中定义的组织名称：%s", paste(bad_tissues, collapse = ", "))
}

expr_mat_raw <- as.matrix(expr_df[, -1, drop = FALSE])
mode(expr_mat_raw) <- "character"
expr_mat_num <- suppressWarnings(apply(expr_mat_raw, 2, as.numeric))
expr_mat_num <- as.matrix(expr_mat_num)

if (!all(is.finite(expr_mat_num) | is.na(expr_mat_num))) {
  stop_with("表达矩阵中存在无法转换为数值的内容，请检查是否混入非数字字符")
}

rownames(expr_mat_num) <- expr_df$Gene_ID
colnames(expr_mat_num) <- sample_tissues

idx <- match(gene_info$Gene_ID, rownames(expr_mat_num))
if (any(is.na(idx))) {
  missing_ids <- gene_info$Gene_ID[is.na(idx)]
  stop_with("以下 Gene_ID 在表达矩阵第一列中找不到对应行：%s", paste(missing_ids, collapse = ", "))
}

expr_sub <- expr_mat_num[idx, , drop = FALSE]
rownames(expr_sub) <- gene_info$Gene_ID
if (anyNA(expr_sub)) expr_sub[is.na(expr_sub)] <- 0

if (isTRUE(USE_ROW_ZSCORE)) {
  plot_mat <- clamp_matrix(row_zscore(expr_sub), ZSCORE_MIN, ZSCORE_MAX)
  col_fun <- circlize::colorRamp2(c(ZSCORE_MIN, 0, ZSCORE_MAX), HEATMAP_COLORS)
  legend_at <- c(ZSCORE_MIN, 0, ZSCORE_MAX)
  legend_title <- HEATMAP_LEGEND_TITLE
} else {
  plot_mat <- expr_sub
  rng <- range(plot_mat, na.rm = TRUE)
  if (rng[1] == rng[2]) rng <- c(rng[1] - 1, rng[2] + 1)
  col_fun <- circlize::colorRamp2(c(rng[1], mean(rng), rng[2]), HEATMAP_COLORS)
  legend_at <- pretty(rng, n = 5)
  legend_title <- "expression"
}

col_split <- factor(colnames(plot_mat), levels = TISSUE_ORDER)
ord_col <- order(col_split)
plot_mat <- plot_mat[, ord_col, drop = FALSE]
col_split <- factor(colnames(plot_mat), levels = TISSUE_ORDER)

row_split <- factor(gene_info$Module, levels = MODULE_ORDER)
row_title_opt <- rep("", length(levels(row_split)))
row_labels <- gene_info$Gene_symbol

left_ha <- NULL
if (isTRUE(SHOW_MODULE_COLOR_BAR)) {
  mod_chr <- as.character(gene_info$Module)
  if (!all(mod_chr %in% names(MODULE_COLORS))) {
    bad <- unique(mod_chr[!mod_chr %in% names(MODULE_COLORS)])
    stop_with("以下 Module 未在 MODULE_COLORS 中配置颜色：%s", paste(bad, collapse = " | "))
  }

  left_ha <- rowAnnotation(
    module = mod_chr,
    col = list(module = MODULE_COLORS),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    border = FALSE,
    width = unit(MODULE_COLOR_BAR_WIDTH_MM, "mm"),
    simple_anno_size = unit(MODULE_COLOR_BAR_WIDTH_MM, "mm")
  )
}

build_ht <- function() {
  ht <- Heatmap(
    plot_mat,
    name = HEATMAP_NAME,
    col = col_fun,

    cluster_rows = FALSE,
    cluster_columns = FALSE,

    row_split = row_split,
    row_gap = if (isTRUE(SHOW_MODULE_SEPARATORS)) unit(ROW_GAP_MM, "mm") else unit(0, "mm"),
    row_title = row_title_opt,

    column_split = col_split,
	column_gap = unit(COLUMN_GAP_MM, "mm"),

    left_annotation = left_ha,

    show_row_names = TRUE,
    show_column_names = FALSE,
    row_labels = row_labels,
    row_names_side = ROW_NAMES_SIDE,
    row_names_gp = gpar(
      fontsize = as.numeric(ROW_NAME_FONT_SIZE %||% 8),
      fontfamily = FONT_FAMILY_TARGET,
      fontface = make_bold_face(ROW_NAME_BOLD)
    ),

    column_title = "%s",
    column_title_side = "top",
    column_title_gp = gpar(
      fontsize = as.numeric(COLUMN_GROUP_TITLE_FONT_SIZE %||% 8),
      fontfamily = FONT_FAMILY_TARGET,
      fontface = make_bold_face(COLUMN_GROUP_TITLE_BOLD)
    ),
    column_title_rot = 0,

    rect_gp = gpar(type = "none"),

    cell_fun = function(j, i, x, y, w, h, fill) {
      border_col <- if (isTRUE(SHOW_CELL_BORDER)) CELL_BORDER_COLOR else NA
      border_lwd <- if (isTRUE(SHOW_CELL_BORDER)) CELL_BORDER_LWD else 0

      grid.roundrect(
        x = x,
        y = y,
        width = w * TILE_SHRINK,
        height = h * TILE_SHRINK,
        r = unit(ROUNDED_TILE_RADIUS_MM, "mm"),
        just = "centre",
        gp = gpar(
          fill = fill,
          col = border_col,
          lwd = border_lwd
        )
      )
    },

    heatmap_legend_param = list(
      title = legend_title,
      at = legend_at,
      title_gp = gpar(
        fontfamily = FONT_FAMILY_TARGET,
        fontsize = as.numeric(LEGEND_TITLE_FONT_SIZE %||% 9)
      ),
      labels_gp = gpar(
        fontfamily = FONT_FAMILY_TARGET,
        fontsize = as.numeric(LEGEND_LABEL_FONT_SIZE %||% 8)
      ),
      legend_direction = "vertical",
      border = NA
    ),

    use_raster = FALSE
  )

  ht
}

draw_one <- function(type, file) {
  if (type == "png") open_png(file, PNG_WIDTH, PNG_HEIGHT, PNG_DPI)
  if (type == "pdf") open_pdf(file, PDF_WIDTH, PDF_HEIGHT)
  on.exit(invisible(dev.off()), add = TRUE)

  ht <- build_ht()

  if (isTRUE(SHOW_MAIN_TITLE)) {
    grid.newpage()
    lay <- grid.layout(
      nrow = 2,
      heights = unit.c(unit(1.2, "lines"), unit(1, "npc") - unit(1.2, "lines"))
    )
    pushViewport(viewport(layout = lay))

    pushViewport(viewport(layout.pos.row = 1))
    grid.text(
      MAIN_TITLE,
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      gp = gpar(
        fontsize = MAIN_TITLE_FONT_SIZE,
        fontfamily = FONT_FAMILY_TARGET,
        fontface = make_bold_face(MAIN_TITLE_BOLD)
      )
    )
    upViewport()

    pushViewport(viewport(layout.pos.row = 2))
    draw(
      ht,
      heatmap_legend_side = HEATMAP_LEGEND_SIDE,
      annotation_legend_side = ANNOTATION_LEGEND_SIDE,
      padding = unit(PLOT_PADDING_MM, "mm"),
      merge_legends = FALSE
    )
    upViewport(2)
  } else {
    draw(
      ht,
      heatmap_legend_side = HEATMAP_LEGEND_SIDE,
      annotation_legend_side = ANNOTATION_LEGEND_SIDE,
      padding = unit(PLOT_PADDING_MM, "mm"),
      merge_legends = FALSE
    )
  }
}

png_file <- file.path(OUTDIR, paste0(OUT_PREFIX, ".png"))
pdf_file <- file.path(OUTDIR, paste0(OUT_PREFIX, ".pdf"))

if (isTRUE(OUTPUT_PNG)) {
  draw_one("png", png_file)
  cat("[OK] PNG:", png_file, "\n")
}

if (isTRUE(OUTPUT_PDF)) {
  draw_one("pdf", pdf_file)
  cat("[OK] PDF:", pdf_file, "\n")
}
