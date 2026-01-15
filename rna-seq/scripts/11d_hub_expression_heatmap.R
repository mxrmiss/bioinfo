#!/usr/bin/env Rscript
# =============================================================================
# 11d_hub_expression_heatmap.R
# Hub genes expression heatmap (Option A: main heatmap no gene labels)
# - Legend: ONE column (z-score on top, group below) -- NO duplicate legends
# - Independent: no config.yaml
# - Inputs: wgcna_objects.rds (datExpr), hub_genes_by_module.tsv, samples.tsv
# - Outputs: PNG/PDF + supplementary TSV (selected hub genes table)
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 顶部参数区（手动改这里） ---------------------------

WGCNA_RDS   <- "results/10_wgcna/rds/wgcna_objects.rds"
HUB_TSV     <- "results/10_wgcna/hub/hub_genes_by_module.tsv"
SAMPLES_TSV <- "data/samples.tsv"
OUTDIR <- "results/10_wgcna/figures"

# 选择哪些模块进入热图（module_color，不带 ME 前缀）
MODULES_FOR_HEATMAP <- c("lightgreen")

# 每个模块取多少个 hub genes（常用：25；两个模块合计约 50）
TOP_N_PER_MODULE <- 50

# 样本排序：是否按 group 分块排在一起
ORDER_SAMPLES_BY_GROUP <- TRUE

# 组织顺序（可选）：如果不想手工指定，设为 NULL；否则按你给的顺序排
TISSUE_ORDER <- c("foot", "gill", "intest", "lips", "mantle", "siphon")
# TISSUE_ORDER <- NULL

# 行/列聚类
CLUSTER_ROWS <- TRUE
CLUSTER_COLUMNS <- FALSE

# 是否显示样本名（默认不显示，避免拥挤）
SHOW_SAMPLE_NAMES <- FALSE

# 方案A：主图不显示基因名（默认 TRUE）
HIDE_GENE_NAMES <- FALSE

# z-score 截断范围（更稳定好看）
ZSCORE_CAP <- 2

# 热图配色（清新一点的蓝-白-珊瑚）
COL_LOW  <- "#3B8BC2"
COL_MID  <- "white"
COL_HIGH <- "#E76F51"

# group 注释配色（清爽马卡龙；不够会自动补）
GROUP_PALETTE <- c(
  "#7FD3C8","#95C8F2","#F6CD96","#F79D93","#A99BEF","#B6E0B6",
  "#A7D9F7","#F6C6E7","#F5B07E","#9FD5CB"
)

# 标题（默认关）
SHOW_TITLE <- FALSE
MAIN_TITLE <- "Hub gene expression (z-score)"

# 图片尺寸
PNG_WIDTH_IN  <- 8
PNG_HEIGHT_IN <- 7
PNG_RES_DPI   <- 300
PDF_WIDTH_IN  <- 10
PDF_HEIGHT_IN <- 7

# 补表输出（TSV）
SUPP_TSV <- "11d_hub_genes.tsv"

# ------------------------------------------------------------------------------

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}
need_pkg(c("data.table", "ComplexHeatmap", "circlize", "grid"))

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

open_png <- function(path, w, h, res) {
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(path, w, h, units = "in", res = res)
  else png(path, width = w * res, height = h * res, res = res, type = "cairo")
}
open_pdf <- function(path, w, h) {
  if (capabilities("cairo")) cairo_pdf(path, width = w, height = h)
  else pdf(path, width = w, height = h, useDingbats = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# 生成足够长的配色（不强依赖额外包）
extend_palette <- function(base_cols, n) {
  base_cols <- as.character(base_cols)
  if (length(base_cols) >= n) return(base_cols[seq_len(n)])
  grDevices::colorRampPalette(base_cols)(n)
}

# =============================================================================
# 主流程
# =============================================================================

safe_mkdir(OUTDIR)

if (!file.exists(WGCNA_RDS)) stop("找不到：", WGCNA_RDS)
if (!file.exists(HUB_TSV)) stop("找不到：", HUB_TSV)
if (!file.exists(SAMPLES_TSV)) stop("找不到：", SAMPLES_TSV)

obj <- readRDS(WGCNA_RDS)
if (is.null(obj$datExpr)) stop("wgcna_objects.rds 内缺少 datExpr。")
datExpr <- obj$datExpr   # samples x genes

hub <- data.table::fread(HUB_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
need_cols <- c("gene_id", "module_color", "kME_primary_module", "GS", "primary_trait")
if (!all(need_cols %in% colnames(hub))) {
  stop("hub_genes_by_module.tsv 缺少列：", paste(setdiff(need_cols, colnames(hub)), collapse = ", "))
}

samples_dt <- data.table::fread(SAMPLES_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!all(c("sample", "group") %in% colnames(samples_dt))) stop("samples.tsv 需要 sample 与 group 列。")

# ---- 1) 选基因：每模块 TOP_N_PER_MODULE ----
hub$abs_score <- abs(hub$kME_primary_module) + abs(hub$GS)

sel_list <- list()
for (mc in MODULES_FOR_HEATMAP) {
  sub <- hub[hub$module_color == mc, , drop = FALSE]
  if (nrow(sub) == 0) next
  sub <- sub[order(sub$abs_score, decreasing = TRUE), , drop = FALSE]
  sel_list[[mc]] <- head(sub, n = min(TOP_N_PER_MODULE, nrow(sub)))
}
sel <- do.call(rbind, sel_list)
if (is.null(sel) || nrow(sel) == 0) stop("未选到任何 hub genes：请检查 MODULES_FOR_HEATMAP 或 hub 表内容。")

# 保持基因顺序：按模块顺序 + score
sel$module_color <- factor(sel$module_color, levels = MODULES_FOR_HEATMAP)
sel <- sel[order(sel$module_color, -sel$abs_score), , drop = FALSE]

genes <- as.character(sel$gene_id)
genes <- genes[!duplicated(genes)]

# datExpr 基因匹配
genes_in <- intersect(genes, colnames(datExpr))
if (length(genes_in) < 5) stop("选中的 hub genes 在 datExpr 中匹配太少（<5）。")

# 若有缺失，给出提示但继续
if (length(setdiff(genes, genes_in)) > 0) {
  message("[WARN] 有部分 hub genes 不在 datExpr 中，将自动忽略：", length(setdiff(genes, genes_in)))
}
genes <- genes_in

# ---- 2) 样本对齐与排序 ----
common_samples <- intersect(rownames(datExpr), samples_dt$sample)
if (length(common_samples) < 6) stop("datExpr 与 samples.tsv 匹配样本太少（<6）。")

samples_dt <- samples_dt[match(common_samples, samples_dt$sample), , drop = FALSE]
datExpr <- datExpr[common_samples, , drop = FALSE]

# 如果指定了 TISSUE_ORDER，就按这个顺序；否则按 samples.tsv 出现顺序
if (is.null(TISSUE_ORDER)) {
  tissue_levels <- unique(as.character(samples_dt$group))
} else {
  # 允许出现额外组织：追加到末尾
  extra <- setdiff(unique(as.character(samples_dt$group)), TISSUE_ORDER)
  tissue_levels <- c(TISSUE_ORDER, extra)
}

if (isTRUE(ORDER_SAMPLES_BY_GROUP)) {
  ord <- order(factor(samples_dt$group, levels = tissue_levels), samples_dt$sample)
} else {
  ord <- seq_len(nrow(samples_dt))
}
samples_dt <- samples_dt[ord, , drop = FALSE]
datExpr <- datExpr[samples_dt$sample, , drop = FALSE]

groups <- as.character(samples_dt$group)
group_levels <- unique(groups)
group_cols <- extend_palette(GROUP_PALETTE, length(group_levels))
col_g <- setNames(group_cols, group_levels)

# ---- 3) 取表达并 z-score（按基因 across samples）----
mat <- datExpr[, genes, drop = FALSE]  # samples x genes
mat <- t(mat)                          # genes x samples
z <- t(scale(t(mat)))
z[is.na(z)] <- 0

# 截断 z-score
cap <- as.numeric(ZSCORE_CAP %||% 2)
z[z >  cap] <-  cap
z[z < -cap] <- -cap

# ---- 4) Legend：只保留一列（z-score 在上，group 在下），杜绝重复 ----
col_fun <- circlize::colorRamp2(c(-cap, 0, cap), c(COL_LOW, COL_MID, COL_HIGH))

lgd_z <- Legend(
  title = "z-score",
  col_fun = col_fun,
  at = c(-cap, 0, cap),
  labels = c(paste0("-", cap), "0", as.character(cap))
)

lgd_g <- Legend(
  title = "group",
  labels = names(col_g),
  legend_gp = gpar(fill = unname(col_g))
)

lgd_onecol <- packLegend(lgd_z, lgd_g, direction = "vertical", gap = unit(3, "mm"))

# 顶部注释（group）——注意：这里直接关掉 annotation 自带 legend
ha <- HeatmapAnnotation(
  group = groups,
  col = list(group = col_g),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# ---- 5) Heatmap：关掉自带 heatmap legend（避免 z-score 复制一份）----
ht <- Heatmap(
  z,
  name = "z-score",
  col = col_fun,
  top_annotation = ha,
  show_row_names = !isTRUE(HIDE_GENE_NAMES),
  show_column_names = isTRUE(SHOW_SAMPLE_NAMES),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = isTRUE(CLUSTER_ROWS),
  cluster_columns = isTRUE(CLUSTER_COLUMNS),
  column_split = factor(groups, levels = group_levels),
  gap = unit(2, "mm"),
  show_heatmap_legend = FALSE
)

# 标题（可选）
if (isTRUE(SHOW_TITLE)) {
  ht <- ht + Heatmap(
    matrix(0, nrow = 1, ncol = 1),
    name = "",
    show_heatmap_legend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE
  )
}

png_out <- file.path(OUTDIR, "11d_hub_expression_heatmap.png")
pdf_out <- file.path(OUTDIR, "11d_hub_expression_heatmap.pdf")

open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     heatmap_legend_list = list(lgd_onecol),
     annotation_legend_list = NULL)
invisible(dev.off())

open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     heatmap_legend_list = list(lgd_onecol),
     annotation_legend_list = NULL)
invisible(dev.off())

cat("[OK]", png_out, "\n")
cat("[OK]", pdf_out, "\n")

# ---- 6) 输出补表 TSV（严格对应热图选中的基因与顺序）----
supp <- sel[sel$gene_id %in% genes, c("gene_id", "module_color", "module_ME", "primary_trait", "GS", "kME_primary_module", "abs_score"), drop = FALSE]
supp$rank_in_selected <- seq_len(nrow(supp))

supp_out <- file.path(OUTDIR, SUPP_TSV)
data.table::fwrite(supp, supp_out, sep = "\t", quote = FALSE)
cat("[OK]", supp_out, "\n")

