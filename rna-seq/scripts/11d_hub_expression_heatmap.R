#!/usr/bin/env Rscript
# =============================================================================
# 11d_hub_expression_heatmap.R
# Hub genes expression heatmap
# - Legend: ONE column (z-score on top, group below) - NO duplicate legends
# - Independent: no config.yaml
# - Inputs: wgcna_objects.rds (datExpr), hub_genes_by_module.tsv OR module lists, samples.tsv
# - Outputs: PNG/PDF + supplementary TSV (selected genes table)
#
# 2026-02-04 改进点（彻底版 + 排版增强）：
# - 支持从 results/10_wgcna/lists/module_<color>.list 读取模块基因（默认）
# - module list 可为 1 列（gene_id）或 2 列（gene_id + gene_name）
# - 纵轴行名规则：list 第2列优先；否则 Swiss-Prot gene_name；否则 locus ID
# - 自动检测重复英文名并追加 -1/-2...（仅重复的才追加）
# - 彻底移除 Unannotated/annotation 标记：无色条、无图例、补表无该字段
# - 新增：可将“仅 locus ID”的基因统一放到最下方（不显示任何标题）
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 顶部参数区（手动改这里） ---------------------------

WGCNA_RDS   <- "results/10_wgcna/rds/wgcna_objects.rds"
HUB_TSV     <- "results/10_wgcna/hub/hub_genes_by_module.tsv"
SAMPLES_TSV <- "data/samples.tsv"
OUTDIR <- "results/10_wgcna/figures"

# 基因来源：推荐 "list"（读 results/10_wgcna/lists 下的模块基因文件）
# 可选："list" / "hub"
GENE_SOURCE <- "list"

# 模块 list 文件所在目录（10_wgcna_modules.R 导出的 lists/）
LIST_DIR <- "results/10_wgcna/lists"

# Swiss-Prot 注释表（用于 gene_name 映射）
SWISSPROT_ANNOT_TSV <- "~/data/anno/swissprot/results/99_reports/2025_04/Sinonovacula_constricta.annot.tsv"
SWISSPROT_ID_COL    <- 1
SWISSPROT_GENE_COL  <- 6

# 选择哪些模块进入热图（module_color，不带 ME 前缀）
MODULES_FOR_HEATMAP <- c("lightgreen")

# 每个模块取多少个基因
# - GENE_SOURCE="hub"：按 hub 表排序取 Top N
# - GENE_SOURCE="list"：按 list 文件出现顺序取前 N；设为 NULL 表示全取
TOP_N_PER_MODULE <- 50
# TOP_N_PER_MODULE <- NULL

# 是否将“只有 locus ID（无 list 第2列且无 Swiss-Prot gene_name）”的基因统一放到最下方
MOVE_LOCUS_ONLY_TO_BOTTOM <- TRUE

# 两块之间的间隔（mm）
ROW_SPLIT_GAP_MM <- 1

# 样本排序：是否按 group 分块排在一起
ORDER_SAMPLES_BY_GROUP <- TRUE

# 组织顺序（可选）：如果不想手工指定，设为 NULL；否则按你给的顺序排
TISSUE_ORDER <- c("foot", "gill", "intest", "lips", "mantle", "siphon")
# TISSUE_ORDER <- NULL

# 行/列聚类
CLUSTER_ROWS <- FALSE
CLUSTER_COLUMNS <- FALSE

# 是否显示样本名（默认不显示，避免拥挤）
SHOW_SAMPLE_NAMES <- FALSE

# 主图是否显示基因名（默认显示）
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
PNG_RES_DPI   <- 600
PDF_WIDTH_IN  <- 10
PDF_HEIGHT_IN <- 8

# 补表输出（TSV）
SUPP_TSV <- "11d_hub_genes.tsv"

# 字体
FONT_FAMILY <- "Arial"

# 热图字体控制（新增）
# - 行名字体：基因名
# - 列名字体：样本名（当 SHOW_SAMPLE_NAMES=TRUE 时才会显示）
# - 图例字体：z-score 与 group 两个 legend 的 title/labels
FONT_SIZE_ROW_NAMES     <- 11
FONT_SIZE_COL_NAMES     <- 8
FONT_SIZE_LEGEND_TITLE  <- 10
FONT_SIZE_LEGEND_LABELS <- 9

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

extend_palette <- function(base_cols, n) {
  base_cols <- as.character(base_cols)
  if (length(base_cols) >= n) return(base_cols[seq_len(n)])
  grDevices::colorRampPalette(base_cols)(n)
}

strip_isoform <- function(x) {
  x <- as.character(x)
  sub("\\.[0-9]+$", "", x)
}

read_swissprot_map <- function(path, id_col = 1, gene_col = 6) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    message("[WARN] Swiss-Prot 注释表不存在或未设置：", path)
    return(list(map = character(), raw = data.frame()))
  }
  dt <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  if (ncol(dt) < max(id_col, gene_col)) {
    stop("Swiss-Prot 注释表列数不足：需要至少 ", max(id_col, gene_col), " 列。")
  }
  ids <- strip_isoform(dt[[id_col]])
  gnm <- as.character(dt[[gene_col]])
  gnm[is.na(gnm)] <- ""
  ids[is.na(ids)] <- ""
  keep <- nzchar(ids)
  ids <- ids[keep]
  gnm <- gnm[keep]
  m <- setNames(gnm, ids)
  list(map = m, raw = dt)
}

read_module_list <- function(path) {
  if (!file.exists(path)) stop("找不到模块基因文件：", path)
  dt <- data.table::fread(path, sep = "\t", header = FALSE, data.table = FALSE, fill = TRUE, quote = "")
  if (ncol(dt) < 1) stop("模块基因文件为空：", path)

  gene_id <- as.character(dt[[1]])
  gene_id <- gene_id[!is.na(gene_id)]
  gene_id <- trimws(gene_id)
  gene_id <- gene_id[nzchar(gene_id)]

  label <- rep("", length(gene_id))
  if (ncol(dt) >= 2) {
    label <- as.character(dt[[2]])
    label[is.na(label)] <- ""
    label <- trimws(label)
  }

  keep <- !duplicated(gene_id)
  data.frame(gene_id = gene_id[keep], list_label = label[keep], stringsAsFactors = FALSE)
}

# 生成 plot_label（不输出/不维护任何 annotation 标记）
make_plot_labels <- function(gene_ids, list_labels, swiss_map) {
  gene_ids <- as.character(gene_ids)
  list_labels <- as.character(list_labels %||% rep("", length(gene_ids)))
  if (length(list_labels) != length(gene_ids)) list_labels <- rep("", length(gene_ids))

  swiss_raw <- swiss_map[gene_ids]
  swiss_raw[is.na(swiss_raw)] <- ""
  swiss_raw <- as.character(swiss_raw)

  base_label <- ifelse(nzchar(list_labels), list_labels,
                       ifelse(nzchar(swiss_raw), swiss_raw, gene_ids))

  # 只对“非 locus ID”的命名做重复后缀
  to_suffix <- (base_label != gene_ids) & nzchar(base_label)
  name_vec <- base_label[to_suffix]

  total <- table(name_vec)
  total_cnt <- as.integer(total[name_vec])

  seen <- new.env(parent = emptyenv())
  idx <- integer(length(name_vec))
  for (i in seq_along(name_vec)) {
    nm <- name_vec[i]
    k <- if (exists(nm, envir = seen, inherits = FALSE)) get(nm, envir = seen) + 1L else 1L
    assign(nm, k, envir = seen)
    idx[i] <- k
  }

  suffixed <- name_vec
  need <- total_cnt > 1
  suffixed[need] <- paste0(name_vec[need], "-", idx[need])

  plot_label <- base_label
  plot_label[to_suffix] <- suffixed

  # 再保险：保证最终 plot_label 唯一
  if (any(duplicated(plot_label))) {
    dup <- plot_label
    seen2 <- new.env(parent = emptyenv())
    idx2 <- integer(length(dup))
    for (i in seq_along(dup)) {
      nm <- dup[i]
      k <- if (exists(nm, envir = seen2, inherits = FALSE)) get(nm, envir = seen2) + 1L else 1L
      assign(nm, k, envir = seen2)
      idx2[i] <- k
    }
    need2 <- duplicated(dup) | duplicated(dup, fromLast = TRUE)
    plot_label[need2] <- paste0(plot_label[need2], "-", idx2[need2])
  }

  data.frame(
    gene_id = gene_ids,
    list_label_raw = list_labels,
    swissprot_gene_name_raw = swiss_raw,
    plot_label = plot_label,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 主流程
# =============================================================================

safe_mkdir(OUTDIR)

if (!file.exists(WGCNA_RDS)) stop("找不到：", WGCNA_RDS)
if (!file.exists(SAMPLES_TSV)) stop("找不到：", SAMPLES_TSV)

obj <- readRDS(WGCNA_RDS)
if (is.null(obj$datExpr)) stop("wgcna_objects.rds 内缺少 datExpr。")
datExpr <- obj$datExpr

samples_dt <- data.table::fread(SAMPLES_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!all(c("sample", "group") %in% colnames(samples_dt))) stop("samples.tsv 需要 sample 与 group 列。")

# ---- 0) 读取 Swiss-Prot 映射 ----
sp <- read_swissprot_map(SWISSPROT_ANNOT_TSV, id_col = SWISSPROT_ID_COL, gene_col = SWISSPROT_GENE_COL)
sp_map <- sp$map

# ---- 1) 选基因：list 模式 or hub 模式 ----
GENE_SOURCE <- tolower(GENE_SOURCE %||% "list")
if (!GENE_SOURCE %in% c("list", "hub")) stop("GENE_SOURCE 只能是 'list' 或 'hub'。")

sel <- NULL

if (GENE_SOURCE == "hub") {
  if (!file.exists(HUB_TSV)) stop("找不到：", HUB_TSV)

  hub <- data.table::fread(HUB_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  need_cols <- c("gene_id", "module_color", "kME_primary_module", "GS", "primary_trait")
  if (!all(need_cols %in% colnames(hub))) {
    stop("hub_genes_by_module.tsv 缺少列：", paste(setdiff(need_cols, colnames(hub)), collapse = ", "))
  }

  hub$abs_score <- abs(hub$kME_primary_module) + abs(hub$GS)

  sel_list <- list()
  for (mc in MODULES_FOR_HEATMAP) {
    sub <- hub[hub$module_color == mc, , drop = FALSE]
    if (nrow(sub) == 0) next
    sub <- sub[order(sub$abs_score, decreasing = TRUE), , drop = FALSE]
    n_take <- TOP_N_PER_MODULE
    if (is.null(n_take)) n_take <- nrow(sub)
    sel_list[[mc]] <- head(sub, n = min(as.integer(n_take), nrow(sub)))
  }

  sel <- do.call(rbind, sel_list)
  if (is.null(sel) || nrow(sel) == 0) stop("未选到任何基因：请检查 MODULES_FOR_HEATMAP 或 hub 表内容。")

  sel$module_color <- factor(sel$module_color, levels = MODULES_FOR_HEATMAP)
  sel <- sel[order(sel$module_color, -sel$abs_score), , drop = FALSE]

  genes <- as.character(sel$gene_id)
  genes <- genes[!duplicated(genes)]

  list_label_raw <- rep("", length(genes))
  names(list_label_raw) <- genes

} else {
  all_rows <- list()
  for (mc in MODULES_FOR_HEATMAP) {
    f <- file.path(LIST_DIR, paste0("module_", mc, ".list"))
    df <- read_module_list(f)
    if (nrow(df) == 0) next

    n_take <- TOP_N_PER_MODULE
    if (!is.null(n_take)) {
      n_take <- as.integer(n_take)
      df <- head(df, n = min(n_take, nrow(df)))
    }

    df$module_color <- mc
    all_rows[[mc]] <- df
  }

  sel <- do.call(rbind, all_rows)
  if (is.null(sel) || nrow(sel) == 0) stop("未从 lists 读取到任何基因：请检查 LIST_DIR 与 MODULES_FOR_HEATMAP。")

  sel$module_color <- factor(sel$module_color, levels = MODULES_FOR_HEATMAP)
  sel <- sel[order(sel$module_color), , drop = FALSE]

  genes <- as.character(sel$gene_id)
  genes <- genes[!duplicated(genes)]

  tmp <- sel[match(genes, sel$gene_id), , drop = FALSE]
  list_label_raw <- as.character(tmp$list_label %||% rep("", length(genes)))
  list_label_raw[is.na(list_label_raw)] <- ""
  names(list_label_raw) <- genes
}

# datExpr 基因匹配
genes_in <- intersect(genes, colnames(datExpr))
if (length(genes_in) < 5) stop("选中的基因在 datExpr 中匹配太少（<5）。")

if (length(setdiff(genes, genes_in)) > 0) {
  message("[WARN] 有部分基因不在 datExpr 中，将自动忽略：", length(setdiff(genes, genes_in)))
}
genes <- genes[genes %in% genes_in]

# ---- 2) 样本对齐与排序 ----
common_samples <- intersect(rownames(datExpr), samples_dt$sample)
if (length(common_samples) < 6) stop("datExpr 与 samples.tsv 匹配样本太少（<6）。")

samples_dt <- samples_dt[match(common_samples, samples_dt$sample), , drop = FALSE]
datExpr <- datExpr[common_samples, , drop = FALSE]

if (is.null(TISSUE_ORDER)) {
  tissue_levels <- unique(as.character(samples_dt$group))
} else {
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

# ---- 3) 取表达并 z-score ----
mat <- datExpr[, genes, drop = FALSE]
mat <- t(mat)
z <- t(scale(t(mat)))
z[is.na(z)] <- 0

cap <- as.numeric(ZSCORE_CAP %||% 2)
z[z >  cap] <-  cap
z[z < -cap] <- -cap

# ---- 4) 生成纵轴标签，并可将“仅 locus ID”的基因整体下移 ----
list_labels_vec <- unname(list_label_raw[rownames(z)])
list_labels_vec[is.na(list_labels_vec)] <- ""
lab_df <- make_plot_labels(rownames(z), list_labels_vec, sp_map)

# 定义“仅 locus ID”：list 第2列为空 且 Swiss-Prot gene_name 也为空
is_locus_only <- (!nzchar(lab_df$list_label_raw)) & (!nzchar(lab_df$swissprot_gene_name_raw))
lab_df$is_locus_only <- is_locus_only

# 如需下移：通过 row_split 强制分成上下两块，并保持两块顺序为：先有名，后 locus-only
row_split <- NULL
if (isTRUE(MOVE_LOCUS_ONLY_TO_BOTTOM)) {
  ord_rows <- order(lab_df$is_locus_only)  # FALSE 在前，TRUE 在后
  z <- z[ord_rows, , drop = FALSE]
  lab_df <- lab_df[ord_rows, , drop = FALSE]
  row_split <- factor(ifelse(lab_df$is_locus_only, "2", "1"), levels = c("1", "2"))
}

rownames(z) <- lab_df$plot_label

# ---- 5) Legend：只保留一列（z-score 在上，group 在下） ----
col_fun <- circlize::colorRamp2(c(-cap, 0, cap), c(COL_LOW, COL_MID, COL_HIGH))

lgd_z <- Legend(
  title = "z-score",
  col_fun = col_fun,
  at = c(-cap, 0, cap),
  labels = c(paste0("-", cap), "0", as.character(cap)),
  title_gp  = gpar(fontfamily = FONT_FAMILY, fontsize = as.numeric(FONT_SIZE_LEGEND_TITLE %||% 9)),
  labels_gp = gpar(fontfamily = FONT_FAMILY, fontsize = as.numeric(FONT_SIZE_LEGEND_LABELS %||% 8))
)

lgd_g <- Legend(
  title = "group",
  labels = names(col_g),
  legend_gp = gpar(fill = unname(col_g)),
  title_gp  = gpar(fontfamily = FONT_FAMILY, fontsize = as.numeric(FONT_SIZE_LEGEND_TITLE %||% 9)),
  labels_gp = gpar(fontfamily = FONT_FAMILY, fontsize = as.numeric(FONT_SIZE_LEGEND_LABELS %||% 8))
)

lgd_onecol <- packLegend(lgd_z, lgd_g, direction = "vertical", gap = unit(3, "mm"))

ha <- HeatmapAnnotation(
  group = groups,
  col = list(group = col_g),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# row_split 不显示任何标题（避免出现“未注释/Unannotated”等文字）
# ComplexHeatmap: row_title 可设为空字符串向量，以达到“无标题”的效果
row_title_opt <- NULL
if (!is.null(row_split)) row_title_opt <- rep("", length(levels(row_split)))

ht <- Heatmap(
  z,
  name = "z-score",
  col = col_fun,
  top_annotation = ha,
  show_row_names = !isTRUE(HIDE_GENE_NAMES),
  show_column_names = isTRUE(SHOW_SAMPLE_NAMES),
  column_names_gp = gpar(fontsize = as.numeric(FONT_SIZE_COL_NAMES %||% 8), fontfamily = FONT_FAMILY),
  row_names_gp = gpar(fontsize = as.numeric(FONT_SIZE_ROW_NAMES %||% 8), fontfamily = FONT_FAMILY),

  # 行拆分：让 locus-only 的基因整体移动到下方
  row_split = row_split,
  row_gap = if (!is.null(row_split)) unit(as.numeric(ROW_SPLIT_GAP_MM %||% 2), "mm") else unit(0, "mm"),
  row_title = row_title_opt,

  # 仍允许每个 slice 内聚类（由 CLUSTER_ROWS 控制）
  cluster_rows = isTRUE(CLUSTER_ROWS),
  cluster_columns = isTRUE(CLUSTER_COLUMNS),

  column_split = factor(groups, levels = group_levels),
  gap = unit(2, "mm"),
  show_heatmap_legend = FALSE
)

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

# ---- 6) 输出补表 TSV（对应热图选中的基因，含最终 plot_label）----
supp_base <- data.frame(
  gene_id = lab_df$gene_id,
  plot_label = lab_df$plot_label,
  swissprot_gene_name_raw = lab_df$swissprot_gene_name_raw,
  list_label_raw = lab_df$list_label_raw,
  stringsAsFactors = FALSE
)

mod_map <- NULL
if (!is.null(sel) && "module_color" %in% colnames(sel)) {
  mod_map <- setNames(as.character(sel$module_color), as.character(sel$gene_id))
}
supp_base$module_color <- if (!is.null(mod_map)) unname(mod_map[supp_base$gene_id]) else NA_character_

if (GENE_SOURCE == "hub") {
  extra_cols <- intersect(
    c("module_ME", "primary_trait", "GS", "kME_primary_module", "abs_score"),
    colnames(sel)
  )
  if (length(extra_cols) > 0) {
    tmp <- sel[match(supp_base$gene_id, sel$gene_id), extra_cols, drop = FALSE]
    supp <- cbind(supp_base, tmp)
  } else {
    supp <- supp_base
  }
} else {
  supp <- supp_base
}

supp$rank_in_plot <- seq_len(nrow(supp))

supp_out <- file.path(OUTDIR, SUPP_TSV)
data.table::fwrite(supp, supp_out, sep = "\t", quote = FALSE)
cat("[OK]", supp_out, "\n")

