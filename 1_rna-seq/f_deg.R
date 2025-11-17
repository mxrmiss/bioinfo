#!/usr/bin/env Rscript
# =========================================================
# F 模块：差异表达 + 顶刊级图件（皇上定制配色）
# 配色（全局固定）：
#   上调：#ECA8A9（暖红粉）  下调：#74AED4（蓝）
#   非显著/背景：#D3E2B7      辅助色：#CFAFD4（紫）、#F7C97E（黄）
# 图件：无灰网格（theme_classic），导出 PNG+PDF
# =========================================================

suppressPackageStartupMessages({
  # library(optparse)  # ← 已移除命令行依赖
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(matrixStats)
  library(systemfonts)   # ← 新增：仅用于字体检测
})

# ---------- 皇上御用配色（统一管理） ----------
COL_UP   <- "#ECA8A9"  # 上调
COL_DOWN <- "#74AED4"  # 下调
COL_NS   <- "#D3E2B7"  # 非显著/背景
COL_AUX1 <- "#CFAFD4"  # 辅助紫
COL_AUX2 <- "#F7C97E"  # 辅助黄
PALETTE5 <- c(COL_UP, COL_DOWN, COL_NS, COL_AUX1, COL_AUX2)

# ---------- 顶刊字体自检与锁定（仅 Arial；缺失即报错退出） ----------
installed_fonts <- unique(systemfonts::system_fonts()$family)
if (!"Arial" %in% installed_fonts) {
  stop(paste0(
    "[F][字体错误] 未检测到 Arial。\n",
    "请先在系统中安装 Arial 字体（Linux 可安装 ttf-mscorefonts-installer，",
    "或将 Arial.ttf 放入 ~/.fonts 并刷新字体缓存）后重试。"
  ))
}
FONT_FAMILY <- "Arial"

# 统一 ggplot2 / ggrepel 字体（不改主题类型、不改 base_size）
theme_set(theme_classic(base_size = 14, base_family = FONT_FAMILY))
update_geom_defaults("text",  list(family = FONT_FAMILY))
update_geom_defaults("label", list(family = FONT_FAMILY))
try(update_geom_defaults("text_repel",  list(family = FONT_FAMILY)), silent = TRUE)
try(update_geom_defaults("label_repel", list(family = FONT_FAMILY)), silent = TRUE)
ggplot2::theme_update(text = element_text(family = FONT_FAMILY))

# —— Cairo PNG 设备包装（供 ggsave 使用）—— 仅新增这段
cairo_png_dev <- function(filename, width, height, dpi = 300, ...) {
  grDevices::png(filename, width = width, height = height, units = "in",
                 res = dpi, type = "cairo", ...)
}

# ---------- 顶部集中参数（原命令行参数已移至此处） ----------
CFG <- list(
  mode        = "txi",                    # "txi" 或 "counts"
  txi_rds     = "",                       # mode="txi" 时必填：tximport RDS 路径
  counts_tsv  = "",                       # mode="counts" 时必填：基因×样本 counts 矩阵（首列基因ID）
  samples     = "samples.tsv",            # 列：sample, group[, batch]
  contrasts   = "contrasts.tsv",          # 列：case, control
  outdeg      = "results/deg",
  outpca      = "results/pca",
  padj        = 0.05,                     # FDR 阈值
  lfc         = 1.0,                      # |log2FC| 阈值
  minCount    = 10,                       # 低表达过滤阈值
  transform   = "vst",                    # "vst" 或 "rlog"
  use_batch   = FALSE                     # 是否启用 batch（samples.tsv 需有 batch 列）
)

dir.create(CFG$outdeg, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$outpca, showWarnings = FALSE, recursive = TRUE)

# ---------- 读表与校验 ----------
samples <- readr::read_tsv(CFG$samples, show_col_types = FALSE)
stopifnot(all(c("sample","group") %in% colnames(samples)))
use_batch <- isTRUE(CFG$use_batch)
if (use_batch && !("batch" %in% colnames(samples))) stop("[F] 已启用 batch，但 samples.tsv 无 batch 列。")

contrasts <- readr::read_tsv(CFG$contrasts, show_col_types = FALSE)
stopifnot(all(c("case","control") %in% colnames(contrasts)))
missing_groups <- setdiff(unique(c(contrasts$case, contrasts$control)), unique(samples$group))
if (length(missing_groups) > 0) stop(paste0("[F] contrasts.tsv 存在未在 samples.tsv 出现的组：", paste(missing_groups, collapse=", ")))

# ---------- 读取表达矩阵 ----------
dds <- NULL
if (CFG$mode == "txi") {
  stopifnot(CFG$txi_rds != "" && file.exists(CFG$txi_rds))
  txi <- readRDS(CFG$txi_rds)
  sample_names <- samples$sample
  stopifnot(all(sample_names %in% colnames(txi$counts)))
  txi$counts <- txi$counts[, sample_names, drop=FALSE]
  colData <- S4Vectors::DataFrame(row.names = sample_names, group = factor(samples$group, levels = unique(samples$group)))
  if (use_batch) colData$batch <- factor(samples$batch)
  design_formula <- if (use_batch) ~ batch + group else ~ group
  dds <- DESeqDataSetFromTximport(txi = txi, colData = colData, design = design_formula)

} else if (CFG$mode == "counts") {
  stopifnot(CFG$counts_tsv != "" && file.exists(CFG$counts_tsv))
  counts_tbl <- readr::read_tsv(CFG$counts_tsv, show_col_types = FALSE)
  stopifnot(colnames(counts_tbl)[1] %in% c("gene","gene_id","id"))
  gene_col <- colnames(counts_tbl)[1]
  counts_mat <- as.matrix(counts_tbl[,-1])
  rownames(counts_mat) <- counts_tbl[[gene_col]]
  sample_names <- samples$sample
  stopifnot(all(sample_names %in% colnames(counts_mat)))
  counts_mat <- counts_mat[, sample_names, drop=FALSE]
  colData <- S4Vectors::DataFrame(row.names = sample_names, group = factor(samples$group, levels = unique(samples$group)))
  if (use_batch) colData$batch <- factor(samples$batch)
  design_formula <- if (use_batch) ~ batch + group else ~ group
  dds <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = colData, design = design_formula)
} else stop("[F] CFG$mode 仅支持 txi 或 counts")

# ---------- 低表达过滤：k = 每组最小生物重复数 ----------
group_sizes <- dplyr::count(samples, group)
k <- min(group_sizes$n)
keep <- rowSums(counts(dds) >= CFG$minCount) >= k
dds <- dds[keep, ]
message(sprintf("[F] 低表达过滤：k=%d, minCount=%d；保留基因 %d 个（占 %.1f%%）",
                k, CFG$minCount, nrow(dds), 100*sum(keep)/length(keep)))

# ---------- 差异分析 ----------
dds <- DESeq(dds)

# =========================================================
# 漂亮图件函数：火山图 / PCA / 样本距离热图
# =========================================================

# 火山图（仅补：legend 的 family、PNG 600dpi）
plot_volcano_pretty <- function(
  df, case_label="Case", ctrl_label="Control",
  lfc_cut=1.0, padj_cut=0.05, x_lim=8, topn_label=10,
  title="Volcano Plot", subtitle=NULL, out_png=NULL, out_pdf=NULL
){
  stopifnot(all(c("gene_id","log2FoldChange","padj") %in% colnames(df)))
  df <- df %>%
    mutate(
      neglog10_p = -log10(padj),
      sig = dplyr::case_when(
        is.na(padj) ~ "NS",
        padj < padj_cut & log2FoldChange >=  lfc_cut ~ "UP",
        padj < padj_cut & log2FoldChange <= -lfc_cut ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  n_up <- sum(df$sig=="UP",na.rm=TRUE); n_down <- sum(df$sig=="DOWN",na.rm=TRUE); n_ns <- sum(df$sig=="NS",na.rm=TRUE)
  lab_df <- df %>% filter(sig!="NS", !is.na(padj)) %>% mutate(score=abs(log2FoldChange) * -log10(padj)) %>% arrange(desc(score)) %>% slice_head(n = topn_label)
  if (is.null(subtitle)) subtitle <- paste0("|log2FC| \u2265 ", lfc_cut, " & padj < ", padj_cut)

  p <- ggplot(df, aes(x = pmin(pmax(log2FoldChange, -x_lim), x_lim), y = neglog10_p)) +
    geom_point(aes(color = sig), size = 1.9, alpha = if_else(df$sig=="NS", 0.35, 0.75), stroke = 0) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.4, color = "grey40") +
    scale_color_manual(
      values = c("UP"=COL_UP, "DOWN"=COL_DOWN, "NS"=COL_NS),
      breaks = c("UP","DOWN","NS"),
      labels = c(paste0("↑ UP (", n_up, ")"), paste0("↓ DOWN (", n_down, ")"), paste0("NS (", n_ns, ")"))
    ) +
    coord_cartesian(xlim = c(-x_lim, x_lim), expand = TRUE) +
    labs(x = "log2(Fold Change)", y = expression(-log[10](adj.~p)),
         color = expression(padj<0.05~"&"~"|"*log[2]*"FC|">=1)) +
    theme_classic(base_size = 14, base_family = FONT_FAMILY) +
    theme(legend.position = "right",
          legend.title = element_text(size = 11, family = FONT_FAMILY),
          legend.text  = element_text(size = 10, family = FONT_FAMILY),
          plot.subtitle= element_blank())
  if (topn_label > 0 && nrow(lab_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      aes(x = pmin(pmax(log2FoldChange, -x_lim), x_lim), y = neglog10_p, label = gene_id, color = sig),
      size = 3, show.legend = FALSE, max.overlaps = 100, box.padding = 0.25, point.padding = 0.2,
      family = FONT_FAMILY
    )
  }
  p <- p +
    annotate("text", x =  x_lim*0.75, y = Inf, vjust = -0.8, label = paste0("Up in ", case_label), size = 3.5, family = FONT_FAMILY) +
    annotate("text", x = -x_lim*0.75, y = Inf, vjust = -0.8, label = paste0("Down vs ", ctrl_label), size = 3.5, family = FONT_FAMILY)
  if (!is.null(out_png)) ggsave(out_png, p, width = 7, height = 6, dpi = 600, device = cairo_png_dev)
  if (!is.null(out_pdf)) ggsave(out_pdf, p, width = 7, height = 6, device = cairo_pdf)
  p
}

# PCA（仅补：legend 的 family、PNG 600dpi）
plot_pca_pretty <- function(vsd, groups, out_png, out_pdf, ntop=500){
  pal <- rep(PALETTE5, length.out = length(unique(groups)))
  rv <- matrixStats::rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]), scale.=FALSE)
  percentVar <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=groups, sample=colnames(vsd))

  p <- ggplot(df, aes(PC1, PC2, color = group, label = sample)) +
    geom_point(size = 3, alpha = 0.95) +
    ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = 50, size = 3, family = FONT_FAMILY) +
    scale_color_manual(values = pal) +
    labs(x = paste0("PC1 (", sprintf("%.1f", percentVar[1]), "%)"),
         y = paste0("PC2 (", sprintf("%.1f", percentVar[2]), "%)")) +
    theme_classic(base_size = 14, base_family = FONT_FAMILY) +
    theme(legend.position = "right",
          legend.title = element_text(size = 11, family = FONT_FAMILY),
          legend.text  = element_text(size = 10, family = FONT_FAMILY))
  ggsave(out_png, p, width = 7, height = 6, dpi = 600, device = cairo_png_dev)
  ggsave(out_pdf, p, width = 7, height = 6, device = cairo_pdf)
  p
}

# 样本距离热图（补：图例字体 family、PNG 600dpi）
plot_heatmap_pretty <- function(vsd, groups, out_png, out_pdf) {
  # 1) 距离矩阵
  d <- dist(t(assay(vsd)))
  mat <- as.matrix(d)

  # 2) 更直觉的渐变：蓝(近) → 白 → 红(远)，用分位数更稳健
  qs <- quantile(mat, c(0, 0.5, 1), na.rm = TRUE)
  col_fun <- circlize::colorRamp2(qs, c("#74AED4", "#FFFFFF", "#ECA8A9"))

  # 3) 组别注释颜色（任意组数自动循环皇上 5 色）
  pal5 <- c("#ECA8A9", "#74AED4", "#D3E2B7", "#CFAFD4", "#F7C97E")
  grp <- factor(groups)
  grp_levels <- levels(grp)
  grp_cols <- setNames(rep(pal5, length.out = length(grp_levels)), grp_levels)

  top_anno  <- HeatmapAnnotation(
    Group = grp,
    col = list(Group = grp_cols),
    annotation_legend_param = list(
      title = "Group",
      title_gp  = grid::gpar(fontfamily = FONT_FAMILY, fontsize = 10),
      labels_gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = 9)
    )
  )

  ht <- Heatmap(
    mat, name = "Distance", col = col_fun,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8, fontfamily = FONT_FAMILY),
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = FONT_FAMILY),
    column_names_rot = 45,
    border = FALSE,
    top_annotation = top_anno,
    heatmap_legend_param = list(
      title = "Distance",
      title_gp  = grid::gpar(fontfamily = FONT_FAMILY, fontsize = 10),
      labels_gp  = grid::gpar(fontfamily = FONT_FAMILY, fontsize = 9)
    )
  )

  # 导出（PNG+PDF）
  grDevices::png(out_png, width = 7, height = 7, units = "in", res = 600, type = "cairo")
  draw(ht, merge_legend = TRUE)
  dev.off()

  grDevices::cairo_pdf(out_pdf, width = 7, height = 7)
  draw(ht, merge_legend = TRUE)
  dev.off()
}

# ---------- 变换（用于 PCA/热图） ----------
vsd <- if (CFG$transform == "rlog") rlog(dds, blind=TRUE) else vst(dds, blind=TRUE)

plot_pca_pretty(
  vsd, groups = colData(dds)$group,
  out_png = file.path(CFG$outpca, "PCA_pretty.png"),
  out_pdf = file.path(CFG$outpca, "PCA_pretty.pdf")
)

plot_heatmap_pretty(
  vsd, groups = colData(dds)$group,
  out_png = file.path(CFG$outpca, "sample_distance_pretty.png"),
  out_pdf = file.path(CFG$outpca, "sample_distance_pretty.pdf")
)

# ---------- 按对比输出 DEG + 漂亮火山图 ----------
summary_list <- list()
for (i in seq_len(nrow(contrasts))) {
  case <- as.character(contrasts$case[i])
  ctrl <- as.character(contrasts$control[i])
  contrast_dir <- file.path(CFG$outdeg, paste0(case, "_vs_", ctrl))
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

  dds_tmp <- dds
  dds_tmp$group <- relevel(dds_tmp$group, ref = ctrl)
  dds_tmp <- nbinomWaldTest(dds_tmp)

  coef_name <- paste0("group_", case, "_vs_", ctrl)
  rn <- resultsNames(dds_tmp)
  if (!coef_name %in% rn) stop(sprintf("[F] 未找到系数 %s；可用：%s", coef_name, paste(rn, collapse=", ")))

  res    <- results(dds_tmp, name = coef_name, alpha = CFG$padj)
  resLFC <- lfcShrink(dds_tmp, coef = coef_name, type = "apeglm")

  res_tbl <- as.data.frame(resLFC) %>% tibble::rownames_to_column("gene_id") %>% arrange(padj)
  sig <- res_tbl %>% filter(!is.na(padj) & padj < CFG$padj & abs(log2FoldChange) >= CFG$lfc)
  up   <- sig %>% filter(log2FoldChange > 0)
  down <- sig %>% filter(log2FoldChange < 0)

  readr::write_tsv(res_tbl, file.path(contrast_dir, "DEG.tsv"))
  readr::write_tsv(up,      file.path(contrast_dir, "DEG_up.tsv"))
  readr::write_tsv(down,    file.path(contrast_dir, "DEG_down.tsv"))
  readr::write_tsv(tibble(gene_id = rownames(dds)), file.path(contrast_dir, "background.tsv"))

  plot_volcano_pretty(
    df = res_tbl,
    case_label = case, ctrl_label = ctrl,
    lfc_cut = CFG$lfc, padj_cut = CFG$padj, x_lim = 8,
    topn_label = 10, title = paste0(case, " vs ", ctrl),
    subtitle = NULL,
    out_png = file.path(contrast_dir, "volcano_pretty.png"),
    out_pdf = file.path(contrast_dir, "volcano_pretty.pdf")
  )

  summary_list[[paste0(case,"_vs_",ctrl)]] <- tibble(
    contrast = paste0(case,"_vs_",ctrl),
    n_total = nrow(res_tbl),
    n_sig = nrow(sig),
    n_up = nrow(up),
    n_down = nrow(down)
  )
}
summary_tbl <- dplyr::bind_rows(summary_list)
readr::write_tsv(summary_tbl, file.path(CFG$outdeg, "summary.tsv"))

message("[F] 完成。结果已写入：")
message(paste0(" - DEG目录：", CFG$outdeg))
message(paste0(" - PCA目录：", CFG$outpca))
