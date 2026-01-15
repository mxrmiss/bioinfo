#!/usr/bin/env Rscript
# =============================================================================
# 11b_ME_bar_point_compact.R
# Module eigengene (ME) across tissues: compact bars + points + error bars
#
# 设计目标：
# - n=3/组依然“诚实”：每个样本点必须显示
# - 同时更紧凑：用条形/柱子表达 summary，叠加点 + 误差线
# - 默认：bar=mean，error=SE（更符合读者直觉；避免 n=3 IQR 的尴尬）
# - 输出：PNG + PDF（ragg/cairo 优先，减少字体设备坑）
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 顶部参数区（手动改这里） ---------------------------

INPUT_ME_TSV <- "results/10_wgcna/modules/module_eigengenes.tsv"
SAMPLES_TSV  <- "data/samples.tsv"
OUTDIR <- "results/10_wgcna/figures"

# 要画哪些模块（module_eigengenes.tsv 的列名）
MODULES_TO_PLOT <- c("MElightgreen")

# 组织顺序（按 samples.tsv 的 group 值；若出现额外组织，会自动追加到末尾）
TISSUE_ORDER <- c("foot", "gill", "intest", "lips", "mantle", "siphon")

# facet：是否允许每个模块自由 y 轴（默认 FALSE，利于模块间对比）
FREE_Y <- TRUE

# Summary：柱子显示什么（推荐 mean）
SUMMARY_STAT <- "mean"   # "mean" / "median"

# Error bar：误差线含义（推荐 SE；也可 SD / RANGE / IQR）
# - mean: "SE" 或 "SD" 更直观
# - median: "RANGE"（min~max）通常比 IQR 更适合 n=3
ERROR_STAT <- "SE"       # "SE" / "SD" / "RANGE" / "IQR"

# 点的样式（每个样本）
POINT_ALPHA <- 0.90
POINT_SIZE  <- 2.6
POINT_WIDTH <- 0.16      # quasirandom 横向抖动宽度

# 柱子样式
BAR_ALPHA <- 0.55
BAR_WIDTH <- 0.40

# 误差线样式
ERRORBAR_WIDTH <- 0.18
ERRORBAR_LWD   <- 0.85

# 是否显示图例（通常不需要）
SHOW_LEGEND <- FALSE

# 是否显示标题（默认关）
SHOW_TITLE <- FALSE
MAIN_TITLE <- "Module eigengene (ME) by tissue"

# 是否添加 y=0 参考线（强烈推荐）
SHOW_ZERO_LINE <- TRUE
ZERO_LINE_LWD  <- 0.6
ZERO_LINE_LTY  <- "dashed"  # "solid" / "dashed"
ZERO_LINE_COL  <- "black"

# 清新配色（按 TISSUE_ORDER 顺序）
PALETTE_TISSUE <- c(
  "#7FD3C8", "#95C8F2", "#F6CD96", "#F79D93", "#A99BEF", "#B6E0B6",
  "#A7D9F7", "#F6C6E7"
)

# 图片尺寸（更紧凑）
PNG_WIDTH_IN  <- 7
PNG_HEIGHT_IN <- 4.2
PNG_RES_DPI   <- 600
PDF_WIDTH_IN  <- 7
PDF_HEIGHT_IN <- 4.2

# 输出文件名（建议保持 11b 前缀）
OUT_PREFIX <- "11b_ME_bar_point_compact"

# 可选：导出每个 Tissue×Module 的 summary 表（TSV）
EXPORT_SUMMARY_TSV <- TRUE

# ------------------------------------------------------------------------------

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}

need_pkg(c("data.table", "ggplot2", "ggbeeswarm"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggbeeswarm)
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

title_case_one <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  w <- strsplit(tolower(x), "\\s+")[[1]]
  w <- paste0(toupper(substring(w, 1, 1)), substring(w, 2))
  paste(w, collapse = " ")
}
title_case_vec <- function(v) vapply(v, title_case_one, FUN.VALUE = character(1))

# --------------------------- Summary / Error 辅助函数 ---------------------------

fun_summary <- function(x, stat = c("mean", "median")) {
  stat <- match.arg(stat)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (stat == "median") return(stats::median(x, na.rm = TRUE))
  return(mean(x, na.rm = TRUE))
}

fun_error_bounds <- function(x, summary_stat = c("mean", "median"),
                            err = c("SE", "SD", "RANGE", "IQR")) {
  summary_stat <- match.arg(summary_stat)
  err <- match.arg(err)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(data.frame(y = NA_real_, ymin = NA_real_, ymax = NA_real_))

  y <- fun_summary(x, summary_stat)

  if (err == "RANGE") {
    return(data.frame(y = y, ymin = min(x, na.rm = TRUE), ymax = max(x, na.rm = TRUE)))
  }

  if (err == "IQR") {
    qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE, type = 7)
    return(data.frame(y = y, ymin = qs[1], ymax = qs[2]))
  }

  # SD / SE
  sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv)) sdv <- 0
  if (err == "SD") {
    return(data.frame(y = y, ymin = y - sdv, ymax = y + sdv))
  }
  se <- sdv / sqrt(length(x))
  return(data.frame(y = y, ymin = y - se, ymax = y + se))
}

# =============================================================================
# 主流程
# =============================================================================

safe_mkdir(OUTDIR)

if (!file.exists(INPUT_ME_TSV)) stop("找不到：", INPUT_ME_TSV)
if (!file.exists(SAMPLES_TSV)) stop("找不到：", SAMPLES_TSV)

me_dt <- data.table::fread(INPUT_ME_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(me_dt)) stop("module_eigengenes.tsv 必须包含 sample 列。")

samples_dt <- data.table::fread(SAMPLES_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!all(c("sample", "group") %in% colnames(samples_dt))) stop("samples.tsv 必须包含 sample 与 group 列。")

miss_me <- setdiff(MODULES_TO_PLOT, colnames(me_dt))
if (length(miss_me) > 0) stop("module_eigengenes.tsv 缺少模块列：", paste(miss_me, collapse = ", "))

common <- intersect(me_dt$sample, samples_dt$sample)
if (length(common) < 6) stop("样本匹配太少（<6），请检查 sample 名是否一致。")

me_dt <- me_dt[match(common, me_dt$sample), , drop = FALSE]
samples_dt <- samples_dt[match(common, samples_dt$sample), , drop = FALSE]

plot_dt <- data.frame(
  sample = me_dt$sample,
  group  = as.character(samples_dt$group),
  me_dt[, MODULES_TO_PLOT, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

long <- data.table::melt(
  data.table::as.data.table(plot_dt),
  id.vars = c("sample", "group"),
  variable.name = "module",
  value.name = "ME"
)
long <- as.data.frame(long)

# 组织顺序：TISSUE_ORDER + 自动追加未知组织
seen_groups <- unique(long$group)
extra_groups <- setdiff(seen_groups, TISSUE_ORDER)
final_order <- c(TISSUE_ORDER, extra_groups)

long$group <- factor(long$group, levels = final_order)
long$Tissue <- factor(title_case_vec(as.character(long$group)), levels = title_case_vec(final_order))
long$module <- factor(as.character(long$module), levels = MODULES_TO_PLOT)

# 配色
tissue_levels <- levels(long$Tissue)
pal <- rep(PALETTE_TISSUE, length.out = length(tissue_levels))
names(pal) <- tissue_levels

# 计算 summary：每个 Tissue×module 一行
split_key <- interaction(long$module, long$Tissue, drop = TRUE)
sp <- split(long, split_key)

sum_list <- lapply(sp, function(df) {
  x <- df$ME
  b <- fun_error_bounds(x, summary_stat = SUMMARY_STAT, err = ERROR_STAT)
  data.frame(
    module = as.character(df$module[1]),
    Tissue = as.character(df$Tissue[1]),
    n = sum(is.finite(x)),
    y = b$y, ymin = b$ymin, ymax = b$ymax,
    stringsAsFactors = FALSE
  )
})
sum_df <- do.call(rbind, sum_list)
sum_df$module <- factor(sum_df$module, levels = MODULES_TO_PLOT)
sum_df$Tissue <- factor(sum_df$Tissue, levels = tissue_levels)

if (isTRUE(EXPORT_SUMMARY_TSV)) {
  out_tsv <- file.path(OUTDIR, paste0(OUT_PREFIX, ".summary.tsv"))
  data.table::fwrite(sum_df, out_tsv, sep = "\t", quote = FALSE)
}

# 主题：更紧凑、更“论文感”
base_family <- "sans"
theme_set(theme_classic(base_size = 13, base_family = base_family))

p <- ggplot() +
  # y=0 参考线
  {if (isTRUE(SHOW_ZERO_LINE)) geom_hline(yintercept = 0, linewidth = ZERO_LINE_LWD,
                                         linetype = ZERO_LINE_LTY, color = ZERO_LINE_COL)} +

  # 柱子（summary）
  geom_col(
    data = sum_df,
    aes(x = Tissue, y = y, fill = Tissue),
    width = BAR_WIDTH,
    alpha = BAR_ALPHA,
    color = NA
  ) +

  # 误差线
  geom_errorbar(
    data = sum_df,
    aes(x = Tissue, ymin = ymin, ymax = ymax),
    width = ERRORBAR_WIDTH,
    linewidth = ERRORBAR_LWD,
    color = "black"
  ) +

  # 每个样本点（quasirandom）
  ggbeeswarm::geom_quasirandom(
    data = long,
    aes(x = Tissue, y = ME, color = Tissue),
    width = POINT_WIDTH,
    groupOnX = TRUE,
    size = POINT_SIZE,
    alpha = POINT_ALPHA
  ) +

  facet_wrap(~ module, nrow = 1, scales = if (isTRUE(FREE_Y)) "free_y" else "fixed") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = NULL, y = "Module eigengene (ME)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks  = element_blank(),
    # 保持干净，但不要“厚重外框”
    axis.line   = element_blank(),
    legend.position = if (isTRUE(SHOW_LEGEND)) "right" else "none",

    # strip 更轻：不画大黑框
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),

    # 更紧凑的边距
    plot.margin = margin(6, 6, 6, 6)
  )

if (isTRUE(SHOW_TITLE)) p <- p + ggtitle(MAIN_TITLE)

png_out <- file.path(OUTDIR, paste0(OUT_PREFIX, ".png"))
pdf_out <- file.path(OUTDIR, paste0(OUT_PREFIX, ".pdf"))

open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI); print(p); invisible(dev.off())
open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN); print(p); invisible(dev.off())

cat("[OK]", png_out, "\n")
cat("[OK]", pdf_out, "\n")
if (isTRUE(EXPORT_SUMMARY_TSV)) cat("[OK]", file.path(OUTDIR, paste0(OUT_PREFIX, ".summary.tsv")), "\n")

