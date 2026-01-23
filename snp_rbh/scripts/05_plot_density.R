#!/usr/bin/env Rscript

# =============================================================================
# 05_plot_density.R
# 功能：基于 step5_stats 产物绘制 3 张图
#   1) chr 柱状图：Chromosome-level divergence
#   2) 1Mb 滑窗分面折线图
#   3) 1Mb 滑窗热图
#
# 重要约束（皇上要求）：
#   - 只支持 step5_stats 输入（不兼容旧 step5/main）
#   - 支持指标开关：METRIC_MODE = "pdistance" 或 "differences"
#   - pdistance 严格 substitutions-only：sum_substitutions / sum_aligned_nuc * 1e6
#   - differences 优先使用 density_per_Mb_aligned（并做一致性检查）
#   - MIN_PAIRS_WINDOW 默认 0：不启用过滤；且图中不显示 min_pairs
#   - 柱状图 y 轴刻度策略：同 03_plot_density_B.R（不手动 breaks，让 ggplot 自动判定）
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ----------------------------- 参数区（皇上只改这里即可） -----------------------------
# 只允许 step5_stats
CHR_TSV <- "results/step5_stats/chr_summary.tsv"
WIN_TSV <- "results/step5_stats/window_1Mb.tsv"

# 指标开关：二选一
#   "pdistance"   : substitutions-only（以 per Mb aligned 的形式展示 = p_distance * 1e6）
#   "differences" : substitutions + indel_bases（优先用 density_per_Mb_aligned）
METRIC_MODE <- "pdistance"

# 滑窗过滤：默认 0=不启用；>0 才启用（仅作用于折线/热图）
MIN_PAIRS_WINDOW <- 0

# 仅绘制主染色体 chr01..chr99（建议保持 chr01..chr19）
DROP_NON_MAIN_CHR <- TRUE
MAIN_CHR_REGEX <- "^chr\\d{2}$"

# 输出
OUTDIR <- "results/step6_plots"
OUT_PREFIX <- "density"

SAVE_PNG <- TRUE
SAVE_PDF <- TRUE
DPI <- 300

# 画布尺寸
W_BAR_IN <- 13
H_BAR_IN <- 5.2

W_FACET_IN <- 18
H_FACET_ROW_IN <- 3.2
FACET_NCOL <- 5

W_HM_IN <- 18
H_HM_IN <- 10

# 折线分面 y 上限：稳健分位
Y_CAP_PCT <- 99.5

# 热图：色阶范围与显示策略
HEATMAP_VMIN_PCT <- 2.0
HEATMAP_VMAX_PCT <- 98.0
APPLY_SQRT_TRANS <- TRUE
HEATMAP_SHOW_MISSING_AS_WHITE <- TRUE
HEATMAP_MAX_XTICKS <- 12

# 热图列对齐：mid_Mb 四舍五入到 1 位作为列键（和你 03_plot_density_B.R 同款）
MID_KEY_DIGITS <- 1

# 颜色
HM_LOW  <- "#D7ECFA"
HM_HIGH <- "#0B4F8A"
BAR_FILL <- "#95C8F2"
BAR_EDGE <- "#3C5488"

# 是否标注 n_pairs（柱上方）
LABEL_N_PAIRS <- TRUE

# 字号（柱状图风格与 03_plot_density_B.R 一致：theme_classic + 默认 breaks）
BAR_GEOM_WIDTH <- 0.75
BAR_THEME_BASE <- 12
BAR_SUBTITLE_SIZE <- 10
BAR_LABEL_SIZE <- 3.2

FACET_BASE_SIZE <- 14
FACET_TITLE_SIZE <- 20
FACET_SUBTITLE_SIZE <- 13
FACET_AXIS_TITLE_SIZE <- 16
FACET_AXIS_TEXT_SIZE <- 12
FACET_STRIP_TEXT_SIZE <- 13

HM_BASE_SIZE <- 14
HM_TITLE_SIZE <- 20
HM_SUBTITLE_SIZE <- 13
HM_AXIS_TITLE_SIZE <- 16
HM_AXIS_TEXT_SIZE <- 12
HM_LEGEND_TITLE_SIZE <- 14
HM_LEGEND_TEXT_SIZE <- 12
# ---------------------------------------------------------------------------

# ----------------------------- 轻量日志输出 -----------------------------
ts_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
msg <- function(...) cat(ts_now(), "[INFO]", ..., "\n")
warn <- function(...) cat(ts_now(), "[WARN]", ..., "\n")
die <- function(...) stop(paste0(...), call. = FALSE)

# ----------------------------- 工具函数 -----------------------------
read_tsv <- function(path) {
  if (!file.exists(path)) die("找不到输入文件：", path)
  df <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0) die("输入文件为空：", path)
  df
}

has_cols <- function(df, cols) all(cols %in% colnames(df))

chr_sort_levels <- function(chrs) {
  chrs <- unique(as.character(chrs))
  key <- function(x) {
    x2 <- trimws(x)
    if (grepl("^chr\\d{2}$", x2)) {
      n <- suppressWarnings(as.integer(sub("^chr", "", x2)))
      if (!is.na(n)) return(sprintf("0_%03d", n))
    }
    if (grepl("^chr\\d+$", x2)) {
      n <- suppressWarnings(as.integer(sub("^chr", "", x2)))
      if (!is.na(n)) return(sprintf("0_%03d", n))
    }
    paste0("1_", x2)
  }
  chrs[order(vapply(chrs, key, character(1)))]
}

save_png <- function(plot_obj, out, w, h, dpi) {
  ok <- FALSE
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = out, width = w, height = h, units = "in", res = dpi)
    print(plot_obj); grDevices::dev.off(); ok <- TRUE
  }
  if (!ok) {
    grDevices::png(filename = out, width = w, height = h, units = "in", res = dpi, type = "cairo")
    print(plot_obj); grDevices::dev.off()
  }
}

save_pdf <- function(plot_obj, out, w, h) {
  grDevices::cairo_pdf(filename = out, width = w, height = h)
  print(plot_obj); grDevices::dev.off()
}

metric_suffix <- function() if (METRIC_MODE == "pdistance") "pdistance" else "differences"

make_bar_subtitle <- function() {
  if (METRIC_MODE == "pdistance") {
    "Metric: substitution differences per Mb aligned (p-distance × 1e6)"
  } else {
    "Metric: differences per Mb aligned (substitutions + indel bases)"
  }
}

make_subtitle_window <- function() {
  base <- "Window: 1Mb"
  base <- if (METRIC_MODE == "pdistance") {
    paste0(base, " | Metric: p-distance (substitutions-only), shown as differences per Mb aligned")
  } else {
    paste0(base, " | Metric: differences per Mb aligned (substitutions + indel bases)")
  }
  if (is.finite(MIN_PAIRS_WINDOW) && MIN_PAIRS_WINDOW > 0) {
    paste0(base, " | min_pairs=", MIN_PAIRS_WINDOW)
  } else {
    base
  }
}

# 数字格式：避免科学计数
lab_num_axis <- label_number(accuracy = 1, big.mark = ",")
lab_num_legend <- label_number(accuracy = 1, big.mark = ",")

# ----------------------------- 指标计算（严格 step5_stats 合同） -----------------------------
calc_chr_metric <- function(df) {
  need <- c("chr", "n_pairs", "sum_aligned_nuc", "sum_substitutions", "sum_indel_bases", "density_per_Mb_aligned")
  if (!has_cols(df, need)) {
    miss <- setdiff(need, colnames(df))
    die("chr_summary.tsv 缺少列：", paste(miss, collapse = ", "))
  }

  aln <- suppressWarnings(as.numeric(df$sum_aligned_nuc))
  sub <- suppressWarnings(as.numeric(df$sum_substitutions))
  ind <- suppressWarnings(as.numeric(df$sum_indel_bases))
  dens <- suppressWarnings(as.numeric(df$density_per_Mb_aligned))

  if (any(!is.finite(aln) | aln <= 0, na.rm = TRUE)) die("chr_summary.tsv 中 sum_aligned_nuc 存在非正值/NA，无法计算")

  y_pd <- (sub / aln) * 1e6
  y_diff_calc <- ((sub + ind) / aln) * 1e6

  if (METRIC_MODE == "pdistance") {
    return(list(y = y_pd,
                y_label = "Divergence (substitution differences per Mb aligned)",
                check_hint = NULL))
  }

  # differences：优先 dens，同时做一致性检查
  ok <- is.finite(dens) & is.finite(y_diff_calc)
  if (any(ok)) {
    rel <- abs(dens[ok] - y_diff_calc[ok]) / pmax(1, abs(y_diff_calc[ok]))
    med_rel <- suppressWarnings(median(rel, na.rm = TRUE))
    if (is.finite(med_rel) && med_rel > 0.01) {
      warn("density_per_Mb_aligned 与 (sub+ind)/aligned*1e6 的中位相对差异约为 ",
           sprintf("%.3f", med_rel),
           "（>1%）。请确认 step5_stats 的计算口径是否一致。")
    }
  }

  return(list(y = dens,
              y_label = "Divergence density (differences per Mb aligned)",
              check_hint = "using density_per_Mb_aligned (checked vs (sub+ind)/aligned*1e6)"))
}

calc_win_metric <- function(df) {
  # 强制需要窗口坐标：win_start/win_end（因为你 03_plot_density_B.R 就是这么用的）
  need_base <- c("chr", "n_pairs", "win_start", "win_end", "sum_aligned_nuc", "sum_substitutions", "sum_indel_bases", "density_per_Mb_aligned")
  if (!has_cols(df, need_base)) {
    miss <- setdiff(need_base, colnames(df))
    die("window_1Mb.tsv 缺少列：", paste(miss, collapse = ", "))
  }

  ws <- suppressWarnings(as.numeric(df$win_start))
  we <- suppressWarnings(as.numeric(df$win_end))
  mid_bp <- (ws + we) / 2
  mid_Mb <- mid_bp / 1e6

  aln <- suppressWarnings(as.numeric(df$sum_aligned_nuc))
  sub <- suppressWarnings(as.numeric(df$sum_substitutions))
  ind <- suppressWarnings(as.numeric(df$sum_indel_bases))
  dens <- suppressWarnings(as.numeric(df$density_per_Mb_aligned))

  if (METRIC_MODE == "pdistance") {
    y <- (sub / aln) * 1e6
    return(list(mid_Mb = mid_Mb, y = y))
  } else {
    return(list(mid_Mb = mid_Mb, y = dens))
  }
}

# ----------------------------- 主流程 -----------------------------
msg("Step05 开始：绘制 divergence density 图（step5_stats 专用）")
msg("METRIC_MODE：", METRIC_MODE)
msg("输入1：", CHR_TSV)
msg("输入2：", WIN_TSV)

if (!(METRIC_MODE %in% c("pdistance", "differences"))) {
  die("METRIC_MODE 只能是 'pdistance' 或 'differences'，你现在是：", METRIC_MODE)
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

chr_raw <- read_tsv(CHR_TSV)
win_raw <- read_tsv(WIN_TSV)

# 主染色体过滤（绘图层做，不改上游数据）
if (DROP_NON_MAIN_CHR) {
  chr_raw <- chr_raw %>% filter(grepl(MAIN_CHR_REGEX, chr))
  win_raw <- win_raw %>% filter(grepl(MAIN_CHR_REGEX, chr))
}

# 指标计算
chr_met <- calc_chr_metric(chr_raw)
win_met <- calc_win_metric(win_raw)

# 组装绘图数据
chr_df <- chr_raw %>%
  mutate(
    chr = as.character(chr),
    n_pairs = suppressWarnings(as.numeric(n_pairs)),
    metric_y = chr_met$y
  ) %>%
  filter(is.finite(metric_y))

win_df <- win_raw %>%
  mutate(
    chr = as.character(chr),
    n_pairs = suppressWarnings(as.numeric(n_pairs)),
    mid_Mb = win_met$mid_Mb,
    metric_y = win_met$y
  ) %>%
  filter(is.finite(mid_Mb) & is.finite(metric_y))

# min_pairs 过滤（仅滑窗）
if (is.finite(MIN_PAIRS_WINDOW) && MIN_PAIRS_WINDOW > 0) {
  before_n <- nrow(win_df)
  win_df <- win_df %>% filter(is.finite(n_pairs) & n_pairs >= MIN_PAIRS_WINDOW)
  after_n <- nrow(win_df)
  msg("应用 min_pairs 过滤：", MIN_PAIRS_WINDOW, "| window 行数：", before_n, "->", after_n)
} else {
  msg("min_pairs 过滤：OFF（MIN_PAIRS_WINDOW=0）")
}

# 染色体顺序
chr_levels <- chr_sort_levels(win_df$chr)
win_df$chr <- factor(win_df$chr, levels = chr_levels)

chr_levels2 <- chr_sort_levels(c(chr_levels, chr_df$chr))
chr_df$chr <- factor(chr_df$chr, levels = chr_levels2) %>% droplevels()
chr_df <- chr_df %>% arrange(chr)

msg("chr_summary 染色体数：", nrow(chr_df))
msg("window 数据行数（过滤后）：", nrow(win_df))

suffix <- metric_suffix()

# =============================================================================
# 图1：染色体级柱状图（刻度策略：不手动设 breaks，跟 03_plot_density_B.R 同款）
# =============================================================================
msg("绘图1/3：Chromosome-level bar")

p1 <- ggplot(chr_df, aes(x = chr, y = metric_y)) +
  geom_col(width = BAR_GEOM_WIDTH, fill = BAR_FILL, color = BAR_EDGE, linewidth = 0.3) +
  theme_classic(base_size = BAR_THEME_BASE) +
  theme(
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.text.x  = element_text(angle = 0, vjust = 0.5),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(size = BAR_SUBTITLE_SIZE)
  ) +
  scale_y_continuous(labels = lab_num_axis) +
  labs(
    title = "Chromosome-level divergence density (one-to-one orthologs)",
    subtitle = make_bar_subtitle(),
    x = "Chromosome",
    y = chr_met$y_label
  )

if (LABEL_N_PAIRS) {
  p1 <- p1 + geom_text(
    aes(label = ifelse(is.na(n_pairs), "", as.integer(n_pairs))),
    vjust = -0.4,
    size = BAR_LABEL_SIZE
  )
}

out1_png <- file.path(OUTDIR, paste0(OUT_PREFIX, ".", suffix, ".chr_bar.png"))
out1_pdf <- file.path(OUTDIR, paste0(OUT_PREFIX, ".", suffix, ".chr_bar.pdf"))

if (SAVE_PNG) { msg("写出 PNG：", out1_png); save_png(p1, out1_png, W_BAR_IN, H_BAR_IN, DPI) }
if (SAVE_PDF) { msg("写出 PDF：", out1_pdf); save_pdf(p1, out1_pdf, W_BAR_IN, H_BAR_IN) }

# =============================================================================
# 图2：滑窗分面折线图
# =============================================================================
msg("绘图2/3：Sliding-window facet line")

all_y <- win_df$metric_y
if (length(all_y) == 0) die("window_1Mb.tsv 没有可用 metric 数据（可能被过滤空了）")

y_cap <- suppressWarnings(as.numeric(quantile(all_y, probs = Y_CAP_PCT / 100, na.rm = TRUE, names = FALSE)))
if (!is.finite(y_cap) || y_cap <= 0) y_cap <- max(all_y, na.rm = TRUE)

n_chr <- length(chr_levels)
nrow_facet <- ceiling(n_chr / FACET_NCOL)
h_facet <- max(4, H_FACET_ROW_IN * nrow_facet)

p2 <- ggplot(win_df, aes(x = mid_Mb, y = metric_y)) +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ chr, ncol = FACET_NCOL, scales = "free_x") +
  coord_cartesian(ylim = c(0, y_cap)) +
  theme_classic(base_size = FACET_BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = FACET_TITLE_SIZE),
    plot.subtitle = element_text(size = FACET_SUBTITLE_SIZE),
    strip.text = element_text(size = FACET_STRIP_TEXT_SIZE),
    axis.title.x = element_text(size = FACET_AXIS_TITLE_SIZE, margin = margin(t = 8)),
    axis.title.y = element_text(size = FACET_AXIS_TITLE_SIZE, margin = margin(r = 10)),
    axis.text.x = element_text(size = FACET_AXIS_TEXT_SIZE),
    axis.text.y = element_text(size = FACET_AXIS_TEXT_SIZE)
  ) +
  scale_y_continuous(labels = lab_num_axis) +
  labs(
    title = "Sliding-window divergence density along chromosomes",
    subtitle = make_subtitle_window(),
    x = "Window midpoint (Mb)",
    y = chr_met$y_label
  )

out2_png <- file.path(OUTDIR, paste0(OUT_PREFIX, ".", suffix, ".window_1Mb.facet_line.png"))
out2_pdf <- file.path(OUTDIR, paste0(OUT_PREFIX, ".", suffix, ".window_1Mb.facet_line.pdf"))

if (SAVE_PNG) { msg("写出 PNG：", out2_png); save_png(p2, out2_png, W_FACET_IN, h_facet, DPI) }
if (SAVE_PDF) { msg("写出 PDF：", out2_pdf); save_pdf(p2, out2_pdf, W_FACET_IN, h_facet) }

# =============================================================================
# 图3：滑窗热图
# =============================================================================
msg("绘图3/3：Sliding-window heatmap")

hm_df <- win_df %>%
  mutate(mid_key = round(mid_Mb, digits = MID_KEY_DIGITS)) %>%
  select(chr, mid_key, metric_y)

all_mid <- sort(unique(hm_df$mid_key))

hm_df_full <- hm_df %>%
  complete(chr = factor(chr_levels, levels = chr_levels),
           mid_key = all_mid) %>%
  arrange(chr, mid_key)

vals <- hm_df_full$metric_y
vals <- vals[is.finite(vals)]
if (length(vals) == 0) die("热图没有任何有效数值（全部 NA）")

vmin <- as.numeric(quantile(vals, probs = HEATMAP_VMIN_PCT / 100, na.rm = TRUE, names = FALSE))
vmax <- as.numeric(quantile(vals, probs = HEATMAP_VMAX_PCT / 100, na.rm = TRUE, names = FALSE))
if (!is.finite(vmin) || !is.finite(vmax) || vmin >= vmax) {
  vmin <- min(vals, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
}

mid_levels <- all_mid
step <- max(1, floor(length(mid_levels) / HEATMAP_MAX_XTICKS))
x_breaks <- mid_levels[seq(1, length(mid_levels), by = step)]

fill_trans <- if (APPLY_SQRT_TRANS) "sqrt" else "identity"

p3 <- ggplot(hm_df_full, aes(x = mid_key, y = chr, fill = metric_y)) +
  geom_tile() +
  theme_classic(base_size = HM_BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = HM_TITLE_SIZE),
    plot.subtitle = element_text(size = HM_SUBTITLE_SIZE),
    axis.title.x = element_text(size = HM_AXIS_TITLE_SIZE, margin = margin(t = 8)),
    axis.title.y = element_text(size = HM_AXIS_TITLE_SIZE, margin = margin(r = 10)),
    axis.text.x = element_text(size = HM_AXIS_TEXT_SIZE),
    axis.text.y = element_text(size = HM_AXIS_TEXT_SIZE),
    legend.title = element_text(size = HM_LEGEND_TITLE_SIZE),
    legend.text  = element_text(size = HM_LEGEND_TEXT_SIZE)
  ) +
  labs(
    title = "Chromosome heatmap of sliding-window divergence density",
    subtitle = make_subtitle_window(),
    x = "Window midpoint (Mb)",
    y = "Chromosome",
    fill = "Metric value"
  ) +
  scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(
    low = HM_LOW,
    high = HM_HIGH,
    na.value = if (HEATMAP_SHOW_MISSING_AS_WHITE) "white" else "grey95",
    limits = c(vmin, vmax),
    oob = squish,
    trans = fill_trans,
    labels = lab_num_legend
  )

out3_png <- file.path(OUTDIR, paste0(OUT_PREFIX, ".", suffix, ".window_1Mb.heatmap.png"))
out3_pdf <- file.path(OUTDIR, paste0(OUT_PREFIX, ".", suffix, ".window_1Mb.heatmap.pdf"))

if (SAVE_PNG) { msg("写出 PNG：", out3_png); save_png(p3, out3_png, W_HM_IN, H_HM_IN, DPI) }
if (SAVE_PDF) { msg("写出 PDF：", out3_pdf); save_pdf(p3, out3_pdf, W_HM_IN, H_HM_IN) }

msg("Step05 完成。输出目录：", OUTDIR)

