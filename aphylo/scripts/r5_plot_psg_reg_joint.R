#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# =============================================================================
# r5_plot_psg_reg_joint.R
# =============================================================================
# 功能：
# 1）读取 PSG_REG_joint.tsv
# 2）绘制 PSG vs REG 四象限散点图
# 3）横轴为 -log10(PSG q-value)
# 4）纵轴为 -log10(REG q-value)
# 5）图中标签由脚本顶部手动指定 gene id -> 英文标签 控制
# 6）只有 joint_class 属于 Both 或 PSG_only 的指定标签才会显示
# 7）基因标签以斜体显示
# 8）修复左下角 x 轴 0 刻度位置问题
# 9）不再通过负坐标边界制造留白，而是让贴边点轻微内缩
# 10）输出 PNG / PDF / plot_data.tsv
# =============================================================================

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(scales)
}))

# ==============================
# 顶部参数区（皇上只需要改这里）
# ==============================

# ---- 输入输出路径 ----
INPUT_JOINT_TSV <- "results/05_sbranch_model/summary/PSG_REG_joint.tsv"
OUTPUT_DIR <- "results/05_sbranch_model/plot"
OUTPUT_PREFIX <- "PSG_REG_quadrant"

# ---- 手动指定要显示的 gene id 及其英文标签 ----
# 写法：
# "gene_id" = "英文标签"
# 如果想显示 gene id 本身，可写成 ""，脚本会自动回退显示 gene id
SELECTED_GENE_LABELS <- c(
  "Sco06g02960.1" = "PHKG1",
  "Sco13g10430.1" = "GRIK3"
)

# ---- 允许显示标签的类别 ----
ALLOWED_LABEL_CLASSES <- c("Both", "PSG_only")

# ---- 阈值设置（两轴都用 FDR）----
PSG_ALPHA <- 0.05
REG_ALPHA <- 0.05

# ---- 自动截断设置 ----
USE_AUTO_LIMIT <- TRUE
AUTO_LIMIT_QUANTILE <- 0.995

# ---- 手动坐标上限（若 USE_AUTO_LIMIT=FALSE，则使用这里）----
X_MAX_MANUAL <- NA
Y_MAX_MANUAL <- NA

# ---- 是否把超出上限的点压到边界内 ----
CLIP_TO_AXIS_MAX <- TRUE

# ---- 超出上限时，距离边框保留一点点空隙 ----
CLIP_INSET_RATIO <- 0.03

# ---- 右/上边界额外留白比例 ----
X_AXIS_PADDING_RATIO <- 0.04
Y_AXIS_PADDING_RATIO <- 0.04

# ---- 点与坐标轴的最小显示距离 ----
# 这两个参数用于“点自身轻微内缩”，而不是移动坐标轴边界
# 建议保持很小，避免左下角点群整体被推开太多
POINT_X_AXIS_MIN_RATIO <- 0.010
POINT_Y_AXIS_MIN_RATIO <- 0.020

# ---- 是否启用点的贴边内缩 ----
APPLY_POINT_AXIS_PADDING <- TRUE

# ---- 仅对小于该阈值的点执行内缩 ----
# 防止把离轴已经很远的点也一起推开
POINT_PAD_APPLY_X_BELOW_RATIO <- 0.10
POINT_PAD_APPLY_Y_BELOW_RATIO <- 0.10

# ---- 绘图尺寸 ----
FIG_WIDTH_IN <- 10
FIG_HEIGHT_IN <- 6
PNG_DPI <- 600

# ---- 点样式 ----
POINT_SIZE <- 2.4
POINT_ALPHA <- 0.85
STROKE_SIZE <- 0.15

# ---- 阈值线样式 ----
THRESHOLD_LINE_COLOR <- "grey55"
THRESHOLD_LINE_TYPE <- "dashed"
THRESHOLD_LINE_SIZE <- 0.45

# ---- 标签样式 ----
LABEL_SIZE <- 3.8
LABEL_FORCE <- 4.2
LABEL_BOX_PADDING <- 0.55
LABEL_POINT_PADDING <- 0.25
LABEL_MAX_OVERLAPS <- Inf
LABEL_MIN_SEGMENT_LENGTH <- 0
LABEL_SEED <- 123

# ---- 标签线样式 ----
LABEL_SEGMENT_SIZE <- 0.45
LABEL_SEGMENT_COLOR <- "grey35"

# ---- 标签背景 ----
LABEL_BG_FILL <- "#FFFFFF"
LABEL_BG_ALPHA <- 0.0

# ---- 标签初始偏移 ----
# PSG_only 更倾向往下，Both 更倾向往上
LABEL_NUDGE_X_DEFAULT <- 0
LABEL_NUDGE_Y_DEFAULT <- 0
LABEL_NUDGE_Y_PSG_ONLY <- -1.2
LABEL_NUDGE_Y_BOTH <- 1.2

# ---- ggrepel 排布方向 ----
# 可选 "both" / "y"
LABEL_REPEL_DIRECTION <- "both"

# ---- 是否为标签提供更宽的活动边界 ----
# 点的坐标轴从 0 开始，但标签仍可允许在图外活动
LABEL_XMIN_RATIO <- -0.02
LABEL_XMAX_RATIO <- 1.02
LABEL_YMIN_RATIO <- -0.12
LABEL_YMAX_RATIO <- 1.02

# ---- 信息块样式 ----
INFO_BLOCK_X_RATIO <- 0.58
INFO_BLOCK_Y_RATIO <- 0.93
INFO_LINE_SPACING_RATIO <- 0.06
INFO_TEXT_SIZE <- 4.2
INFO_LABEL_FILL_ALPHA <- 0.78

# ---- 颜色设置 ----
COLOR_NONE <- "#BFBFBF"
COLOR_PSG_ONLY <- "#00A087"
COLOR_REG_ONLY <- "#4DBBD5"
COLOR_BOTH <- "#E64B35"

# ---- 主题与字体 ----
BASE_SIZE <- 15
BASE_FAMILY <- "Arial"

# ---- 边框样式 ----
BORDER_COLOR <- "black"
BORDER_SIZE <- 0.5

# ---- 输出 plot_data 是否保留被截断前原值 ----
KEEP_ORIGINAL_COORDS <- TRUE


# ==============================
# 基础工具函数
# ==============================

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

stop_if_not_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少文件：%s", path), call. = FALSE)
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

calc_axis_max <- function(x, q = 0.995, fallback = 5) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(fallback)
  xmax <- as.numeric(stats::quantile(x, probs = q, na.rm = TRUE))
  xmax <- ceiling(xmax * 10) / 10
  if (!is.finite(xmax) || xmax <= 0) xmax <- fallback
  xmax
}

clip_to_max_with_inset <- function(x, xmax, inset_ratio = 0.03) {
  if (!is.finite(xmax) || xmax <= 0) {
    return(x)
  }
  inset_ratio <- max(0, min(inset_ratio, 0.2))
  inset_max <- xmax * (1 - inset_ratio)
  out <- x
  over <- !is.na(x) & x > xmax
  out[over] <- inset_max
  out
}

choose_device_png <- function(file, width, height, dpi) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = file, width = width, height = height, units = "in", res = dpi)
    return(invisible(TRUE))
  }

  caps <- capabilities()
  if (!is.null(caps) && "cairo" %in% names(caps) && isTRUE(caps["cairo"])) {
    grDevices::png(filename = file, width = width, height = height, units = "in", res = dpi, type = "cairo")
    return(invisible(TRUE))
  }

  grDevices::png(filename = file, width = width, height = height, units = "in", res = dpi)
  invisible(TRUE)
}

choose_device_pdf <- function(file, width, height, family = "sans") {
  caps <- capabilities()
  if (!is.null(caps) && "cairo" %in% names(caps) && isTRUE(caps["cairo"])) {
    grDevices::cairo_pdf(filename = file, width = width, height = height, family = family)
    return(invisible(TRUE))
  }

  grDevices::pdf(file = file, width = width, height = height, family = family, useDingbats = FALSE)
  invisible(TRUE)
}

calc_minus_log10 <- function(x, zero_cap = 300) {
  x <- safe_num(x)
  out <- rep(NA_real_, length(x))
  out[is.finite(x) & x > 0] <- -log10(x[is.finite(x) & x > 0])
  out[is.finite(x) & x <= 0] <- zero_cap
  out
}

make_info_block_df <- function(xmin, xmax, ymin, ymax, counts_df) {
  none_n <- counts_df$n[counts_df$class == "None"]
  psg_n  <- counts_df$n[counts_df$class == "PSG_only"]
  reg_n  <- counts_df$n[counts_df$class == "REG_only"]
  both_n <- counts_df$n[counts_df$class == "Both"]

  none_n <- ifelse(length(none_n) == 0, 0, none_n)
  psg_n  <- ifelse(length(psg_n) == 0, 0, psg_n)
  reg_n  <- ifelse(length(reg_n) == 0, 0, reg_n)
  both_n <- ifelse(length(both_n) == 0, 0, both_n)

  xr <- xmax - xmin
  yr <- ymax - ymin

  x0 <- xmin + xr * INFO_BLOCK_X_RATIO
  y0 <- ymin + yr * INFO_BLOCK_Y_RATIO
  dy <- yr * INFO_LINE_SPACING_RATIO

  tibble::tibble(
    x = c(x0, x0, x0, x0),
    y = c(y0, y0 - dy, y0 - 2 * dy, y0 - 3 * dy),
    label = c(
      sprintf("Both (%d)", both_n),
      sprintf("Positive selection only (%d)", psg_n),
      sprintf("Rapid evolution only (%d)", reg_n),
      sprintf("None (%d)", none_n)
    ),
    color = c(COLOR_BOTH, COLOR_PSG_ONLY, COLOR_REG_ONLY, COLOR_NONE)
  )
}

split_fg_ids <- function(fg_gene_ids) {
  ids <- unlist(strsplit(ifelse(is.na(fg_gene_ids), "", fg_gene_ids), ";", fixed = TRUE))
  ids <- trimws(ids)
  ids[nzchar(ids)]
}

has_selected_gene <- function(fg_gene_ids, selected_ids) {
  if (length(selected_ids) == 0) return(FALSE)
  ids <- split_fg_ids(fg_gene_ids)
  any(ids %in% selected_ids)
}

pick_label_from_selected_map <- function(fg_gene_ids, gene_label, selected_map) {
  ids <- split_fg_ids(fg_gene_ids)
  hit <- ids[ids %in% names(selected_map)]

  if (length(hit) > 0) {
    gid <- hit[1]
    lab <- unname(selected_map[gid])
    if (!is.na(lab) && nzchar(trimws(lab))) {
      return(trimws(lab))
    } else {
      return(gid)
    }
  }

  if (!is.na(gene_label) && nzchar(trimws(gene_label))) {
    return(trimws(gene_label))
  }

  return("")
}

make_italic_parse_label <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("'", "\\\\'", x)
  sprintf("italic('%s')", x)
}

make_axis_breaks <- function(min_val, max_val, n = 6, keep_zero = TRUE) {
  brks <- pretty(c(min_val, max_val), n = n)
  brks <- brks[is.finite(brks)]
  brks <- brks[brks >= min_val & brks <= max_val]
  if (keep_zero) {
    brks <- sort(unique(c(0, brks[brks > 0])))
  }
  brks
}

apply_axis_padding_to_points <- function(x, axis_min_pad, apply_below) {
  out <- x
  idx <- !is.na(out) & is.finite(out) & out < axis_min_pad & out <= apply_below
  out[idx] <- axis_min_pad
  out
}


# ==============================
# 读取数据
# ==============================

msg("开始读取输入文件")
stop_if_not_file(INPUT_JOINT_TSV)

joint_df <- readr::read_tsv(
  INPUT_JOINT_TSV,
  col_types = cols(.default = col_character())
)

if (nrow(joint_df) == 0) {
  stop("PSG_REG_joint.tsv 为空，无法绘图。", call. = FALSE)
}

required_cols <- c(
  "OG", "foreground", "gene_label", "foreground_gene_ids", "joint_class",
  "psg_q", "reg_q"
)

missing_cols <- setdiff(required_cols, colnames(joint_df))
if (length(missing_cols) > 0) {
  stop(sprintf("联合表缺少关键列：%s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

joint_df <- joint_df %>%
  mutate(
    psg_p = safe_num(psg_p),
    psg_q = safe_num(psg_q),
    reg_p = safe_num(reg_p),
    reg_q = safe_num(reg_q),
    x_psg_fdr = calc_minus_log10(psg_q),
    y_reg_fdr = calc_minus_log10(reg_q),
    joint_class = factor(joint_class, levels = c("None", "REG_only", "PSG_only", "Both"))
  )

plot_df <- joint_df %>%
  filter(is.finite(x_psg_fdr) & is.finite(y_reg_fdr)) %>%
  mutate(
    x_raw = x_psg_fdr,
    y_raw = y_reg_fdr
  )

if (nrow(plot_df) == 0) {
  stop("没有可用于绘图的有效数值行。", call. = FALSE)
}

msg("读取完成：有效绘图点数 = ", nrow(plot_df))


# ==============================
# 计算坐标范围与阈值
# ==============================

x_thr <- -log10(PSG_ALPHA)
y_thr <- -log10(REG_ALPHA)

if (USE_AUTO_LIMIT) {
  x_base_max <- calc_axis_max(plot_df$x_psg_fdr, q = AUTO_LIMIT_QUANTILE, fallback = max(3, x_thr + 1))
  y_base_max <- calc_axis_max(plot_df$y_reg_fdr, q = AUTO_LIMIT_QUANTILE, fallback = max(3, y_thr + 1))
} else {
  x_base_max <- ifelse(is.na(X_MAX_MANUAL), max(plot_df$x_psg_fdr, na.rm = TRUE), X_MAX_MANUAL)
  y_base_max <- ifelse(is.na(Y_MAX_MANUAL), max(plot_df$y_reg_fdr, na.rm = TRUE), Y_MAX_MANUAL)
}

x_base_max <- max(x_base_max, x_thr + 0.5)
y_base_max <- max(y_base_max, y_thr + 0.5)

if (CLIP_TO_AXIS_MAX) {
  plot_df <- plot_df %>%
    mutate(
      x_plot_prepad = clip_to_max_with_inset(x_psg_fdr, x_base_max, inset_ratio = CLIP_INSET_RATIO),
      y_plot_prepad = clip_to_max_with_inset(y_reg_fdr, y_base_max, inset_ratio = CLIP_INSET_RATIO),
      clipped_x = ifelse(x_psg_fdr > x_base_max, TRUE, FALSE),
      clipped_y = ifelse(y_reg_fdr > y_base_max, TRUE, FALSE)
    )
} else {
  plot_df <- plot_df %>%
    mutate(
      x_plot_prepad = x_psg_fdr,
      y_plot_prepad = y_reg_fdr,
      clipped_x = FALSE,
      clipped_y = FALSE
    )
}

# 坐标轴边界固定从 0 开始，避免 0 刻度位置漂移
x_min <- 0
y_min <- 0
x_max <- x_base_max * (1 + X_AXIS_PADDING_RATIO)
y_max <- y_base_max * (1 + Y_AXIS_PADDING_RATIO)

# 点贴边内缩逻辑
point_x_min_pad <- x_max * POINT_X_AXIS_MIN_RATIO
point_y_min_pad <- y_max * POINT_Y_AXIS_MIN_RATIO
point_x_apply_below <- x_max * POINT_PAD_APPLY_X_BELOW_RATIO
point_y_apply_below <- y_max * POINT_PAD_APPLY_Y_BELOW_RATIO

if (APPLY_POINT_AXIS_PADDING) {
  plot_df <- plot_df %>%
    mutate(
      x_plot = apply_axis_padding_to_points(
        x = x_plot_prepad,
        axis_min_pad = point_x_min_pad,
        apply_below = point_x_apply_below
      ),
      y_plot = apply_axis_padding_to_points(
        x = y_plot_prepad,
        axis_min_pad = point_y_min_pad,
        apply_below = point_y_apply_below
      )
    )
} else {
  plot_df <- plot_df %>%
    mutate(
      x_plot = x_plot_prepad,
      y_plot = y_plot_prepad
    )
}

# 标签活动区
label_xlim <- c(x_max * LABEL_XMIN_RATIO, x_max * LABEL_XMAX_RATIO)
label_ylim <- c(y_max * LABEL_YMIN_RATIO, y_max * LABEL_YMAX_RATIO)

# 明确刻度，保证 0 刻度正常显示
x_breaks <- make_axis_breaks(0, x_max, n = 6, keep_zero = TRUE)
y_breaks <- make_axis_breaks(0, y_max, n = 6, keep_zero = TRUE)

msg("坐标范围：x_min = ", round(x_min, 4), " ; x_max = ", round(x_max, 4),
    " ; y_min = ", round(y_min, 4), " ; y_max = ", round(y_max, 4))
msg("点内缩参数：point_x_min_pad = ", round(point_x_min_pad, 4),
    " ; point_y_min_pad = ", round(point_y_min_pad, 4))
msg("标签活动范围：label_xmin = ", round(label_xlim[1], 4), " ; label_xmax = ", round(label_xlim[2], 4),
    " ; label_ymin = ", round(label_ylim[1], 4), " ; label_ymax = ", round(label_ylim[2], 4))


# ==============================
# 四类计数
# ==============================

counts_df <- plot_df %>%
  count(joint_class, name = "n") %>%
  mutate(class = as.character(joint_class)) %>%
  select(class, n)

for (cls in c("None", "PSG_only", "REG_only", "Both")) {
  if (!cls %in% counts_df$class) {
    counts_df <- bind_rows(counts_df, tibble(class = cls, n = 0))
  }
}
counts_df <- counts_df %>% arrange(match(class, c("None", "PSG_only", "REG_only", "Both")))

msg("四类计数：")
for (i in seq_len(nrow(counts_df))) {
  msg("  ", counts_df$class[i], " = ", counts_df$n[i])
}


# ==============================
# 手动标签数据
# ==============================

selected_map <- SELECTED_GENE_LABELS
if (is.null(selected_map)) {
  selected_map <- character(0)
}
selected_map <- selected_map[nzchar(names(selected_map))]
selected_ids <- unique(names(selected_map))

if (length(selected_ids) == 0) {
  msg("未指定 SELECTED_GENE_LABELS，本图将不显示任何基因标签")
  label_df <- plot_df[0, , drop = FALSE]
} else {
  label_df <- plot_df %>%
    rowwise() %>%
    mutate(
      has_selected = has_selected_gene(foreground_gene_ids, selected_ids),
      manual_label = pick_label_from_selected_map(foreground_gene_ids, gene_label, selected_map)
    ) %>%
    ungroup() %>%
    filter(
      has_selected,
      as.character(joint_class) %in% ALLOWED_LABEL_CLASSES
    ) %>%
    distinct(OG, foreground, manual_label, .keep_all = TRUE) %>%
    mutate(
      manual_label_parse = make_italic_parse_label(manual_label),
      nudge_x = LABEL_NUDGE_X_DEFAULT,
      nudge_y = dplyr::case_when(
        as.character(joint_class) == "PSG_only" ~ LABEL_NUDGE_Y_PSG_ONLY,
        as.character(joint_class) == "Both" ~ LABEL_NUDGE_Y_BOTH,
        TRUE ~ LABEL_NUDGE_Y_DEFAULT
      )
    )

  msg("手动指定 gene_id 数量 = ", length(selected_ids))
  msg("最终实际显示标签数量 = ", nrow(label_df))
}


# ==============================
# 四象限说明块
# ==============================

info_df <- make_info_block_df(
  xmin = x_min,
  xmax = x_max,
  ymin = y_min,
  ymax = y_max,
  counts_df = counts_df
)


# ==============================
# 配色
# ==============================

class_colors <- c(
  "None" = COLOR_NONE,
  "PSG_only" = COLOR_PSG_ONLY,
  "REG_only" = COLOR_REG_ONLY,
  "Both" = COLOR_BOTH
)


# ==============================
# 构图
# ==============================

msg("开始构图")

p <- ggplot(plot_df, aes(x = x_plot, y = y_plot)) +
  geom_point(
    data = plot_df %>% filter(joint_class == "None"),
    aes(color = joint_class),
    size = POINT_SIZE,
    alpha = POINT_ALPHA * 0.65,
    stroke = STROKE_SIZE
  ) +
  geom_point(
    data = plot_df %>% filter(joint_class == "REG_only"),
    aes(color = joint_class),
    size = POINT_SIZE,
    alpha = POINT_ALPHA,
    stroke = STROKE_SIZE
  ) +
  geom_point(
    data = plot_df %>% filter(joint_class == "PSG_only"),
    aes(color = joint_class),
    size = POINT_SIZE,
    alpha = POINT_ALPHA,
    stroke = STROKE_SIZE
  ) +
  geom_point(
    data = plot_df %>% filter(joint_class == "Both"),
    aes(color = joint_class),
    size = POINT_SIZE,
    alpha = POINT_ALPHA,
    stroke = STROKE_SIZE
  ) +
  geom_vline(
    xintercept = x_thr,
    linetype = THRESHOLD_LINE_TYPE,
    color = THRESHOLD_LINE_COLOR,
    linewidth = THRESHOLD_LINE_SIZE
  ) +
  geom_hline(
    yintercept = y_thr,
    linetype = THRESHOLD_LINE_TYPE,
    color = THRESHOLD_LINE_COLOR,
    linewidth = THRESHOLD_LINE_SIZE
  ) +
  scale_color_manual(values = class_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = x_breaks,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    x = expression(-log[10] * "(FDR of positive selection)"),
    y = expression(-log[10] * "(FDR of rapid evolution)")
  ) +
  theme_classic(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "black", size = BASE_SIZE + 1),
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.border = element_rect(color = BORDER_COLOR, fill = NA, linewidth = BORDER_SIZE),
    plot.margin = margin(12, 18, 28, 12)
  )

p <- p +
  geom_label(
    data = info_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    linewidth = 0,
    fill = scales::alpha("white", INFO_LABEL_FILL_ALPHA),
    color = info_df$color,
    size = INFO_TEXT_SIZE,
    family = BASE_FAMILY,
    fontface = "plain",
    label.padding = grid::unit(0.18, "lines")
  )

if (nrow(label_df) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(
          x = x_plot,
          y = y_plot,
          label = manual_label_parse,
          color = joint_class
        ),
        parse = TRUE,
        nudge_x = label_df$nudge_x,
        nudge_y = label_df$nudge_y,
        size = LABEL_SIZE,
        family = BASE_FAMILY,
        force = LABEL_FORCE,
        box.padding = LABEL_BOX_PADDING,
        point.padding = LABEL_POINT_PADDING,
        max.overlaps = LABEL_MAX_OVERLAPS,
        min.segment.length = LABEL_MIN_SEGMENT_LENGTH,
        segment.size = LABEL_SEGMENT_SIZE,
        segment.color = LABEL_SEGMENT_COLOR,
        seed = LABEL_SEED,
        bg.color = scales::alpha(LABEL_BG_FILL, LABEL_BG_ALPHA),
        direction = LABEL_REPEL_DIRECTION,
        xlim = label_xlim,
        ylim = label_ylim,
        show.legend = FALSE
      )
  } else {
    msg("未检测到 ggrepel，标签将使用普通文本，可能发生重叠")
    p <- p +
      geom_text(
        data = label_df,
        aes(
          x = x_plot,
          y = y_plot,
          label = manual_label_parse,
          color = joint_class
        ),
        parse = TRUE,
        size = LABEL_SIZE,
        family = BASE_FAMILY,
        vjust = -0.6,
        show.legend = FALSE
      )
  }
}


# ==============================
# 导出 plot_data
# ==============================

ensure_dir(OUTPUT_DIR)

plot_data_out <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_plot_data.tsv"))

need_export_cols <- c(
  "OG", "foreground", "gene_label", "foreground_gene_ids",
  "psg_p", "psg_q", "reg_p", "reg_q",
  "omega_foreground", "omega_background", "delta_omega",
  "psg_sig", "reg_sig", "joint_class",
  "x_raw", "y_raw",
  "x_plot_prepad", "y_plot_prepad",
  "x_plot", "y_plot",
  "clipped_x", "clipped_y"
)

present_export_cols <- intersect(need_export_cols, colnames(plot_df))

plot_export_df <- plot_df %>%
  select(all_of(present_export_cols))

if (!KEEP_ORIGINAL_COORDS) {
  drop_cols <- intersect(c("x_raw", "y_raw"), colnames(plot_export_df))
  if (length(drop_cols) > 0) {
    plot_export_df <- plot_export_df %>% select(-all_of(drop_cols))
  }
}

readr::write_tsv(plot_export_df, plot_data_out)
msg("已写出 plot_data：", plot_data_out)


# ==============================
# 导出图片
# ==============================

png_out <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, ".png"))
pdf_out <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, ".pdf"))

msg("开始导出 PNG：", png_out)
choose_device_png(
  file = png_out,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  dpi = PNG_DPI
)
print(p)
dev.off()

msg("开始导出 PDF：", pdf_out)
choose_device_pdf(
  file = pdf_out,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  family = BASE_FAMILY
)
print(p)
dev.off()

msg("绘图完成")
msg("PNG：", png_out)
msg("PDF：", pdf_out)
msg("PLOT_DATA：", plot_data_out)
