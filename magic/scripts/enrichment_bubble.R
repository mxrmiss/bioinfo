#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: viz_enrichment_bubble_v5.1.R (Viridis-only Color Maps)
# Description:
#   - [配色] 仅使用 viridis 家族连续色板：
#       option = "viridis"  紫 → 蓝 → 绿 → 黄  (默认)
#       option = "magma"    深紫 → 红 → 橙黄
#       option = "plasma"   紫红 → 橙黄
#       option = "inferno"  黑紫 → 红橙 → 黄
#       option = "cividis"  蓝 → 黄（色盲友好、稳重）
#   - [逻辑] 保持原有宽度算法、筛选与网格逻辑不变
#   - [说明] 顶部只需修改 COLOR_MAP_OPTION 即可切换配色
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tools)
  library(scales)
  library(viridis)  # 连续色板：viridis/magma/plasma/inferno/cividis
})

# ==============================================================================
# 1. 用户配置区域 (USER CONFIGURATION)
#    >>> 皇上只需要在这里改参数 <<<
# ==============================================================================

# ------------------------------------------------------------------------------
# [A] 尺寸控制 (DIMENSIONS)
# ------------------------------------------------------------------------------

# TRUE  = [锁定模式] 强制使用 FIXED_HEIGHT_VAL (推荐，保证气泡纵向不被拉伸)
# FALSE = [自动模式] 根据气泡数量自动伸缩高度
USE_FIXED_HEIGHT    <- TRUE

# [参数] 锁定模式下的高度 (单位: inch)
FIXED_HEIGHT_VAL    <- 9.5

# [参数] 锁定模式下单面图固定宽度 (单位: inch)
FIXED_WIDTH_VAL     <- 15

# [参数] 单面图的宽高比 (Aspect Ratio)
# 例：1.9 代表宽度 ≈ 高度 / 1.9，保持气泡区域视觉统一
SINGLE_ASPECT_RATIO <- 2.0

# 自动高度参数 (仅在 USE_FIXED_HEIGHT = FALSE 时生效)
ROW_HEIGHT_INCH     <- 0.35    # 每一行气泡的高度
BASE_HEIGHT_INCH    <- 1.5     # 基础边距

# 分面图专用尺寸 (仅在数据包含多个 Ontology 时生效)
FACET_WIDTH_INCH    <- 12.0
FACET_HEIGHT_INCH   <- 10.0

# 文字排版：通路名称超过多少字符自动换行
MANUAL_WRAP_WIDTH   <- 40

# ------------------------------------------------------------------------------
# [B] 外观与配色 (AESTHETICS & COLORS)
# ------------------------------------------------------------------------------

# 是否显示图片主标题？ (TRUE=显示, FALSE=隐藏/投稿模式)
SHOW_PLOT_TITLE     <- TRUE

# 是否显示网格线？ (TRUE=显示, FALSE=白底)
SHOW_GRID           <- TRUE

# 网格线颜色开关：FALSE = 很浅灰，TRUE = 稍深一点
USE_DARK_GRID       <- FALSE
GRID_COLOR_LIGHT    <- "grey92"
GRID_COLOR_DARK     <- "grey80"

# 颜色标尺模式
# TRUE  = [动态范围] 根据当前数据的 P 值极值自动调整颜色深浅 (对比度最强)
# FALSE = [固定范围] 使用 GLOBAL_P_MIN/MAX
USE_DYNAMIC_COLOR_SCALE <- TRUE
GLOBAL_P_MIN <- 1e-10
GLOBAL_P_MAX <- 0.05

# --- [唯一配色开关：viridis 家族] ---
# 可选值：
#   "viridis"  紫 → 蓝 → 绿 → 黄  (默认 & 推荐)
#   "magma"    深紫 → 红 → 橙黄
#   "plasma"   紫红 → 橙黄
#   "inferno"  黑紫 → 红橙 → 黄
#   "cividis"  蓝 → 黄（色盲友好）
COLOR_MAP_OPTION <- "viridis"

# viridis 颜色采样数量（越大渐变越平滑，一般 64~256 均可）
VIRIDIS_N <- 256

# 字体设置
FONT_FAMILY  <- "Arial"     # 字体 (Arial/Helvetica)
BASE_SIZE    <- 22         # 基础字号

# ------------------------------------------------------------------------------
# [C] 路径设置 (PATHS)
# ------------------------------------------------------------------------------
IN_DIR        <- "input"
SCRIPT_PREFIX <- "bubble_enrichment" # 输出子文件夹名称
OUT_DIR       <- file.path("output", SCRIPT_PREFIX)
TEMPLATE_DIR  <- "templates"

# ==============================================================================
# 2. 辅助函数 (HELPER FUNCTIONS) - 无需修改
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
}

parse_ratio_hybrid <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") return(0)
    if (is.numeric(val)) return(val)
    val_str <- as.character(val)
    if (grepl("/", val_str)) {
      tryCatch({ eval(parse(text = val_str)) }, error = function(e) 0)
    } else {
      as.numeric(val_str)
    }
  })
}

scientific_p_formatter <- function(x) {
  formatC(x, format = "e", digits = 1)
}

generate_fixed_breaks <- function(min_val, max_val, n) {
  if (min_val <= 0) min_val <- 1e-300
  log_seq <- seq(log10(max_val), log10(min_val), length.out = n)
  10^log_seq
}

force_integer_breaks <- function(n = 5) {
  function(x) {
    rng <- range(x, na.rm = TRUE)
    seq_breaks <- seq(rng[1], rng[2], length.out = n)
    unique(round(seq_breaks))
  }
}

# 校验 viridis option，避免用户误输
sanitize_viridis_option <- function(opt) {
  valid_opts <- c("viridis", "magma", "plasma", "inferno", "cividis")
  if (!(opt %in% valid_opts)) {
    warning(sprintf(
      "[WARN] COLOR_MAP_OPTION = '%s' 非法，已自动回退为 'viridis'\n",
      opt
    ))
    return("viridis")
  }
  opt
}

# ==============================================================================
# 3. 主程序循环 (MAIN LOOP)
# ==============================================================================

init_template_system()
files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("[ERROR] input 文件夹下未找到 .tsv 文件。")

# 校正颜色选项
COLOR_MAP_OPTION <- sanitize_viridis_option(COLOR_MAP_OPTION)

for (f_path in files) {

  fname  <- basename(f_path)
  f_base <- file_path_sans_ext(fname)
  cat(sprintf(">> 正在处理文件: %s\n", fname))

  df_raw <- tryCatch({ read_tsv(f_path, show_col_types = FALSE) }, error = function(e) NULL)
  if (is.null(df_raw) || nrow(df_raw) == 0) next

  # Cleaning
  df_plot <- df_raw %>%
    mutate(
      p_adjust          = as.numeric(p_adjust),
      gene_count        = as.numeric(gene_count),
      GeneRatioNumeric  = parse_ratio_hybrid(gene_ratio),
      TermWrapped       = str_wrap(term_name, width = MANUAL_WRAP_WIDTH)
    ) %>%
    filter(!is.na(p_adjust), GeneRatioNumeric > 0)

  if (nrow(df_plot) == 0) next

  # Sorting & Slicing
  unique_ontologies <- if ("ontology" %in% colnames(df_plot)) unique(df_plot$ontology) else c()
  is_multi_ontology <- length(unique_ontologies) > 1
  target_n <- if (is_multi_ontology) 10 else 20

  df_ready <- df_plot %>%
    group_by(if (is_multi_ontology) ontology else NULL) %>%
    arrange(p_adjust) %>%
    slice_head(n = target_n) %>%
    ungroup() %>%
    arrange(GeneRatioNumeric) %>%
    mutate(
      TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)),
      negLog10FDR = -log10(p_adjust)
    )

  final_count <- nrow(df_ready)

  # =====================================================
  # [尺寸计算] 保持气泡区统一宽高逻辑
  # =====================================================

  if (is_multi_ontology) {
    final_width_inch  <- FACET_WIDTH_INCH
    final_height_inch <- FACET_HEIGHT_INCH
  } else {
    if (USE_FIXED_HEIGHT) {
      final_height_inch <- FIXED_HEIGHT_VAL
      final_width_inch  <- FIXED_WIDTH_VAL
    } else {
      final_height_inch <- (final_count * ROW_HEIGHT_INCH) + BASE_HEIGHT_INCH
      final_height_inch <- max(final_height_inch, 3.0)

      grid_width_inch <- final_height_inch / SINGLE_ASPECT_RATIO

      max_char_len    <- max(nchar(as.character(df_ready$TermWrapped)), na.rm = TRUE)
      text_width_inch <- max_char_len * 0.09 + 0.5

      legend_width_inch <- 2.0

      final_width_inch <- text_width_inch + grid_width_inch + legend_width_inch
    }
  }

  cat(sprintf("   [尺寸] %.2f x %.2f inch | 模式: %s | 颜色: %s\n",
              final_width_inch, final_height_inch,
              if (is_multi_ontology) "分面" else "单面",
              COLOR_MAP_OPTION))

  # --- Color scale breaks & limits ---
  if (USE_DYNAMIC_COLOR_SCALE) {
    current_p_min <- min(df_ready$p_adjust, na.rm = TRUE)
    current_p_max <- max(df_ready$p_adjust, na.rm = TRUE)
    if (current_p_min == current_p_max) {
      current_p_max <- if (current_p_max == 0) 0.01 else current_p_max * 10
    }
    plot_breaks <- generate_fixed_breaks(current_p_min, current_p_max, 5)
    plot_limits <- c(current_p_min, current_p_max)
  } else {
    plot_breaks <- generate_fixed_breaks(GLOBAL_P_MIN, GLOBAL_P_MAX, 5)
    plot_limits <- c(GLOBAL_P_MIN, GLOBAL_P_MAX)
  }

  # 使用 viridis 家族生成连续色板
  color_vec <- viridis(
    n       = VIRIDIS_N,
    option  = COLOR_MAP_OPTION,
    direction = -1  # -1: 让“最显著的小 p 值”对应明亮端 or 深色端，可按喜好调
  )

  # --- Plotting ---
  p <- ggplot(df_ready, aes(x = GeneRatioNumeric, y = TermWrapped)) +

    geom_point(aes(size = gene_count, color = negLog10FDR), alpha = 0.9) +

    scale_color_gradientn(
      colours = color_vec,
      breaks  = -log10(plot_breaks),
      limits  = -log10(rev(plot_limits)),
      labels  = label_number(accuracy = 0.1),
      oob     = scales::squish,
      name    = "-log10(FDR)\n"
    ) +

    scale_size_continuous(
      range  = c(3, 8),
      breaks = force_integer_breaks(4),
      name   = "Count\n"
    ) +

    labs(
      x     = "GeneRatio",
      y     = NULL,
      title = if (SHOW_PLOT_TITLE) gsub("_", " ", f_base) else NULL
    ) +

    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY)

  # 网格线样式
  grid_col <- if (USE_DARK_GRID) GRID_COLOR_DARK else GRID_COLOR_LIGHT

  p <- p +
    theme(
      aspect.ratio = if (!is_multi_ontology) SINGLE_ASPECT_RATIO else NULL,

      panel.grid.major.x = if (SHOW_GRID) element_line(colour = grid_col, linewidth = 0.3) else element_blank(),
      panel.grid.major.y = if (SHOW_GRID) element_line(colour = grid_col, linewidth = 0.3) else element_blank(),
      panel.grid.minor.x = if (SHOW_GRID) element_line(colour = grid_col, linewidth = 0.25) else element_blank(),
      panel.grid.minor.y = if (SHOW_GRID) element_line(colour = grid_col, linewidth = 0.25) else element_blank(),

      axis.text  = element_text(colour = "black", size = BASE_SIZE),
      axis.title = element_text(size = BASE_SIZE, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = BASE_SIZE + 3),

      legend.key.height = unit(0.8, "cm")
    )

  if (is_multi_ontology) {
    p <- p +
      facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
      theme(strip.background = element_rect(fill = "grey95", color = NA))
  }

  ggsave(
    file.path(OUT_DIR, paste0(f_base, ".bubble.pdf")),
    p,
    width  = final_width_inch,
    height = final_height_inch,
    device = cairo_pdf
  )
  ggsave(
    file.path(OUT_DIR, paste0(f_base, ".bubble.png")),
    p,
    width  = final_width_inch,
    height = final_height_inch,
    dpi    = 600,
    bg     = "white"
  )
  cat("   [成功] 图片已保存。\n")
}

cat("所有任务处理完成。\n")
