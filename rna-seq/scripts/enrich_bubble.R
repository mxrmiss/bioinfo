#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: enrichment_bubble.R
# Description:
#   基于原始 magic/scripts/enrichment_bubble.R 做最小修改：
#   1) 保留原有核心绘图逻辑、尺寸逻辑、配色逻辑
#   2) 修改输入读取逻辑：
#        - auto   : 自动扫描 results/08_enrich 下各 contrast 目录中的标准富集结果
#        - manual : 手动指定某个 contrast 目录或某个具体富集结果文件
#   3) 修改输出目录为：
#        results/plots_enrich/bubble
#   4) 不使用 match_font()/match_fonts()，避免 systemfonts 弃用警告
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tools)
  library(scales)
  library(viridis)
})

# ==============================================================================
# 1. 用户配置区域 (USER CONFIGURATION)
#    >>> 皇上只需要在这里改参数 <<<
# ==============================================================================

# ------------------------------------------------------------------------------
# [A] 运行模式
# ------------------------------------------------------------------------------
# auto   = 自动扫描 results/08_enrich 下所有 contrast 目录
# manual = 手动指定某个 contrast 目录，或某个具体文件
RUN_MODE <- "auto"

# 手动模式输入：
# 1）指定目录，例如：
#    "results/08_enrich/B_vs_M"
# 2）指定文件，例如：
#    "results/08_enrich/B_vs_M/GO_BP_by_term_up.tsv"
# RUN_MODE = "manual" 时生效
MANUAL_INPUT <- ""

# 项目根目录：
# 留空时自动检测（推荐）
# 若自动检测失败，可手工指定，例如：
# PROJECT_ROOT <- "/home/xxx/3project/10_sc_rna"
PROJECT_ROOT <- ""

# ------------------------------------------------------------------------------
# [B] 输入输出路径
# ------------------------------------------------------------------------------
ENRICH_DIR_REL <- "results/08_enrich"
OUT_DIR_REL    <- "results/plots_enrich/bubble"

# ------------------------------------------------------------------------------
# [C] 自动扫描时的筛选
# ------------------------------------------------------------------------------
# 可选：
# "all"        = 全部标准富集结果
# "go_by_term" = 仅 GO_BP/CC/MF_by_term
# "go_sig"     = 仅 GO_sig_up/down
# "kegg"       = 仅 KEGG_by_term_up/down
TARGET_TYPE <- "all"

# 可选：
# "all" / "up" / "down"
DIRECTION_FILTER <- "all"

# 若只想处理部分 contrast，可填写：
# CONTRAST_FILTER <- c("B_vs_M", "J_vs_E")
# 不限制则保持 character(0)
CONTRAST_FILTER <- character(0)

# ------------------------------------------------------------------------------
# [D] 尺寸控制 (DIMENSIONS)
# ------------------------------------------------------------------------------

# TRUE  = [锁定模式] 强制使用 FIXED_HEIGHT_VAL
# FALSE = [自动模式] 根据气泡数量自动伸缩高度
USE_FIXED_HEIGHT    <- TRUE

# [参数] 锁定模式下的高度 (单位: inch)
FIXED_HEIGHT_VAL    <- 9.5

# [参数] 锁定模式下单面图固定宽度 (单位: inch)
FIXED_WIDTH_VAL     <- 15

# [参数] 单面图的宽高比 (Aspect Ratio)
SINGLE_ASPECT_RATIO <- 2.0

# 自动高度参数 (仅在 USE_FIXED_HEIGHT = FALSE 时生效)
ROW_HEIGHT_INCH     <- 0.35
BASE_HEIGHT_INCH    <- 1.5

# 分面图专用尺寸 (仅在数据包含多个 Ontology 时生效)
FACET_WIDTH_INCH    <- 12.0
FACET_HEIGHT_INCH   <- 10.0

# 文字排版：通路名称超过多少字符自动换行
MANUAL_WRAP_WIDTH   <- 60

# ------------------------------------------------------------------------------
# [E] 外观与配色 (AESTHETICS & COLORS)
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
# TRUE  = [动态范围] 根据当前数据的 P 值极值自动调整颜色深浅
# FALSE = [固定范围] 使用 GLOBAL_P_MIN/MAX
USE_DYNAMIC_COLOR_SCALE <- TRUE
GLOBAL_P_MIN <- 1e-10
GLOBAL_P_MAX <- 0.05

# viridis 家族配色
# 可选：
#   "viridis" / "magma" / "plasma" / "inferno" / "cividis"
COLOR_MAP_OPTION <- "viridis"

# viridis 颜色采样数量
VIRIDIS_N <- 256

# 字体设置
# 不再做 systemfonts 探测，直接使用你原脚本的方式，避免弃用警告
FONT_FAMILY  <- "Arial"
BASE_SIZE    <- 14

# ------------------------------------------------------------------------------
# [F] 图片输出
# ------------------------------------------------------------------------------
WRITE_PNG <- TRUE
WRITE_PDF <- TRUE
PNG_DPI   <- 600

# ==============================================================================
# 2. 辅助函数 (HELPER FUNCTIONS)
# ==============================================================================

normalize_path2 <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0(...), "\n")
}

detect_project_root <- function(user_root = "", enrich_rel = "results/08_enrich") {
  if (!is.null(user_root) && nzchar(user_root)) {
    return(normalize_path2(user_root))
  }

  wd <- normalize_path2(getwd())
  if (dir.exists(file.path(wd, enrich_rel))) {
    return(wd)
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_dir  <- dirname(normalize_path2(script_path))
    candidate_root <- normalize_path2(file.path(script_dir, ".."))
    if (dir.exists(file.path(candidate_root, enrich_rel))) {
      return(candidate_root)
    }
  }

  wd
}

parse_ratio_hybrid <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") return(0)
    if (is.numeric(val)) return(val)
    val_str <- as.character(val)
    if (grepl("/", val_str)) {
      tryCatch({ eval(parse(text = val_str)) }, error = function(e) 0)
    } else {
      suppressWarnings(as.numeric(val_str))
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

sanitize_viridis_option <- function(opt) {
  valid_opts <- c("viridis", "magma", "plasma", "inferno", "cividis")
  if (!(opt %in% valid_opts)) {
    warning(sprintf(
      "[WARN] COLOR_MAP_OPTION = '%s' 非法，已自动回退为 'viridis'",
      opt
    ))
    return("viridis")
  }
  opt
}

get_standard_filenames <- function(target_type = "all", direction_filter = "all") {
  all_files <- c(
    "GO_BP_by_term_up.tsv",
    "GO_CC_by_term_up.tsv",
    "GO_MF_by_term_up.tsv",
    "GO_BP_by_term_down.tsv",
    "GO_CC_by_term_down.tsv",
    "GO_MF_by_term_down.tsv",
    "GO_sig_up.tsv",
    "GO_sig_down.tsv",
    "KEGG_by_term_up.tsv",
    "KEGG_by_term_down.tsv"
  )

  if (target_type == "go_by_term") {
    all_files <- all_files[grepl("^GO_(BP|CC|MF)_by_term_(up|down)\\.tsv$", all_files)]
  } else if (target_type == "go_sig") {
    all_files <- all_files[grepl("^GO_sig_(up|down)\\.tsv$", all_files)]
  } else if (target_type == "kegg") {
    all_files <- all_files[grepl("^KEGG_by_term_(up|down)\\.tsv$", all_files)]
  }

  if (direction_filter == "up") {
    all_files <- all_files[grepl("_up\\.tsv$", all_files)]
  } else if (direction_filter == "down") {
    all_files <- all_files[grepl("_down\\.tsv$", all_files)]
  }

  all_files
}

is_valid_contrast_dir <- function(path) {
  if (!dir.exists(path)) return(FALSE)
  b <- basename(path)
  if (b %in% c("background", "inputs")) return(FALSE)
  TRUE
}

discover_auto_files <- function(enrich_dir, target_type, direction_filter, contrast_filter) {
  target_files <- get_standard_filenames(target_type, direction_filter)

  subdirs <- list.dirs(enrich_dir, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[vapply(subdirs, is_valid_contrast_dir, logical(1))]

  if (length(contrast_filter) > 0) {
    subdirs <- subdirs[basename(subdirs) %in% contrast_filter]
  }

  out <- character(0)
  for (dir_i in subdirs) {
    files_i <- file.path(dir_i, target_files)
    files_i <- files_i[file.exists(files_i)]
    out <- c(out, files_i)
  }

  unique(normalize_path2(out))
}

discover_manual_files <- function(manual_input, root_dir, target_type, direction_filter) {
  if (!nzchar(manual_input)) {
    stop("RUN_MODE='manual' 时，MANUAL_INPUT 不能为空。", call. = FALSE)
  }

  target_files <- get_standard_filenames(target_type, direction_filter)

  manual_path <- manual_input
  if (!grepl("^/", manual_path)) {
    manual_path <- file.path(root_dir, manual_path)
  }
  manual_path <- normalize_path2(manual_path)

  if (file.exists(manual_path) && !dir.exists(manual_path)) {
    if (!(basename(manual_path) %in% target_files)) {
      stop("手动指定的文件不是标准富集结果文件：", manual_path, call. = FALSE)
    }
    return(manual_path)
  }

  if (dir.exists(manual_path)) {
    files_i <- file.path(manual_path, target_files)
    files_i <- files_i[file.exists(files_i)]
    return(normalize_path2(files_i))
  }

  stop("MANUAL_INPUT 不存在：", manual_input, call. = FALSE)
}

make_output_prefix <- function(file_path) {
  contrast_name <- basename(dirname(file_path))
  file_base <- file_path_sans_ext(basename(file_path))
  paste0(contrast_name, "__", file_base)
}

# ==============================================================================
# 3. 主程序
# ==============================================================================

COLOR_MAP_OPTION <- sanitize_viridis_option(COLOR_MAP_OPTION)

PROJECT_ROOT_FINAL <- detect_project_root(PROJECT_ROOT, ENRICH_DIR_REL)
ENRICH_DIR <- normalize_path2(file.path(PROJECT_ROOT_FINAL, ENRICH_DIR_REL))
OUT_DIR    <- normalize_path2(file.path(PROJECT_ROOT_FINAL, OUT_DIR_REL))

if (!dir.exists(ENRICH_DIR)) {
  stop("[ERROR] 找不到富集结果目录：", ENRICH_DIR, call. = FALSE)
}
ensure_dir(OUT_DIR)

log_msg("项目根目录：", PROJECT_ROOT_FINAL)
log_msg("富集输入目录：", ENRICH_DIR)
log_msg("图片输出目录：", OUT_DIR)
log_msg("运行模式：", RUN_MODE)

files <- character(0)

if (RUN_MODE == "auto") {
  files <- discover_auto_files(
    enrich_dir = ENRICH_DIR,
    target_type = TARGET_TYPE,
    direction_filter = DIRECTION_FILTER,
    contrast_filter = CONTRAST_FILTER
  )
} else if (RUN_MODE == "manual") {
  files <- discover_manual_files(
    manual_input = MANUAL_INPUT,
    root_dir = PROJECT_ROOT_FINAL,
    target_type = TARGET_TYPE,
    direction_filter = DIRECTION_FILTER
  )
} else {
  stop("RUN_MODE 只能是 'auto' 或 'manual'。", call. = FALSE)
}

if (length(files) == 0) {
  stop("[ERROR] 未找到可绘图的标准富集结果文件。", call. = FALSE)
}

log_msg("共发现待绘图文件数：", length(files))

# ==============================================================================
# 4. 主程序循环 (MAIN LOOP)
#    >>> 下面核心绘图逻辑基本沿用原始脚本 <<<
# ==============================================================================

for (f_path in files) {

  fname  <- basename(f_path)
  f_base <- file_path_sans_ext(fname)
  out_prefix <- make_output_prefix(f_path)

  cat(sprintf(">> 正在处理文件: %s\n", f_path))

  df_raw <- tryCatch({
    read_tsv(f_path, show_col_types = FALSE)
  }, error = function(e) NULL)

  if (is.null(df_raw) || nrow(df_raw) == 0) {
    cat("   [跳过] 文件为空或读取失败。\n")
    next
  }

  # Cleaning
  df_plot <- df_raw %>%
    mutate(
      p_adjust         = as.numeric(p_adjust),
      gene_count       = as.numeric(gene_count),
      GeneRatioNumeric = parse_ratio_hybrid(gene_ratio),
      TermWrapped      = str_wrap(term_name, width = MANUAL_WRAP_WIDTH)
    ) %>%
    filter(!is.na(p_adjust), GeneRatioNumeric > 0)

  if (nrow(df_plot) == 0) {
    cat("   [跳过] 过滤后无有效数据。\n")
    next
  }

  # Sorting & Slicing
  unique_ontologies <- if ("ontology" %in% colnames(df_plot)) unique(df_plot$ontology) else c()
  is_multi_ontology <- length(unique_ontologies) > 1
  target_n <- if (is_multi_ontology) 10 else 20

  if (is_multi_ontology) {
    df_ready <- df_plot %>%
      group_by(ontology) %>%
      arrange(p_adjust) %>%
      slice_head(n = target_n) %>%
      ungroup() %>%
      arrange(GeneRatioNumeric) %>%
      mutate(
        TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)),
        negLog10FDR = -log10(p_adjust)
      )
  } else {
    df_ready <- df_plot %>%
      arrange(p_adjust) %>%
      slice_head(n = target_n) %>%
      arrange(GeneRatioNumeric) %>%
      mutate(
        TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)),
        negLog10FDR = -log10(p_adjust)
      )
  }

  final_count <- nrow(df_ready)
  if (final_count == 0) {
    cat("   [跳过] 取前 N 后无可绘制数据。\n")
    next
  }

  # =====================================================
  # [尺寸计算] 沿用原脚本逻辑
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
    n = VIRIDIS_N,
    option = COLOR_MAP_OPTION,
    direction = -1
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
      title = if (SHOW_PLOT_TITLE) gsub("_", " ", out_prefix) else NULL
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

  if (WRITE_PDF) {
    ggsave(
      file.path(OUT_DIR, paste0(out_prefix, ".bubble.pdf")),
      p,
      width  = final_width_inch,
      height = final_height_inch,
      device = cairo_pdf
    )
  }

  if (WRITE_PNG) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggsave(
        file.path(OUT_DIR, paste0(out_prefix, ".bubble.png")),
        p,
        width  = final_width_inch,
        height = final_height_inch,
        dpi    = PNG_DPI,
        bg     = "white",
        device = ragg::agg_png
      )
    } else {
      ggsave(
        file.path(OUT_DIR, paste0(out_prefix, ".bubble.png")),
        p,
        width  = final_width_inch,
        height = final_height_inch,
        dpi    = PNG_DPI,
        bg     = "white"
      )
    }
  }

  cat("   [成功] 图片已保存。\n")
}

cat("所有任务处理完成。\n")
