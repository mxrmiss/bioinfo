#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: viz_enrichment_bubble_v5.1.R (Final Color Restored)
# Description: 
#   - [配色] 已恢复为你最初代码中的配色 (Blue -> Orange -> Red)
#   - [逻辑] 保持 V5.0 的完美宽度算法 (气泡区大小固定，背景随文字伸缩)
#   - [配置] 全中文详细注释
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tools)
  library(scales)
})

# ==============================================================================
# 1. 用户配置区域 (USER CONFIGURATION)
#    >>> 请在此处修改参数 <<<
# ==============================================================================

# ------------------------------------------------------------------------------
# [A] 尺寸控制 (DIMENSIONS)
# ------------------------------------------------------------------------------

# 1. 高度设置 (推荐锁定)
# --------------------
# TRUE  = [锁定模式] 强制使用下方 FIXED_HEIGHT_VAL 的高度 (推荐，保证气泡纵向不被拉伸)
# FALSE = [自动模式] 根据气泡数量自动伸缩高度
USE_FIXED_HEIGHT    <- TRUE

# [参数] 锁定模式下的高度 (单位: inch)
# 默认 9.5 是你最习惯的标准高度
FIXED_HEIGHT_VAL    <- 9.5     

# [参数] 单面图的宽高比 (Aspect Ratio)
# 这个参数决定了气泡绘图区的胖瘦。1.7 代表宽度是高度的 1.7 倍。
# 无论文字有多长，气泡区域永远保持这个比例，确保视觉统一。
SINGLE_ASPECT_RATIO <- 1.7     

# 2. 自动高度参数 (仅在 USE_FIXED_HEIGHT = FALSE 时生效)
# --------------------
ROW_HEIGHT_INCH     <- 0.35    # 每一行气泡的高度
BASE_HEIGHT_INCH    <- 1.5     # 基础边距

# 3. 分面图专用尺寸 (仅在数据包含多个 Ontology 时生效)
# --------------------
FACET_WIDTH_INCH    <- 12.0
FACET_HEIGHT_INCH   <- 10.0

# 4. 文字排版
# --------------------
MANUAL_WRAP_WIDTH   <- 60      # 通路名称超过多少个字符自动换行

# ------------------------------------------------------------------------------
# [B] 外观与配色 (AESTHETICS & COLORS)
# ------------------------------------------------------------------------------

# 是否显示图片主标题？ (TRUE=显示, FALSE=隐藏/投稿模式)
SHOW_PLOT_TITLE     <- TRUE

# 是否显示垂直网格线？ (TRUE=显示, FALSE=白底)
SHOW_GRID           <- FALSE   

# 颜色标尺模式
# TRUE  = [动态范围] 根据当前数据的 P值 极值自动调整颜色深浅 (对比度最强)
# FALSE = [固定范围] 强制使用下方 GLOBAL_P_MIN/MAX
USE_DYNAMIC_COLOR_SCALE <- TRUE
GLOBAL_P_MIN <- 1e-10
GLOBAL_P_MAX <- 0.05

# --- [配色方案] ---
# 这里恢复了你一开始提供的配色代码
# Low (不显著) -> Mid (过渡) -> High (极显著)
COLOR_LOW    <- "#4DBBD5"  # Cyan (蓝青色)
COLOR_MID    <- "#F39B7F"  # Light Orange (浅橙色)
COLOR_HIGH   <- "#E64B35"  # Maroon (深红/栗色)

# --- [字体设置] ---
FONT_FAMILY  <- "sans"     # 字体 (Arial/Helvetica)
BASE_SIZE    <- 14         # 基础字号

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
  return(10^log_seq)
}

force_integer_breaks <- function(n = 5) {
  function(x) {
    rng <- range(x, na.rm = TRUE)
    seq_breaks <- seq(rng[1], rng[2], length.out = n)
    unique(round(seq_breaks))
  }
}

# ==============================================================================
# 3. 主程序循环 (MAIN LOOP)
# ==============================================================================

init_template_system()
files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("[ERROR] input 文件夹下未找到 .tsv 文件。")

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
  target_n <- if(is_multi_ontology) 10 else 20
  
  df_ready <- df_plot %>%
    group_by(if(is_multi_ontology) ontology else NULL) %>%
    arrange(p_adjust) %>%
    slice_head(n = target_n) %>%
    ungroup() %>%
    arrange(GeneRatioNumeric) %>%
    mutate(TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)))
  
  final_count <- nrow(df_ready)
  
  # =====================================================
  # [尺寸计算] 你的“黄金逻辑”已回归
  # 逻辑：Grid大小固定 + Text长短自适应 = 完美统一的气泡区
  # =====================================================
  
  if (is_multi_ontology) {
    # 分面图维持固定大尺寸
    final_width_inch  <- FACET_WIDTH_INCH
    final_height_inch <- FACET_HEIGHT_INCH
  } else {
    # 单面图逻辑
    
    # 1. 确定高度
    if (USE_FIXED_HEIGHT) {
      final_height_inch <- FIXED_HEIGHT_VAL
    } else {
      final_height_inch <- (final_count * ROW_HEIGHT_INCH) + BASE_HEIGHT_INCH
      final_height_inch <- max(final_height_inch, 3.0) 
    }
    
    # 2. 计算气泡绘图区的宽度 (Grid Width)
    # 只要高度不变，Ratio不变，Grid Width 就永远不变。
    # 这保证了所有图片里的气泡区域看起来是一模一样的。
    grid_width_inch <- final_height_inch / SINGLE_ASPECT_RATIO
    
    # 3. 计算文字区宽度 (Text Width)
    # 这里的动态计算是为了让画布适应文字，而不是改变气泡区
    max_char_len    <- max(nchar(as.character(df_ready$TermWrapped)), na.rm = TRUE)
    text_width_inch <- max_char_len * 0.09 + 0.5
    
    # 4. 图例宽度
    legend_width_inch <- 2.0
    
    # 5. 最终总宽度
    final_width_inch <- text_width_inch + grid_width_inch + legend_width_inch
  }

  cat(sprintf("   [尺寸] %.2f x %.2f inch (Grid宽: %.2f) | 模式: %s\n", 
              final_width_inch, final_height_inch, 
              if(!is_multi_ontology) final_height_inch/SINGLE_ASPECT_RATIO else 0,
              if(is_multi_ontology) "分面" else "单面"))

  # --- Color ---
  if (USE_DYNAMIC_COLOR_SCALE) {
    current_p_min <- min(df_ready$p_adjust, na.rm = TRUE)
    current_p_max <- max(df_ready$p_adjust, na.rm = TRUE)
    if (current_p_min == current_p_max) current_p_max <- if(current_p_max == 0) 0.01 else current_p_max * 10
    plot_breaks <- generate_fixed_breaks(current_p_min, current_p_max, 5)
    plot_limits <- c(current_p_min, current_p_max)
  } else {
    plot_breaks <- generate_fixed_breaks(GLOBAL_P_MIN, GLOBAL_P_MAX, 5)
    plot_limits <- c(GLOBAL_P_MIN, GLOBAL_P_MAX)
  }
  
  # --- Plotting ---
  p <- ggplot(df_ready, aes(x = GeneRatioNumeric, y = TermWrapped)) +
    
    geom_point(aes(size = gene_count, color = p_adjust), alpha = 0.9) +
    
    scale_color_gradientn(
      colours = c(COLOR_HIGH, COLOR_MID, COLOR_LOW), # 恢复了 Red -> Orange -> Blue
      trans   = "log10", breaks = plot_breaks, limits = plot_limits,
      labels  = scientific_p_formatter, oob = scales::squish, name = "P.adjust\n"
    ) +
    
    scale_size_continuous(range = c(3, 8), breaks = force_integer_breaks(4), name = "Count\n") +
    
    labs(x = "GeneRatio", y = NULL, title = if(SHOW_PLOT_TITLE) gsub("_", " ", f_base) else NULL) +
    
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      # [核心] 锁定气泡区比例
      aspect.ratio = if (!is_multi_ontology) SINGLE_ASPECT_RATIO else NULL,
      
      panel.grid.major.y = if (SHOW_GRID) element_line(colour = "grey92") else element_blank(),
      panel.grid.major.x = element_blank(), 
      panel.grid.minor = element_blank(),
      
      axis.text  = element_text(colour = "black", size = BASE_SIZE),
      axis.title = element_text(size = BASE_SIZE, face="bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = BASE_SIZE + 3),
      
      legend.key.height = unit(0.8, "cm")
    )
  
  if (is_multi_ontology) {
    p <- p + facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
      theme(strip.background = element_rect(fill = "grey95", color = NA))
  }
  
  ggsave(file.path(OUT_DIR, paste0(f_base, ".bubble.pdf")), p, width=final_width_inch, height=final_height_inch, device=cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(f_base, ".bubble.png")), p, width=final_width_inch, height=final_height_inch, dpi=300, bg="white")
  cat("   [成功] 图片已保存。\n")
}
cat("所有任务处理完成。\n")
