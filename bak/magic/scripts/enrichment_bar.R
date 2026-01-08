#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: viz_enrichment_bar_v4.2.R (True Original Logic)
# Description: 
#   - [核心修复] 还原了最初代码的“宽高比约束”算法 (Aspect Ratio)
#   - [效果] 绘图区(Grid)形状恒定，不再被强行拉伸，画布宽度随文字自动适应
#   - [配色] Nature Comm 蓝紫冷色调
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
# ==============================================================================

# ------------------------------------------------------------------------------
# [A] 尺寸控制 (回归原始逻辑)
# ------------------------------------------------------------------------------
# 1. 高度控制
# 强烈建议 TRUE (锁死)，这也是你最初代码的默认行为
USE_FIXED_HEIGHT    <- TRUE
FIXED_HEIGHT_VAL    <- 9.5     # [原版标准] 
SINGLE_ASPECT_RATIO <- 2.0     # [原版标准] 宽高比 (决定了绘图区的胖瘦)

# 2. 自动高度 (仅在 USE_FIXED_HEIGHT = FALSE 时生效)
ROW_HEIGHT_INCH     <- 0.35
BASE_HEIGHT_INCH    <- 1.5

# 3. 分面图专用尺寸
FACET_WIDTH_INCH    <- 12.0
FACET_HEIGHT_INCH   <- 10.0

# 4. 文字换行
MANUAL_WRAP_WIDTH   <- 60

# ------------------------------------------------------------------------------
# [B] 外观与配色
# ------------------------------------------------------------------------------
SHOW_PLOT_TITLE   <- TRUE
SHOW_GRID         <- FALSE 

USE_DYNAMIC_COLOR_SCALE <- TRUE 
GLOBAL_P_MIN <- 1e-10
GLOBAL_P_MAX <- 0.05

# --- Nature Comm Style (Blue-Purple) ---
# 适合 "Down-regulated"
COLOR_LOW    <- "#BDD7E7"  # 浅灰蓝
COLOR_MID    <- "#6A51A3"  # 中紫
COLOR_HIGH   <- "#08306B"  # 深靛蓝

# 如果你想用回红蓝配色，取消下面注释即可：
# COLOR_LOW    <- "#4DBBD5"; COLOR_MID <- "#F39B7F"; COLOR_HIGH <- "#E64B35"

FONT_FAMILY  <- "sans"
BASE_SIZE    <- 14     

# ------------------------------------------------------------------------------
# [C] 路径
# ------------------------------------------------------------------------------
IN_DIR        <- "input"
SCRIPT_PREFIX <- "bar_enrichment"
OUT_DIR       <- file.path("output", SCRIPT_PREFIX)
TEMPLATE_DIR  <- "templates"

# ==============================================================================
# 2. 辅助函数
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

# ==============================================================================
# 3. 主程序循环
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
  
  # Sorting
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
  # [逻辑回归] 1. 动态柱宽 (Anti-bloat Logic)
  # =====================================================
  # 恢复你代码中那个精妙的公式：通路少时柱子自动变细
  dynamic_bar_width <- min(0.7, 0.7 * (final_count / 15))
  
  # =====================================================
  # [逻辑回归] 2. 尺寸计算 (Grid Ratio Logic)
  # =====================================================
  if (is_multi_ontology) {
    # 分面图维持固定大尺寸
    final_width_inch  <- FACET_WIDTH_INCH
    final_height_inch <- FACET_HEIGHT_INCH
    dynamic_bar_width <- 0.7 # 分面图强制柱宽正常
  } else {
    # [单面图]
    
    # A. 确定高度
    if (USE_FIXED_HEIGHT) {
      final_height_inch <- FIXED_HEIGHT_VAL
    } else {
      final_height_inch <- (final_count * ROW_HEIGHT_INCH) + BASE_HEIGHT_INCH
      final_height_inch <- max(final_height_inch, 3.0) 
    }
    
    # B. 确定绘图区(Grid)宽度 -> 由高度反推，保持比例!
    # 这步是之前漏掉的关键，保证图看起来永远是 9.5 : 5.6 的长方形
    grid_width_inch <- final_height_inch / SINGLE_ASPECT_RATIO
    
    # C. 确定文字宽度 -> 随文字长短延伸
    max_char_len    <- max(nchar(as.character(df_ready$TermWrapped)), na.rm = TRUE)
    text_width_inch <- max_char_len * 0.09 + 0.5
    
    # D. 图例宽度
    legend_width_inch <- 2.0
    
    # E. 总宽度
    final_width_inch <- text_width_inch + grid_width_inch + legend_width_inch
  }
  
  cat(sprintf("   [尺寸] Total: %.2f x %.2f inch | Grid宽: %.2f | 柱宽系数: %.2f\n", 
              final_width_inch, final_height_inch, 
              if(!is_multi_ontology) final_height_inch/SINGLE_ASPECT_RATIO else 0,
              dynamic_bar_width))

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
    
    # 柱子
    geom_col(aes(fill = p_adjust), width = dynamic_bar_width) +
    
    # 标签
    geom_text(aes(label = gene_count), 
              hjust = -0.5, 
              size = 4, 
              color = "black") +
    
    scale_fill_gradientn(
      colours = c(COLOR_HIGH, COLOR_MID, COLOR_LOW), 
      trans   = "log10", breaks = plot_breaks, limits = plot_limits,
      labels  = scientific_p_formatter, oob = scales::squish, name = "P.adjust\n"
    ) +
    
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    
    labs(x = "GeneRatio", y = NULL, title = if(SHOW_PLOT_TITLE) gsub("_", " ", f_base) else NULL) +
    
    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      # [回归] 确保 Grid 比例固定，不会被拉伸变形
      aspect.ratio = if (!is_multi_ontology) SINGLE_ASPECT_RATIO else NULL,
      
      panel.grid.major.x = if (SHOW_GRID) element_line(colour = "white", linewidth = 0.5) else element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = if (SHOW_GRID) element_rect(fill = "#F0F0F0") else element_rect(fill = "white"),
      
      axis.text  = element_text(colour = "black", size = BASE_SIZE),
      axis.title = element_text(size = BASE_SIZE, face="bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = BASE_SIZE + 3),
      
      legend.key.height = unit(0.8, "cm")
    )
  
  if (is_multi_ontology) {
    p <- p + facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
      theme(strip.background = element_rect(fill = "grey95", color = NA))
  }
  
  ggsave(file.path(OUT_DIR, paste0(f_base, ".bar.pdf")), p, width=final_width_inch, height=final_height_inch, device=cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(f_base, ".bar.png")), p, width=final_width_inch, height=final_height_inch, dpi=300, bg="white")
  cat("   [成功] 图片已保存。\n")
}
cat("所有任务处理完成。\n")
