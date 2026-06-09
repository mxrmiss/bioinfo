#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tools)
  library(scales)
  library(grid)
})

# ==============================================================================
# 1. 用户配置区域
# ==============================================================================

# ------------------------------------------------------------------------------
# [A] 尺寸控制
# ------------------------------------------------------------------------------

USE_FIXED_HEIGHT    <- TRUE

FIXED_HEIGHT_VAL    <- 11.5
FIXED_WIDTH_VAL     <- 13

SINGLE_ASPECT_RATIO <- 2.2

ROW_HEIGHT_INCH     <- 0.35
BASE_HEIGHT_INCH    <- 1.5

FACET_WIDTH_INCH    <- 12.0
FACET_HEIGHT_INCH   <- 10.0

MANUAL_WRAP_WIDTH   <- 40

# ------------------------------------------------------------------------------
# [B] 外观设置
# ------------------------------------------------------------------------------

SHOW_PLOT_TITLE     <- FALSE

SHOW_GRID           <- FALSE
USE_DARK_GRID       <- FALSE
GRID_COLOR_LIGHT    <- "grey92"
GRID_COLOR_DARK     <- "grey80"

FONT_FAMILY         <- "Arial"
BASE_SIZE           <- 22.8

POINT_ALPHA         <- 0.90
POINT_SIZE_MIN      <- 3.5
POINT_SIZE_MAX      <- 8.5

# ------------------------------------------------------------------------------
# [C] 横坐标设置
# ------------------------------------------------------------------------------

USE_FIXED_X_LIMITS  <- FALSE
X_LIMIT_MIN         <- 0.00
X_LIMIT_MAX         <- 0.10
X_BREAKS            <- c(0.02, 0.04, 0.06, 0.08)

# ------------------------------------------------------------------------------
# [D] FDR 颜色设置
# ------------------------------------------------------------------------------

USE_DYNAMIC_COLOR_SCALE <- TRUE

GLOBAL_P_MIN <- 1e-10
GLOBAL_P_MAX <- 0.05

CUSTOM_FDR_COLORS <- c(
  "#C8E8C8",
  "#68CDD2",
  "#4DBBD5",
  "#3C5488"
)

COLOR_N <- 256

# ------------------------------------------------------------------------------
# [E] 路径设置
# ------------------------------------------------------------------------------

IN_DIR        <- "input"
SCRIPT_PREFIX <- "bubble_enrichment"
OUT_DIR       <- file.path("output", SCRIPT_PREFIX)
TEMPLATE_DIR  <- "templates"

# ==============================================================================
# 2. 辅助函数
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) {
    dir.create(TEMPLATE_DIR, recursive = TRUE)
  }

  if (!dir.exists(OUT_DIR)) {
    dir.create(OUT_DIR, recursive = TRUE)
  }
}

parse_ratio_hybrid <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") {
      return(0)
    }

    if (is.numeric(val)) {
      return(as.numeric(val))
    }

    val_str <- as.character(val)

    if (grepl("/", val_str, fixed = TRUE)) {
      parts <- strsplit(val_str, "/", fixed = TRUE)[[1]]

      if (length(parts) != 2) {
        return(0)
      }

      numerator <- suppressWarnings(as.numeric(parts[1]))
      denominator <- suppressWarnings(as.numeric(parts[2]))

      if (is.na(numerator) || is.na(denominator) || denominator == 0) {
        return(0)
      }

      return(numerator / denominator)
    }

    suppressWarnings(as.numeric(val_str))
  })
}

generate_fixed_breaks <- function(min_val, max_val, n) {
  if (min_val <= 0) {
    min_val <- 1e-300
  }

  if (max_val <= 0) {
    max_val <- min_val * 10
  }

  if (min_val == max_val) {
    min_val <- min_val / 10
    max_val <- max_val * 10
  }

  log_seq <- seq(log10(max_val), log10(min_val), length.out = n)
  10^log_seq
}

force_integer_breaks <- function(n = 5) {
  function(x) {
    rng <- range(x, na.rm = TRUE)

    if (!all(is.finite(rng))) {
      return(NULL)
    }

    if (rng[1] == rng[2]) {
      return(round(rng[1]))
    }

    seq_breaks <- seq(rng[1], rng[2], length.out = n)
    unique(round(seq_breaks))
  }
}

safe_neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  p[p <= 0 | is.na(p)] <- NA_real_
  -log10(p)
}

make_color_vector <- function(colors, n = 256) {
  colorRampPalette(colors)(n)
}

# ==============================================================================
# 3. 主程序
# ==============================================================================

init_template_system()

files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("[ERROR] input 文件夹下未找到 .tsv 文件。")
}

color_vec <- make_color_vector(CUSTOM_FDR_COLORS, COLOR_N)

for (f_path in files) {

  fname  <- basename(f_path)
  f_base <- file_path_sans_ext(fname)

  cat(sprintf(">> 正在处理文件: %s\n", fname))

  df_raw <- tryCatch(
    read_tsv(f_path, show_col_types = FALSE),
    error = function(e) NULL
  )

  if (is.null(df_raw) || nrow(df_raw) == 0) {
    warning(sprintf("[WARN] 文件为空或读取失败: %s，已跳过。", fname))
    next
  }

  required_cols <- c("p_adjust", "gene_count", "gene_ratio", "term_name")
  missing_cols <- setdiff(required_cols, colnames(df_raw))

  if (length(missing_cols) > 0) {
    warning(sprintf(
      "[WARN] 文件 %s 缺少必要列: %s，已跳过。",
      fname,
      paste(missing_cols, collapse = ", ")
    ))
    next
  }

  df_plot <- df_raw %>%
    mutate(
      p_adjust         = as.numeric(p_adjust),
      gene_count       = as.numeric(gene_count),
      GeneRatioNumeric = parse_ratio_hybrid(gene_ratio),
      TermWrapped      = str_wrap(term_name, width = MANUAL_WRAP_WIDTH),
      negLog10FDR      = safe_neglog10(p_adjust)
    ) %>%
    filter(
      !is.na(p_adjust),
      p_adjust > 0,
      !is.na(gene_count),
      gene_count > 0,
      GeneRatioNumeric > 0,
      !is.na(negLog10FDR)
    )

  if (nrow(df_plot) == 0) {
    warning(sprintf("[WARN] 文件 %s 清洗后无可绘制数据，已跳过。", fname))
    next
  }

  unique_ontologies <- if ("ontology" %in% colnames(df_plot)) {
    unique(df_plot$ontology[!is.na(df_plot$ontology)])
  } else {
    character(0)
  }

  is_multi_ontology <- length(unique_ontologies) > 1
  target_n <- if (is_multi_ontology) 10 else 20

  if (is_multi_ontology) {
    df_ready <- df_plot %>%
      group_by(ontology) %>%
      arrange(p_adjust, .by_group = TRUE) %>%
      slice_head(n = target_n) %>%
      ungroup() %>%
      arrange(GeneRatioNumeric) %>%
      mutate(
        TermWrapped = factor(TermWrapped, levels = unique(TermWrapped))
      )
  } else {
    df_ready <- df_plot %>%
      arrange(p_adjust) %>%
      slice_head(n = target_n) %>%
      arrange(GeneRatioNumeric) %>%
      mutate(
        TermWrapped = factor(TermWrapped, levels = unique(TermWrapped))
      )
  }

  final_count <- nrow(df_ready)

  if (is_multi_ontology) {
    final_width_inch  <- FACET_WIDTH_INCH
    final_height_inch <- FACET_HEIGHT_INCH
  } else {
    if (USE_FIXED_HEIGHT) {
      final_height_inch <- FIXED_HEIGHT_VAL
      final_width_inch  <- FIXED_WIDTH_VAL
    } else {
      final_height_inch <- final_count * ROW_HEIGHT_INCH + BASE_HEIGHT_INCH
      final_height_inch <- max(final_height_inch, 3.0)

      grid_width_inch <- final_height_inch / SINGLE_ASPECT_RATIO

      max_char_len <- max(nchar(as.character(df_ready$TermWrapped)), na.rm = TRUE)
      text_width_inch <- max_char_len * 0.09 + 0.5

      legend_width_inch <- 2.0

      final_width_inch <- text_width_inch + grid_width_inch + legend_width_inch
    }
  }

  cat(sprintf(
    "   [尺寸] %.2f x %.2f inch | 模式: %s | 条目数: %d\n",
    final_width_inch,
    final_height_inch,
    if (is_multi_ontology) "分面" else "单面",
    final_count
  ))

  if (USE_DYNAMIC_COLOR_SCALE) {
    current_p_min <- min(df_ready$p_adjust, na.rm = TRUE)
    current_p_max <- max(df_ready$p_adjust, na.rm = TRUE)

    if (current_p_min == current_p_max) {
      current_p_min <- current_p_min / 10
      current_p_max <- current_p_max * 10
    }

    plot_breaks <- generate_fixed_breaks(current_p_min, current_p_max, 5)
    plot_limits <- range(safe_neglog10(c(current_p_min, current_p_max)), na.rm = TRUE)
  } else {
    plot_breaks <- generate_fixed_breaks(GLOBAL_P_MIN, GLOBAL_P_MAX, 5)
    plot_limits <- range(safe_neglog10(c(GLOBAL_P_MIN, GLOBAL_P_MAX)), na.rm = TRUE)
  }

  p <- ggplot(df_ready, aes(x = GeneRatioNumeric, y = TermWrapped)) +
    geom_point(
      aes(size = gene_count, color = negLog10FDR),
      alpha = POINT_ALPHA
    ) +
    scale_color_gradientn(
      colours = color_vec,
      breaks  = sort(safe_neglog10(plot_breaks)),
      limits  = plot_limits,
      labels  = label_number(accuracy = 0.1),
      oob     = scales::squish,
      name    = "-log10(FDR)\n"
    ) +
    scale_size_continuous(
      range  = c(POINT_SIZE_MIN, POINT_SIZE_MAX),
      breaks = force_integer_breaks(4),
      name   = "Count\n"
    ) +
    labs(
      x     = "GeneRatio",
      y     = NULL,
      title = if (SHOW_PLOT_TITLE) gsub("_", " ", f_base) else NULL
    ) +
    theme_bw(
      base_size = BASE_SIZE,
      base_family = FONT_FAMILY
    )

  if (USE_FIXED_X_LIMITS) {
    p <- p +
      scale_x_continuous(
        limits = c(X_LIMIT_MIN, X_LIMIT_MAX),
        breaks = X_BREAKS,
        labels = label_number(accuracy = 0.01),
        expand = expansion(mult = c(0.02, 0.05))
      )
  } else {
    p <- p +
      scale_x_continuous(
        labels = label_number(accuracy = 0.01),
        expand = expansion(mult = c(0.03, 0.08))
      )
  }

  grid_col <- if (USE_DARK_GRID) GRID_COLOR_DARK else GRID_COLOR_LIGHT

  p <- p +
    theme(
      aspect.ratio = if (!is_multi_ontology) SINGLE_ASPECT_RATIO else NULL,

      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),

      panel.grid.major.x = if (SHOW_GRID) {
        element_line(colour = grid_col, linewidth = 0.3)
      } else {
        element_blank()
      },
      panel.grid.major.y = if (SHOW_GRID) {
        element_line(colour = grid_col, linewidth = 0.3)
      } else {
        element_blank()
      },
      panel.grid.minor.x = if (SHOW_GRID) {
        element_line(colour = grid_col, linewidth = 0.25)
      } else {
        element_blank()
      },
      panel.grid.minor.y = if (SHOW_GRID) {
        element_line(colour = grid_col, linewidth = 0.25)
      } else {
        element_blank()
      },

      axis.text = element_text(
        colour = "black",
        size = BASE_SIZE
      ),
      axis.title = element_text(
        size = BASE_SIZE,
        face = "bold",
        colour = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = BASE_SIZE + 3,
        colour = "black"
      ),

      legend.title = element_text(
        size = BASE_SIZE - 1,
        face = "plain",
        colour = "black"
      ),
      legend.text = element_text(
        size = BASE_SIZE - 1,
		face = "plain",
        colour = "black"
      ),
      legend.key = element_blank(),
      legend.key.height = unit(0.8, "cm"),
      legend.box.spacing = unit(0.3, "cm"),

      plot.margin = margin(t = 8, r = 10, b = 8, l = 8)
    )

  if (is_multi_ontology) {
    p <- p +
      facet_grid(
        ontology ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      theme(
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(
          size = BASE_SIZE,
          face = "bold",
          colour = "black"
        )
      )
  }

  ggsave(
    file.path(OUT_DIR, paste0(f_base, ".bubble.pdf")),
    p,
    width  = final_width_inch,
    height = final_height_inch,
    device = cairo_pdf,
    bg     = "white"
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
