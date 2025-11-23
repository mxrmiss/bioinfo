#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: viz_enrichment_bubble.R (v2.5 Font Control Edition)
# Description:
#   - [Feat] Toggle Plot Title ON/OFF for publication.
#   - [Feat] Auto-removes underscores from titles (e.g. "GO_BP" -> "GO BP").
#   - [Feat] Centralized font size control for journal submission.
#   - [System] Template Guard & Dual Layout Standard included.
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
# 1. USER CONFIGURATION
# ==============================================================================

# --- [NEW] Title Control ---
SHOW_PLOT_TITLE     <- TRUE   # TRUE = Show Title; FALSE = No Title (Clean for Paper)

# --- A. Single Plot Standards (Dynamic Width) ---
SINGLE_HEIGHT_INCH  <- 9.5
SINGLE_ASPECT_RATIO <- 1.3
TOP_N_SINGLE        <- 20

# --- B. Facet Plot Standards (Fixed Dimensions) ---
FACET_WIDTH_INCH    <- 12.0
FACET_HEIGHT_INCH   <- 10.0
TOP_N_FACET         <- 10

# --- C. Global Legend Standards ---
GLOBAL_P_MIN <- 1e-10
GLOBAL_P_MAX <- 0.05
LEGEND_TICKS <- 5
COUNT_TICKS  <- 5

# --- D. Aesthetics & Sorting ---
SORT_BY      <- "GeneRatio"
COLOR_LOW    <- "blue"
COLOR_MID    <- "orange"
COLOR_HIGH   <- "red"
FONT_FAMILY  <- "sans"

# =======================
# [NEW] Font Size Control
# =======================
# 建议值（投稿友好）：
# - 坐标刻度：12–14
# - 坐标标题：14–16
# - 图例文字：11–12
# - 图例标题：12–13
# - 分面标题：13–14
# - 主标题：16–18
BASE_SIZE         <- 14  # 全局基础字号（theme_bw 的 base_size）
AXIS_TEXT_SIZE    <- 13
AXIS_TITLE_SIZE   <- 14
LEGEND_TEXT_SIZE  <- 13
LEGEND_TITLE_SIZE <- 14
STRIP_TEXT_SIZE   <- 13
TITLE_SIZE        <- 17

WRAP_WIDTH   <- 60
SHOW_GRID    <- TRUE

# Paths
IN_DIR       <- "input"

# 输出到脚本同名子文件夹
SCRIPT_PREFIX <- "bubble_enrichment"
OUT_DIR       <- file.path("output", SCRIPT_PREFIX)

TEMPLATE_DIR <- "templates"
TEMPLATE_FILE<- "enrichment_input_template.tsv"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) {
    dir.create(TEMPLATE_DIR, recursive = TRUE)
  }
  tmpl_path <- file.path(TEMPLATE_DIR, TEMPLATE_FILE)
  if (!file.exists(tmpl_path)) {
    example_data <- data.frame(
      term_name  = c("Example Term (Fraction Format)", "Example Term (Decimal Format)"),
      gene_ratio = c("25/100", "0.05"),
      p_adjust   = c("1.2e-10", "0.004"),
      gene_count = c(25, 10),
      ontology   = c("BP", "NA"),
      stringsAsFactors = FALSE
    )
    write_tsv(example_data, tmpl_path)
  }
}

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

parse_ratio_hybrid <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") return(0)
    if (is.numeric(val)) return(val)
    val_str <- as.character(val)
    if (grepl("/", val_str)) {
      tryCatch({ eval(parse(text = val_str)) }, error = function(e) 0)
    } else {
      res <- as.numeric(val_str)
      if (is.na(res)) 0 else res
    }
  })
}

scientific_p_formatter <- function(x) {
  formatC(x, format = "e", digits = 1)
}

generate_fixed_breaks <- function(min_val, max_val, n) {
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
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic Bubble Plotter v2.5 (Font Control)\n")
cat("================================================================\n")

init_template_system()

files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("[ERROR] No .tsv files found in magic/input/")
}

cat(sprintf("Detected %d input files. Processing...\n", length(files)))
cat(sprintf("Show Title: %s\n", SHOW_PLOT_TITLE))
cat("----------------------------------------------------------------\n")

for (f_path in files) {

  fname  <- basename(f_path)
  f_base <- file_path_sans_ext(fname)
  cat(sprintf(">> File: %s\n", fname))

  # --- Step 1: Read ---
  df_raw <- tryCatch({
    read_tsv(f_path, show_col_types = FALSE, progress = FALSE)
  }, error = function(e) NULL)

  if (is.null(df_raw) || nrow(df_raw) == 0) {
    cat("   [WARNING] Empty or unreadable file. Skipping.\n"); next
  }

  # --- Step 2: Header Check ---
  required_cols <- c("term_name", "gene_ratio", "p_adjust", "gene_count")
  missing_cols  <- setdiff(required_cols, colnames(df_raw))

  if (length(missing_cols) > 0) {
    cat(sprintf("   [ERROR] INVALID HEADERS! Missing: %s\n", paste(missing_cols, collapse=", ")))
    next
  }

  # --- Step 3: Clean ---
  df_plot <- df_raw %>%
    mutate(
      p_adjust        = as.numeric(p_adjust),
      gene_count      = as.numeric(gene_count),
      GeneRatioNumeric= parse_ratio_hybrid(gene_ratio),
      TermWrapped     = str_wrap(term_name, width = WRAP_WIDTH)
    ) %>%
    filter(!is.na(p_adjust), !is.na(GeneRatioNumeric))

  df_valid <- df_plot %>% filter(GeneRatioNumeric > 0)

  if (nrow(df_valid) == 0) {
    cat("   [WARNING] No valid data rows. Skipping.\n"); next
  }

  # --- Step 4: Mode & Sorting ---
  unique_ontologies <- if ("ontology" %in% colnames(df_valid)) unique(df_valid$ontology) else c()
  is_multi_ontology <- length(unique_ontologies) > 1

  if (is_multi_ontology) {
    mode_label <- "Facet"
    target_n   <- TOP_N_FACET
    df_ready <- df_valid %>%
      group_by(ontology) %>%
      arrange(p_adjust) %>%
      slice_head(n = target_n) %>%
      ungroup()
  } else {
    mode_label <- "Single"
    target_n   <- TOP_N_SINGLE
    df_ready <- df_valid %>%
      arrange(p_adjust) %>%
      slice_head(n = target_n)
  }

  if (SORT_BY == "GeneRatio") {
    df_ready <- df_ready %>%
      arrange(GeneRatioNumeric) %>%
      mutate(TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)))
  } else if (SORT_BY == "p_adjust") {
    df_ready <- df_ready %>%
      arrange(desc(p_adjust)) %>%
      mutate(TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)))
  } else {
    df_ready <- df_ready %>%
      arrange(gene_count) %>%
      mutate(TermWrapped = factor(TermWrapped, levels = unique(TermWrapped)))
  }

  final_count <- nrow(df_ready)
  cat(sprintf("   [INFO] Mode: %s | Plotting: %d terms\n", mode_label, final_count))

  # --- Step 5: Width Calculation ---
  if (is_multi_ontology) {
    final_width_inch  <- FACET_WIDTH_INCH
    final_height_inch <- FACET_HEIGHT_INCH
  } else {
    max_char_len     <- max(nchar(as.character(df_ready$TermWrapped)), na.rm = TRUE)
    grid_width_inch  <- SINGLE_HEIGHT_INCH / SINGLE_ASPECT_RATIO
    text_width_inch  <- max_char_len * 0.09 + 0.5
    legend_width_inch<- 2.0
    final_width_inch <- text_width_inch + grid_width_inch + legend_width_inch
    final_height_inch<- SINGLE_HEIGHT_INCH
  }

  # --- Step 6: Plotting ---
  fixed_breaks <- generate_fixed_breaks(GLOBAL_P_MIN, GLOBAL_P_MAX, LEGEND_TICKS)
  legend_h_cm  <- 3.5

  # Title Logic
  if (SHOW_PLOT_TITLE) {
    plot_title_text <- gsub("_", " ", f_base)
  } else {
    plot_title_text <- NULL
  }

  p <- ggplot(df_ready, aes(x = GeneRatioNumeric, y = TermWrapped)) +
    geom_point(aes(size = gene_count, color = p_adjust)) +

    scale_color_gradientn(
      colours = c(COLOR_HIGH, COLOR_MID, COLOR_LOW),
      trans   = "log10",
      breaks  = fixed_breaks,
      labels  = scientific_p_formatter,
      limits  = c(GLOBAL_P_MIN, GLOBAL_P_MAX),
      oob     = scales::squish,
      name    = "P.adjust\n"
    ) +

    scale_size_continuous(
      range  = c(3, 8),
      breaks = force_integer_breaks(n = COUNT_TICKS),
      name   = "Gene Count\n"
    ) +

    guides(
      color = guide_colorbar(order = 1, barheight = unit(legend_h_cm, "cm")),
      size  = guide_legend(order = 2)
    ) +

    labs(
      x     = "GeneRatio",
      y     = NULL,
      title = plot_title_text
    ) +

    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      aspect.ratio = if (!is_multi_ontology) SINGLE_ASPECT_RATIO else NULL,

      panel.grid.major = if (SHOW_GRID) element_line(colour = "grey92") else element_blank(),
      panel.grid.minor = element_blank(),

      axis.line  = element_line(colour = "black"),
      axis.text  = element_text(colour = "black", size = AXIS_TEXT_SIZE),
      axis.title = element_text(size = AXIS_TITLE_SIZE),

      axis.text.y = element_text(lineheight = 0.8, size = AXIS_TEXT_SIZE),

      legend.position = "right",
      legend.text     = element_text(size = LEGEND_TEXT_SIZE),
      legend.title    = element_text(size = LEGEND_TITLE_SIZE),

      plot.title = element_text(
        hjust = 0.5,
        face  = "bold",
        size  = TITLE_SIZE
      )
    )

  if (is_multi_ontology) {
    p <- p +
      facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
      theme(
        strip.background = element_rect(fill = "#E5E5E5", color = NA),
        strip.text       = element_text(face = "bold", size = STRIP_TEXT_SIZE)
      )
  }

  # --- Step 7: Save ---
  tryCatch({
    ggsave(
      file.path(OUT_DIR, paste0(f_base, ".bubble.pdf")),
      plot   = p,
      width  = final_width_inch,
      height = final_height_inch,
      device = cairo_pdf
    )

    ggsave(
      file.path(OUT_DIR, paste0(f_base, ".bubble.png")),
      plot   = p,
      width  = final_width_inch,
      height = final_height_inch,
      dpi    = 300,
      bg     = "white"
    )

    cat("   [SUCCESS] Saved.\n")
  }, error = function(e) {
    cat(sprintf("   [ERROR] Save failed: %s\n", e$message))
  })
}

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")

