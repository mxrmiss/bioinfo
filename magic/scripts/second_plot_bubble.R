#!/usr/bin/env Rscript

# =============================================================================
# second_plot_bubble.R
# Purpose:
#   Read GO BP and KEGG enrichment tables (4 species; 7 terms each),
#   build a long-format table, and draw two-panel bubble plots:
#     Panel A: GO BP
#     Panel B: KEGG
#   Mapping:
#     dot size  = gene_ratio
#     dot color = -log10(FDR) where FDR = p_adjust
#
# Requirements (implemented):
#   - Run from magic/ with relative paths
#   - Inputs in magic/input
#   - Outputs in magic/output/<OUT_SUBDIR> (no timestamp)
#   - Script location: magic/scripts/second_plot_bubble.R
#   - Species order: Sc Sr Tg Rp
#   - Y axis shows term English names only (no IDs)
#   - No underscores in any visible plot text
#   - viridis family palette
#   - Grid toggle
#   - User-settable font family
#   - User-settable axis font sizes (x and y separately)
#   - Streaming console output + log file
# =============================================================================

# =========================
# User settings (edit only here)
# =========================

INPUT_DIR  <- "input2"
OUTPUT_DIR <- "output"

# Output subfolder name (no timestamp)
OUT_SUBDIR <- "bubble"

# Species order (fixed)
SPECIES_ORDER <- c("Sc", "Sr", "Tg", "Rp")

# Input file templates
BP_FILE_TMPL   <- "{sp}.GO_BP_by_term_test.tsv"
KEGG_FILE_TMPL <- "{sp}.KEGG_by_term_test.tsv"

# Keep lists (exactly 7 + 7)
BP_KEEP <- c(
  "GO:0099565",
  "GO:0060078",
  "GO:0032970",
  "GO:0008344",
  "GO:0006936",
  "GO:0045214",
  "GO:0042692"
)

KEGG_KEEP <- c(
  "ko04080",
  "ko04260",
  "ko04261",
  "ko04814",
  "ko04713",
  "ko04371",
  "ko05410"
)

# Font family ("" => system default sans)
FONT_FAMILY <- "Arial"

# Grid toggle
GRID_ON <- FALSE

# Axis text sizes (publication-ready; adjust as needed)
AXIS_TEXT_X_SIZE <- 17
AXIS_TEXT_Y_SIZE <- 18

# Strip/title and legend text sizes
STRIP_TITLE_SIZE <- 17
PANEL_TITLE_SIZE <- 18
LEGEND_TITLE_SIZE <- 18
LEGEND_TEXT_SIZE  <- 17

# Output figures
OUT_PNG <- TRUE
OUT_PDF <- TRUE
PNG_DPI <- 600

# Figure size (inches)
FIG_WIDTH  <- 9.5
FIG_HEIGHT <- 8

# Bubble size range
POINT_SIZE_MIN <- 2.2
POINT_SIZE_MAX <- 10

# viridis option: viridis, magma, inferno, plasma, cividis, turbo
VIRIDIS_OPTION <- "viridis"

# Cap for -log10(FDR) to avoid colorbar dominated by extreme values (NA to disable)
NEGLOG10_FDR_CAP <- 20

# Overwrite existing output folder if it exists
OVERWRITE_OUTPUT_DIR <- TRUE

# =========================
# Libraries
# =========================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

# =========================
# Logging helpers (streaming + file)
# =========================

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

rm_dir <- function(p) {
  if (dir.exists(p)) unlink(p, recursive = TRUE, force = TRUE)
}

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_init <- function(log_path) {
  file(log_path, open = "wt")
}

log_msg <- function(con, ...) {
  msg <- paste0("[", now_str(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
  writeLines(msg, con = con)
  flush(con)
}

stop_with_log <- function(con, ...) {
  log_msg(con, "ERROR: ", paste0(..., collapse = ""))
  close(con)
  stop(paste0(..., collapse = ""), call. = FALSE)
}

# =========================
# IO utilities
# =========================

build_path <- function(dir, tmpl, sp) {
  fn <- str_replace_all(tmpl, fixed("{sp}"), sp)
  file.path(dir, fn)
}

# =========================
# Read and normalize
# =========================

read_bp_one <- function(sp, fp, keep_ids, con) {
  log_msg(con, "Read GO BP: ", fp)
  if (!file.exists(fp)) stop_with_log(con, "Missing file: ", fp)

  df <- suppressWarnings(read_tsv(fp, show_col_types = FALSE, progress = FALSE))
  need_cols <- c("term_id", "term_name", "gene_ratio", "p_adjust")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0) stop_with_log(con, "GO BP file missing columns: ", paste(miss, collapse = ", "))

  df2 <- df %>%
    filter(.data$term_id %in% keep_ids) %>%
    transmute(
      panel = "GO BP",
      species = sp,
      id = .data$term_id,
      term_name = .data$term_name,
      gene_ratio = as.numeric(.data$gene_ratio),
      p_adjust = as.numeric(.data$p_adjust),
      gene_count = if ("gene_count" %in% colnames(df)) as.numeric(.data$gene_count) else NA_real_
    )

  log_msg(con, "  Rows kept: ", nrow(df2))
  df2
}

read_kegg_one <- function(sp, fp, keep_ids, con) {
  log_msg(con, "Read KEGG: ", fp)
  if (!file.exists(fp)) stop_with_log(con, "Missing file: ", fp)

  df <- suppressWarnings(read_tsv(fp, show_col_types = FALSE, progress = FALSE))
  need_cols <- c("pathway_id", "term_name", "gene_ratio", "p_adjust")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0) stop_with_log(con, "KEGG file missing columns: ", paste(miss, collapse = ", "))

  df2 <- df %>%
    filter(.data$pathway_id %in% keep_ids) %>%
    transmute(
      panel = "KEGG",
      species = sp,
      id = .data$pathway_id,
      term_name = .data$term_name,
      gene_ratio = as.numeric(.data$gene_ratio),
      p_adjust = as.numeric(.data$p_adjust),
      gene_count = if ("gene_count" %in% colnames(df)) as.numeric(.data$gene_count) else NA_real_
    )

  log_msg(con, "  Rows kept: ", nrow(df2))
  df2
}

qc_and_finalize <- function(df, keep_bp, keep_kegg, con) {
  log_msg(con, "QC: checking completeness and numeric columns")

  df$species <- factor(df$species, levels = SPECIES_ORDER)
  df$panel <- factor(df$panel, levels = c("GO BP", "KEGG"))

  # Replace p_adjust == 0 to avoid -Inf in -log10
  n0 <- sum(!is.na(df$p_adjust) & df$p_adjust == 0)
  if (n0 > 0) {
    log_msg(con, "  Found p_adjust == 0 rows: ", n0, " -> set to 1e-300")
    df$p_adjust[df$p_adjust == 0] <- 1e-300
  }

  if (any(is.na(df$gene_ratio))) stop_with_log(con, "gene_ratio has NA after numeric conversion")
  if (any(is.na(df$p_adjust))) stop_with_log(con, "p_adjust has NA after numeric conversion")
  if (any(df$p_adjust <= 0 | df$p_adjust > 1)) stop_with_log(con, "p_adjust out of (0,1] detected")

  df <- df %>%
    mutate(neglog10_fdr = -log10(.data$p_adjust))

  if (!is.na(NEGLOG10_FDR_CAP)) {
    df <- df %>% mutate(neglog10_fdr = pmin(.data$neglog10_fdr, NEGLOG10_FDR_CAP))
  }

  # Completeness checks: each term should exist and have 4 species rows
  check_block <- function(panel_name, keep_ids, df_block) {
    missing_terms <- setdiff(keep_ids, unique(df_block$id))
    if (length(missing_terms) > 0) {
      stop_with_log(con, panel_name, " missing terms: ", paste(missing_terms, collapse = ", "))
    }
    bad <- df_block %>% count(id) %>% filter(n != length(SPECIES_ORDER))
    if (nrow(bad) > 0) {
      stop_with_log(con, panel_name, " has terms with incomplete species rows: ",
                    paste0(bad$id, "(", bad$n, ")", collapse = ", "))
    }
  }

  check_block("GO BP", keep_bp, df %>% filter(panel == "GO BP"))
  check_block("KEGG", keep_kegg, df %>% filter(panel == "KEGG"))

  # Remove underscores from visible text fields
  df <- df %>%
    mutate(
      term_name = str_replace_all(.data$term_name, "_", " "),
      panel = str_replace_all(as.character(.data$panel), "_", " ")
    ) %>%
    mutate(panel = factor(panel, levels = c("GO BP", "KEGG")))

  df
}

# =========================
# Plot helpers (publication-ready theme)
# =========================

make_theme_pub <- function(grid_on, font_family) {
  base_family <- ifelse(is.null(font_family) || font_family == "", "sans", font_family)

  th <- theme_classic(base_family = base_family) +
    theme(
      plot.title = element_text(size = PANEL_TITLE_SIZE, face = "bold", hjust = 0),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = AXIS_TEXT_X_SIZE, color = "black"),
      axis.text.y  = element_text(size = AXIS_TEXT_Y_SIZE, color = "black"),
      axis.ticks   = element_line(color = "black", linewidth = 0.4),
      legend.title = element_text(size = LEGEND_TITLE_SIZE),
      legend.text  = element_text(size = LEGEND_TEXT_SIZE),
      legend.key.height = unit(0.9, "lines"),
      legend.key.width  = unit(1.4, "lines"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = STRIP_TITLE_SIZE, face = "bold"),
      panel.spacing = unit(0.8, "lines")
    )

  if (grid_on) {
    th <- th + theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2)
    )
  } else {
    th <- th + theme(panel.grid = element_blank())
  }
  th
}

make_bubble_plot <- function(df, panel_name, keep_ids, con) {
  d <- df %>% filter(panel == panel_name)

  # Fix y-order by keep list order
  d <- d %>% mutate(id = factor(id, levels = keep_ids))

  # Map id -> term_name (should be unique per id)
  y_map <- d %>% select(id, term_name) %>% distinct() %>% arrange(id)
  if (nrow(y_map) != length(keep_ids)) {
    stop_with_log(con, "Panel ", panel_name, " term_name mapping not 1-to-1 with keep list")
  }

  # y-axis uses term_name only (no IDs)
  y_levels <- y_map$term_name
  d <- d %>% mutate(term_name = factor(term_name, levels = y_levels))

  log_msg(con, "Plot panel: ", panel_name, " terms=", length(y_levels), " rows=", nrow(d))

  ggplot(d, aes(x = species, y = term_name)) +
    geom_point(aes(size = gene_ratio, color = neglog10_fdr), alpha = 0.95) +
    scale_color_viridis_c(
      option = VIRIDIS_OPTION,
      name = "-log10 FDR"
    ) +
    scale_size_continuous(
      name = "Gene ratio",
      range = c(POINT_SIZE_MIN, POINT_SIZE_MAX)
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    ) +
    labs(title = panel_name) +
    make_theme_pub(GRID_ON, FONT_FAMILY)
}

# =========================
# Main
# =========================

ensure_dir(OUTPUT_DIR)
OUTDIR <- file.path(OUTPUT_DIR, OUT_SUBDIR)

if (OVERWRITE_OUTPUT_DIR) {
  if (dir.exists(OUTDIR)) rm_dir(OUTDIR)
}
ensure_dir(OUTDIR)

LOG_PATH <- file.path(OUTDIR, "second_plot_bubble.log")
log_con <- log_init(LOG_PATH)

log_msg(log_con, "Start second_plot_bubble.R")
log_msg(log_con, "Working dir: ", getwd())
log_msg(log_con, "Input dir: ", INPUT_DIR)
log_msg(log_con, "Output dir: ", OUTDIR)
log_msg(log_con, "Species order: ", paste(SPECIES_ORDER, collapse = " "))

bp_list <- list()
kegg_list <- list()

for (sp in SPECIES_ORDER) {
  bp_fp <- build_path(INPUT_DIR, BP_FILE_TMPL, sp)
  kg_fp <- build_path(INPUT_DIR, KEGG_FILE_TMPL, sp)
  bp_list[[sp]]   <- read_bp_one(sp, bp_fp, BP_KEEP, log_con)
  kegg_list[[sp]] <- read_kegg_one(sp, kg_fp, KEGG_KEEP, log_con)
}

df_all <- bind_rows(bind_rows(bp_list), bind_rows(kegg_list))
log_msg(log_con, "Merge done. Total rows: ", nrow(df_all))

df_all <- qc_and_finalize(df_all, BP_KEEP, KEGG_KEEP, log_con)

# Write reproducible plot input table
plot_input_path <- file.path(OUTDIR, "plot_input.tsv")
write_tsv(df_all, plot_input_path)
log_msg(log_con, "Write plot input: ", plot_input_path)

# Plot panels
p_bp   <- make_bubble_plot(df_all, "GO BP", BP_KEEP, log_con)
p_kegg <- make_bubble_plot(df_all, "KEGG", KEGG_KEEP, log_con)

# Combine (vertical stack)
p_all <- p_bp / p_kegg + plot_layout(heights = c(1, 1))

# Save outputs
if (OUT_PNG) {
  out_png <- file.path(OUTDIR, "bubble_GO_BP_KEGG.png")
  ggsave(out_png, p_all, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = PNG_DPI)
  log_msg(log_con, "Write PNG: ", out_png)
}
if (OUT_PDF) {
  out_pdf <- file.path(OUTDIR, "bubble_GO_BP_KEGG.pdf")
  ggsave(out_pdf, p_all, width = FIG_WIDTH, height = FIG_HEIGHT, device = cairo_pdf)
  log_msg(log_con, "Write PDF: ", out_pdf)
}

log_msg(log_con, "Done.")
close(log_con)

