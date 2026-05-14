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
# =============================================================================

# =========================
# User settings
# =========================

INPUT_DIR  <- "input2"
OUTPUT_DIR <- "output"

# 输出子文件夹名称
OUT_SUBDIR <- "bubble"

# 物种顺序
SPECIES_ORDER <- c("Sc", "Sr", "Tg", "Rp")

# 输入文件模板
BP_FILE_TMPL   <- "{sp}.GO_BP_by_term_test.tsv"
KEGG_FILE_TMPL <- "{sp}.KEGG_by_term_test.tsv"

# GO BP 保留条目
BP_KEEP <- c(
  "GO:0099565",
  "GO:0060078",
  "GO:0032970",
  "GO:0008344",
  "GO:0006936",
  "GO:0045214",
  "GO:0042692"
)

# KEGG 保留条目
KEGG_KEEP <- c(
  "ko04080",
  "ko04260",
  "ko04261",
  "ko04814",
  "ko04713",
  "ko04371",
  "ko05410"
)

# 字体设置，空字符串表示系统默认 sans
FONT_FAMILY <- "Arial"

# 是否显示网格线
GRID_ON <- FALSE

# 坐标轴文字字号
AXIS_TEXT_X_SIZE <- 17
AXIS_TEXT_Y_SIZE <- 18

# 标题、strip 和图例字号
STRIP_TITLE_SIZE  <- 17
PANEL_TITLE_SIZE  <- 18
LEGEND_TITLE_SIZE <- AXIS_TEXT_Y_SIZE * 0.8
LEGEND_TEXT_SIZE  <- 15

# 输出设置
OUT_PNG <- TRUE
OUT_PDF <- TRUE
PNG_DPI <- 600

# 图像尺寸
FIG_WIDTH  <- 9.5
FIG_HEIGHT <- 11

# 气泡大小范围
POINT_SIZE_MIN <- 2.2
POINT_SIZE_MAX <- 10

# 御用连续配色：用于 -log10 FDR
CUSTOM_FDR_COLORS <- c(
  "#95C8F2",
  "#A99BEF",
  "#F3A6C9",
  "#F79D93",
  "#F5B07E"
)

# -log10(FDR) 上限，避免极端值支配颜色条；NA 表示不截断
NEGLOG10_FDR_CAP <- 20

# 是否覆盖已有输出文件夹
OVERWRITE_OUTPUT_DIR <- TRUE

# =========================
# Libraries
# =========================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

# =========================
# Logging helpers
# =========================

ensure_dir <- function(p) {
  if (!dir.exists(p)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

rm_dir <- function(p) {
  if (dir.exists(p)) {
    unlink(p, recursive = TRUE, force = TRUE)
  }
}

now_str <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

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

  if (!file.exists(fp)) {
    stop_with_log(con, "Missing file: ", fp)
  }

  df <- suppressWarnings(read_tsv(fp, show_col_types = FALSE, progress = FALSE))

  need_cols <- c("term_id", "term_name", "gene_ratio", "p_adjust")
  miss <- setdiff(need_cols, colnames(df))

  if (length(miss) > 0) {
    stop_with_log(con, "GO BP file missing columns: ", paste(miss, collapse = ", "))
  }

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

  if (!file.exists(fp)) {
    stop_with_log(con, "Missing file: ", fp)
  }

  df <- suppressWarnings(read_tsv(fp, show_col_types = FALSE, progress = FALSE))

  need_cols <- c("pathway_id", "term_name", "gene_ratio", "p_adjust")
  miss <- setdiff(need_cols, colnames(df))

  if (length(miss) > 0) {
    stop_with_log(con, "KEGG file missing columns: ", paste(miss, collapse = ", "))
  }

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

  # p_adjust 为 0 时无法计算 -log10，因此替换成极小正数
  n0 <- sum(!is.na(df$p_adjust) & df$p_adjust == 0)

  if (n0 > 0) {
    log_msg(con, "  Found p_adjust == 0 rows: ", n0, " -> set to 1e-300")
    df$p_adjust[df$p_adjust == 0] <- 1e-300
  }

  if (any(is.na(df$gene_ratio))) {
    stop_with_log(con, "gene_ratio has NA after numeric conversion")
  }

  if (any(is.na(df$p_adjust))) {
    stop_with_log(con, "p_adjust has NA after numeric conversion")
  }

  if (any(df$p_adjust <= 0 | df$p_adjust > 1)) {
    stop_with_log(con, "p_adjust out of (0,1] detected")
  }

  df <- df %>%
    mutate(neglog10_fdr = -log10(.data$p_adjust))

  if (!is.na(NEGLOG10_FDR_CAP)) {
    df <- df %>%
      mutate(neglog10_fdr = pmin(.data$neglog10_fdr, NEGLOG10_FDR_CAP))
  }

  check_block <- function(panel_name, keep_ids, df_block) {
    missing_terms <- setdiff(keep_ids, unique(df_block$id))

    if (length(missing_terms) > 0) {
      stop_with_log(con, panel_name, " missing terms: ", paste(missing_terms, collapse = ", "))
    }

    bad <- df_block %>%
      count(id) %>%
      filter(n != length(SPECIES_ORDER))

    if (nrow(bad) > 0) {
      stop_with_log(
        con,
        panel_name,
        " has terms with incomplete species rows: ",
        paste0(bad$id, "(", bad$n, ")", collapse = ", ")
      )
    }
  }

  check_block("GO BP", keep_bp, df %>% filter(panel == "GO BP"))
  check_block("KEGG", keep_kegg, df %>% filter(panel == "KEGG"))

  # 可见文本去掉下划线
  df <- df %>%
    mutate(
      term_name = str_replace_all(.data$term_name, "_", " "),
      panel = str_replace_all(as.character(.data$panel), "_", " ")
    ) %>%
    mutate(panel = factor(panel, levels = c("GO BP", "KEGG")))

  df
}

# =========================
# Scale helpers
# =========================

make_gene_ratio_breaks <- function(x, n = 4) {
  rng <- range(x, na.rm = TRUE)

  if (!all(is.finite(rng))) {
    return(NULL)
  }

  if (rng[1] == rng[2]) {
    return(rng[1])
  }

  raw_breaks <- pretty(rng, n = n)
  raw_breaks <- raw_breaks[raw_breaks >= rng[1] & raw_breaks <= rng[2]]

  if (length(raw_breaks) == 0) {
    raw_breaks <- seq(rng[1], rng[2], length.out = n)
  }

  raw_breaks
}

make_color_breaks <- function(x, n = 5) {
  rng <- range(x, na.rm = TRUE)

  if (!all(is.finite(rng))) {
    return(NULL)
  }

  if (rng[1] == rng[2]) {
    return(rng[1])
  }

  pretty(rng, n = n)
}

get_global_scales <- function(df) {
  size_limits <- range(df$gene_ratio, na.rm = TRUE)
  color_limits <- range(df$neglog10_fdr, na.rm = TRUE)

  size_breaks <- make_gene_ratio_breaks(df$gene_ratio, n = 4)
  color_breaks <- make_color_breaks(df$neglog10_fdr, n = 5)

  list(
    size_limits = size_limits,
    color_limits = color_limits,
    size_breaks = size_breaks,
    color_breaks = color_breaks
  )
}

# =========================
# Plot helpers
# =========================

make_theme_pub <- function(grid_on, font_family) {
  base_family <- ifelse(is.null(font_family) || font_family == "", "sans", font_family)

  th <- theme_classic(base_family = base_family) +
    theme(
      plot.title = element_text(
        size = PANEL_TITLE_SIZE,
        face = "bold",
        hjust = 0
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = AXIS_TEXT_X_SIZE,
        color = "black"
      ),
      axis.text.y = element_text(
        size = AXIS_TEXT_Y_SIZE,
        color = "black"
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.4
      ),
      axis.line = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.35
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      legend.title = element_text(
        size = LEGEND_TITLE_SIZE,
        color = "black"
      ),
      legend.text = element_text(
        size = LEGEND_TEXT_SIZE,
        color = "black"
      ),
      legend.key.height = unit(0.9, "lines"),
      legend.key.width = unit(1.4, "lines"),
      strip.background = element_rect(
        fill = "grey95",
        color = NA
      ),
      strip.text = element_text(
        size = STRIP_TITLE_SIZE,
        face = "bold"
      ),
      panel.spacing = unit(0.8, "lines")
    )

  if (grid_on) {
    th <- th +
      theme(
        panel.grid.major = element_line(
          color = "grey85",
          linewidth = 0.3
        ),
        panel.grid.minor = element_line(
          color = "grey92",
          linewidth = 0.2
        )
      )
  } else {
    th <- th +
      theme(panel.grid = element_blank())
  }

  th
}

make_bubble_plot <- function(df, panel_name, keep_ids, global_scales, con) {
  d <- df %>%
    filter(panel == panel_name)

  # y 轴顺序按 keep list 固定
  d <- d %>%
    mutate(id = factor(id, levels = keep_ids))

  y_map <- d %>%
    select(id, term_name) %>%
    distinct() %>%
    arrange(id)

  if (nrow(y_map) != length(keep_ids)) {
    stop_with_log(con, "Panel ", panel_name, " term_name mapping not 1-to-1 with keep list")
  }

  y_levels <- y_map$term_name

  d <- d %>%
    mutate(term_name = factor(term_name, levels = y_levels))

  log_msg(con, "Plot panel: ", panel_name, " terms=", length(y_levels), " rows=", nrow(d))

  ggplot(d, aes(x = species, y = term_name)) +
    geom_point(
      aes(size = gene_ratio, color = neglog10_fdr),
      alpha = 0.95
    ) +
    scale_color_gradientn(
      colours = CUSTOM_FDR_COLORS,
      limits = global_scales$color_limits,
      breaks = global_scales$color_breaks,
      labels = label_number(accuracy = 1),
      name = "-log10 FDR",
      oob = scales::squish
    ) +
    scale_size_continuous(
      name = "Gene ratio",
      range = c(POINT_SIZE_MIN, POINT_SIZE_MAX),
      limits = global_scales$size_limits,
      breaks = global_scales$size_breaks,
      labels = label_number(accuracy = 0.01)
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title.position = "top"
      ),
      color = guide_colorbar(
        order = 2,
        title.position = "top",
        barwidth = unit(1.26, "lines"),
        barheight = unit(5.0, "lines")
      )
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
  if (dir.exists(OUTDIR)) {
    rm_dir(OUTDIR)
  }
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

  bp_list[[sp]] <- read_bp_one(
    sp = sp,
    fp = bp_fp,
    keep_ids = BP_KEEP,
    con = log_con
  )

  kegg_list[[sp]] <- read_kegg_one(
    sp = sp,
    fp = kg_fp,
    keep_ids = KEGG_KEEP,
    con = log_con
  )
}

df_all <- bind_rows(
  bind_rows(bp_list),
  bind_rows(kegg_list)
)

log_msg(log_con, "Merge done. Total rows: ", nrow(df_all))

df_all <- qc_and_finalize(
  df = df_all,
  keep_bp = BP_KEEP,
  keep_kegg = KEGG_KEEP,
  con = log_con
)

plot_input_path <- file.path(OUTDIR, "plot_input.tsv")
write_tsv(df_all, plot_input_path)
log_msg(log_con, "Write plot input: ", plot_input_path)

global_scales <- get_global_scales(df_all)

log_msg(
  log_con,
  "Global gene_ratio limits: ",
  paste(round(global_scales$size_limits, 4), collapse = " - ")
)

log_msg(
  log_con,
  "Global -log10 FDR limits: ",
  paste(round(global_scales$color_limits, 4), collapse = " - ")
)

p_bp <- make_bubble_plot(
  df = df_all,
  panel_name = "GO BP",
  keep_ids = BP_KEEP,
  global_scales = global_scales,
  con = log_con
)

p_kegg <- make_bubble_plot(
  df = df_all,
  panel_name = "KEGG",
  keep_ids = KEGG_KEEP,
  global_scales = global_scales,
  con = log_con
)

# 上下组合，并收集为一份共用图例
p_all <- p_bp / p_kegg +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

if (OUT_PNG) {
  out_png <- file.path(OUTDIR, "bubble_GO_BP_KEGG.png")

  ggsave(
    filename = out_png,
    plot = p_all,
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    dpi = PNG_DPI,
    bg = "white"
  )

  log_msg(log_con, "Write PNG: ", out_png)
}

if (OUT_PDF) {
  out_pdf <- file.path(OUTDIR, "bubble_GO_BP_KEGG.pdf")

  ggsave(
    filename = out_pdf,
    plot = p_all,
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    device = cairo_pdf,
    bg = "white"
  )

  log_msg(log_con, "Write PDF: ", out_pdf)
}

log_msg(log_con, "Done.")
close(log_con)

