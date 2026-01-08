#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: venn.R (v1.6 Layout Polish)
# Description: 
#   - [Visual] Threshold info moved to Bottom Caption (Single line).
#   - [Config] Toggles for Title and Caption (Legend).
#   - [Config] Custom Title support.
#   - [Visual] Pastel colors, Borderless, Clean style.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tools)
  
  if (!require("ggvenn", quietly = TRUE)) {
    stop("[ERROR] 'ggvenn' is missing. Install via: install.packages('ggvenn')")
  }
  library(ggvenn)
  
  has_complex_heatmap <- require("ComplexHeatmap", quietly = TRUE)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Input Control ---
TARGET_FILES <- c()

# --- Thresholds ---
FC_CUTOFF   <- 1.0    
P_CUTOFF    <- 0.05   

# --- Display Control ---
SHOW_PLOT_TITLE   <- TRUE       # Toggle Top Title
SHOW_CAPTION      <- TRUE       # Toggle Bottom Threshold Info
CUSTOM_PLOT_TITLE <- ""         # Override title string

# --- Aesthetics ---
FONT_FAMILY  <- "sans"
EDGE_SIZE    <- 0.3             
# Fresh Pastel Palette
SET_COLORS   <- c("#8DD3C7", "#FB8072", "#BEBADA", "#80B1D3") 
PLOT_MARGIN  <- 20              

#Paths
IN_DIR       <- "input"
OUT_DIR      <- file.path("output", "venn")
TEMPLATE_DIR <- "templates"
TEMPLATE_FILE<- "volcano_input_template.tsv" 

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  tmpl_path <- file.path(TEMPLATE_DIR, TEMPLATE_FILE)
  if (!file.exists(tmpl_path)) {
    example_data <- data.frame(
      gene_id  = c("Gene_A", "Gene_B", "Gene_C"),
      log2fc   = c(2.5, -3.1, 0.2),
      p_adjust = c(0.0001, 1.2e-10, 0.8),
      stringsAsFactors = FALSE
    )
    write_tsv(example_data, tmpl_path)
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic Venn/UpSet Plotter v1.6 (Layout Polish)\n")
cat("================================================================\n")

init_template_system()

# --- Step 1: Identify Files ---
files_to_process <- c()

if (length(TARGET_FILES) > 0) {
  for (f in TARGET_FILES) {
    f_path <- file.path(IN_DIR, f)
    if (file.exists(f_path)) files_to_process <- c(files_to_process, f_path)
    else cat(sprintf("   [WARN] Specified file not found: %s\n", f))
  }
} else {
  cat("Mode: Auto-scan 'input/' for DEG tables...\n")
  all_files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
  files_to_process <- all_files[!grepl("tpm|count|sample|meta|template", basename(all_files), ignore.case = TRUE)]
}

if (length(files_to_process) < 2) stop("[ERROR] Need at least 2 DEG files.")

cat(sprintf(">> Found %d files. Extracting significant genes...\n", length(files_to_process)))

# --- Step 2: Load & Filter Data ---
gene_sets <- list()

for (f_path in files_to_process) {
  fname <- basename(f_path)
  set_name <- file_path_sans_ext(fname)
  set_name <- gsub("_DEG_all", "", set_name)
  set_name <- gsub("_vs_D0", "", set_name) 
  
  df <- tryCatch(read_tsv(f_path, show_col_types = FALSE), error = function(e) NULL)
  
  if (is.null(df)) next
  if (!all(c("gene_id", "log2fc", "p_adjust") %in% colnames(df))) {
    cat(sprintf("   [SKIP] Invalid headers in %s\n", fname))
    next
  }
  
  sig_genes <- df %>% 
    filter(p_adjust < P_CUTOFF & abs(log2fc) > FC_CUTOFF) %>%
    pull(gene_id)
  
  cat(sprintf("   + %s: %d sig genes\n", set_name, length(sig_genes)))
  
  if (length(sig_genes) > 0) {
    gene_sets[[set_name]] <- sig_genes
  }
}

n_sets <- length(gene_sets)
if (n_sets < 2) stop("[ERROR] Not enough valid gene sets (<2).")

# --- Step 3: Plotting Logic ---

# [FIX] Single Line Caption at Bottom
caption_text <- NULL
if (SHOW_CAPTION) {
  caption_text <- sprintf("Thresholds: P.adjust < %s & |log2FC| > %s", P_CUTOFF, FC_CUTOFF)
}

# Title Logic
final_title <- NULL
if (SHOW_PLOT_TITLE) {
  if (CUSTOM_PLOT_TITLE != "") {
    final_title <- CUSTOM_PLOT_TITLE
  } else {
    final_title <- "Venn Analysis"
  }
}

if (n_sets <= 4) {
  # >>> VENN DIAGRAM MODE <<<
  cat(">> Mode: Pastel Venn Diagram (<= 4 sets)\n")
  
  p <- ggvenn(
    gene_sets, 
    fill_color = SET_COLORS[1:n_sets],
    fill_alpha = 0.5,
    stroke_color = "white",
    stroke_size = 0.5,
    set_name_size = 6,
    text_size = 4
  ) +
  # [FIX] Title Top, Caption Bottom
  labs(title = final_title, caption = caption_text) +
  theme(
    text = element_text(family = FONT_FAMILY),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    # [FIX] Center align caption at bottom
    plot.caption = element_text(hjust = 0.5, size = 12, color = "grey30", margin = margin(t = 10)),
    plot.margin = margin(t = PLOT_MARGIN, r = PLOT_MARGIN, b = PLOT_MARGIN, l = PLOT_MARGIN, unit = "pt"),
    legend.position = "bottom" # Just in case
  )
  
  out_name <- paste0("Venn_", paste(names(gene_sets), collapse="_"))
  if (nchar(out_name) > 50) out_name <- "Venn_Combined_Sets"
  
  tryCatch({
    ggsave(file.path(OUT_DIR, paste0(out_name, ".pdf")), p, width = 8, height = 8, device = cairo_pdf)
    ggsave(file.path(OUT_DIR, paste0(out_name, ".png")), p, width = 8, height = 8, dpi = 300, bg = "white")
    cat("   [SUCCESS] Venn Diagram Saved.\n")
  }, error = function(e) cat(sprintf("   [ERROR] Save failed: %s\n", e$message)))
  
} else {
  # >>> UPSET PLOT MODE <<<
  cat(">> Mode: UpSet Plot (> 4 sets)\n")
  if (!has_complex_heatmap) stop("[ERROR] >4 sets requires 'ComplexHeatmap'.")
  
  m <- make_comb_mat(gene_sets)
  
  # For UpSet, we append caption to title because it handles layout differently
  upset_title <- final_title
  if (!is.null(caption_text)) {
    upset_title <- paste0(upset_title, "\n", caption_text)
  }
  
  out_name <- "UpSet_Combined_Sets"
  
  tryCatch({
    pdf(file.path(OUT_DIR, paste0(out_name, ".pdf")), width = 10, height = 6)
    draw(UpSet(m, column_title = upset_title))
    dev.off()
    
    png(file.path(OUT_DIR, paste0(out_name, ".png")), width = 3000, height = 1800, res = 300)
    draw(UpSet(m, column_title = upset_title))
    dev.off()
    cat("   [SUCCESS] UpSet Plot Saved.\n")
  }, error = function(e) { if (dev.cur() > 1) dev.off(); cat(sprintf("   [ERROR] Save failed: %s\n", e$message)) })
}

# --- Step 4: Export Intersections ---
cat(">> Exporting Intersection Lists...\n")
core_genes <- Reduce(intersect, gene_sets)
if (length(core_genes) > 0) {
  writeLines(core_genes, file.path(OUT_DIR, "Genes_Intersection_ALL.txt"))
}
for (set_name in names(gene_sets)) {
  others <- gene_sets[names(gene_sets) != set_name]
  unique_genes <- setdiff(gene_sets[[set_name]], unlist(others))
  if (length(unique_genes) > 0) {
    writeLines(unique_genes, file.path(OUT_DIR, paste0("Genes_Unique_", set_name, ".txt")))
  }
}

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
