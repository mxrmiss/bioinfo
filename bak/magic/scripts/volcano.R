#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: volcano.R (v1.5 Final Polish)
# Description: 
#   - [Legend] Title removed "Regulation", split thresholds into 2 lines.
#   - [Legend] Counts format simplified: "Up (100)" (removed 'n=').
#   - [Config] Target File support & Strict Mode enabled.
#   - [Visual] Gradient Coloring, Hybrid Labeling, Symmetrical Axis.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tools)
  library(scales)
  if (!require("ggrepel", quietly = TRUE)) {
    stop("[ERROR] Package 'ggrepel' is missing. Please install it via install.packages('ggrepel').")
  }
  library(ggrepel)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Input Control ---
# Leave empty "" to scan all files, or specify "DEG_all.tsv"
TARGET_FILE <- ""  

# --- Thresholds ---
FC_CUTOFF   <- 1.0    # log2FC threshold
P_CUTOFF    <- 0.05   # P-value threshold

# --- Labeling (Hybrid) ---
SHOW_LABELS  <- TRUE  
TOP_N_LABELS <- 10    

# --- Aesthetics ---
COLOR_UP     <- "#CC0000"  
COLOR_DOWN   <- "#0000CC"  
COLOR_NS     <- "grey70"   
POINT_SIZE   <- 2.0        

# --- Canvas ---
FIXED_WIDTH  <- 8.0
FIXED_HEIGHT <- 8.0   
SHOW_TITLE   <- TRUE  

#Paths
IN_DIR       <- "input"
OUT_DIR      <- file.path("output", "volcano") 
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
      gene_id  = c("Gene_A", "Gene_B", "Gene_C", "Gene_D"),
      log2fc   = c(2.5, -3.1, 0.2, 1.5),
      p_adjust = c(0.0001, 1.2e-10, 0.8, 0.03),
      stringsAsFactors = FALSE
    )
    write_tsv(example_data, tmpl_path)
    cat(sprintf("   [SYSTEM] Generated standard template: %s\n", tmpl_path))
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic Volcano Plotter v1.5 (Clean Legend)\n")
cat("================================================================\n")

init_template_system()

# --- Step 0: File Selection ---
files_to_process <- c()

if (TARGET_FILE != "") {
  target_path <- file.path(IN_DIR, TARGET_FILE)
  if (!file.exists(target_path)) {
    stop(sprintf("[FATAL ERROR] Specified file not found: %s", target_path))
  }
  files_to_process <- c(target_path)
  cat(sprintf("Mode: Single File Target -> %s\n", TARGET_FILE))
} else {
  cat("Mode: Auto-scan 'input/' directory...\n")
  files_to_process <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
  files_to_process <- files_to_process[!grepl("template", basename(files_to_process), ignore.case = TRUE)]
  
  if (length(files_to_process) == 0) {
    stop("[ERROR] No .tsv files found in magic/input/")
  }
}

cat(sprintf("Thresholds: |log2FC| > %.1f & P < %.2f\n", FC_CUTOFF, P_CUTOFF))
cat("----------------------------------------------------------------\n")

for (f_path in files_to_process) {
  
  fname <- basename(f_path)
  f_base <- file_path_sans_ext(fname)
  cat(sprintf(">> Processing: %s\n", fname))
  
  # --- Step 1: Read ---
  df_raw <- tryCatch({
    read_tsv(f_path, show_col_types = FALSE, progress = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(df_raw) || nrow(df_raw) == 0) { 
    cat("   [ERROR] File is empty or unreadable. Skipping.\n"); next 
  }
  
  # --- Step 2: Header Check ---
  required_cols <- c("gene_id", "log2fc", "p_adjust")
  missing_cols <- setdiff(required_cols, colnames(df_raw))
  
  if (length(missing_cols) > 0) {
    cat(sprintf("   [FATAL ERROR] Invalid Headers! Missing: %s\n", paste(missing_cols, collapse=", ")))
    cat("   ------------------------------------------------------------\n")
    cat("   [SOLUTION] Fix headers to match template: magic/templates/\n")
    next 
  }
  
  # --- Step 3: Processing ---
  df_plot <- df_raw %>%
    filter(!is.na(log2fc), !is.na(p_adjust)) %>%
    mutate(
      logP = -log10(p_adjust),
      Group = case_when(
        p_adjust < P_CUTOFF & log2fc >= FC_CUTOFF  ~ "Up",
        p_adjust < P_CUTOFF & log2fc <= -FC_CUTOFF ~ "Down",
        TRUE ~ "Normal"
      ),
      RankScore = logP * abs(log2fc)
    )
  
  # --- Step 4: Labels & Title Construction ---
  n_up   <- sum(df_plot$Group == "Up")
  n_down <- sum(df_plot$Group == "Down")
  n_ns   <- sum(df_plot$Group == "Normal")
  
  # [MODIFIED] Simplified labels (removed "n=")
  label_up   <- sprintf("Up (%d)", n_up)
  label_down <- sprintf("Down (%d)", n_down)
  label_ns   <- sprintf("Normal (%d)", n_ns)
  
  cat(sprintf("   [INFO] %s | %s\n", label_up, label_down))
  
  df_plot <- df_plot %>%
    mutate(RegLabel = case_when(
      Group == "Up"   ~ label_up,
      Group == "Down" ~ label_down,
      TRUE            ~ label_ns
    )) %>%
    mutate(RegLabel = factor(RegLabel, levels = c(label_down, label_ns, label_up)))

  # [MODIFIED] Legend Title: No "Regulation", 2 lines, Full Names
  legend_title_text <- sprintf("P.adjust < %s\n|log2(FoldChange)| > %s", P_CUTOFF, FC_CUTOFF)

  # --- Step 5: Label Selection ---
  label_data <- NULL
  if (SHOW_LABELS) {
    top_up <- df_plot %>% filter(Group == "Up") %>% arrange(desc(RankScore)) %>% head(TOP_N_LABELS)
    top_down <- df_plot %>% filter(Group == "Down") %>% arrange(desc(RankScore)) %>% head(TOP_N_LABELS)
    label_data <- bind_rows(top_up, top_down)
  }

  # --- Step 6: Symmetrical Axis ---
  max_fc <- max(abs(df_plot$log2fc), na.rm = TRUE)
  x_lim <- max(max_fc * 1.1, 1)

  # --- Step 7: Plotting ---
  plot_title <- if(SHOW_TITLE) gsub("_", " ", f_base) else NULL
  
  my_colors <- c()
  my_colors[label_up]   <- COLOR_UP
  my_colors[label_down] <- COLOR_DOWN
  my_colors[label_ns]   <- COLOR_NS

  p <- ggplot(df_plot, aes(x = log2fc, y = logP)) +
    
    # Points
    geom_point(aes(color = RegLabel, alpha = logP), size = POINT_SIZE) +
    
    scale_alpha_continuous(range = c(0.2, 1.0), guide = "none") +
    scale_color_manual(values = my_colors, name = legend_title_text) +
    
    geom_hline(yintercept = -log10(P_CUTOFF), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed", color = "black", alpha = 0.5) +
    
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    
    geom_text_repel(
      data = label_data,
      aes(label = gene_id),
      size = 3.5, box.padding = 0.5, max.overlaps = 20, show.legend = FALSE, color = "black"
    ) +
    
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
    
    # Axis Labels match Legend Title
    labs(x = "log2(FoldChange)", y = "-log10(P.adjust)", title = plot_title) +
    
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # --- Step 8: Save ---
  tryCatch({
    ggsave(file.path(OUT_DIR, paste0(f_base, ".volcano.pdf")), p, width = FIXED_WIDTH, height = FIXED_HEIGHT, device = cairo_pdf)
    ggsave(file.path(OUT_DIR, paste0(f_base, ".volcano.png")), p, width = FIXED_WIDTH, height = FIXED_HEIGHT, dpi = 300, bg = "white")
    cat("   [SUCCESS] Saved.\n")
  }, error = function(e) cat(sprintf("   [ERROR] Save failed: %s\n", e$message)))
}

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
