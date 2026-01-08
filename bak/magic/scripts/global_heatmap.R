#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: global_heatmap.R (v1.1 Custom Title & Legend Order)
# Description: 
#   - [Config] Custom Title Control (Custom > Default > None).
#   - [Visual] Legend Order: Z-score (Top) -> Group (Bottom).
#   - [Logic] Selects Top N variable genes -> Bidirectional Clustering.
#   - [System] Auto-template & Magic directory standard.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(circlize)
  library(grid)
  library(tools)
  library(grDevices)
  if (!require("ComplexHeatmap", quietly = TRUE)) {
    stop("[ERROR] 'ComplexHeatmap' is missing. Install via: BiocManager::install('ComplexHeatmap')")
  }
  library(ComplexHeatmap)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Filtering ---
# How many genes to keep? (Ranking by Variance)
TOP_N_VAR_GENES <- 2000 

# --- Title Control [NEW] ---
SHOW_PLOT_TITLE   <- TRUE       # Master Switch
# Leave empty "" to use default title; Type string to override.
CUSTOM_PLOT_TITLE <- ""         

# --- Aesthetics ---
SHOW_ROW_NAMES  <- FALSE    # Hide gene names (too many to read)
SHOW_COL_NAMES  <- TRUE     # Show sample names
CLUSTER_SAMPLES <- TRUE     # Show Dendrogram

# --- Colors ---
COLOR_LOW    <- "navy"      
COLOR_MID    <- "white"     
COLOR_HIGH   <- "firebrick" 
FONT_FAMILY  <- "sans"      

# --- Paths ---
IN_DIR       <- "input"
OUT_DIR      <- file.path("output", "global_heatmap") # Updated folder name
TEMPLATE_DIR <- "templates"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  
  t1 <- file.path(TEMPLATE_DIR, "heatmap_matrix_template.tsv")
  if (!file.exists(t1)) {
    write_tsv(data.frame(
      gene_id = c("GeneA", "GeneB", "GeneC"),
      S1 = c(10, 5, 1), S2 = c(11, 6, 2), S3 = c(2, 1, 10), S4 = c(3, 2, 11)
    ), t1)
  }
  
  t2 <- file.path(TEMPLATE_DIR, "heatmap_samples_template.tsv")
  if (!file.exists(t2)) {
    write_tsv(data.frame(sample = c("S1","S2","S3","S4"), group = c("Ctrl","Ctrl","Trt","Trt")), t2)
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

generate_smart_palette <- function(groups) {
  n <- length(unique(groups))
  base_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", 
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  if (n <= length(base_colors)) {
    cols <- base_colors[1:n]
  } else {
    cols <- hcl.colors(n, palette = "Spectral")
  }
  names(cols) <- sort(unique(groups))
  return(cols)
}

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic Global Heatmap v1.1 (Custom Title)\n")
cat("================================================================\n")

init_template_system()

# --- Step 1: Load Data ---
all_files <- list.files(IN_DIR, full.names = TRUE)

# Find Matrix
mat_file <- all_files[grepl("tpm|count|matrix", basename(all_files), ignore.case = TRUE)][1]
if (is.na(mat_file)) stop("[ERROR] Expression Matrix (TPM) not found.")

# Find Metadata
meta_file <- file.path(IN_DIR, "samples.tsv")
if (!file.exists(meta_file)) meta_file <- all_files[grepl("sample|meta", basename(all_files), ignore.case = TRUE)][1]
if (is.na(meta_file)) stop("[ERROR] Samples Metadata not found.")

cat(sprintf(">> Matrix: %s\n>> Metadata: %s\n", basename(mat_file), basename(meta_file)))

df_mat_raw <- read_tsv(mat_file, show_col_types = FALSE)
df_meta <- read_tsv(meta_file, show_col_types = FALSE)

cnames <- colnames(df_meta)
samp_col <- cnames[grep("sample", cnames, ignore.case = TRUE)[1]]
grp_col  <- cnames[grep("group|condition", cnames, ignore.case = TRUE)[1]]
if (is.na(samp_col) || is.na(grp_col)) stop("[ERROR] Metadata missing 'sample' or 'group' column.")

sample_info <- df_meta %>% select(sample = all_of(samp_col), group = all_of(grp_col))
valid_samples <- sample_info$sample

gene_col <- colnames(df_mat_raw)[1]
df_mat_clean <- df_mat_raw %>% select(all_of(gene_col), any_of(valid_samples))
mat_data <- as.matrix(df_mat_clean[,-1])
rownames(mat_data) <- df_mat_clean[[1]]

# --- Step 2: Select Top Variable Genes ---
cat(">> Calculating Variance...\n")
gene_vars <- apply(mat_data, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(TOP_N_VAR_GENES, length(gene_vars))]

cat(sprintf("   Selected Top %d genes.\n", length(top_genes)))

sub_mat <- mat_data[top_genes, ]
z_mat <- t(scale(t(sub_mat)))
z_mat[is.na(z_mat)] <- 0
z_mat[z_mat > 2] <- 2; z_mat[z_mat < -2] <- -2

# --- Step 3: Plotting ---
groups <- sample_info$group
group_cols <- generate_smart_palette(groups)

col_ha <- HeatmapAnnotation(
  Group = groups,
  col = list(Group = group_cols),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(0.3, "cm")
)

col_fun <- colorRamp2(c(-2, 0, 2), c(COLOR_LOW, COLOR_MID, COLOR_HIGH))

# Title Logic
final_title <- ""
if (SHOW_PLOT_TITLE) {
  if (CUSTOM_PLOT_TITLE != "") {
    final_title <- CUSTOM_PLOT_TITLE
  } else {
    final_title <- paste0("Global Heatmap (Top ", length(top_genes), " Var. Genes)")
  }
}

ht <- Heatmap(
  z_mat,
  name = "Z-score",
  col = col_fun,
  
  cluster_rows = TRUE,    
  show_row_dend = FALSE,  
  show_row_names = SHOW_ROW_NAMES,
  
  cluster_columns = CLUSTER_SAMPLES, 
  show_column_dend = TRUE,           
  show_column_names = SHOW_COL_NAMES,
  column_names_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY),
  
  top_annotation = col_ha,
  column_title = final_title,
  column_title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT_FAMILY),
  border = TRUE,
  use_raster = TRUE
)

# --- Step 4: Save ---
w_inch <- 8
h_inch <- 10

out_pdf <- file.path(OUT_DIR, "Global_Clustering_Heatmap.pdf")
out_png <- file.path(OUT_DIR, "Global_Clustering_Heatmap.png")

tryCatch({
  pdf(out_pdf, width = w_inch, height = h_inch)
  # [FIX] merge_legend = TRUE forces stacking: Heatmap Legend (Z-score) on Top, Annotation (Group) below.
  draw(ht, merge_legend = TRUE) 
  dev.off()
  
  png(out_png, width = w_inch*300, height = h_inch*300, res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  
  cat("   [SUCCESS] Global Heatmap Saved.\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat(sprintf("   [ERROR] Plot failed: %s\n", e$message))
})

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
