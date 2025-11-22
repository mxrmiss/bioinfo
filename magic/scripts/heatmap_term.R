#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: heatmap_term.R (v1.4 Custom Term Edition)
# Description: 
#   - [Feat] 'TARGET_TERMS': Specify exact terms to plot (Overrides Top N).
#   - [Fix] Legend Order: Z-score (Top) -> Group (Bottom).
#   - [System] Smart Color Engine & Template Guard maintained.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
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

# --- Input Control ---
# Leave empty "" to auto-scan all enrichment files.
TARGET_FILE <- ""

# --- [NEW] Target Specific Terms ---
# Option A: Leave empty c() to use TOP_N_TERMS logic (Auto-select).
# Option B: Fill in specific term names to plot ONLY these (e.g. c("DNA replication")).
TARGET_TERMS <- c() 

# --- Heatmap Logic ---
TOP_N_TERMS        <- 5    # Used only if TARGET_TERMS is empty
MAX_GENES_PER_TERM <- 5     # Max genes per pathway
SHOW_GENE_NAMES    <- TRUE  

# --- Title Control ---
SHOW_PLOT_TITLE   <- TRUE   
CUSTOM_PLOT_TITLE <- ""     

# --- Aesthetics ---
COLOR_LOW    <- "navy"      
COLOR_MID    <- "white"     
COLOR_HIGH   <- "firebrick" 
FONT_FAMILY  <- "sans"      

# --- Paths ---
IN_DIR       <- "input"
OUT_DIR      <- file.path("output", "heatmap_term")
TEMPLATE_DIR <- "templates"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  
  t1 <- file.path(TEMPLATE_DIR, "heatmap_enrich_template.tsv")
  if (!file.exists(t1)) {
    write_tsv(data.frame(
      term_name = c("Cell cycle", "DNA replication"),
      p_adjust = c(1e-5, 1e-4),
      gene_ids = c("GeneA;GeneB", "GeneC;GeneA") 
    ), t1)
  }
  
  t2 <- file.path(TEMPLATE_DIR, "heatmap_matrix_template.tsv")
  if (!file.exists(t2)) {
    write_tsv(data.frame(
      gene_id = c("GeneA", "GeneB", "GeneC"),
      S1 = c(10, 5, 1), S2 = c(11, 6, 2), S3 = c(2, 1, 10), S4 = c(3, 2, 11)
    ), t2)
  }
  
  t3 <- file.path(TEMPLATE_DIR, "heatmap_samples_template.tsv")
  if (!file.exists(t3)) {
    write_tsv(data.frame(sample = c("S1","S2","S3","S4"), group = c("Ctrl","Ctrl","Trt","Trt")), t3)
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
cat(" Magic Pathway Heatmap v1.4 (Custom Term)\n")
cat("================================================================\n")

init_template_system()

# --- Step 1: Load Matrix & Metadata ---
all_files <- list.files(IN_DIR, full.names = TRUE)

# Find Matrix
mat_file <- all_files[grepl("tpm|count|matrix", basename(all_files), ignore.case = TRUE)][1]
if (is.na(mat_file)) stop("[ERROR] No expression matrix found.")

# Find Metadata
meta_file <- file.path(IN_DIR, "samples.tsv")
if (!file.exists(meta_file)) meta_file <- all_files[grepl("sample|meta", basename(all_files), ignore.case = TRUE)][1]
if (is.na(meta_file)) stop("[ERROR] No sample metadata found.")

cat(sprintf(">> Matrix: %s\n>> Metadata: %s\n", basename(mat_file), basename(meta_file)))

df_mat_raw <- read_tsv(mat_file, show_col_types = FALSE)
df_meta <- read_tsv(meta_file, show_col_types = FALSE)

cnames <- colnames(df_meta)
samp_col <- cnames[grep("sample", cnames, ignore.case = TRUE)[1]]
grp_col  <- cnames[grep("group|condition", cnames, ignore.case = TRUE)[1]]
if (is.na(samp_col) || is.na(grp_col)) stop("[ERROR] Metadata issue.")

sample_info <- df_meta %>% 
  select(sample = all_of(samp_col), group = all_of(grp_col)) %>%
  arrange(group) 
valid_samples <- sample_info$sample

mat_gene_col <- colnames(df_mat_raw)[1]
df_mat_clean <- df_mat_raw %>% select(all_of(mat_gene_col), any_of(valid_samples))
mat_data <- as.matrix(df_mat_clean[,-1])
rownames(mat_data) <- df_mat_clean[[1]]

gene_vars <- apply(mat_data, 1, var)

# --- Step 2: Locate Enrichment Files ---
enrich_files <- c()
if (TARGET_FILE != "") {
  t_path <- file.path(IN_DIR, TARGET_FILE)
  if (file.exists(t_path)) enrich_files <- c(t_path)
  else stop(sprintf("[ERROR] Target file not found: %s", t_path))
} else {
  candidates <- all_files[!grepl("tpm|count|matrix|sample|meta|template", basename(all_files), ignore.case = TRUE)]
  enrich_files <- candidates
}

if (length(enrich_files) == 0) stop("[ERROR] No enrichment files found.")

# --- Step 3: Iterate & Plot ---
for (f_path in enrich_files) {
  
  fname <- basename(f_path)
  f_base <- file_path_sans_ext(fname)
  cat(sprintf(">> Processing: %s\n", fname))
  
  df_enrich <- tryCatch(read_tsv(f_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df_enrich) || nrow(df_enrich) == 0) next
  
  if (!all(c("term_name", "p_adjust", "gene_ids") %in% colnames(df_enrich))) {
    cat("   [SKIP] Invalid headers.\n"); next
  }
  
  # 3.1 Select Terms (Custom vs Auto)
  if (length(TARGET_TERMS) > 0) {
    # [NEW] Filter by specific terms
    target_terms_data <- df_enrich %>% filter(term_name %in% TARGET_TERMS)
    if (nrow(target_terms_data) == 0) {
      cat("   [SKIP] None of the TARGET_TERMS found in this file.\n")
      next
    }
    # Keep user order if possible, or sort by P
    plot_data <- target_terms_data
    cat(sprintf("   [INFO] Found %d / %d specified terms.\n", nrow(plot_data), length(TARGET_TERMS)))
  } else {
    # Auto Top N
    plot_data <- df_enrich %>% arrange(p_adjust) %>% head(TOP_N_TERMS)
  }
  
  # 3.2 Expand & Filter
  expanded_data <- plot_data %>%
    select(term_name, gene_ids) %>%
    mutate(gene_ids = str_split(gene_ids, "[;/]")) %>%
    unnest(gene_ids) %>%
    mutate(gene_ids = str_trim(gene_ids)) %>%
    filter(gene_ids %in% rownames(mat_data)) %>%
    mutate(variance = gene_vars[gene_ids]) %>%
    group_by(term_name) %>%
    arrange(desc(variance)) %>%
    slice_head(n = MAX_GENES_PER_TERM) %>%
    ungroup()
  
  if (nrow(expanded_data) == 0) { 
    cat("   [SKIP] No matching genes found.\n"); next 
  }
  
  # 3.3 Prepare Matrix
  plot_genes <- expanded_data$gene_ids
  plot_terms <- expanded_data$term_name
  
  sub_mat <- mat_data[plot_genes, , drop=FALSE]
  sub_mat <- sub_mat[, sample_info$sample, drop=FALSE]
  
  z_mat <- t(scale(t(sub_mat)))
  z_mat[is.na(z_mat)] <- 0
  z_mat[z_mat > 2] <- 2
  z_mat[z_mat < -2] <- -2
  
  # 3.4 Annotation
  groups <- sample_info$group
  group_cols <- generate_smart_palette(groups) 
  
  col_ha <- HeatmapAnnotation(
    Group = groups,
    col = list(Group = group_cols),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.3, "cm")
  )
  
  # 3.5 Plot
  final_title <- ""
  if (SHOW_PLOT_TITLE) {
    if (CUSTOM_PLOT_TITLE != "") {
      final_title <- CUSTOM_PLOT_TITLE
    } else {
      final_title <- paste("Heatmap:", gsub("_", " ", f_base))
    }
  }
  
  col_fun <- colorRamp2(c(-2, 0, 2), c(COLOR_LOW, COLOR_MID, COLOR_HIGH))
  
  ht <- Heatmap(
    z_mat,
    name = "Z-score",
    col = col_fun,
    
    cluster_rows = FALSE,
    row_split = factor(plot_terms, levels = unique(plot_terms)), 
    
    row_title_rot = 0,    
    row_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = FONT_FAMILY),
    show_row_names = SHOW_GENE_NAMES,
    row_names_gp = gpar(fontsize = 8, fontface = "plain", fontfamily = FONT_FAMILY),
    
    cluster_columns = FALSE, 
    column_split = sample_info$group, 
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY),
    
    top_annotation = col_ha,
    column_title = final_title,
    border = TRUE
  )
  
  # 3.6 Save
  h_inch <- max(4, nrow(z_mat) * 0.15 + 2) 
  w_inch <- max(8, ncol(z_mat) * 0.3 + 4) 
  
  out_pdf <- file.path(OUT_DIR, paste0("heatmap_", f_base, ".pdf"))
  out_png <- file.path(OUT_DIR, paste0("heatmap_", f_base, ".png"))
  
  tryCatch({
    pdf(out_pdf, width = w_inch, height = h_inch)
    # [FIX] merge_legend = TRUE puts Z-score legend on top of Group legend
    draw(ht, merge_legend = TRUE) 
    dev.off()
    
    png(out_png, width = w_inch*300, height = h_inch*300, res = 300)
    draw(ht, merge_legend = TRUE)
    dev.off()
    
    cat(sprintf("   [SUCCESS] Saved Heatmap (%d genes)\n", nrow(z_mat)))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("   [ERROR] Plot failed: %s\n", e$message))
  })
}

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
