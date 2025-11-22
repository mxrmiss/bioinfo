#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: kmeans_trend.R (v1.0 Time-Course Trend)
# Description: 
#   - [Input] 1. Expression Matrix (TPM).
#             2. Sample Metadata (Time points).
#             3. MULTIPLE DEG Tables (for filtering significant genes).
#   - [Logic] 1. Filter DEGs by P.adjust & FC -> Take UNION.
#             2. Calculate Mean TPM per Time Point.
#             3. Z-score Transform -> K-means Clustering.
#   - [Output] Faceted Trend Plot, Cluster Tables, Gene Lists per Cluster.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(reshape2) # for melting data
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Time Order (CRITICAL) ---
# Define the chronological order of your groups.
# Must match the 'group' column in samples.tsv.
TIME_ORDER <- c("D0", "D1", "D8", "D15")

# --- DEG Filtering ---
# Only genes passing these thresholds in AT LEAST ONE comparison will be used.
FC_CUTOFF   <- 1.0    # |log2FC| threshold
P_CUTOFF    <- 0.05   # P.adjust threshold

# --- Clustering ---
K_CLUSTERS  <- 6      # Number of patterns to find (e.g., 6 or 9)
SEED_NUM    <- 123    # For reproducibility

# --- Input Control ---
# Pattern to identify DEG files in input/ (e.g. "vs", "DEG", "diff")
# Leave "" to scan all TSVs except matrix/samples.
DEG_FILE_PATTERN <- "" 

# --- Aesthetics ---
LINE_COLOR   <- "firebrick"  # Color of the centroid line
LINE_SIZE    <- 1.2          # Thickness of centroid
BG_ALPHA     <- 0.1          # Transparency of background genes (grey lines)
FONT_FAMILY  <- "sans"
BASE_SIZE    <- 12

# --- Paths ---
IN_DIR       <- "input"
OUT_DIR      <- file.path("output", "kmeans_trend")
LIST_DIR     <- file.path(OUT_DIR, "cluster_genes_lists") # Sub-folder for lists
TEMPLATE_DIR <- "templates"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  
  # 1. Matrix Template
  t1 <- file.path(TEMPLATE_DIR, "kmeans_matrix_template.tsv")
  if (!file.exists(t1)) {
    write_tsv(data.frame(
      gene_id = c("GeneA", "GeneB", "GeneC"),
      D0_1 = c(10, 100, 1), D0_2 = c(12, 90, 2), 
      D15_1 = c(50, 10, 50), D15_2 = c(48, 12, 45)
    ), t1)
  }
  
  # 2. Metadata Template
  t2 <- file.path(TEMPLATE_DIR, "kmeans_samples_template.tsv")
  if (!file.exists(t2)) {
    write_tsv(data.frame(sample = c("D0_1", "D0_2", "D15_1", "D15_2"), group = c("D0", "D0", "D15", "D15")), t2)
  }
  
  # 3. DEG Template
  t3 <- file.path(TEMPLATE_DIR, "kmeans_deg_template.tsv")
  if (!file.exists(t3)) {
    write_tsv(data.frame(gene_id = "GeneA", log2fc = 2.5, p_adjust = 0.001), t3)
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(LIST_DIR)) dir.create(LIST_DIR, recursive = TRUE)

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic K-means Trend Plotter v1.0\n")
cat("================================================================\n")

init_template_system()

# --- Step 1: Load Metadata & Matrix ---
all_files <- list.files(IN_DIR, full.names = TRUE)

# Find Matrix & Metadata
mat_file  <- all_files[grepl("tpm|count|matrix", basename(all_files), ignore.case = TRUE)][1]
meta_file <- file.path(IN_DIR, "samples.tsv")
if (!file.exists(meta_file)) meta_file <- all_files[grepl("sample|meta", basename(all_files), ignore.case = TRUE)][1]

if (is.na(mat_file)) stop("[ERROR] Expression Matrix not found.")
if (is.na(meta_file)) stop("[ERROR] Samples Metadata not found.")

cat(sprintf(">> Matrix: %s\n>> Metadata: %s\n", basename(mat_file), basename(meta_file)))

# Load
df_meta <- read_tsv(meta_file, show_col_types = FALSE)
cnames <- colnames(df_meta)
samp_col <- cnames[grep("sample", cnames, ignore.case = TRUE)[1]]
grp_col  <- cnames[grep("group|condition", cnames, ignore.case = TRUE)[1]]

# Filter Metadata to keep only groups in TIME_ORDER
sample_info <- df_meta %>% 
  select(sample = all_of(samp_col), group = all_of(grp_col)) %>%
  filter(group %in% TIME_ORDER)

if (nrow(sample_info) == 0) stop("[ERROR] No samples matched the groups in TIME_ORDER. Check your config.")

# Load Matrix & Subset
df_mat_raw <- read_tsv(mat_file, show_col_types = FALSE)
gene_col <- colnames(df_mat_raw)[1]
df_mat_clean <- df_mat_raw %>% select(gene_id = all_of(gene_col), any_of(sample_info$sample))

# --- Step 2: Scan DEG Files & Get Union of Genes ---
deg_files <- all_files[!grepl("tpm|count|matrix|sample|meta|template", basename(all_files), ignore.case = TRUE)]
if (DEG_FILE_PATTERN != "") deg_files <- deg_files[grepl(DEG_FILE_PATTERN, basename(deg_files))]

if (length(deg_files) == 0) stop("[ERROR] No DEG files found to define the gene set.")

cat(sprintf(">> Scanning %d DEG files for significant genes...\n", length(deg_files)))
sig_genes_union <- c()

for (f in deg_files) {
  df <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df)) next
  
  # Check headers
  if (!all(c("gene_id", "log2fc", "p_adjust") %in% colnames(df))) {
    cat(sprintf("   [SKIP] Invalid headers: %s\n", basename(f)))
    next
  }
  
  # Filter
  sig_df <- df %>% filter(p_adjust < P_CUTOFF & abs(log2fc) > FC_CUTOFF)
  count_n <- nrow(sig_df)
  
  if (count_n > 0) {
    cat(sprintf("   + %s: %d sig genes\n", basename(f), count_n))
    sig_genes_union <- c(sig_genes_union, sig_df$gene_id)
  }
}

sig_genes_union <- unique(sig_genes_union)
target_genes <- intersect(sig_genes_union, df_mat_clean$gene_id)
cat(sprintf(">> Final Union: %d genes found in expression matrix.\n", length(target_genes)))

if (length(target_genes) < 10) stop("[ERROR] Too few significant genes (<10) for clustering.")

# --- Step 3: Calculate Group Means & Z-score ---
cat(">> Calculating Group Means & Z-scores...\n")

# Filter matrix to target genes
sub_mat <- df_mat_clean %>% filter(gene_id %in% target_genes)

# Reshape to Long format for averaging
long_mat <- sub_mat %>%
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "tpm") %>%
  left_join(sample_info, by = "sample") %>%
  group_by(gene_id, group) %>%
  summarise(mean_tpm = mean(tpm), .groups = "drop")

# Pivot back to Wide (Genes x Groups) for Scaling
wide_mean_mat <- long_mat %>%
  pivot_wider(names_from = group, values_from = mean_tpm) %>%
  column_to_rownames("gene_id")

# Reorder columns by TIME_ORDER
wide_mean_mat <- wide_mean_mat[, TIME_ORDER]

# Z-score Transform (Row-wise)
# t(scale(t(x))) scales across columns (Time points) for each row (Gene)
z_mat <- t(scale(t(as.matrix(wide_mean_mat))))
z_mat <- na.omit(z_mat) # Remove genes with 0 variance (all zeros)

# --- Step 4: K-means Clustering ---
cat(sprintf(">> Running K-means (k=%d)...\n", K_CLUSTERS))
set.seed(SEED_NUM)
km_res <- kmeans(z_mat, centers = K_CLUSTERS, iter.max = 100, nstart = 25)

# Attach Cluster ID
gene_clusters <- data.frame(gene_id = rownames(z_mat), cluster = paste0("Cluster ", km_res$cluster))

# Prepare Plot Data
plot_data <- as.data.frame(z_mat) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_clusters, by = "gene_id") %>%
  pivot_longer(cols = all_of(TIME_ORDER), names_to = "Time", values_to = "Z_Score")

# Enforce Time Order Factor
plot_data$Time <- factor(plot_data$Time, levels = TIME_ORDER)

# Calculate Centroids (Mean Z-score per Cluster per Time)
centroids <- plot_data %>%
  group_by(cluster, Time) %>%
  summarise(Z_Score = mean(Z_Score), .groups = "drop")

# --- Step 5: Plotting ---
p <- ggplot() +
  # 1. Background Lines (All genes) - Grey & Transparent
  geom_line(data = plot_data, aes(x = Time, y = Z_Score, group = gene_id), 
            color = "grey80", alpha = BG_ALPHA) +
  
  # 2. Centroid Lines (Mean trend) - Colored & Thick
  geom_line(data = centroids, aes(x = Time, y = Z_Score, group = cluster, color = cluster), 
            linewidth = LINE_SIZE) +
  
  # 3. Facet by Cluster
  facet_wrap(~cluster, scales = "free_y") +
  
  # Aesthetics
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time Points", y = "Scaled Expression (Z-score)", 
       title = sprintf("Expression Trends (K=%d)", K_CLUSTERS),
       subtitle = sprintf("Based on %d significant genes (Union)", length(target_genes))) +
  theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  theme(
    legend.position = "none", # Color already shown by facet labels
    strip.background = element_rect(fill = "#E5E5E5"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# --- Step 6: Output ---
tryCatch({
  # Save Plot
  ggsave(file.path(OUT_DIR, "kmeans_trend_plot.pdf"), p, width = 10, height = 7, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, "kmeans_trend_plot.png"), p, width = 10, height = 7, dpi = 300, bg = "white")
  
  # Save Cluster Table
  write_tsv(gene_clusters, file.path(OUT_DIR, "gene_clusters.tsv"))
  
  # Save Individual Lists for Enrichment
  cat(">> Saving individual cluster lists...\n")
  clusters <- unique(gene_clusters$cluster)
  for (cl in clusters) {
    # Clean filename (Cluster 1 -> Cluster_1)
    safe_name <- gsub(" ", "_", cl)
    genes_in_cl <- gene_clusters %>% filter(cluster == cl) %>% pull(gene_id)
    
    out_txt <- file.path(LIST_DIR, paste0(safe_name, ".txt"))
    writeLines(genes_in_cl, out_txt)
  }
  
  cat(sprintf("   [SUCCESS] Saved Plot + Table + %d Gene Lists.\n", length(clusters)))
}, error = function(e) cat(sprintf("   [ERROR] Save failed: %s\n", e$message)))

cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
