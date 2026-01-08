#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: pca.R (v3.9 Bugfix Edition)
# Description: 
#   - [Fix] Restored missing 'SHOW_HULL' configuration variable.
#   - [Feat] Custom Titles, Biplot Arrows, Convex Hulls, Boxed 3D.
#   - [System] Input Control & Smart Color Engine.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tools)
  library(ggrepel)
  library(grDevices) 
  if (!require("scatterplot3d", quietly = TRUE)) {
    stop("[ERROR] Package 'scatterplot3d' is missing. Please install: install.packages('scatterplot3d')")
  }
  library(scatterplot3d)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Input Control ---
# Leave empty "" to auto-scan, or specify e.g. "gene_tpm.tsv"
TARGET_FILE <- ""

# --- Title Control ---
SHOW_PLOT_TITLE   <- TRUE       # Master Switch
# Custom Title: Leave empty "" to use filename; Type string to override.
CUSTOM_PLOT_TITLE <- ""         

# --- Biplot Settings (Gene Arrows) ---
SHOW_ARROWS     <- FALSE       # Draw variable arrows?
TOP_N_LOADINGS  <- 5          # How many top genes to show

# --- Label Control ---
SHOW_SAMPLE_LABELS <- FALSE    # Show sample names?

# --- Aesthetics ---
SHOW_HULL       <- TRUE       # [RESTORED] Draw Convex Hulls?
SHOW_ELLIPSE    <- FALSE      # (Keep FALSE if using Hull)
POINT_SIZE      <- 3.5        
FONT_FAMILY     <- "sans"     
BASE_SIZE       <- 12

# --- Canvas (2D) ---
FIXED_WIDTH_2D  <- 8.0
FIXED_HEIGHT_2D <- 7.5        

# --- Output ---
IN_DIR          <- "input"
OUT_DIR         <- file.path("output", "pca")
TEMPLATE_DIR    <- "templates"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

init_template_system <- function() {
  if (!dir.exists(TEMPLATE_DIR)) dir.create(TEMPLATE_DIR, recursive = TRUE)
  
  mat_tmpl <- file.path(TEMPLATE_DIR, "pca_expression_template.tsv")
  if (!file.exists(mat_tmpl)) {
    ex_mat <- data.frame(
      gene_id = c("GeneA", "GeneB", "GeneC", "GeneD"),
      Control_1 = c(10, 0, 500, 20),
      Control_2 = c(11, 0, 480, 22),
      Treat_1   = c(2, 100, 10, 50),
      Treat_2   = c(3, 95, 12, 45)
    )
    write_tsv(ex_mat, mat_tmpl)
  }
  
  meta_tmpl <- file.path(TEMPLATE_DIR, "pca_metadata_template.tsv")
  if (!file.exists(meta_tmpl)) {
    ex_meta <- data.frame(
      sample = c("Control_1", "Control_2", "Treat_1", "Treat_2"),
      group  = c("Control", "Control", "Treat", "Treat")
    )
    write_tsv(ex_meta, meta_tmpl)
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

generate_smart_palette <- function(groups) {
  n <- length(unique(groups))
  base_colors <- c(
  "#F79D93", "#95C8F2", "#F6CD96", "#9FD5CB", "#F6C6E7",
  "#A99BEF", "#F5B07E", "#A7D9F7", "#F3A6C9", "#B6E0B6",
  "#FBE3A8", "#8FE3F2", "#F4B1A8", "#8FB1F2", "#F0D07A",
  "#7FD3C8", "#F7B6D2", "#C3A3F5", "#D6E7B5", "#9E90E6"
)
  if (n <= length(base_colors)) return(base_colors[1:n])
  else return(hcl.colors(n, palette = "Spectral"))
}

generate_shapes <- function(groups) {
  n <- length(unique(groups))
  base_shapes <- c(19, 17, 15, 18, 3, 4, 8, 10, 16, 0, 1, 2, 5, 6) 
  if (n <= length(base_shapes)) return(base_shapes[1:n])
  else return(rep(19, n)) 
}

get_hull_data <- function(df) {
  df %>% group_by(group) %>% slice(chull(PC1, PC2))
}

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================

cat("================================================================\n")
cat(" Magic PCA Plotter v3.9 (Bugfix Edition)\n")
cat("================================================================\n")

init_template_system()

# --- Step 0: File Selection ---
target_path <- NULL
if (TARGET_FILE != "") {
  t_path <- file.path(IN_DIR, TARGET_FILE)
  if (file.exists(t_path)) target_path <- t_path
  else stop(sprintf("[ERROR] Specified file not found: %s", t_path))
} else {
  cat("Mode: Auto-scan 'input/' for expression matrices...\n")
  all_files <- list.files(IN_DIR, pattern = "\\.tsv$", full.names = TRUE)
  candidates <- all_files[!grepl("sample|meta|template", basename(all_files), ignore.case = TRUE)]
  priority <- candidates[grepl("tpm|count|matrix", basename(candidates), ignore.case = TRUE)]
  if (length(priority) > 0) target_path <- priority[1]
  else if (length(candidates) > 0) target_path <- candidates[1]
  
  if (is.null(target_path)) stop("[ERROR] No valid expression matrix found.")
}

f_base <- file_path_sans_ext(basename(target_path))
cat(sprintf(">> Processing: %s\n", basename(target_path)))

# --- Step 1: Load Data ---
df_mat <- read_tsv(target_path, show_col_types = FALSE)
mat_data <- df_mat %>% select(-1) %>% as.matrix()
rownames(mat_data) <- df_mat[[1]]
mat_data <- mat_data[rowSums(mat_data) > 0, ] 
log_mat <- log2(mat_data + 1)

# Metadata
meta_file <- file.path(IN_DIR, "samples.tsv")
if (!file.exists(meta_file)) meta_file <- file.path(IN_DIR, "pca_metadata.tsv")

sample_info <- NULL
if (file.exists(meta_file)) {
  cat(sprintf("   [INFO] Metadata found: %s\n", basename(meta_file)))
  raw_meta <- read_tsv(meta_file, show_col_types = FALSE)
  cnames <- colnames(raw_meta)
  samp_col <- cnames[grep("sample", cnames, ignore.case = TRUE)[1]]
  grp_col  <- cnames[grep("group|condition", cnames, ignore.case = TRUE)[1]]
  if (!is.na(samp_col) && !is.na(grp_col)) {
    sample_info <- raw_meta %>% select(sample = all_of(samp_col), group = all_of(grp_col))
  }
}

# --- Step 2: PCA Calculation ---
pca_input <- t(log_mat)
pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE) 

var_exp <- round(summary(pca_res)$importance[2, ] * 100, 2)
pc1_lab <- sprintf("PC1 (%s%%)", var_exp[1])
pc2_lab <- sprintf("PC2 (%s%%)", var_exp[2])
pc3_lab <- sprintf("PC3 (%s%%)", var_exp[3])

plot_df <- as.data.frame(pca_res$x)
plot_df$sample <- rownames(plot_df)

if (!is.null(sample_info)) {
  plot_df <- left_join(plot_df, sample_info, by = "sample")
  plot_df$group[is.na(plot_df$group)] <- "Unknown"
} else {
  cat("   [WARN] No metadata. Auto-inferring groups.\n")
  plot_df$group <- str_remove(plot_df$sample, "[_\\d]+$")
}
plot_df$group <- as.factor(plot_df$group)

# --- Step 3: Aesthetics ---
final_palette <- generate_smart_palette(plot_df$group)
names(final_palette) <- levels(plot_df$group)
if(length(levels(plot_df$group)) > length(final_palette)) final_palette <- generate_smart_palette(plot_df$group)

final_shapes <- generate_shapes(plot_df$group)
names(final_shapes) <- levels(plot_df$group)

# --- Step 4: Arrow Calculation ---
loadings <- as.data.frame(pca_res$rotation)
loadings$gene <- rownames(loadings)
loadings$magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)
top_loadings <- loadings %>% arrange(desc(magnitude)) %>% head(TOP_N_LOADINGS)

r <- min(max(plot_df$PC1)-min(plot_df$PC1), max(plot_df$PC2)-min(plot_df$PC2))
r_arrow <- max(top_loadings$magnitude)
top_loadings$PC1_scaled <- top_loadings$PC1 / r_arrow * (r/2) * 0.8
top_loadings$PC2_scaled <- top_loadings$PC2 / r_arrow * (r/2) * 0.8

# --- Prepare Title ---
final_title <- NULL
if (SHOW_PLOT_TITLE) {
  if (CUSTOM_PLOT_TITLE != "") {
    final_title <- CUSTOM_PLOT_TITLE
  } else {
    clean_name <- gsub("_", " ", f_base)
    final_title <- paste("PCA:", clean_name)
  }
}

# --- Step 5: Plotting 2D ---
limit_val <- max(abs(c(plot_df$PC1, plot_df$PC2))) * 1.2 
limit_range <- c(-limit_val, limit_val)
hull_df <- get_hull_data(plot_df)

p2d <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +

  # Convex Hulls
  {
    if (SHOW_HULL && min(table(plot_df$group)) >= 3) {
      geom_polygon(data = hull_df, aes(x = PC1, y = PC2, fill = group), 
                   alpha = 0.15, show.legend = FALSE)
    }
  } +
  
  # Arrows
  {
    if (SHOW_ARROWS) {
      geom_segment(data = top_loadings, 
                   aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
                   arrow = arrow(length = unit(0.2, "cm")), 
                   color = "grey40", alpha = 0.8)
    }
  } +
  
  # Arrow Labels
  {
    if (SHOW_ARROWS) {
      geom_text_repel(data = top_loadings, 
                      aes(x = PC1_scaled, y = PC2_scaled, label = gene),
                      color = "black", size = 3, fontface = "italic",
                      box.padding = 0.5, segment.color = "grey60")
    }
  } +

  # Points
  geom_point(data = plot_df, aes(x = PC1, y = PC2, color = group, shape = group), 
             size = POINT_SIZE) +
  
  # Sample Labels
  {
    if (SHOW_SAMPLE_LABELS) {
      geom_text_repel(data = plot_df, aes(x = PC1, y = PC2, label = sample, color = group), 
                      size = 3.5, show.legend = FALSE, max.overlaps = 20, 
                      box.padding = 0.6)
    }
  } +
  
  coord_fixed(ratio = 1, xlim = limit_range, ylim = limit_range) +
  scale_color_manual(values = final_palette) +
  scale_fill_manual(values = final_palette) +
  scale_shape_manual(values = final_shapes) +
  labs(x = pc1_lab, y = pc2_lab, title = final_title) +
  theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  theme(
    panel.grid = element_line(color = "grey92"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
  )

tryCatch({
  ggsave(file.path(OUT_DIR, paste0("pca_2d_", f_base, ".pdf")), p2d, width = FIXED_WIDTH_2D, height = FIXED_HEIGHT_2D, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0("pca_2d_", f_base, ".png")), p2d, width = FIXED_WIDTH_2D, height = FIXED_HEIGHT_2D, dpi = 300, bg = "white")
  cat("   [SUCCESS] 2D PCA Saved.\n")
}, error = function(e) cat(sprintf("   [ERROR] 2D Save failed: %s\n", e$message)))

# --- Step 6: Plotting 3D ---
cat("   [CALC] Generating 3D Plot...\n")
pdf_3d_path <- file.path(OUT_DIR, paste0("pca_3d_", f_base, ".pdf"))
png_3d_path <- file.path(OUT_DIR, paste0("pca_3d_", f_base, ".png"))

group_idx <- as.numeric(plot_df$group)
colors_3d <- final_palette[group_idx]
shapes_3d <- final_shapes[group_idx] 

all_coords <- c(plot_df$PC1, plot_df$PC2, plot_df$PC3)
max_range <- max(abs(all_coords)) * 1.1
lim_cube <- c(-max_range, max_range)

draw_3d <- function() {
  s3d <- scatterplot3d(
    x = plot_df$PC1, y = plot_df$PC2, z = plot_df$PC3,
    color = colors_3d, pch = shapes_3d, cex.symbols = 1.5, 
    box = TRUE, 
    xlim = lim_cube, ylim = lim_cube, zlim = lim_cube,
    main = if(SHOW_PLOT_TITLE) paste("3D", final_title) else NULL,
    xlab = pc1_lab, ylab = pc2_lab, zlab = pc3_lab,
    angle = 45, grid = TRUE, mar = c(5, 3, 4, 3) + 0.1
  )
  
  legend("topleft", legend = levels(plot_df$group),
         col = final_palette, pch = final_shapes, 
         inset = 0.05, bty = "n", title = "Group", cex = 0.8)
  
  if (SHOW_SAMPLE_LABELS) {
    text_coords <- s3d$xyz.convert(plot_df$PC1, plot_df$PC2, plot_df$PC3)
    text(text_coords$x, text_coords$y, labels = plot_df$sample, 
         cex = 0.7, pos = 4, offset = 0.5, col = colors_3d)
  }
}

pdf(pdf_3d_path, width = 8, height = 8)
draw_3d()
garbage <- dev.off()

png(png_3d_path, width = 2400, height = 2400, res = 300)
draw_3d()
garbage <- dev.off()

cat("   [SUCCESS] 3D PCA Saved.\n")
cat("----------------------------------------------------------------\n")
cat("All processing finished.\n")
