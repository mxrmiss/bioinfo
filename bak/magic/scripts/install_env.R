#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: install_env.R (v2.0 Complete Dependency Manager)
# Description: 
#   - Automates installation of ALL R packages required by Magic Plotting Kit.
#   - Handles CRAN packages.
#   - Handles Bioconductor packages (ComplexHeatmap).
#   - Checks system-level capabilities (Cairo).
# ==============================================================================

# --- 1. CRAN Packages List ---
cran_packages <- c(
  "ggplot2",      # Plotting engine
  "dplyr",        # Data manipulation
  "readr",        # File reading
  "stringr",      # String processing
  "scales",       # Axis/Legend scaling
  "tibble",       # Data frames
  "tidyr",        # Data reshaping
  "ggrepel",      # Non-overlapping labels (Volcano/PCA)
  "scatterplot3d",# 3D PCA
  "ggvenn",       # Venn diagrams
  "magick",       # High-quality rasterization for Heatmaps
  "reshape2"      # Legacy data reshaping support
)

# --- 2. Install CRAN Packages ---
cat(">> Checking CRAN packages...\n")
installed_pkgs <- installed.packages()[,"Package"]
new_pkgs <- cran_packages[!(cran_packages %in% installed_pkgs)]

if(length(new_pkgs)) {
  cat(sprintf("   Installing %d missing CRAN packages: %s\n", length(new_pkgs), paste(new_pkgs, collapse=", ")))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
} else {
  cat("   All CRAN packages are already installed.\n")
}

# --- 3. Install Bioconductor Packages ---
cat(">> Checking Bioconductor packages (ComplexHeatmap)...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!require("ComplexHeatmap", quietly = TRUE)) {
  cat("   Installing ComplexHeatmap...\n")
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
} else {
  cat("   ComplexHeatmap is already installed.\n")
}

# --- 4. System Capabilities Check ---
cat(">> Checking System Capabilities...\n")

# Check Cairo (Critical for PDF output on Linux)
if (capabilities("cairo")) {
  cat("   [OK] Cairo graphics is supported (High-quality PDF ready).\n")
} else {
  cat("   [WARNING] Cairo graphics is NOT supported.\n")
  cat("             PDF output might fail or look pixelated on Linux.\n")
  cat("             Fix: Install 'libcairo2-dev' (Ubuntu) or 'cairo-devel' (CentOS).\n")
}

cat("\n================================================================\n")
cat(" Environment Setup Complete! You are ready to use Magic Plotting Kit.\n")
cat("================================================================\n")
