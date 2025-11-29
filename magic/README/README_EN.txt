================================================================================
MAGIC PLOTTING KIT (v2.5 Final)
Top-tier Transcriptome Data Visualization Toolkit
================================================================================

This toolkit contains 9 standardized R scripts designed to generate 
publication-quality (Nature/Science/Cell standard) figures. 
It covers Quality Control, Differential Expression, Enrichment, and Time-course 
Trend analysis.

All scripts feature:
1. Template Guard (Auto-validation of inputs).
2. Publication Aesthetics (Sans-serif fonts, Standard palettes).
3. Engineering Workflow (Auto-archiving outputs).

================================================================================
1. DIRECTORY STRUCTURE
================================================================================

Please strictly maintain the following structure:

magic/
├── input/          # [USER ACTION] Place your .tsv data files here.
├── output/         # [AUTO] Results are automatically saved here.
├── scripts/        # [CORE] Contains 9 .R scripts. DO NOT modify logic inside.
├── templates/      # [AUTO] Contains header templates (Reference only).
└── README_EN.txt   # This documentation.

================================================================================
2. QUICK START
================================================================================

Step 1: Install Environment (Run Once)
--------------------------------------
$ Rscript scripts/install_env.R

Step 2: Prepare Data
--------------------------------------
Place your analysis results (e.g., gene_tpm.tsv, DEG_all.tsv, samples.tsv) 
into the 'input/' folder. Check 'templates/' if you are unsure about headers.

Step 3: Run Scripts (THE IRON RULE)
--------------------------------------
ALWAYS run from the 'magic/' root directory.

[Correct] $ Rscript scripts/pca.R
[Wrong]   $ cd scripts; Rscript pca.R

================================================================================
3. SCRIPT INVENTORY (The Magic 9)
================================================================================

--- Group A: QC & Overview ---

1. pca.R
   - Function: Draws 2D (Biplot) and 3D PCA plots.
   - Input: gene_tpm.tsv, samples.tsv.
   - Highlight: 2D with arrows/hulls; 3D enclosed box style (No drop lines).

2. global_heatmap.R
   - Function: Heatmap of Top variable genes.
   - Use Case: Supplementary Figure for sample reproducibility check.

--- Group B: Differential Expression ---

3. volcano.R
   - Function: Volcano plot for DEGs.
   - Input: DEG table (gene_id, log2fc, p_adjust).
   - Highlight: Dual-gradient coloring, hybrid labeling, symmetrical axis.

4. venn.R
   - Function: Venn Diagram (<5 sets) or UpSet Plot (>4 sets).
   - Input: Multiple DEG tables.
   - Highlight: Auto-switches mode based on file count. Exports intersection lists.

--- Group C: Enrichment (Mechanism) ---

5. enrichment_bubble.R
   - Function: Bubble plot for GO/KEGG.
   - Highlight: Smooth S-curve sorting, fixed aspect ratio for sparse data.

6. enrichment_bar.R
   - Function: Bar plot for GO/KEGG.
   - Highlight: Inverted pyramid sorting, integer count labels.

7. heatmap_term.R
   - Function: Pathway-specific gene heatmap.
   - Input: Enrichment table + TPM Matrix + Metadata.
   - Use Case: Visualizes core genes for specific pathways (e.g., "Cell Cycle").
   - Config: Use 'TARGET_TERMS' to specify exact pathways to plot.

--- Group D: Time-Course Trends ---

8. kmeans_trend.R
   - Function: K-means clustering of expression trends.
   - Input: TPM Matrix + Metadata + DEG tables (Union).
   - Use Case: Identifies major temporal patterns (e.g., Early vs Late response).
   - Note: Set 'TIME_ORDER' in the script to match your experimental design.

--- Infrastructure ---

9. install_env.R
   - Function: Auto-installs CRAN & Bioconductor dependencies.

================================================================================
4. TROUBLESHOOTING
================================================================================

Q: Error "Invalid Headers"?
A: Check the templates generated in 'magic/templates/'. Your input columns must 
   match exactly (case-sensitive).

Q: Error "failed to load cairo DLL"?
A: Your Linux system misses graphics libraries. Contact admin to install:
   - Ubuntu: libcairo2-dev
   - CentOS: cairo-devel

Q: 3D PCA not rotating?
A: The script produces static PDF (vector) for publication. For interactive 
   rotation, use Python/Plotly (not included).

================================================================================
