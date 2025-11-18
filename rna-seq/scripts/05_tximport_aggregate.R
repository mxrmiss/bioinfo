#!/usr/bin/env Rscript
# 05_tximport_aggregate.R —— 按契约从 quant.sf 聚合到“基因层”
# 只读 config.yaml；不接收命令行参数；不使用 ignoreTxVersion。

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(tximport)
  library(readr)
})

CONFIG_PATH <- "config.yaml"

merge_list <- function(user, defaults) {
  if (is.null(user)) return(defaults)
  for (nm in names(defaults)) {
    if (is.list(defaults[[nm]]) && !is.null(user[[nm]]) && is.list(user[[nm]])) {
      defaults[[nm]] <- merge_list(user[[nm]], defaults[[nm]])
    } else if (!is.null(user[[nm]])) {
      defaults[[nm]] <- user[[nm]]
    }
  }
  defaults
}

defaults <- list(
  dirs = list(
    maps   = "results/03_maps",
    matrix = "results/05_matrix"
  ),
  tximport = list(
    counts_from_abundance = "no",
    matrix_stats = list(report_libsize_quantiles = c(0.25, 0.5, 0.75))
  )
)

cfg <- read_yaml(CONFIG_PATH); cfg <- merge_list(cfg, defaults)

matrix_dir <- cfg$dirs$matrix
maps_dir   <- cfg$dirs$maps

meta_file  <- file.path(matrix_dir, "tximport_meta.tsv")
tx2gene_fp <- file.path(maps_dir, "tx2gene.clean.tsv")

if (!file.exists(meta_file))  stop(sprintf("未找到 tximport_meta.tsv：%s", meta_file))
if (!file.exists(tx2gene_fp)) stop(sprintf("未找到 tx2gene.clean.tsv：%s", tx2gene_fp))

cat("========== 05-R — tximport 聚合 ==========\n")
cat(sprintf("[Info] meta_file     = %s\n", meta_file))
cat(sprintf("[Info] tx2gene       = %s\n", tx2gene_fp))
cat(sprintf("[Info] countsFromAbundance = %s\n", cfg$tximport$counts_from_abundance))

# 读 meta 与 tx2gene（列名按契约）
meta <- fread(meta_file, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("sample","group","quant_sf") %in% colnames(meta))) {
  stop("tximport_meta.tsv 必须包含列：sample, group, quant_sf")
}
tx2g <- fread(tx2gene_fp, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("transcript_id","gene_id") %in% colnames(tx2g))) {
  stop("tx2gene.clean.tsv 必须包含列：transcript_id, gene_id")
}
colnames(tx2g)[match(c("transcript_id","gene_id"), colnames(tx2g))] <- c("tx","gene")

files <- meta$quant_sf; names(files) <- meta$sample

cat(sprintf("[Run] tximport(...) with %d samples, %d tx→gene mappings\n",
            length(files), nrow(tx2g)))

# 契约：不使用 ignoreTxVersion；版本一致由 02/03 已保证
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2g,
                countsFromAbundance = cfg$tximport$counts_from_abundance,
                dropInfReps = TRUE)

# 输出目录契约
counts_dir <- file.path(matrix_dir, "counts")
tpms_dir   <- file.path(matrix_dir, "tpms")
dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tpms_dir,   recursive = TRUE, showWarnings = FALSE)

# gene_counts.tsv
counts <- as.data.frame(txi$counts)
counts <- cbind(gene_id = rownames(counts), counts)
write_tsv(counts, file.path(counts_dir, "gene_counts.tsv"))

# gene_tpm.tsv
tpms <- as.data.frame(txi$abundance)
tpms <- cbind(gene_id = rownames(tpms), tpms)
write_tsv(tpms, file.path(tpms_dir, "gene_tpm.tsv"))

cat(sprintf("[Out] %s\n", file.path(counts_dir, "gene_counts.tsv")))
cat(sprintf("[Out] %s\n", file.path(tpms_dir,   "gene_tpm.tsv")))

# 矩阵体检
cnt_mat <- txi$counts
lib_sizes <- colSums(cnt_mat, na.rm = TRUE)
qs <- quantile(lib_sizes,
               probs = as.numeric(cfg$tximport$matrix_stats$report_libsize_quantiles),
               na.rm = TRUE, names = FALSE)
zero_gene_ratio <- if (nrow(cnt_mat) > 0) {
  mean(rowSums(cnt_mat > 0, na.rm = TRUE) == 0)
} else 0

stats <- data.frame(
  n_genes = nrow(cnt_mat),
  n_samples = ncol(cnt_mat),
  zero_gene_ratio = round(as.numeric(zero_gene_ratio), 6),
  library_size_q1 = as.numeric(qs[1]),
  library_size_median = as.numeric(qs[2]),
  library_size_q3 = as.numeric(qs[3])
)
write_tsv(stats, file.path(matrix_dir, "matrix_stats.tsv"))
cat(sprintf("[Out] %s\n", file.path(matrix_dir, "matrix_stats.tsv")))
cat("========== 05-R 完成 ==========\n")
