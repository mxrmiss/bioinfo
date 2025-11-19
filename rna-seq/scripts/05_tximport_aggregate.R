#!/usr/bin/env Rscript
# 05_tximport_aggregate.R —— 按契约从 quant.sf 聚合到“基因层”（vNext 版）
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
    if (is.list(defaults[[nm]])) {
      defaults[[nm]] <- merge_list(user[[nm]], defaults[[nm]])
    } else if (!is.null(user[[nm]])) {
      defaults[[nm]] <- user[[nm]]
    }
  }
  defaults
}

defaults <- list(
  dirs = list(
    matrix = "results/05_matrix",
    maps   = "results/03_maps"
  ),
  tximport = list(
    counts_from_abundance = "no",
    matrix_stats = list(report_libsize_quantiles = c(0.25, 0.5, 0.75))
  )
)

cfg <- read_yaml(CONFIG_PATH)
cfg <- merge_list(cfg, defaults)

matrix_dir <- cfg$dirs$matrix
maps_dir   <- cfg$dirs$maps

meta_file  <- file.path(matrix_dir, "tximport_meta.tsv")
tx2gene_fp <- file.path(maps_dir, "tx2gene.clean.tsv")

if (!file.exists(meta_file))  stop(sprintf("未找到 tximport_meta.tsv：%s", meta_file))
if (!file.exists(tx2gene_fp)) stop(sprintf("未找到 tx2gene.clean.tsv：%s", tx2gene_fp))

dir.create(matrix_dir, showWarnings = FALSE, recursive = TRUE)

message("[05] 读取 tximport_meta.tsv: ", meta_file)
meta <- fread(meta_file)

if (!all(c("sample_id", "quant_file") %in% names(meta))) {
  stop("[05] tximport_meta.tsv 需要包含列：sample_id, quant_file")
}

# ========== 预检：quant.sf:Name 与 tx2gene.clean:transcript_id 交集比例 ==========
message("[05] 预检 quant.sf:Name 与 tx2gene.clean:transcript_id 的对齐情况...")

tx2gene_dt <- fread(tx2gene_fp)
if (!all(c("transcript_id", "gene_id") %in% names(tx2gene_dt))) {
  stop("[05] tx2gene.clean.tsv 必须包含列：transcript_id, gene_id")
}
tx_ids <- unique(tx2gene_dt$transcript_id)

# 收集所有 quant.sf 的 Name 列（Name 为 tx 层主键）
all_names <- character(0)
for (qf in meta$quant_file) {
  if (!file.exists(qf)) stop(sprintf("[05] quant 文件不存在：%s", qf))
  # 只读 Name 列，节省内存
  dt <- fread(qf, select = "Name")
  all_names <- c(all_names, dt$Name)
}
all_names <- unique(all_names)

intersect_ids <- intersect(all_names, tx_ids)
ratio <- length(intersect_ids) / length(all_names)

message(sprintf("[05] 交集统计：|Name ∩ transcript_id| / |Name| = %.4f (%d / %d)",
                ratio, length(intersect_ids), length(all_names)))

if (ratio < 0.95) {
  stop(sprintf(
    "[05] 预检失败：quant.sf:Name 与 tx2gene.clean:transcript_id 的交集比例仅为 %.4f (<0.95)。\n请检查 02/03 的 ID 口径或 annotations.id_cleanup 配置。",
    ratio
  ))
}

# ========== 正式 tximport ==========
files <- meta$quant_file
names(files) <- meta$sample_id

message("[05] 使用 tximport 聚合为基因层...")
tx2gene <- as.matrix(tx2gene_dt[, .(transcript_id, gene_id)])
counts_from_abundance <- cfg$tximport$counts_from_abundance

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = counts_from_abundance
)

# 输出 gene_counts.tsv / gene_tpm.tsv
gene_counts_fp <- file.path(matrix_dir, "counts", "gene_counts.tsv")
gene_tpm_fp    <- file.path(matrix_dir, "tpms",   "gene_tpm.tsv")
dir.create(dirname(gene_counts_fp), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(gene_tpm_fp),    showWarnings = FALSE, recursive = TRUE)

message("[05] 写出基因层 counts: ", gene_counts_fp)
write_tsv(as.data.frame(cbind(gene_id = rownames(txi$counts), as.data.frame(txi$counts))), gene_counts_fp)

message("[05] 写出基因层 TPM: ", gene_tpm_fp)
write_tsv(as.data.frame(cbind(gene_id = rownames(txi$abundance), as.data.frame(txi$abundance))), gene_tpm_fp)

# ========== 生成 matrix_stats.tsv（vNext 约定表头） ==========
matrix_stats_fp <- file.path(matrix_dir, "matrix_stats.tsv")
message("[05] 生成矩阵体检表: ", matrix_stats_fp)

counts_mat <- txi$counts
n_samples <- ncol(counts_mat)
n_genes   <- nrow(counts_mat)

# 每个基因的“零表达比例”（以 sample 为单位），取其中位数
zero_fraction_per_gene <- rowSums(counts_mat == 0) / n_samples
zero_fraction_median   <- median(zero_fraction_per_gene)

# 每个样本的 library size
libsize <- colSums(counts_mat)
qs <- cfg$tximport$matrix_stats$report_libsize_quantiles
if (length(qs) < 3) qs <- c(0.25, 0.5, 0.75)
q_vals <- as.numeric(quantile(libsize, probs = qs, names = FALSE))

# 这里约定：第一个作为 Q1，第二个作为 median，第三个作为 Q3
libsize_q1     <- q_vals[1]
libsize_median <- q_vals[2]
libsize_q3     <- q_vals[3]

stats_dt <- data.table(
  n_samples = n_samples,
  n_genes   = n_genes,
  zero_fraction_median = zero_fraction_median,
  libsize_median = libsize_median,
  libsize_q1     = libsize_q1,
  libsize_q3     = libsize_q3
)

fwrite(stats_dt, matrix_stats_fp, sep = "\t")

message("[05] 完成 tximport 聚合与矩阵体检。")