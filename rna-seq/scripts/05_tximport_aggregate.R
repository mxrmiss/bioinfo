#!/usr/bin/env Rscript

# 05_tximport_aggregate.R —— 基于 tximport 的表达矩阵聚合（vNext 版）
# 职责：
#   1) 读取 config.yaml，定位 tximport_meta.tsv 与 tx2gene.clean.tsv
#   2) 调用 tximport 聚合 salmon quant.sf → gene-level counts/TPM
#   3) 输出：
#        - results/05_matrix/counts/gene_counts.tsv
#        - results/05_matrix/tpms/gene_tpm.tsv
#        - results/05_matrix/matrix_stats.tsv
#   4) 兼容不同的 tximport_meta 列名：
#        - 样本列：sample 或 sample_id
#        - 定量文件列：quant_file 或 quant_sf
#
# 约定：
#   - 不从命令行传参，仅从项目根目录下 config.yaml 读取配置。
#   - 所有输出为 UTF-8 / LF / TSV，缺失记为 NA。
#   - 不修改任何 ID，只使用 03 产生的 tx2gene.clean.tsv。

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(tximport)
})

#------------------------- 小工具函数 -------------------------#

log_msg <- function(level = "INFO", ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", level, "] ", paste0(..., collapse = ""))
  cat(ts, msg, "\n")
}

stop_with <- function(...) {
  log_msg("ERROR", ...)
  stop(paste0(...), call. = FALSE)
}

# 从 list 安全取值（带默认）
get_or <- function(lst, ..., default = NULL) {
  cur <- lst
  keys <- list(...)
  for (k in keys) {
    if (is.null(cur) || is.null(cur[[k]])) {
      return(default)
    }
    cur <- cur[[k]]
  }
  cur
}

#------------------------- 读取配置 -------------------------#

config_file <- "config.yaml"
if (!file.exists(config_file)) {
  stop_with("找不到配置文件：", config_file)
}

log_msg("INFO", "使用配置文件：", normalizePath(config_file))
cfg <- yaml::read_yaml(config_file)

dir_matrix <- get_or(cfg, "dirs", "matrix", default = "results/05_matrix")
dir_maps   <- get_or(cfg, "dirs", "maps",   default = "results/03_maps")

dir_counts <- file.path(dir_matrix, "counts")
dir_tpms   <- file.path(dir_matrix, "tpms")

if (!dir.exists(dir_matrix)) dir.create(dir_matrix, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_counts)) dir.create(dir_counts, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_tpms))   dir.create(dir_tpms,   recursive = TRUE, showWarnings = FALSE)

meta_path    <- file.path(dir_matrix, "tximport_meta.tsv")
tx2gene_path <- file.path(dir_maps,   "tx2gene.clean.tsv")
matrix_stats_path <- file.path(dir_matrix, "matrix_stats.tsv")
counts_out_path   <- file.path(dir_counts, "gene_counts.tsv")
tpm_out_path      <- file.path(dir_tpms,   "gene_tpm.tsv")

log_msg("INFO", "dirs.matrix = ", normalizePath(dir_matrix))
log_msg("INFO", "dirs.maps   = ", normalizePath(dir_maps))
log_msg("INFO", "tximport_meta.tsv = ", meta_path)
log_msg("INFO", "tx2gene.clean.tsv = ", tx2gene_path)

#------------------------- 读取 tximport_meta.tsv -------------------------#

if (!file.exists(meta_path)) {
  stop_with("找不到 tximport_meta.tsv：", meta_path,
            " 请先运行 05_matrix_aggregate.py 生成该文件。")
}

meta <- data.table::fread(meta_path)
if (nrow(meta) == 0L) {
  stop_with("tximport_meta.tsv 为空，无法进行聚合。")
}

# 兼容不同列名：sample / sample_id
sample_col <- NULL
if ("sample" %in% names(meta)) {
  sample_col <- "sample"
} else if ("sample_id" %in% names(meta)) {
  sample_col <- "sample_id"
}

# 兼容不同列名：quant_file / quant_sf
quant_col <- NULL
if ("quant_file" %in% names(meta)) {
  quant_col <- "quant_file"
} else if ("quant_sf" %in% names(meta)) {
  quant_col <- "quant_sf"
}

missing_cols <- character()
if (is.null(sample_col)) missing_cols <- c(missing_cols, "sample 或 sample_id")
if (is.null(quant_col))  missing_cols <- c(missing_cols, "quant_file 或 quant_sf")

if (length(missing_cols) > 0L) {
  stop_with(
    "tximport_meta.tsv 缺少必需列：",
    paste(missing_cols, collapse = ", "),
    "；当前列为：", paste(names(meta), collapse = ", ")
  )
}

log_msg("INFO", "检测到样本列：", sample_col,
        "；定量文件列：", quant_col)

# 构建文件向量
files_vec <- meta[[quant_col]]
names(files_vec) <- meta[[sample_col]]

# 检查文件是否存在
missing_files <- files_vec[!file.exists(files_vec)]
if (length(missing_files) > 0L) {
  show_n <- min(5L, length(missing_files))
  stop_with(
    "以下 quant 文件不存在（仅展示前 ", show_n, " 个）：\n  ",
    paste(head(missing_files, show_n), collapse = "\n  "),
    "\n请检查 04_salmon_quant.py 的输出路径与 tximport_meta.tsv 是否一致。"
  )
}

log_msg("INFO", "共 ", length(files_vec), " 个样本将进入 tximport。")

#------------------------- 读取 tx2gene.clean.tsv -------------------------#

if (!file.exists(tx2gene_path)) {
  stop_with("找不到 tx2gene.clean.tsv：", tx2gene_path,
            " 请先运行 03_build_tx2gene_map.py。")
}

tx2gene <- data.table::fread(tx2gene_path)
# 只保留 transcript_id / gene_id 两列（按 vNext 约定）
expected_tx_cols <- c("transcript_id", "gene_id")
missing_tx_cols <- setdiff(expected_tx_cols, names(tx2gene))
if (length(missing_tx_cols) > 0L) {
  stop_with("tx2gene.clean.tsv 需要包含列：",
            paste(expected_tx_cols, collapse = ", "),
            "；缺少：", paste(missing_tx_cols, collapse = ", "))
}

tx2gene_use <- as.data.frame(tx2gene[, ..expected_tx_cols])

log_msg("INFO", "tx2gene.clean.tsv 行数 = ", nrow(tx2gene_use))

#------------------------- 调用 tximport 聚合 -------------------------#

counts_from_abundance <- get_or(cfg, "tximport", "counts_from_abundance",
                                default = "no")

log_msg("INFO", "tximport.counts_from_abundance = ", counts_from_abundance)

txi <- tximport::tximport(
  files   = files_vec,
  type    = "salmon",
  tx2gene = tx2gene_use,
  countsFromAbundance = counts_from_abundance
)

counts <- txi$counts
tpm    <- txi$abundance  # salmon + tximport：abundance 默认是 TPM

log_msg("INFO", "tximport 完成：",
        "基因数 = ", nrow(counts),
        "，样本数 = ", ncol(counts))

#------------------------- 写出 counts / TPM -------------------------#

write_tsv_matrix <- function(mat, path, id_col = "gene_id") {
  dt <- as.data.table(mat, keep.rownames = id_col)
  data.table::fwrite(dt, file = path, sep = "\t", quote = FALSE, na = "NA")
}

write_tsv_matrix(counts, counts_out_path)
write_tsv_matrix(tpm,    tpm_out_path)

log_msg("INFO", "已写出 counts：", counts_out_path)
log_msg("INFO", "已写出 TPM：",    tpm_out_path)

#------------------------- 生成 matrix_stats.tsv -------------------------#

# 库深（每个样本的总 counts）
lib_sizes <- colSums(counts, na.rm = TRUE)

# 零表达比例（每个样本零计数的基因比例）
zero_frac <- colMeans(counts == 0, na.rm = TRUE)

# 全局统计
n_samples <- ncol(counts)
n_genes   <- nrow(counts)

qfun <- function(x, probs) {
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

lib_q <- qfun(lib_sizes, c(0.25, 0.5, 0.75))
zero_median <- stats::median(zero_frac, na.rm = TRUE)

stats_dt <- data.table(
  n_samples            = n_samples,
  n_genes              = n_genes,
  zero_fraction_median = zero_median,
  libsize_q1           = lib_q[1],
  libsize_median       = lib_q[2],
  libsize_q3           = lib_q[3]
)

data.table::fwrite(stats_dt, file = matrix_stats_path,
                   sep = "\t", quote = FALSE, na = "NA")

log_msg("INFO", "已写出 matrix_stats：", matrix_stats_path)
log_msg("INFO", "===== 05 — 表达矩阵聚合完成 =====")