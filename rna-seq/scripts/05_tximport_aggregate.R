#!/usr/bin/env Rscript
# 05_tximport_aggregate.R — 汇总 Salmon 结果为基因层矩阵（无绘图，路径与告警修复）

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(tibble)
  library(yaml)
  library(stringr)
  library(rlang)   # 为 sym()/!! 选择列，避免 tidyselect 告警
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -------- 健壮版 getp --------
getp <- function(lst, ..., default=NULL) {
  ref <- lst
  keys <- list(...)
  for (k in keys) {
    if (is.null(ref)) return(default)
    if (!(is.list(ref) || is.environment(ref))) return(default)
    if (is.null(k) || !nzchar(as.character(k))) return(default)
    v <- ref[[as.character(k)]]
    if (is.null(v)) return(default)
    ref <- v
  }
  ref
}

normalize_path <- function(x, fallback) {
  if (is.null(x) || length(x) == 0) return(as.character(fallback))
  x <- as.character(x)[1]
  if (!nzchar(x)) return(as.character(fallback))
  x
}

# -------- 读取配置 --------
cfg_path <- "config.yaml"
if (!file.exists(cfg_path)) stop("缺少 config.yaml（请在项目根运行）")
cfg <- tryCatch(yaml::read_yaml(cfg_path), error=function(e) NULL)
if (is.null(cfg) || !(is.list(cfg) || is.environment(cfg))) {
  stop("config.yaml 解析失败：不是有效 YAML 映射")
}

quant_dir   <- normalize_path(getp(cfg, "paths","quant_dir"),   "results/quant")
matrix_dir  <- normalize_path(getp(cfg, "paths","matrix_dir"),  "results/matrix")
logs_dir    <- normalize_path(getp(cfg, "paths","logs_dir"),    "logs")
samples_tsv <- normalize_path(getp(cfg, "tables","samples"),    "data/samples.tsv")
tx2gene_fp  <- normalize_path(getp(cfg, "paths","tx2gene"),     "ref/tx2gene.geneMap.tsv")

message("[INFO] Resolved paths:")
message("  - quant_dir  = ", quant_dir)
message("  - matrix_dir = ", matrix_dir)
message("  - logs_dir   = ", logs_dir)
message("  - samples    = ", samples_tsv)
message("  - tx2gene    = ", tx2gene_fp)

# -------- 创建输出目录 --------
ok_m <- try(dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
if (inherits(ok_m, "try-error")) {
  stop("[ERR] 无法创建 matrix_dir：", matrix_dir)
}
ok_meta <- try(dir.create(file.path(matrix_dir, "meta"), recursive = TRUE, showWarnings = FALSE), silent = TRUE)
if (inherits(ok_meta, "try-error")) {
  stop("[ERR] 无法创建 meta 目录：", file.path(matrix_dir, "meta"))
}

# -------- 基础健诊 --------
missing <- c()
if (!dir.exists(quant_dir))    missing <- c(missing, paste0("paths.quant_dir 不存在: ", quant_dir))
if (!file.exists(samples_tsv)) missing <- c(missing, paste0("tables.samples 缺失: ", samples_tsv))
if (!file.exists(tx2gene_fp))  missing <- c(missing, paste0("paths.tx2gene 缺失: ", tx2gene_fp))
if (length(missing) > 0) {
  stop("[ERR] 配置或输入缺失：\n - ", paste(missing, collapse = "\n - "),
       "\n请检查 config.yaml 与相对路径（以项目根为基准）。")
}

# -------- 样本与 quant.sf 列表 --------
samples <- readr::read_tsv(samples_tsv, show_col_types = FALSE)
if (!nrow(samples)) stop("[ERR] samples.tsv 为空：", samples_tsv)
names_lower <- setNames(names(samples), tolower(names(samples)))
col_sample <- names_lower[["sample"]] %||% names_lower[["samples"]] %||% names_lower[["id"]] %||% names_lower[["sid"]]
if (is.null(col_sample)) stop("[ERR] samples.tsv 需要包含列：sample")
sample_ids <- samples[[col_sample]] %>% as.character() %>% trimws()
sample_ids <- sample_ids[nzchar(sample_ids)]
if (!length(sample_ids)) stop("[ERR] 未在样本表中解析到有效样本名")

files <- file.path(quant_dir, sample_ids, "quant.sf")
names(files) <- sample_ids
missing_sf <- files[!file.exists(files)]
if (length(missing_sf) > 0) {
  msg <- paste0("以下样本缺少 quant.sf：\n", paste0(" - ", names(missing_sf), ": ", missing_sf, collapse="\n"))
  stop(msg)
}

# -------- 读取 tx2gene（03 清洗版）--------
tx2gene <- readr::read_tsv(tx2gene_fp, show_col_types = FALSE)
if (!nrow(tx2gene)) stop("[ERR] tx2gene 为空：", tx2gene_fp)
names_lower2 <- tolower(names(tx2gene))
if (!all(c("tx","gene") %in% names_lower2)) {
  stop("[ERR] tx2gene 应包含列名 TX 与 GENE（不区分大小写）。实际列：", paste(names(tx2gene), collapse=", "))
}
col_tx   <- names(tx2gene)[match("tx",   names_lower2)]
col_gene <- names(tx2gene)[match("gene", names_lower2)]
tx2gene <- tx2gene %>%
  dplyr::select(TX = !!sym(col_tx), GENE = !!sym(col_gene)) %>%  # 改为 rlang::sym，消除 tidyselect 告警
  distinct()

# -------- 聚合（tximport）--------
message("reading in files with read_tsv")
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE,
                ignoreAfterBar  = TRUE,
                countsFromAbundance = "lengthScaledTPM")

counts <- as.data.frame(txi$counts)    %>% rownames_to_column("gene_id")
abund  <- as.data.frame(txi$abundance) %>% rownames_to_column("gene_id")

# -------- 写出矩阵与元信息 --------
out_counts <- file.path(matrix_dir, "gene_counts.tsv")
out_tpm    <- file.path(matrix_dir, "gene_tpm.tsv")
readr::write_tsv(counts, out_counts)
readr::write_tsv(abund,  out_tpm)

sample_info <- tibble(sample = sample_ids, file = unname(files))
readr::write_tsv(sample_info, file.path(matrix_dir,"meta","sample_info.tsv"))

# 修正：把可表达基因写入 meta/ 目录（与 Python 外壳期望一致）
expr_genes <- abund %>%
  mutate(any_tpm = if_any(-gene_id, ~ .x > 0)) %>%
  filter(any_tpm) %>%
  select(gene_id)
readr::write_tsv(expr_genes, file.path(matrix_dir,"meta","expressed_genes.tsv"))

message("汇总完成 ✅")

