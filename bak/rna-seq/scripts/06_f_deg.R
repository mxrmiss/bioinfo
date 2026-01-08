#!/usr/bin/env Rscript
# 06_f_deg.R —— 差异表达主程（严格契约）
# 只读 config.yaml；不接命令行参数；输出契约文件

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(matrixStats)
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
  data = list(
    samples_tsv   = "data/samples.tsv",
    contrasts_tsv = "data/contrasts.tsv"
  ),
  dirs = list(
    qc     = "results/01_qc",
    matrix = "results/05_matrix",
    deg    = "results/06_deg"
  ),
  deg = list(
    lfc = 1.0,
    fdr = 0.05,
    use_apeglm = TRUE,
    independent_filter = TRUE,
    allow_batch = TRUE
  )
)

cfg <- read_yaml(CONFIG_PATH); cfg <- merge_list(cfg, defaults)

samples_fp   <- cfg$data$samples_tsv
contrasts_fp <- cfg$data$contrasts_tsv
matrix_dir   <- cfg$dirs$matrix
deg_root     <- cfg$dirs$deg
qc_dir       <- cfg$dirs$qc

counts_fp <- file.path(matrix_dir, "counts", "gene_counts.tsv")
tpms_fp   <- file.path(matrix_dir, "tpms",   "gene_tpm.tsv")

cat("========== 06-R — DEG ==========\n")
cat(sprintf("[Info] samples   = %s\n", samples_fp))
cat(sprintf("[Info] contrasts = %s\n", contrasts_fp))
cat(sprintf("[Info] counts    = %s\n", counts_fp))
cat(sprintf("[Info] tpms      = %s\n", tpms_fp))

# 读入
samples <- fread(samples_fp, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("sample","group") %in% colnames(samples))) {
  stop("samples.tsv 必须包含列：sample, group")
}
contrasts <- fread(contrasts_fp, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("contrast","case","control") %in% colnames(contrasts))) {
  stop("contrasts.tsv 必须包含列：contrast, case, control")
}
counts <- fread(counts_fp, sep = "\t", header = TRUE, data.table = FALSE)
tpms   <- fread(tpms_fp,   sep = "\t", header = TRUE, data.table = FALSE)

rownames_assign <- function(df, idcol) {
  rn <- df[[idcol]]; df[[idcol]] <- NULL; rn <- make.unique(as.character(rn))
  rownames(df) <- rn; df
}
counts <- rownames_assign(counts, "gene_id")
tpms   <- rownames_assign(tpms,   "gene_id")

# 剔除 rejects.tsv 样本
rej_fp <- file.path(qc_dir, "rejects.tsv")
if (file.exists(rej_fp)) {
  rej <- fread(rej_fp, sep = "\t", header = TRUE, data.table = FALSE)
  if ("sample" %in% colnames(rej)) {
    samples <- subset(samples, !(samples$sample %in% rej$sample))
    cat(sprintf("[Info] 剔除 rejects 样本：%d 个\n", nrow(rej)))
  }
}

# DE 主循环
for (i in seq_len(nrow(contrasts))) {
  label   <- contrasts$contrast[i]
  grp1    <- contrasts$case[i]
  grp0    <- contrasts$control[i]

  cat(sprintf("\n[Run] %s : %s vs %s\n", label, grp1, grp0))

  keep_samp <- samples[samples$group %in% c(grp1, grp0), , drop = FALSE]
  if (nrow(keep_samp) < 2) {
    warning(sprintf("对比 %s 有效样本数不足；跳过", label)); next
  }

  # 提取计数矩阵
  common_cols <- intersect(colnames(counts), keep_samp$sample)
  if (length(common_cols) < 2) {
    warning(sprintf("对比 %s 可用列不足；跳过", label)); next
  }
  cnt <- counts[, common_cols, drop = FALSE]
  coldata <- data.frame(row.names = common_cols,
                        group = factor(keep_samp$group[match(common_cols, keep_samp$sample)],
                                       levels = c(grp0, grp1)))

  # 批次可选
  use_batch <- FALSE
  if (isTRUE(cfg$deg$allow_batch) && "batch" %in% colnames(samples)) {
    bvals <- samples$batch[match(common_cols, samples$sample)]
    if (length(unique(bvals[!is.na(bvals)])) >= 2) {
      coldata$batch <- factor(bvals)
      use_batch <- TRUE
    }
  }

  # 设计式与记录
  design_txt <- file.path(deg_root, label, "design.txt")
  dir.create(file.path(deg_root, label), showWarnings = FALSE, recursive = TRUE)
  if (use_batch) {
    design_formula <- ~ batch + group
    writeLines(c(
      sprintf("design : ~ batch + group"),
      sprintf("case   : %s", grp1),
      sprintf("control: %s", grp0),
      sprintf("allow_batch: TRUE (本对比已纳入 batch)")
    ), con = design_txt)
  } else {
    design_formula <- ~ group
    writeLines(c(
      sprintf("design : ~ group"),
      sprintf("case   : %s", grp1),
      sprintf("control: %s", grp0),
      sprintf("allow_batch: %s（未纳入）", as.character(cfg$deg$allow_batch))
    ), con = design_txt)
  }

  # 构建 DESeq2
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(cnt)),
                                colData   = coldata,
                                design    = design_formula)

  # 运行 DESeq2
  dds <- DESeq(dds)

  # 结果与收缩
  res <- results(dds, contrast = c("group", grp1, grp0),
                 alpha = cfg$deg$fdr, independentFiltering = cfg$deg$independent_filter)
  if (isTRUE(cfg$deg$use_apeglm)) {
    resLFC <- lfcShrink(dds, coef = paste0("group_", grp1, "_vs_", grp0),
                        type = "apeglm")
    res$log2FoldChange <- resLFC$log2FoldChange
    res$lfcSE <- resLFC$lfcSE
  }

  # 整理表头（契约）
  out <- data.frame(
    gene_id   = rownames(res),
    log2fc    = as.numeric(res$log2FoldChange),
    lfc_se    = as.numeric(res$lfcSE),
    stat      = as.numeric(res$stat),
    p_value   = as.numeric(res$pvalue),
    p_adjust  = as.numeric(res$padj),
    base_mean = as.numeric(res$baseMean),
    stringsAsFactors = FALSE
  )
  # 三套表
  deg_dir <- file.path(deg_root, label)
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(out, file.path(deg_dir, "DEG_all.tsv"))
  up   <- subset(out,  !is.na(p_adjust) & p_adjust <= cfg$deg$fdr & log2fc >= cfg$deg$lfc)
  down <- subset(out, !is.na(p_adjust) & p_adjust <= cfg$deg$fdr & log2fc <= -cfg$deg$lfc)
  write_tsv(up,   file.path(deg_dir, "DEG_up.tsv"))
  write_tsv(down, file.path(deg_dir, "DEG_down.tsv"))

  # 诊断：Top-var 与 RLE 范围
  vsd <- vst(dds, blind = TRUE)
  vmat <- assay(vsd)
  gvar <- rowVars(vmat)
  top_idx <- order(gvar, decreasing = TRUE)[seq_len(min(100, length(gvar)))]
  writeLines(rownames(vmat)[top_idx], con = file.path(deg_dir, "varTop100.list"))

  # RLE 范围
  med <- rowMedians(vmat)
  rle <- sweep(vmat, 1, med, FUN = "-")
  rle_range <- data.frame(
    min = min(rle, na.rm = TRUE),
    max = max(rle, na.rm = TRUE),
    iqr = IQR(as.numeric(rle), na.rm = TRUE)
  )
  write_tsv(rle_range, file.path(deg_dir, "rle_range.tsv"))

  cat(sprintf("[Out] %s/{DEG_all.tsv, DEG_up.tsv, DEG_down.tsv, varTop100.list, rle_range.tsv}\n",
              deg_dir))
}

cat("========== 06-R 完成 ==========\n")
