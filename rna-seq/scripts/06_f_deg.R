#!/usr/bin/env Rscript

# 06_f_deg.R —— DESeq2 实施（严格契约·只读 config.yaml）
# 约定要点：
#  1) 只读 config.yaml / 标准路径，不接受命令行/环境变量。
#  2) contrasts.tsv 列严格为：contrast, case, control。
#  3) 计数矩阵严格为：results/05_matrix/counts/gene_counts.tsv（无回退）。
#  4) 剔除 results/01_qc/rejects.tsv 中样本（若存在）。
#  5) 产物（每对比）固定 6 件套；DEG 三表列名统一 snake_case。
#  6) 批次：deg.allow_batch=true 且可估计 → ~ batch + group；否则回退 ~ group，并在 design.txt 写明原因。

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(matrixStats)
  library(readr)
})

CONFIG_PATH <- "config.yaml"

# ---------- 工具：递归合并默认值 ----------
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
    deg    = "results/06_deg",
    logs   = "logs"
  ),
  deg = list(
    lfc = 1.0,
    fdr = 0.05,
    use_apeglm = TRUE,
    independent_filter = TRUE,
    allow_batch = TRUE
  )
)

cfg <- read_yaml(CONFIG_PATH)
cfg <- merge_list(cfg, defaults)

# ---------- 路径 ----------
samples_fp   <- cfg$data$samples_tsv
contrasts_fp <- cfg$data$contrasts_tsv
counts_fp    <- file.path(cfg$dirs$matrix, "counts", "gene_counts.tsv")
qc_rejects   <- file.path(cfg$dirs$qc, "rejects.tsv")
deg_root     <- cfg$dirs$deg

# ---------- 必需文件检查 ----------
if (!file.exists(samples_fp))   stop(sprintf("未找到 samples.tsv：%s", samples_fp))
if (!file.exists(contrasts_fp)) stop(sprintf("未找到 contrasts.tsv：%s", contrasts_fp))
if (!file.exists(counts_fp))    stop(sprintf("未找到 gene_counts.tsv：%s（请先完成 05）", counts_fp))

# ---------- 读取数据 ----------
cat("========== 06-R — DESeq2 ==========\n")
cat(sprintf("[Info] samples   = %s\n", samples_fp))
cat(sprintf("[Info] contrasts = %s\n", contrasts_fp))
cat(sprintf("[Info] counts    = %s\n", counts_fp))
cat(sprintf("[Info] qc.rejects= %s\n", qc_rejects))

samples <- fread(samples_fp, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("sample","group") %in% colnames(samples))) {
  stop("samples.tsv 必含列：sample, group")
}
has_batch <- "batch" %in% colnames(samples)

contr <- fread(contrasts_fp, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("contrast","case","control") %in% colnames(contr))) {
  stop("contrasts.tsv 必含列：contrast, case, control")
}

cnt <- fread(counts_fp, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (!"gene_id" %in% colnames(cnt)) stop("gene_counts.tsv 必含列：gene_id")
rownames(cnt) <- cnt$gene_id
cnt$gene_id <- NULL

# 剔除 rejects.tsv 中样本（若存在）
rej <- character(0)
if (file.exists(qc_rejects)) {
  rej_tab <- tryCatch(fread(qc_rejects, sep = "\t", header = TRUE, data.table = FALSE),
                      error = function(e) NULL)
  if (!is.null(rej_tab) && "sample" %in% colnames(rej_tab)) {
    rej <- unique(as.character(rej_tab$sample))
  }
}
if (length(rej) > 0) {
  keep_samp <- setdiff(colnames(cnt), rej)
  cnt <- cnt[, keep_samp, drop = FALSE]
  samples <- subset(samples, !(sample %in% rej))
  cat(sprintf("[Info] 剔除低质样本（来自 rejects.tsv）：%s\n",
              paste(rej, collapse = ",")))
}

# 与计数矩阵对齐（交集）
common <- intersect(colnames(cnt), samples$sample)
if (length(common) < 2) stop("有效样本数不足：与计数矩阵的交集 < 2")
cnt <- cnt[, common, drop = FALSE]
samples <- samples[samples$sample %in% common, , drop = FALSE]
row.names(samples) <- samples$sample

# ---------- 配置参数 ----------
deg_lfc  <- as.numeric(cfg$deg$lfc)
deg_fdr  <- as.numeric(cfg$deg$fdr)
use_apeg <- isTRUE(cfg$deg$use_apeglm)
use_iflt <- isTRUE(cfg$deg$independent_filter)
allow_bt <- isTRUE(cfg$deg$allow_batch)

cat(sprintf("[Parm] lfc=%.3f  fdr=%.3f  use_apeglm=%s  independent_filter=%s  allow_batch=%s\n",
            deg_lfc, deg_fdr, as.character(use_apeg), as.character(use_iflt), as.character(allow_bt)))

# ---------- 每个对比逐一计算 ----------
for (i in seq_len(nrow(contr))) {
  lb   <- as.character(contr$contrast[i])
  case <- as.character(contr$case[i])
  ctrl <- as.character(contr$control[i])

  outdir <- file.path(deg_root, lb)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  cat(sprintf("\n[Run] contrast=%s  case=%s  control=%s\n", lb, case, ctrl))

  # 只保留 case/control 两组样本
  samp_sub <- subset(samples, group %in% c(case, ctrl))
  sids <- samp_sub$sample
  if (length(sids) < 2) {
    cat("[Warn] 对比有效样本数不足，跳过\n")
    next
  }
  cnt_sub <- cnt[, sids, drop = FALSE]

  # 因子水平（control 作为基线）
  samp_sub$group <- factor(samp_sub$group, levels = c(ctrl, case))
  if (has_batch) {
    # 若存在 batch 列，转为因子；否则后续根据 allow_batch 决策是否纳入
    samp_sub$batch <- factor(samp_sub$batch)
  }

  # 设计式与满秩检查
  reason <- "—"
  use_batch <- FALSE
  design_formula <- ~ group

  if (allow_bt && has_batch) {
    # 只有 batch 至少两个水平且每水平至少有样本时才考虑
    if (nlevels(droplevels(samp_sub$batch)) >= 2) {
      mm <- model.matrix(~ batch + group, data = samp_sub)
      full_rank <- qr(mm)$rank == ncol(mm)
      if (full_rank) {
        design_formula <- ~ batch + group
        use_batch <- TRUE
      } else {
        reason <- "设计矩阵不可估（秩亏），回退为 ~ group"
      }
    } else {
      reason <- "batch 单一或无效，回退为 ~ group"
    }
  } else if (!allow_bt && has_batch) {
    reason <- "按配置未纳入 batch（allow_batch=false）"
  }

  # 构造 DESeq2 对象并运行
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(cnt_sub)),
                                colData   = samp_sub,
                                design    = design_formula)
  # 低表达过滤（温和）：保留总计数>=10的基因
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  # 结果与 LFC 收缩
  res <- results(dds,
                 contrast = c("group", case, ctrl),
                 alpha = deg_fdr,
                 independentFiltering = use_iflt)

  if (use_apeg) {
    # 使用 apeglm 收缩
    resL <- lfcShrink(dds,
                      contrast = c("group", case, ctrl),
                      type = "apeglm",
                      res = res)
    # 用收缩后的 lfc 和 SE 回填
    res$log2FoldChange <- resL$log2FoldChange
    res$lfcSE <- resL$lfcSE
  }

  # 组装表格并统一列名（snake_case）
  dt <- as.data.frame(res)
  dt$gene_id <- rownames(dt)
  dt <- dt[, c("gene_id", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]
  colnames(dt) <- c("gene_id", "log2fc", "lfc_se", "stat", "p_value", "p_adjust", "base_mean")

  # 写 All
  all_fp <- file.path(outdir, "DEG_all.tsv")
  write_tsv(dt, all_fp)

  # 写 Up / Down（按阈值与 FDR）
  up_dt   <- subset(dt, !is.na(p_adjust) & p_adjust <= deg_fdr & log2fc >=  deg_lfc)
  down_dt <- subset(dt, !is.na(p_adjust) & p_adjust <= deg_fdr & log2fc <= -deg_lfc)
  write_tsv(up_dt,   file.path(outdir, "DEG_up.tsv"))
  write_tsv(down_dt, file.path(outdir, "DEG_down.tsv"))

  cat(sprintf("[Out] %s\n[Out] %s\n[Out] %s\n",
              all_fp,
              file.path(outdir, "DEG_up.tsv"),
              file.path(outdir, "DEG_down.tsv")))

  # 诊断：varTop100（基于归一化计数的方差）
  ncnt <- counts(dds, normalized = TRUE)
  vars <- rowVars(as.matrix(ncnt))
  ord  <- order(vars, decreasing = TRUE)
  top  <- rownames(ncnt)[head(ord, 100)]
  writeLines(top, con = file.path(outdir, "varTop100.list"))
  cat(sprintf("[Out] %s\n", file.path(outdir, "varTop100.list")))

  # 诊断：RLE 范围（log2(norm+1) 减基因中位数）
  logn <- log2(ncnt + 1)
  med  <- matrixStats::rowMedians(as.matrix(logn))
  rle  <- sweep(logn, 1, med, FUN = "-")
  rle_vals <- as.numeric(as.matrix(rle))
  rle_min <- min(rle_vals, na.rm = TRUE)
  rle_max <- max(rle_vals, na.rm = TRUE)
  rle_iqr <- IQR(rle_vals, na.rm = TRUE)
  rle_tab <- data.frame(rle_min = rle_min, rle_max = rle_max, rle_iqr = rle_iqr)
  write_tsv(rle_tab, file.path(outdir, "rle_range.tsv"))
  cat(sprintf("[Out] %s\n", file.path(outdir, "rle_range.tsv")))

  # design.txt（把口径与原因写清楚）
  des_lines <- c(
    sprintf("contrast: %s  (case=%s vs control=%s)", lb, case, ctrl),
    sprintf("design: %s", deparse(design_formula)),
    sprintf("batch_included: %s", ifelse(use_batch, "true", "false")),
    sprintf("reason_if_no_batch: %s", reason),
    sprintf("thresholds: |log2fc| >= %.3f  , FDR <= %.3f", deg_lfc, deg_fdr),
    sprintf("n_samples: %d  (case=%d, control=%d)",
            nrow(samp_sub), sum(samp_sub$group == case), sum(samp_sub$group == ctrl)),
    sprintf("samples: %s", paste(rownames(samp_sub), collapse = ",")),
    sprintf("rejected_samples_excluded: %s",
            ifelse(length(rej) == 0, "none", paste(rej, collapse = ",")))
  )
  writeLines(des_lines, con = file.path(outdir, "design.txt"))
  cat(sprintf("[Out] %s\n", file.path(outdir, "design.txt")))
}

cat("========== 06-R 完成 ==========\n")

