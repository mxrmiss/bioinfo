#!/usr/bin/env Rscript
# =============================================================================
# 10_wgcna_modules.R —— RNA-seq WGCNA 模块分析（计算版·不画图·顶刊稳妥口径）
#
# 核心特性：
# - 只做“计算 + 导出原材料”，不生成任何图片（避免字体/设备坑）
# - config.yaml 参数优先级最高；脚本顶部仅提供默认值
# - hub genes 计算不使用 signedKME（避免某些环境下 bicor 内部 eval(use=...) 报错）
# - 线程数可在 config.yaml 控制（resources.threads.wgcna > wgcna.threads > 默认）
#
# 2026-01-09 增量改进（不改 config.yaml）：
# - module_trait_cor/p/padj 的 TSV 导出：第一列补上模块/ME 名（module），便于检查和作图
# - 同时额外输出旧格式“纯数值矩阵”备份：*.matrix_only.tsv，避免下游兼容性问题
# - 相关性计算前强制对齐 traits_mat 与 MEs 的样本顺序，防止隐式错配
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# 0) 脚本顶部默认参数（可改；但 config.yaml 若提供则覆盖）
# =============================================================================

CONFIG_YAML <- "config.yaml"

WGCNA_DEFAULT <- list(
  enable = TRUE,
  method = "vst",  # "vst" / "voom"

  input_counts = "results/05_matrix/counts/gene_counts.tsv",
  samples_tsv  = "data/samples.tsv",

  outdir = "results/10_wgcna",

  save_rds = TRUE,
  random_seed = 20260108,

  # 线程控制（如果 config.yaml 未提供 resources.threads.wgcna / wgcna.threads，则使用此默认值）
  threads = 8,

  filter = list(
    min_count = 10,
    min_samples = 6,        # ✅ 修复点：这里必须是正常整数
    var_metric = "mad",      # "mad" / "variance"
    keep_top_fraction = 0.5
    # keep_top_n: 8000  # 可选：若设置则优先于 keep_top_fraction
  ),

  vst = list(
    blind = TRUE
  ),
  voom = list(
    enable_tmm = TRUE
  ),

  network = list(
    networkType = "signed",
    TOMType = "signed",
    corType = "bicor",                 # "bicor" / "pearson"
    maxPOutliers = 0.05,               # bicor 鲁棒参数
    useObs = "pairwise.complete.obs"   # 只能是 "all.obs" 或 "pairwise.complete.obs"
  ),

  power = list(
    candidates = c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20),
    scale_free_R2_target = 0.85,
    min_mean_connectivity = 20
    # fixed: 12  # 可选：固定 power，跳过自动选择
  ),

  modules = list(
    minModuleSize = 30,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25,
    reassignThreshold = 0,
    numericLabels = TRUE
  ),

  traits = list(
    sample_column = "sample",
    group_column  = "group",
    build_onehot = TRUE,
    binary_traits = list(
      list(
        name = "target_vs_rest",
        positive_groups = c("foot"),
        negative_groups = c()
      )
    ),
    extra_columns_as_traits = c()
  ),

  export = list(
    gene_module_table = TRUE,
    module_lists = TRUE,
    hub_genes = list(
      enable = TRUE,
      top_n_per_module = 50,
      rank_by = c("kME", "GS"),
      primary_trait = ""   # 建议在 config.yaml 显式指定，例如 "target_vs_rest"
    )
  )
)

# =============================================================================
# 1) 小工具：日志、目录、YAML、递归合并
# =============================================================================

safe_mkdir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_init <- function(log_file) {
  safe_mkdir(dirname(log_file))
  con <- file(log_file, open = "wt", encoding = "UTF-8")
  list(con = con, file = log_file)
}

log_close <- function(logger) {
  try(close(logger$con), silent = TRUE)
}

log_msg <- function(logger, level = "INFO", msg = "") {
  line <- sprintf("[%s] [%s] %s", now_str(), level, msg)
  cat(line, "\n")
  if (!is.null(logger) && !is.null(logger$con)) {
    writeLines(line, logger$con)
    flush(logger$con)
  }
}

need_pkg <- function(pkgs, logger = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    log_msg(logger, "ERROR",
            paste0("缺少 R 包：", paste(missing, collapse = ", "),
                   "。请先安装后再运行。"))
    stop("Missing packages: ", paste(missing, collapse = ", "))
  }
}

read_yaml_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read config.yaml.")
  }
  yaml::read_yaml(path)
}

write_yaml_safe <- function(x, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to write yaml output.")
  }
  safe_mkdir(dirname(path))
  yaml::write_yaml(x, path)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

deep_merge <- function(base, override) {
  if (is.null(override)) return(base)
  if (!is.list(base) || !is.list(override)) return(override)
  out <- base
  for (nm in names(override)) {
    if (!nm %in% names(out)) {
      out[[nm]] <- override[[nm]]
    } else {
      out[[nm]] <- deep_merge(out[[nm]], override[[nm]])
    }
  }
  out
}

# -----------------------------------------------------------------------------
# ✅ 增量工具：把矩阵写成“首列含 rowname”的 TSV，同时可选写一份旧格式 matrix-only
# -----------------------------------------------------------------------------
write_matrix_with_rownames <- function(mat, out_tsv, rowname_col = "module", out_matrix_only = NULL) {
  # mat: matrix with rownames and colnames
  if (is.null(rownames(mat))) {
    stop("write_matrix_with_rownames: matrix has no rownames, cannot export module names.")
  }
  df <- data.frame(
    tmp_rowname = rownames(mat),
    as.data.frame(mat, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(df)[1] <- rowname_col
  data.table::fwrite(df, out_tsv, sep = "\t", quote = FALSE)

  if (!is.null(out_matrix_only)) {
    # 旧格式：仅数值矩阵（无行名列），用于兼容历史下游脚本
    data.table::fwrite(as.data.frame(mat, check.names = FALSE), out_matrix_only, sep = "\t", quote = FALSE)
  }
}

# =============================================================================
# 2) 读取 config.yaml，并应用“config 优先级最高”的覆盖策略
# =============================================================================

cfg_all <- read_yaml_safe(CONFIG_YAML)

cfg_dirs_wgcna <- NULL
if (!is.null(cfg_all) && !is.null(cfg_all$dirs) && !is.null(cfg_all$dirs$wgcna)) {
  cfg_dirs_wgcna <- cfg_all$dirs$wgcna
}

cfg_wgcna <- NULL
if (!is.null(cfg_all) && !is.null(cfg_all$wgcna)) cfg_wgcna <- cfg_all$wgcna

W <- deep_merge(WGCNA_DEFAULT, cfg_wgcna)

# outdir 优先级：wgcna.outdir > dirs.wgcna > 默认
if (!is.null(cfg_wgcna) && !is.null(cfg_wgcna$outdir) && nzchar(cfg_wgcna$outdir)) {
  W$outdir <- cfg_wgcna$outdir
} else if (!is.null(cfg_dirs_wgcna) && nzchar(cfg_dirs_wgcna)) {
  W$outdir <- cfg_dirs_wgcna
}

# samples_tsv：若 wgcna.samples_tsv 未设，则继承 data.samples_tsv
if (!is.null(cfg_all) && !is.null(cfg_all$data) && !is.null(cfg_all$data$samples_tsv)) {
  if (is.null(cfg_wgcna) || is.null(cfg_wgcna$samples_tsv)) {
    W$samples_tsv <- cfg_all$data$samples_tsv
  }
}

# 线程优先级：resources.threads.wgcna > wgcna.threads > 默认
cfg_threads_wgcna <- NULL
if (!is.null(cfg_all) && !is.null(cfg_all$resources) && !is.null(cfg_all$resources$threads) &&
    !is.null(cfg_all$resources$threads$wgcna)) {
  cfg_threads_wgcna <- cfg_all$resources$threads$wgcna
}
if (!is.null(cfg_threads_wgcna)) {
  W$threads <- as.integer(cfg_threads_wgcna)
} else if (!is.null(cfg_wgcna) && !is.null(cfg_wgcna$threads)) {
  W$threads <- as.integer(cfg_wgcna$threads)
}
if (is.na(W$threads) || W$threads < 1) W$threads <- 1L

# =============================================================================
# 3) 初始化输出目录与日志
# =============================================================================

OUTDIR <- W$outdir
safe_mkdir(OUTDIR)

SUBDIRS <- list(
  logs    = file.path(OUTDIR, "logs"),
  traits  = file.path(OUTDIR, "traits"),
  qc_raw  = file.path(OUTDIR, "QC_raw"),
  power   = file.path(OUTDIR, "power"),
  modules = file.path(OUTDIR, "modules"),
  assoc   = file.path(OUTDIR, "assoc"),
  lists   = file.path(OUTDIR, "lists"),
  hub     = file.path(OUTDIR, "hub"),
  rds     = file.path(OUTDIR, "rds")
)
for (d in SUBDIRS) safe_mkdir(d)

logger <- log_init(file.path(SUBDIRS$logs, "10_wgcna_modules.run.log"))
on.exit(log_close(logger), add = TRUE)

log_msg(logger, "INFO", paste0("Config path: ", CONFIG_YAML))
log_msg(logger, "INFO", paste0("Outdir: ", OUTDIR))

if (!isTRUE(W$enable)) {
  log_msg(logger, "WARN", "wgcna.enable=false，脚本退出（未执行任何分析）。")
  quit(save = "no", status = 0)
}

# =============================================================================
# 4) 依赖包加载 + 线程限制（WGCNA + BLAS/OMP）
# =============================================================================

need_pkg(c("data.table", "yaml", "WGCNA", "DESeq2", "SummarizedExperiment", "parallel"), logger)
if (tolower(W$method) == "voom") {
  need_pkg(c("limma", "edgeR"), logger)
}

Sys.setenv(
  OMP_NUM_THREADS = as.character(W$threads),
  OPENBLAS_NUM_THREADS = as.character(W$threads),
  MKL_NUM_THREADS = as.character(W$threads),
  VECLIB_MAXIMUM_THREADS = as.character(W$threads),
  NUMEXPR_NUM_THREADS = as.character(W$threads)
)

suppressPackageStartupMessages({
  library(WGCNA)
})

WGCNA::allowWGCNAThreads(nThreads = W$threads)
log_msg(logger, "INFO", paste0("Threads (from config): ", W$threads))

set.seed(as.integer(W$random_seed))

# =============================================================================
# 5) 读入 samples.tsv 并构建 traits（完全通用）
# =============================================================================

log_msg(logger, "INFO", paste0("Reading samples: ", W$samples_tsv))
if (!file.exists(W$samples_tsv)) {
  log_msg(logger, "ERROR", paste0("samples.tsv 不存在：", W$samples_tsv))
  stop("Missing samples.tsv")
}
samples_dt <- data.table::fread(W$samples_tsv, sep = "\t", header = TRUE, data.table = FALSE)

sc <- W$traits$sample_column %||% "sample"
gc <- W$traits$group_column  %||% "group"
if (!sc %in% colnames(samples_dt)) stop("samples.tsv missing sample_column: ", sc)
if (!gc %in% colnames(samples_dt)) stop("samples.tsv missing group_column: ", gc)

samples_dt[[sc]] <- as.character(samples_dt[[sc]])
samples_dt[[gc]] <- as.character(samples_dt[[gc]])

trait_list <- list()

if (isTRUE(W$traits$build_onehot)) {
  groups <- sort(unique(samples_dt[[gc]]))
  onehot <- sapply(groups, function(g) as.integer(samples_dt[[gc]] == g))
  colnames(onehot) <- paste0("group_", groups)
  trait_list[["onehot_group"]] <- onehot
}

binary_traits <- W$traits$binary_traits %||% list()
all_groups <- sort(unique(samples_dt[[gc]]))

for (bt in binary_traits) {
  if (is.null(bt$name) || !nzchar(bt$name)) stop("binary_traits: each trait must have a non-empty name")
  pos <- as.character(bt$positive_groups %||% character())
  neg <- as.character(bt$negative_groups %||% character())

  miss_pos <- setdiff(pos, all_groups)
  miss_neg <- setdiff(neg, all_groups)
  if (length(miss_pos) > 0) stop("binary_traits: positive_groups not found: ", paste(miss_pos, collapse = ", "))
  if (length(miss_neg) > 0) stop("binary_traits: negative_groups not found: ", paste(miss_neg, collapse = ", "))

  if (length(neg) == 0) neg <- setdiff(all_groups, pos)
  inter <- intersect(pos, neg)
  if (length(inter) > 0) stop("binary_traits: positive/negative overlap: ", paste(inter, collapse = ", "))

  vec <- rep(NA_integer_, nrow(samples_dt))
  vec[samples_dt[[gc]] %in% pos] <- 1L
  vec[samples_dt[[gc]] %in% neg] <- 0L
  if (any(is.na(vec))) stop("binary_traits produced NA values; please check group definitions.")

  trait_list[[paste0("binary__", bt$name)]] <- matrix(vec, ncol = 1, dimnames = list(NULL, bt$name))
}

extra_cols <- W$traits$extra_columns_as_traits %||% character()
if (length(extra_cols) > 0) {
  for (cn in extra_cols) {
    if (!cn %in% colnames(samples_dt)) stop("extra_columns_as_traits not found: ", cn)
    v <- samples_dt[[cn]]
    if (is.numeric(v)) {
      trait_list[[paste0("extra__", cn)]] <- matrix(v, ncol = 1, dimnames = list(NULL, cn))
    } else {
      v <- as.character(v)
      lv <- sort(unique(v))
      oh <- sapply(lv, function(z) as.integer(v == z))
      colnames(oh) <- paste0(cn, "_", lv)
      trait_list[[paste0("extra_onehot__", cn)]] <- oh
    }
  }
}

traits_mat <- NULL
if (length(trait_list) > 0) {
  traits_mat <- do.call(cbind, trait_list)
  traits_mat <- as.data.frame(traits_mat)
  rownames(traits_mat) <- samples_dt[[sc]]
  data.table::fwrite(
    data.frame(sample = rownames(traits_mat), traits_mat, check.names = FALSE),
    file.path(SUBDIRS$traits, "trait_matrix.tsv"),
    sep = "\t",
    quote = FALSE
  )
  log_msg(logger, "INFO", paste0("Traits built: ", ncol(traits_mat), " columns"))
} else {
  log_msg(logger, "WARN", "未构建 traits：后续将跳过模块-性状关联与 hub genes。")
}

# =============================================================================
# 6) 读入 gene_counts.tsv，匹配样本，过滤低表达
# =============================================================================

log_msg(logger, "INFO", paste0("Reading counts: ", W$input_counts))
if (!file.exists(W$input_counts)) {
  log_msg(logger, "ERROR", paste0("counts 文件不存在：", W$input_counts))
  stop("Missing counts file")
}

counts_dt <- data.table::fread(W$input_counts, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (ncol(counts_dt) < 3) stop("Counts table too small (need gene_id + >=2 samples).")

gene_col <- colnames(counts_dt)[1]
gene_ids <- as.character(counts_dt[[gene_col]])
counts_mat <- as.matrix(counts_dt[, -1, drop = FALSE])
rownames(counts_mat) <- gene_ids
suppressWarnings(storage.mode(counts_mat) <- "numeric")
if (anyNA(counts_mat)) {
  log_msg(logger, "WARN", "counts 矩阵含 NA，已将 NA 视作 0。")
  counts_mat[is.na(counts_mat)] <- 0
}

sample_names_counts <- colnames(counts_mat)
sample_names_meta   <- samples_dt[[sc]]
common_samples <- intersect(sample_names_meta, sample_names_counts)
if (length(common_samples) < 6) stop("Too few matched samples between counts and samples.tsv.")

common_samples <- sample_names_meta[sample_names_meta %in% common_samples]
counts_mat <- counts_mat[, common_samples, drop = FALSE]
samples_dt <- samples_dt[match(common_samples, samples_dt[[sc]]), , drop = FALSE]
if (!is.null(traits_mat)) traits_mat <- traits_mat[common_samples, , drop = FALSE]

log_msg(logger, "INFO", paste0("Matched samples: ", length(common_samples)))

min_count   <- as.numeric(W$filter$min_count %||% 10)
min_samples <- as.integer(W$filter$min_samples %||% max(1, floor(length(common_samples) / 3)))

keep_expr <- rowSums(counts_mat >= min_count) >= min_samples
counts_f1 <- counts_mat[keep_expr, , drop = FALSE]

log_msg(logger, "INFO", sprintf("Low-expression filter: kept %d / %d genes (min_count=%s in >=%s samples)",
                                nrow(counts_f1), nrow(counts_mat), min_count, min_samples))

# =============================================================================
# 7) 变换：VST（推荐）或 voom
# =============================================================================

method <- tolower(W$method %||% "vst")

datExpr <- NULL   # samples x genes
expr_mat <- NULL  # genes x samples

if (method == "vst") {
  log_msg(logger, "INFO", "Transform method: DESeq2 VST")

  colData <- data.frame(row.names = common_samples)
  colData[[gc]] <- factor(samples_dt[[gc]])

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_f1),
    colData   = colData,
    design    = ~ 1
  )
  dds <- DESeq2::estimateSizeFactors(dds)

  blind <- isTRUE(W$vst$blind %||% TRUE)
  vsd <- DESeq2::vst(dds, blind = blind)

  expr_mat <- SummarizedExperiment::assay(vsd)
  datExpr  <- t(expr_mat)

} else if (method == "voom") {
  log_msg(logger, "INFO", "Transform method: limma-voom")

  dge <- edgeR::DGEList(counts = counts_f1)
  if (isTRUE(W$voom$enable_tmm %||% TRUE)) {
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
  }

  design <- matrix(1, nrow = length(common_samples), ncol = 1)
  v <- limma::voom(dge, design = design, plot = FALSE)
  expr_mat <- v$E
  datExpr  <- t(expr_mat)

} else {
  stop("Unknown wgcna.method: ", method, " (use 'vst' or 'voom')")
}

data.table::fwrite(
  data.frame(sample = rownames(datExpr), datExpr, check.names = FALSE),
  file.path(SUBDIRS$modules, "datExpr.transformed.samples_x_genes.tsv"),
  sep = "\t",
  quote = FALSE
)

# =============================================================================
# 8) 二次过滤：低方差（在变换后矩阵上做）
# =============================================================================

var_metric <- tolower(W$filter$var_metric %||% "mad")
keep_top_n <- W$filter$keep_top_n %||% NULL
keep_frac  <- as.numeric(W$filter$keep_top_fraction %||% 0.5)

gene_stats <- NULL
if (var_metric == "mad") {
  gene_stats <- apply(datExpr, 2, mad, na.rm = TRUE)
} else if (var_metric %in% c("variance", "var")) {
  gene_stats <- apply(datExpr, 2, var, na.rm = TRUE)
} else {
  stop("Unknown filter.var_metric: ", var_metric)
}

ord <- order(gene_stats, decreasing = TRUE, na.last = NA)
if (!is.null(keep_top_n)) {
  keep_top_n <- as.integer(keep_top_n)
  k <- min(length(ord), keep_top_n)
} else {
  keep_frac <- max(min(keep_frac, 1), 0.01)
  k <- max(1, floor(length(ord) * keep_frac))
}
keep_genes <- names(gene_stats)[ord[seq_len(k)]]
datExpr <- datExpr[, keep_genes, drop = FALSE]

log_msg(logger, "INFO", sprintf("Low-variance filter (%s): kept %d genes (top_n=%s, top_fraction=%s)",
                                var_metric, ncol(datExpr),
                                ifelse(is.null(keep_top_n), "NA", keep_top_n),
                                ifelse(is.null(keep_top_n), keep_frac, "NA")))

# =============================================================================
# 9) WGCNA QC：好样本好基因 + 保存样本聚类树原材料（不画图）
# =============================================================================

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  badS <- sum(!gsg$goodSamples)
  badG <- sum(!gsg$goodGenes)
  log_msg(logger, "WARN", sprintf("goodSamplesGenes: badSamples=%d, badGenes=%d，将自动剔除。", badS, badG))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  keep_samp <- rownames(datExpr)
  samples_dt <- samples_dt[match(keep_samp, samples_dt[[sc]]), , drop = FALSE]
  if (!is.null(traits_mat)) traits_mat <- traits_mat[keep_samp, , drop = FALSE]
}

sampleTree <- hclust(dist(datExpr), method = "average")
saveRDS(sampleTree, file.path(SUBDIRS$qc_raw, "sampleTree.rds"))

# =============================================================================
# 10) 选择 soft-threshold power（保存 fitIndices 原材料）
# =============================================================================

corType <- tolower(W$network$corType %||% "bicor")
useObs  <- as.character(W$network$useObs %||% "pairwise.complete.obs")
maxPOut <- as.numeric(W$network$maxPOutliers %||% 0.05)
netType <- W$network$networkType %||% "signed"

if (!useObs %in% c("all.obs", "pairwise.complete.obs")) {
  stop("wgcna.network.useObs must be 'all.obs' or 'pairwise.complete.obs', got: ", useObs)
}

power_fixed <- W$power$fixed %||% NULL
powers <- as.numeric(W$power$candidates %||% c(1:10, 12, 14, 16, 18, 20))

sft <- NULL
chosen_power <- NULL

if (!is.null(power_fixed)) {
  chosen_power <- as.numeric(power_fixed)
  log_msg(logger, "INFO", paste0("Power fixed by config: ", chosen_power))
} else {
  log_msg(logger, "INFO", "Selecting soft-threshold power by pickSoftThreshold...")

  sft <- WGCNA::pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = netType,
    corFnc = ifelse(corType == "bicor", "bicor", "cor"),
    corOptions = if (corType == "bicor") list(use = useObs, maxPOutliers = maxPOut) else list(use = useObs),
    verbose = 5
  )

  fit <- sft$fitIndices
  data.table::fwrite(fit, file.path(SUBDIRS$power, "pickSoftThreshold.fitIndices.tsv"), sep = "\t", quote = FALSE)
  saveRDS(sft, file.path(SUBDIRS$power, "pickSoftThreshold.sft.rds"))

  r2_target <- as.numeric(W$power$scale_free_R2_target %||% 0.85)
  min_k     <- as.numeric(W$power$min_mean_connectivity %||% 20)

  if (!all(c("Power", "SFT.R.sq", "mean.k.") %in% colnames(fit))) {
    log_msg(logger, "WARN", "pickSoftThreshold 输出列名与预期不一致，将使用最大 R2 选择 power。")
    chosen_power <- fit$Power[which.max(fit$SFT.R.sq)]
  } else {
    ok <- which(fit$SFT.R.sq >= r2_target & fit$mean.k. >= min_k)
    if (length(ok) > 0) {
      chosen_power <- fit$Power[min(ok)]
    } else {
      log_msg(logger, "WARN", "未找到同时满足 R2 与 mean connectivity 的 power，将采用折中选择。")
      top <- order(fit$SFT.R.sq, decreasing = TRUE)[1:min(5, nrow(fit))]
      chosen_power <- fit$Power[top[which.max(fit$mean.k.[top])]]
    }
  }
}

log_msg(logger, "INFO", paste0("Chosen power: ", chosen_power))

# =============================================================================
# 11) 构建模块：blockwiseModules
# =============================================================================

log_msg(logger, "INFO", "Running WGCNA blockwiseModules...")

TOMType <- W$network$TOMType %||% "signed"
minModuleSize <- as.integer(W$modules$minModuleSize %||% 30)
deepSplit <- as.integer(W$modules$deepSplit %||% 2)
pamRespectsDendro <- isTRUE(W$modules$pamRespectsDendro %||% FALSE)
mergeCutHeight <- as.numeric(W$modules$mergeCutHeight %||% 0.25)
reassignThreshold <- as.numeric(W$modules$reassignThreshold %||% 0)
numericLabels <- isTRUE(W$modules$numericLabels %||% TRUE)

net <- WGCNA::blockwiseModules(
  datExpr,
  power = chosen_power,
  TOMType = TOMType,
  networkType = netType,
  corType = corType,
  maxPOutliers = ifelse(corType == "bicor", maxPOut, 1),
  minModuleSize = minModuleSize,
  deepSplit = deepSplit,
  pamRespectsDendro = pamRespectsDendro,
  mergeCutHeight = mergeCutHeight,
  reassignThreshold = reassignThreshold,
  numericLabels = numericLabels,
  saveTOMs = FALSE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- WGCNA::labels2colors(moduleLabels)

log_msg(logger, "INFO", paste0("Modules found (including grey): ", length(unique(moduleColors))))

# =============================================================================
# 12) 导出模块结果：基因-模块表、模块基因列表
# =============================================================================

gene_ids_final <- colnames(datExpr)
gene_module_df <- data.frame(
  gene_id = gene_ids_final,
  module_label = moduleLabels,
  module_color = moduleColors,
  stringsAsFactors = FALSE
)

if (isTRUE(W$export$gene_module_table)) {
  data.table::fwrite(gene_module_df, file.path(SUBDIRS$modules, "gene_module_assignments.tsv"), sep = "\t", quote = FALSE)
}

if (isTRUE(W$export$module_lists)) {
  unique_mods <- sort(unique(moduleColors))
  unique_mods <- unique_mods[unique_mods != "grey"]
  for (mc in unique_mods) {
    genes <- gene_module_df$gene_id[gene_module_df$module_color == mc]
    out_list <- file.path(SUBDIRS$lists, paste0("module_", mc, ".list"))
    writeLines(genes, out_list, useBytes = TRUE)
  }
}

saveRDS(net$dendrograms, file.path(SUBDIRS$modules, "gene_dendrograms.rds"))
saveRDS(net$blockGenes,  file.path(SUBDIRS$modules, "blockGenes.rds"))

# =============================================================================
# 13) MEs 与模块-性状关联（若 traits 存在）
# =============================================================================

MEs0 <- WGCNA::moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs  <- WGCNA::orderMEs(MEs0)

data.table::fwrite(
  data.frame(sample = rownames(MEs), MEs, check.names = FALSE),
  file.path(SUBDIRS$modules, "module_eigengenes.tsv"),
  sep = "\t",
  quote = FALSE
)

if (!is.null(traits_mat) && ncol(traits_mat) > 0) {
  log_msg(logger, "INFO", "Computing module-trait correlations...")

  # ✅ 增量稳妥：强制对齐 traits_mat 与 MEs 的样本顺序，避免隐式错配
  if (is.null(rownames(traits_mat))) {
    stop("traits_mat has no rownames; expected sample IDs as rownames.")
  }
  miss_in_traits <- setdiff(rownames(MEs), rownames(traits_mat))
  if (length(miss_in_traits) > 0) {
    stop("traits_mat is missing samples present in MEs: ", paste(miss_in_traits, collapse = ", "))
  }
  traits_mat_aligned <- traits_mat[rownames(MEs), , drop = FALSE]

  mt_cor <- WGCNA::cor(MEs, traits_mat_aligned, use = useObs)
  mt_p   <- WGCNA::corPvalueStudent(mt_cor, nSamples = nrow(datExpr))
  mt_padj <- matrix(p.adjust(as.vector(mt_p), method = "BH"),
                    nrow = nrow(mt_p), ncol = ncol(mt_p),
                    dimnames = dimnames(mt_p))

  # ✅ 改进：TSV 第一列写模块名（module=ME name），并额外输出旧格式 matrix-only 备份
  write_matrix_with_rownames(
    mt_cor,
    file.path(SUBDIRS$assoc, "module_trait_cor.tsv"),
    rowname_col = "module",
    out_matrix_only = file.path(SUBDIRS$assoc, "module_trait_cor.matrix_only.tsv")
  )
  write_matrix_with_rownames(
    mt_p,
    file.path(SUBDIRS$assoc, "module_trait_p.tsv"),
    rowname_col = "module",
    out_matrix_only = file.path(SUBDIRS$assoc, "module_trait_p.matrix_only.tsv")
  )
  write_matrix_with_rownames(
    mt_padj,
    file.path(SUBDIRS$assoc, "module_trait_padj.tsv"),
    rowname_col = "module",
    out_matrix_only = file.path(SUBDIRS$assoc, "module_trait_padj.matrix_only.tsv")
  )

  saveRDS(mt_cor,  file.path(SUBDIRS$assoc, "module_trait_cor.rds"))
  saveRDS(mt_p,    file.path(SUBDIRS$assoc, "module_trait_p.rds"))
  saveRDS(mt_padj, file.path(SUBDIRS$assoc, "module_trait_padj.rds"))
} else {
  log_msg(logger, "WARN", "traits 为空：跳过模块-性状关联与 hub genes。")
}

# =============================================================================
# 14) Hub genes：GS + kME（关键：不使用 signedKME，直接 bicor/cor 计算，避开 eval 雷区）
# =============================================================================

safe_cor <- function(x, y, corType, useObs, maxPOut) {
  if (corType == "bicor") {
    return(WGCNA::bicor(x, y, use = useObs, maxPOutliers = maxPOut))
  } else {
    return(stats::cor(x, y, use = useObs, method = "pearson"))
  }
}

if (!is.null(traits_mat) && ncol(traits_mat) > 0 && isTRUE(W$export$hub_genes$enable %||% TRUE)) {

  trait_names <- colnames(traits_mat)

  pick_primary_trait <- function() {
    cfg_primary <- W$export$hub_genes$primary_trait %||% ""
    cfg_primary <- as.character(cfg_primary)

    if (nzchar(cfg_primary)) {
      if (cfg_primary %in% trait_names) return(cfg_primary)
      stop("Configured primary_trait not found in trait matrix: ", cfg_primary)
    }

    if ("target_vs_rest" %in% trait_names) return("target_vs_rest")
    return(trait_names[1])
  }

  primary_trait_name <- pick_primary_trait()
  log_msg(logger, "INFO", paste0("Primary trait for hub ranking: ", primary_trait_name))

  trait_vec <- as.numeric(traits_mat[[primary_trait_name]])
  names(trait_vec) <- rownames(traits_mat)

  trait_mat <- matrix(trait_vec, ncol = 1)
  rownames(trait_mat) <- rownames(datExpr)
  colnames(trait_mat) <- primary_trait_name

  GS_mat <- safe_cor(datExpr, trait_mat, corType = corType, useObs = useObs, maxPOut = maxPOut)
  GS <- as.numeric(GS_mat[, 1])
  names(GS) <- colnames(datExpr)

  GS_p <- WGCNA::corPvalueStudent(GS, nSamples = nrow(datExpr))

  kME_mat <- safe_cor(datExpr, MEs, corType = corType, useObs = useObs, maxPOut = maxPOut)
  kME_df <- as.data.frame(kME_mat, check.names = FALSE)
  colnames(kME_df) <- paste0("kME_", colnames(kME_df))

  hub_base <- gene_module_df
  hub_base$GS <- as.numeric(GS[hub_base$gene_id])
  hub_base$GS_p <- as.numeric(GS_p[hub_base$gene_id])

  hub_df <- cbind(hub_base, kME_df[hub_base$gene_id, , drop = FALSE])

  data.table::fwrite(hub_df, file.path(SUBDIRS$hub, "gene_level_GS_kME.tsv"), sep = "\t", quote = FALSE)

  top_n <- as.integer(W$export$hub_genes$top_n_per_module %||% 50)
  rank_by <- W$export$hub_genes$rank_by %||% c("kME", "GS")

  map_color_to_ME <- function(color) paste0("ME", color)
  get_kME_col_for_module <- function(me_name) paste0("kME_", me_name)

  mods <- sort(unique(hub_df$module_color))
  mods <- mods[mods != "grey"]

  hub_list_all <- list()

  for (mc in mods) {
    sub <- hub_df[hub_df$module_color == mc, , drop = FALSE]
    me_col <- get_kME_col_for_module(map_color_to_ME(mc))
    if (!me_col %in% colnames(sub)) next

    if (length(rank_by) >= 1 && rank_by[1] == "kME") {
      sub <- sub[order(abs(sub[[me_col]]), decreasing = TRUE), , drop = FALSE]
    } else if (length(rank_by) >= 1 && rank_by[1] == "GS") {
      sub <- sub[order(abs(sub$GS), decreasing = TRUE), , drop = FALSE]
    }

    if (length(rank_by) >= 2) {
      if (rank_by[2] == "GS") {
        sub <- sub[order(abs(sub[[me_col]]), abs(sub$GS), decreasing = TRUE), , drop = FALSE]
      } else if (rank_by[2] == "kME") {
        sub <- sub[order(abs(sub$GS), abs(sub[[me_col]]), decreasing = TRUE), , drop = FALSE]
      }
    }

    topk <- head(sub, n = min(top_n, nrow(sub)))
    topk$kME_primary_module <- topk[[me_col]]
    topk$primary_trait <- primary_trait_name
    topk$module_ME <- map_color_to_ME(mc)
    hub_list_all[[mc]] <- topk
  }

  hub_out <- do.call(rbind, hub_list_all)
  if (!is.null(hub_out) && nrow(hub_out) > 0) {
    data.table::fwrite(hub_out, file.path(SUBDIRS$hub, "hub_genes_by_module.tsv"), sep = "\t", quote = FALSE)
  } else {
    log_msg(logger, "WARN", "未生成 hub genes 表（可能模块过少或 traits 定义较弱）。")
  }
}

# =============================================================================
# 15) 保存复现参数与核心对象
# =============================================================================

used_params <- list(
  config_yaml = CONFIG_YAML,
  wgcna = W,
  matched_samples = rownames(datExpr),
  n_samples = nrow(datExpr),
  n_genes_final = ncol(datExpr),
  chosen_power = chosen_power,
  method = method,
  threads = W$threads
)
write_yaml_safe(used_params, file.path(OUTDIR, "run_params.used.yaml"))

if (isTRUE(W$save_rds)) {
  saveRDS(
    list(
      W = W,
      datExpr = datExpr,
      expr_mat = expr_mat,
      sampleTree = sampleTree,
      sft = sft,
      chosen_power = chosen_power,
      moduleLabels = moduleLabels,
      moduleColors = moduleColors,
      MEs = MEs,
      traits = traits_mat,
      samples = samples_dt,
      net = net
    ),
    file.path(SUBDIRS$rds, "wgcna_objects.rds")
  )
}

log_msg(logger, "INFO", "DONE: WGCNA (compute-only) finished successfully.")

