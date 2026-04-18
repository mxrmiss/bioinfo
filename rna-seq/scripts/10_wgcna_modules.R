#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

CONFIG_YAML <- "config.yaml"

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

`%||%` <- function(a, b) if (!is.null(a)) a else b

read_yaml_safe <- function(path) {
  if (!file.exists(path)) stop("找不到配置文件：", path)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("需要先安装 R 包 yaml")
  }
  yaml::read_yaml(path)
}

write_yaml_safe <- function(x, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("需要先安装 R 包 yaml")
  }
  safe_mkdir(dirname(path))
  yaml::write_yaml(x, path)
}

as_char_vec <- function(x) {
  if (is.null(x)) return(character())
  if (length(x) == 0) return(character())
  out <- unlist(x, use.names = FALSE)
  out <- as.character(out)
  out <- out[nzchar(out)]
  unique(out)
}

as_named_numeric <- function(x, arg_name = "mapping") {
  if (is.null(x) || length(x) == 0) {
    stop(arg_name, " is empty.")
  }
  if (!is.null(names(x)) && any(nzchar(names(x)))) {
    out <- suppressWarnings(as.numeric(unlist(x, use.names = TRUE)))
    nms <- names(unlist(x, use.names = TRUE))
    if (any(is.na(out))) stop(arg_name, " contains non-numeric values.")
    names(out) <- nms
    return(out)
  }
  stop(arg_name, " must be a named vector or named list.")
}

write_matrix_with_rownames <- function(mat, out_tsv, rowname_col = "module", out_matrix_only = NULL) {
  if (is.null(rownames(mat))) {
    stop("矩阵没有行名，无法导出。")
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
    data.table::fwrite(as.data.frame(mat, check.names = FALSE), out_matrix_only, sep = "\t", quote = FALSE)
  }
}

safe_cor <- function(x, y, corType, useObs, maxPOut) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (tolower(corType) == "bicor") {
    return(WGCNA::bicor(x, y, use = useObs, maxPOutliers = maxPOut))
  } else {
    return(stats::cor(x, y, use = useObs, method = "pearson"))
  }
}


pick_primary_trait_auto <- function(mt_cor, mt_padj, trait_meta_df, cfg_sel, logger = NULL) {
  if (is.null(mt_cor) || is.null(mt_padj) || nrow(mt_cor) == 0 || ncol(mt_cor) == 0) return(NA_character_)
  trait_names <- colnames(mt_cor)
  max_abs_cor <- apply(abs(mt_cor), 2, max, na.rm = TRUE)
  min_padj <- apply(mt_padj, 2, min, na.rm = TRUE)
  meta <- data.frame(trait = trait_names, stringsAsFactors = FALSE)
  if (!is.null(trait_meta_df) && nrow(trait_meta_df) > 0) {
    meta <- merge(meta, trait_meta_df, by = 'trait', all.x = TRUE, sort = FALSE)
    meta <- meta[match(trait_names, meta$trait), , drop = FALSE]
  }
  meta$max_abs_cor <- as.numeric(max_abs_cor[meta$trait])
  meta$min_padj <- as.numeric(min_padj[meta$trait])
  meta$pass <- !is.na(meta$min_padj) & meta$min_padj <= as.numeric(cfg_sel$padj_cutoff %||% 0.05)

  excl <- as_char_vec(cfg_sel$exclude_prefixes %||% character())
  if (length(excl) > 0) {
    keep <- rep(TRUE, nrow(meta))
    for (px in excl) keep <- keep & !startsWith(meta$trait, px)
    meta <- meta[keep, , drop = FALSE]
  }
  if (nrow(meta) == 0) return(NA_character_)

  prefer_cont <- isTRUE(cfg_sel$prefer_continuous %||% TRUE)
  if (prefer_cont && 'trait_kind' %in% colnames(meta)) {
    meta$kind_rank <- ifelse(meta$trait_kind %in% c('numeric', 'ordered'), 1L,
                             ifelse(meta$trait_kind %in% c('binary'), 2L, 3L))
  } else {
    meta$kind_rank <- 1L
  }

  method <- tolower(as.character(cfg_sel$method %||% 'best_significant'))
  sig <- meta[meta$pass, , drop = FALSE]
  if (method == 'first') return(meta$trait[1])
  if (method == 'largest_abs_cor') {
    cand <- if (nrow(sig) > 0) sig else meta
    cand <- cand[order(cand$kind_rank, -cand$max_abs_cor, cand$min_padj), , drop = FALSE]
    return(cand$trait[1])
  }
  cand <- if (nrow(sig) > 0) sig else meta
  cand <- cand[order(cand$kind_rank, cand$min_padj, -cand$max_abs_cor), , drop = FALSE]
  if (!is.null(logger) && nrow(cand) > 0) {
    log_msg(logger, 'INFO', paste0('Auto primary trait selected by rule: ', cand$trait[1]))
  }
  cand$trait[1]
}

get_group_levels <- function(samples_dt, group_col, ordered_group_levels = character()) {
  og <- as_char_vec(ordered_group_levels)
  present_groups <- unique(as.character(samples_dt[[group_col]]))

  if (length(og) == 0) {
    return(present_groups)
  }

  miss <- setdiff(og, present_groups)
  if (length(miss) > 0) {
    stop("ordered_group_levels 中包含 samples.tsv 里不存在的 group：", paste(miss, collapse = ", "))
  }

  extra <- setdiff(present_groups, og)
  if (length(extra) > 0) {
    stop("samples.tsv 中有 group 不在 ordered_group_levels 里：", paste(extra, collapse = ", "))
  }

  og
}

require_cfg <- function(x, name) {
  if (is.null(x)) stop("配置缺失：", name)
  x
}


normalize_logical <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  if (is.logical(x)) return(isTRUE(x[1]))
  if (is.numeric(x)) return(isTRUE(as.logical(x[1])))
  x1 <- tolower(trimws(as.character(x[1])))
  if (x1 %in% c('true', 't', '1', 'yes', 'y')) return(TRUE)
  if (x1 %in% c('false', 'f', '0', 'no', 'n')) return(FALSE)
  default
}

normalize_trait_kind <- function(x, default = 'numeric') {
  x <- tolower(as.character(x %||% default))
  if (!x %in% c('numeric', 'categorical', 'onehot')) default else x
}

normalize_trait_items <- function(items) {
  if (is.null(items) || length(items) == 0) return(list())
  out <- list()
  idx <- 0L
  for (it in items) {
    idx <- idx + 1L
    if (is.null(it)) next
    nm <- as.character(it$name %||% '')
    if (!nzchar(nm) && !is.null(it$column)) nm <- as.character(it$column)
    if (!nzchar(nm) && !is.null(it$trait_name)) nm <- as.character(it$trait_name)
    if (!nzchar(nm)) nm <- paste0('trait_', idx)
    it$name <- nm
    out[[length(out) + 1L]] <- it
  }
  out
}

normalize_traits_config <- function(W) {
  # 只保留一套通用 trait 框架；仍允许读取少量旧字段作为别名，避免老配置立刻失效
  W$traits$mode <- 'generic'

  if (is.null(W$traits$build_group_onehot)) {
    W$traits$build_group_onehot <- W$traits$build_onehot %||% TRUE
  }
  if (is.null(W$traits$group_onehot_prefix)) {
    W$traits$group_onehot_prefix <- W$traits$onehot_prefix %||% 'group_'
  }

  if (is.null(W$traits$trait_columns) || length(W$traits$trait_columns) == 0) {
    W$traits$trait_columns <- list()
  }
  if (is.null(W$traits$binary_traits)) W$traits$binary_traits <- list()
  if (is.null(W$traits$score_traits)) W$traits$score_traits <- W$traits$group_score_traits %||% list()
  if (is.null(W$traits$sample_numeric_traits)) W$traits$sample_numeric_traits <- list()
  if (is.null(W$traits$sample_categorical_traits)) W$traits$sample_categorical_traits <- list()
  if (is.null(W$traits$sample_trait_columns)) W$traits$sample_trait_columns <- W$traits$extra_columns_as_traits %||% character()
  if (is.null(W$traits$group_order_trait)) W$traits$group_order_trait <- W$traits$auto_group_order_trait %||% list(enable = FALSE, name = 'group_order')

  # 从直接列名列表自动补成 trait_columns
  if (length(as_char_vec(W$traits$sample_trait_columns)) > 0) {
    for (cn in as_char_vec(W$traits$sample_trait_columns)) {
      W$traits$trait_columns[[length(W$traits$trait_columns) + 1L]] <- list(column = cn, name = cn, kind = 'auto')
    }
  }

  # 向后兼容：若老配置里仍残留扩展模式字段，则自动折叠进通用框架，但不再要求/暴露 mode
  if (!is.null(W$traits$target_groups) && length(as_char_vec(W$traits$target_groups)) > 0) {
    target_trait_name <- as.character(W$traits$target_trait_name %||% 'target_vs_rest')
    W$traits$binary_traits[[length(W$traits$binary_traits) + 1L]] <- list(
      name = target_trait_name,
      positive_groups = as_char_vec(W$traits$target_groups),
      negative_groups = c()
    )
    if (!nzchar(as.character(W$export$hub_genes$primary_trait %||% ''))) {
      W$export$hub_genes$primary_trait <- target_trait_name
    }
  }

  if (!is.null(W$traits$developmental_binary_traits) && length(W$traits$developmental_binary_traits) > 0) {
    W$traits$binary_traits <- c(W$traits$binary_traits, W$traits$developmental_binary_traits)
  }
  if (!is.null(W$traits$developmental_score_traits) && length(W$traits$developmental_score_traits) > 0) {
    W$traits$score_traits <- c(W$traits$score_traits, W$traits$developmental_score_traits)
  }
  if (!is.null(W$traits$ordered_group_levels) && length(as_char_vec(W$traits$ordered_group_levels)) > 0) {
    if (is.null(W$traits$group_order_trait) || !isTRUE(W$traits$group_order_trait$enable)) {
      stage_name <- as.character(W$traits$stage_order_trait_name %||% W$traits$group_order_trait$name %||% 'group_order')
      W$traits$group_order_trait <- list(enable = TRUE, name = stage_name)
    }
  }

  W$traits$binary_traits <- normalize_trait_items(W$traits$binary_traits)
  W$traits$score_traits <- normalize_trait_items(W$traits$score_traits)
  W$traits$trait_columns <- normalize_trait_items(W$traits$trait_columns)
  W$traits$sample_numeric_traits <- normalize_trait_items(W$traits$sample_numeric_traits)
  W$traits$sample_categorical_traits <- normalize_trait_items(W$traits$sample_categorical_traits)

  if (is.null(W$traits$redundancy)) W$traits$redundancy <- list()
  if (is.null(W$traits$redundancy$enable)) W$traits$redundancy$enable <- TRUE
  if (is.null(W$traits$redundancy$cor_threshold)) W$traits$redundancy$cor_threshold <- 0.85
  if (is.null(W$traits$redundancy$p_adjust_method)) W$traits$redundancy$p_adjust_method <- 'BH'

  if (is.null(W$export$hub_genes$primary_trait_selection)) W$export$hub_genes$primary_trait_selection <- list()
  if (is.null(W$export$hub_genes$primary_trait_selection$method)) W$export$hub_genes$primary_trait_selection$method <- 'best_significant'
  if (is.null(W$export$hub_genes$primary_trait_selection$padj_cutoff)) W$export$hub_genes$primary_trait_selection$padj_cutoff <- 0.05
  if (is.null(W$export$hub_genes$primary_trait_selection$prefer_continuous)) W$export$hub_genes$primary_trait_selection$prefer_continuous <- TRUE
  if (is.null(W$export$hub_genes$primary_trait_selection$exclude_prefixes)) W$export$hub_genes$primary_trait_selection$exclude_prefixes <- character()

  W
}

cfg_all <- read_yaml_safe(CONFIG_YAML)
cfg_wgcna <- require_cfg(cfg_all$wgcna, "wgcna")

W <- cfg_wgcna

if (!is.null(cfg_all$data) && !is.null(cfg_all$data$samples_tsv) && is.null(W$samples_tsv)) {
  W$samples_tsv <- cfg_all$data$samples_tsv
}

if (!is.null(cfg_all$dirs) && !is.null(cfg_all$dirs$wgcna) && is.null(W$outdir)) {
  W$outdir <- cfg_all$dirs$wgcna
}

if (!is.null(cfg_all$resources) && !is.null(cfg_all$resources$threads) &&
    !is.null(cfg_all$resources$threads$wgcna)) {
  W$threads <- cfg_all$resources$threads$wgcna
}

W <- normalize_traits_config(W)

require_cfg(W$enable, "wgcna.enable")
require_cfg(W$method, "wgcna.method")
require_cfg(W$input_counts, "wgcna.input_counts")
require_cfg(W$samples_tsv, "wgcna.samples_tsv")
require_cfg(W$outdir, "wgcna.outdir")
require_cfg(W$save_rds, "wgcna.save_rds")
require_cfg(W$random_seed, "wgcna.random_seed")
require_cfg(W$threads, "wgcna.threads 或 resources.threads.wgcna")

require_cfg(W$filter, "wgcna.filter")
require_cfg(W$filter$min_count, "wgcna.filter.min_count")
require_cfg(W$filter$min_samples, "wgcna.filter.min_samples")
require_cfg(W$filter$var_metric, "wgcna.filter.var_metric")
if (is.null(W$filter$keep_top_fraction) && is.null(W$filter$keep_top_n)) {
  stop("必须提供 wgcna.filter.keep_top_fraction 或 wgcna.filter.keep_top_n")
}

if (tolower(W$method) == "vst") {
  require_cfg(W$vst, "wgcna.vst")
  require_cfg(W$vst$blind, "wgcna.vst.blind")
}
if (tolower(W$method) == "voom") {
  require_cfg(W$voom, "wgcna.voom")
  require_cfg(W$voom$enable_tmm, "wgcna.voom.enable_tmm")
  if (is.null(W$voom$plot_mean_var)) W$voom$plot_mean_var <- FALSE
}

require_cfg(W$network, "wgcna.network")
require_cfg(W$network$networkType, "wgcna.network.networkType")
require_cfg(W$network$TOMType, "wgcna.network.TOMType")
require_cfg(W$network$corType, "wgcna.network.corType")
require_cfg(W$network$maxPOutliers, "wgcna.network.maxPOutliers")
require_cfg(W$network$useObs, "wgcna.network.useObs")

require_cfg(W$power, "wgcna.power")
require_cfg(W$power$candidates, "wgcna.power.candidates")
require_cfg(W$power$scale_free_R2_target, "wgcna.power.scale_free_R2_target")
require_cfg(W$power$min_mean_connectivity, "wgcna.power.min_mean_connectivity")

require_cfg(W$modules, "wgcna.modules")
require_cfg(W$modules$minModuleSize, "wgcna.modules.minModuleSize")
require_cfg(W$modules$deepSplit, "wgcna.modules.deepSplit")
require_cfg(W$modules$pamRespectsDendro, "wgcna.modules.pamRespectsDendro")
require_cfg(W$modules$mergeCutHeight, "wgcna.modules.mergeCutHeight")
require_cfg(W$modules$reassignThreshold, "wgcna.modules.reassignThreshold")
require_cfg(W$modules$numericLabels, "wgcna.modules.numericLabels")

require_cfg(W$traits, "wgcna.traits")
require_cfg(W$traits$sample_column, "wgcna.traits.sample_column")
require_cfg(W$traits$group_column, "wgcna.traits.group_column")
require_cfg(W$traits$build_group_onehot, 'wgcna.traits.build_group_onehot')
require_cfg(W$traits$group_onehot_prefix, 'wgcna.traits.group_onehot_prefix')
if (is.null(W$traits$trait_columns)) W$traits$trait_columns <- list()
if (is.null(W$traits$sample_trait_columns)) W$traits$sample_trait_columns <- character()

require_cfg(W$export, "wgcna.export")
require_cfg(W$export$gene_module_table, "wgcna.export.gene_module_table")
require_cfg(W$export$module_lists, "wgcna.export.module_lists")
require_cfg(W$export$hub_genes, "wgcna.export.hub_genes")
require_cfg(W$export$hub_genes$enable, "wgcna.export.hub_genes.enable")
require_cfg(W$export$hub_genes$top_n_per_module, "wgcna.export.hub_genes.top_n_per_module")
require_cfg(W$export$hub_genes$rank_by, "wgcna.export.hub_genes.rank_by")
if (is.null(W$export$hub_genes$primary_trait_candidates)) {
  W$export$hub_genes$primary_trait_candidates <- character()
}
if (is.null(W$export$hub_genes$primary_trait)) {
  W$export$hub_genes$primary_trait <- ""
}

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
log_msg(logger, "INFO", "Traits mode: generic (unified trait framework)")
log_msg(logger, "INFO", paste0("Threads: ", W$threads))

if (!isTRUE(W$enable)) {
  log_msg(logger, "WARN", "wgcna.enable=false，脚本退出。")
  quit(save = "no", status = 0)
}

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

WGCNA::allowWGCNAThreads(nThreads = as.integer(W$threads))
set.seed(as.integer(W$random_seed))

log_msg(logger, "INFO", paste0("Reading samples: ", W$samples_tsv))
if (!file.exists(W$samples_tsv)) {
  stop("samples.tsv 不存在：", W$samples_tsv)
}
samples_dt <- data.table::fread(W$samples_tsv, sep = "\t", header = TRUE, data.table = FALSE)

sc <- W$traits$sample_column
gc <- W$traits$group_column
if (!sc %in% colnames(samples_dt)) stop("samples.tsv 缺少列：", sc)
if (!gc %in% colnames(samples_dt)) stop("samples.tsv 缺少列：", gc)

samples_dt[[sc]] <- as.character(samples_dt[[sc]])
samples_dt[[gc]] <- as.character(samples_dt[[gc]])

group_levels <- get_group_levels(
  samples_dt = samples_dt,
  group_col = gc,
  ordered_group_levels = W$traits$ordered_group_levels %||% character()
)


trait_list <- list()
trait_meta <- list()

register_trait_block <- function(key, mat, source, trait_kind, group_column = NA_character_, details = NA_character_) {
  if (is.null(mat) || ncol(mat) == 0) return(NULL)
  mat <- as.matrix(mat)
  rownames(mat) <- samples_dt[[sc]]
  trait_list[[key]] <<- mat
  meta <- data.frame(
    trait = colnames(mat),
    source = source,
    trait_kind = trait_kind,
    group_column = group_column,
    details = details,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  trait_meta[[length(trait_meta) + 1L]] <<- meta
  invisible(NULL)
}

if (isTRUE(W$traits$build_group_onehot)) {
  onehot_prefix <- as.character(W$traits$group_onehot_prefix)
  onehot <- sapply(group_levels, function(g) as.integer(samples_dt[[gc]] == g))
  colnames(onehot) <- paste0(onehot_prefix, group_levels)
  rownames(onehot) <- samples_dt[[sc]]
  register_trait_block('onehot_group', onehot, source = 'group_onehot', trait_kind = 'categorical', group_column = gc)
}

binary_traits <- W$traits$binary_traits %||% list()
all_groups <- unique(as.character(samples_dt[[gc]]))

if (length(binary_traits) > 0) {
  for (bt in binary_traits) {
    trait_name <- as.character(require_cfg(bt$name, 'binary_traits.name'))
    pos <- as_char_vec(require_cfg(bt$positive_groups, paste0('binary_traits[', trait_name, '].positive_groups')))
    neg <- as_char_vec(bt$negative_groups %||% character())

    miss_pos <- setdiff(pos, all_groups)
    miss_neg <- setdiff(neg, all_groups)
    if (length(miss_pos) > 0) stop('binary_traits 正组不存在：', paste(miss_pos, collapse = ', '))
    if (length(miss_neg) > 0) stop('binary_traits 负组不存在：', paste(miss_neg, collapse = ', '))

    if (length(pos) == 0) stop('binary_traits 的 positive_groups 不能为空：', trait_name)
    if (length(neg) == 0) neg <- setdiff(all_groups, pos)

    inter <- intersect(pos, neg)
    if (length(inter) > 0) stop('binary_traits 正负组重叠：', paste(inter, collapse = ', '))

    vec <- rep(NA_integer_, nrow(samples_dt))
    vec[samples_dt[[gc]] %in% pos] <- 1L
    vec[samples_dt[[gc]] %in% neg] <- 0L
    if (any(is.na(vec))) stop('binary_traits 生成了 NA：', trait_name)

    mat <- matrix(vec, ncol = 1, dimnames = list(samples_dt[[sc]], trait_name))
    register_trait_block(paste0('binary__', trait_name), mat,
                         source = 'binary_trait', trait_kind = 'binary',
                         group_column = gc,
                         details = paste0('positive=', paste(pos, collapse = ','), ';negative=', paste(neg, collapse = ',')))
  }
}

trait_columns <- W$traits$trait_columns %||% list()
if (length(trait_columns) > 0) {
  for (tc in trait_columns) {
    col_in <- as.character(require_cfg(tc$column %||% tc$name, 'trait_columns.column'))
    trait_name <- as.character(tc$name %||% col_in)
    kind <- normalize_trait_kind(tc$kind %||% 'numeric', default = 'numeric')
    if (!col_in %in% colnames(samples_dt)) stop('trait_columns 指定的列不存在于 samples.tsv：', col_in)
    v <- samples_dt[[col_in]]

    if (kind == 'numeric' || FALSE) {
      suppressWarnings(v_num <- as.numeric(v))
      if (all(is.na(v_num))) stop('trait_columns 列无法转为 numeric：', col_in)
      mat <- matrix(v_num, ncol = 1, dimnames = list(samples_dt[[sc]], trait_name))
      register_trait_block(paste0('trait_column__', trait_name), mat,
                           source = 'sample_column', trait_kind = 'numeric', group_column = col_in)
    } else {
      v_chr <- as.character(v)
      lv <- unique(v_chr)
      oh <- sapply(lv, function(z) as.integer(v_chr == z))
      colnames(oh) <- paste0(trait_name, '_', lv)
      rownames(oh) <- samples_dt[[sc]]
      register_trait_block(paste0('trait_column__', trait_name), oh,
                           source = 'sample_column', trait_kind = 'categorical', group_column = col_in)
    }
  }
}

ag <- W$traits$group_order_trait %||% list(enable = FALSE, name = 'group_order')
if (isTRUE(ag$enable %||% FALSE)) {
  trait_name <- as.character(require_cfg(ag$name, 'wgcna.traits.group_order_trait.name'))
  idx_map <- seq_along(group_levels)
  names(idx_map) <- group_levels
  vec <- as.numeric(idx_map[samples_dt[[gc]]])
  mat <- matrix(vec, ncol = 1, dimnames = list(samples_dt[[sc]], trait_name))
  register_trait_block(paste0('group_order__', trait_name), mat,
                       source = 'group_order', trait_kind = 'ordered', group_column = gc)
}

score_traits <- W$traits$score_traits %||% list()
if (length(score_traits) > 0) {
  for (gst in score_traits) {
    trait_name <- as.character(require_cfg(gst$name, 'score_traits.name'))
    mapping <- as_named_numeric(require_cfg(gst$mapping, paste0('score_traits[', trait_name, ']$mapping')),
                                arg_name = paste0('score_traits[', trait_name, ']$mapping'))

    miss_map <- setdiff(all_groups, names(mapping))
    extra_map <- setdiff(names(mapping), all_groups)
    if (length(miss_map) > 0) stop('score_traits 缺少 group 映射：', paste(miss_map, collapse = ', '))
    if (length(extra_map) > 0) stop('score_traits 含未知 group：', paste(extra_map, collapse = ', '))

    vec <- as.numeric(mapping[samples_dt[[gc]]])
    if (any(is.na(vec))) stop('score_traits 生成了 NA：', trait_name)

    mat <- matrix(vec, ncol = 1, dimnames = list(samples_dt[[sc]], trait_name))
    register_trait_block(paste0('score__', trait_name), mat,
                         source = 'score_trait', trait_kind = 'numeric', group_column = gc)
  }
}

traits_mat <- NULL
trait_meta_df <- data.frame()
if (length(trait_list) > 0) {
  traits_mat <- do.call(cbind, trait_list)
  traits_mat <- as.data.frame(traits_mat, check.names = FALSE)
  rownames(traits_mat) <- samples_dt[[sc]]

  trait_meta_df <- do.call(rbind, trait_meta)
  trait_meta_df <- trait_meta_df[match(colnames(traits_mat), trait_meta_df$trait), , drop = FALSE]

  data.table::fwrite(
    data.frame(sample = rownames(traits_mat), traits_mat, check.names = FALSE),
    file.path(SUBDIRS$traits, 'trait_matrix.tsv'),
    sep = '\t',
    quote = FALSE
  )
  data.table::fwrite(trait_meta_df, file.path(SUBDIRS$traits, 'trait_metadata.tsv'), sep = '\t', quote = FALSE)
  log_msg(logger, 'INFO', paste0('Traits built: ', ncol(traits_mat), ' columns'))
} else {
  log_msg(logger, 'WARN', '未构建 traits，后续将跳过模块-性状关联与 hub genes。')
}

log_msg(logger, "INFO", paste0("Reading counts: ", W$input_counts))
if (!file.exists(W$input_counts)) {
  stop("counts 文件不存在：", W$input_counts)
}

counts_dt <- data.table::fread(W$input_counts, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
if (ncol(counts_dt) < 3) stop("counts 表至少需要 gene_id + 2 个样本列。")

gene_col <- colnames(counts_dt)[1]
gene_ids <- as.character(counts_dt[[gene_col]])
counts_mat <- as.matrix(counts_dt[, -1, drop = FALSE])
rownames(counts_mat) <- gene_ids
suppressWarnings(storage.mode(counts_mat) <- "numeric")

if (anyNA(counts_mat)) {
  log_msg(logger, "WARN", "counts 含 NA，已按 0 处理。")
  counts_mat[is.na(counts_mat)] <- 0
}

sample_names_counts <- colnames(counts_mat)
sample_names_meta   <- samples_dt[[sc]]
common_samples <- intersect(sample_names_meta, sample_names_counts)
if (length(common_samples) < 6) stop("counts 与 samples.tsv 匹配到的样本数过少。")

common_samples <- sample_names_meta[sample_names_meta %in% common_samples]
counts_mat <- counts_mat[, common_samples, drop = FALSE]
samples_dt <- samples_dt[match(common_samples, samples_dt[[sc]]), , drop = FALSE]
if (!is.null(traits_mat)) traits_mat <- traits_mat[common_samples, , drop = FALSE]

log_msg(logger, "INFO", paste0("Matched samples: ", length(common_samples)))

min_count   <- as.numeric(W$filter$min_count)
min_samples <- as.integer(W$filter$min_samples)

keep_expr <- rowSums(counts_mat >= min_count) >= min_samples
counts_f1 <- counts_mat[keep_expr, , drop = FALSE]

log_msg(
  logger, "INFO",
  sprintf("Low-expression filter: kept %d / %d genes (min_count=%s in >=%s samples)",
          nrow(counts_f1), nrow(counts_mat), min_count, min_samples)
)

method <- tolower(W$method)

datExpr <- NULL
expr_mat <- NULL

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
  vsd <- DESeq2::vst(dds, blind = isTRUE(W$vst$blind))

  expr_mat <- SummarizedExperiment::assay(vsd)
  datExpr  <- t(expr_mat)

} else if (method == "voom") {
  log_msg(logger, "INFO", "Transform method: limma-voom")

  dge <- edgeR::DGEList(counts = counts_f1)
  if (isTRUE(W$voom$enable_tmm)) {
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
  }

  design <- matrix(1, nrow = length(common_samples), ncol = 1)
  v <- limma::voom(dge, design = design, plot = isTRUE(W$voom$plot_mean_var))
  expr_mat <- v$E
  datExpr  <- t(expr_mat)

} else {
  stop("wgcna.method 只能是 vst 或 voom")
}

data.table::fwrite(
  data.frame(sample = rownames(datExpr), datExpr, check.names = FALSE),
  file.path(SUBDIRS$modules, "datExpr.transformed.samples_x_genes.tsv"),
  sep = "\t",
  quote = FALSE
)

var_metric <- tolower(W$filter$var_metric)
keep_top_n <- W$filter$keep_top_n %||% NULL
keep_frac  <- W$filter$keep_top_fraction %||% NULL

if (var_metric == "mad") {
  gene_stats <- apply(datExpr, 2, mad, na.rm = TRUE)
} else if (var_metric %in% c("variance", "var")) {
  gene_stats <- apply(datExpr, 2, var, na.rm = TRUE)
} else {
  stop("wgcna.filter.var_metric 只能是 mad 或 variance")
}

ord <- order(gene_stats, decreasing = TRUE, na.last = NA)
if (!is.null(keep_top_n)) {
  keep_top_n <- as.integer(keep_top_n)
  k <- min(length(ord), keep_top_n)
} else {
  keep_frac <- as.numeric(keep_frac)
  keep_frac <- max(min(keep_frac, 1), 0.01)
  k <- max(1, floor(length(ord) * keep_frac))
}
keep_genes <- names(gene_stats)[ord[seq_len(k)]]
datExpr <- datExpr[, keep_genes, drop = FALSE]

log_msg(
  logger, "INFO",
  sprintf("Low-variance filter (%s): kept %d genes", var_metric, ncol(datExpr))
)

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

corType <- tolower(W$network$corType)
useObs  <- as.character(W$network$useObs)
maxPOut <- as.numeric(W$network$maxPOutliers)
netType <- W$network$networkType

if (!useObs %in% c("all.obs", "pairwise.complete.obs")) {
  stop("wgcna.network.useObs 必须是 all.obs 或 pairwise.complete.obs")
}

power_fixed <- W$power$fixed %||% NULL
powers <- as.numeric(W$power$candidates)

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

  r2_target <- as.numeric(W$power$scale_free_R2_target)
  min_k     <- as.numeric(W$power$min_mean_connectivity)

  if (!all(c("Power", "SFT.R.sq", "mean.k.") %in% colnames(fit))) {
    chosen_power <- fit$Power[which.max(fit$SFT.R.sq)]
  } else {
    ok <- which(fit$SFT.R.sq >= r2_target & fit$mean.k. >= min_k)
    if (length(ok) > 0) {
      chosen_power <- fit$Power[min(ok)]
    } else {
      top <- order(fit$SFT.R.sq, decreasing = TRUE)[1:min(5, nrow(fit))]
      chosen_power <- fit$Power[top[which.max(fit$mean.k.[top])]]
    }
  }
}

log_msg(logger, "INFO", paste0("Chosen power: ", chosen_power))
log_msg(logger, "INFO", "Running WGCNA blockwiseModules...")

net <- WGCNA::blockwiseModules(
  datExpr,
  power = chosen_power,
  TOMType = W$network$TOMType,
  networkType = netType,
  corType = corType,
  maxPOutliers = ifelse(corType == "bicor", maxPOut, 1),
  minModuleSize = as.integer(W$modules$minModuleSize),
  deepSplit = as.integer(W$modules$deepSplit),
  pamRespectsDendro = isTRUE(W$modules$pamRespectsDendro),
  mergeCutHeight = as.numeric(W$modules$mergeCutHeight),
  reassignThreshold = as.numeric(W$modules$reassignThreshold),
  numericLabels = isTRUE(W$modules$numericLabels),
  saveTOMs = FALSE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- WGCNA::labels2colors(moduleLabels)

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
    writeLines(genes, file.path(SUBDIRS$lists, paste0("module_", mc, ".list")), useBytes = TRUE)
  }
}

saveRDS(net$dendrograms, file.path(SUBDIRS$modules, "gene_dendrograms.rds"))
saveRDS(net$blockGenes,  file.path(SUBDIRS$modules, "blockGenes.rds"))

MEs0 <- WGCNA::moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs  <- WGCNA::orderMEs(MEs0)

data.table::fwrite(
  data.frame(sample = rownames(MEs), MEs, check.names = FALSE),
  file.path(SUBDIRS$modules, "module_eigengenes.tsv"),
  sep = "\t",
  quote = FALSE
)

mt_cor <- NULL
mt_p <- NULL
mt_padj <- NULL
auto_primary_trait_name <- NA_character_

if (!is.null(traits_mat) && ncol(traits_mat) > 0) {
  if (is.null(rownames(traits_mat))) stop("traits_mat 缺少样本行名。")

  miss_in_traits <- setdiff(rownames(MEs), rownames(traits_mat))
  if (length(miss_in_traits) > 0) {
    stop("traits_mat 缺少这些样本：", paste(miss_in_traits, collapse = ", "))
  }

  traits_mat_aligned <- traits_mat[rownames(MEs), , drop = FALSE]

  mt_cor <- safe_cor(MEs, traits_mat_aligned, corType = corType, useObs = useObs, maxPOut = maxPOut)
  mt_p   <- WGCNA::corPvalueStudent(mt_cor, nSamples = nrow(datExpr))
  mt_padj <- matrix(
    p.adjust(as.vector(mt_p), method = "BH"),
    nrow = nrow(mt_p), ncol = ncol(mt_p),
    dimnames = dimnames(mt_p)
  )

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

  if (isTRUE(W$traits$redundancy$enable %||% TRUE) && ncol(traits_mat_aligned) >= 2) {
    tt_cor <- safe_cor(traits_mat_aligned, traits_mat_aligned, corType = corType, useObs = useObs, maxPOut = maxPOut)
    tt_p <- WGCNA::corPvalueStudent(tt_cor, nSamples = nrow(traits_mat_aligned))
    diag(tt_p) <- 0
    tt_padj <- matrix(
      p.adjust(as.vector(tt_p), method = as.character(W$traits$redundancy$p_adjust_method %||% 'BH')),
      nrow = nrow(tt_p), ncol = ncol(tt_p), dimnames = dimnames(tt_p)
    )
    write_matrix_with_rownames(tt_cor, file.path(SUBDIRS$assoc, 'trait_trait_cor.tsv'), rowname_col = 'trait')
    write_matrix_with_rownames(tt_padj, file.path(SUBDIRS$assoc, 'trait_trait_padj.tsv'), rowname_col = 'trait')

    thr <- as.numeric(W$traits$redundancy$cor_threshold %||% 0.85)
    upper_idx <- which(upper.tri(tt_cor) & abs(tt_cor) >= thr, arr.ind = TRUE)
    if (nrow(upper_idx) > 0) {
      redun_df <- data.frame(
        trait_a = rownames(tt_cor)[upper_idx[, 1]],
        trait_b = colnames(tt_cor)[upper_idx[, 2]],
        correlation = tt_cor[upper_idx],
        padj = tt_padj[upper_idx],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      redun_df <- data.frame(trait_a = character(), trait_b = character(), correlation = numeric(), padj = numeric(), check.names = FALSE)
    }
    data.table::fwrite(redun_df, file.path(SUBDIRS$assoc, 'trait_redundancy_pairs.tsv'), sep = '\t', quote = FALSE)
  }

  auto_primary_trait_name <- pick_primary_trait_auto(mt_cor, mt_padj, trait_meta_df,
                                                     W$export$hub_genes$primary_trait_selection,
                                                     logger = logger)
}

if (!is.null(traits_mat) && ncol(traits_mat) > 0 && isTRUE(W$export$hub_genes$enable)) {
  if (is.null(rownames(traits_mat))) stop("traits_mat 缺少样本行名。")

  miss_in_traits_hub <- setdiff(rownames(datExpr), rownames(traits_mat))
  if (length(miss_in_traits_hub) > 0) {
    stop("traits_mat 缺少这些 datExpr 样本：", paste(miss_in_traits_hub, collapse = ", "))
  }
  traits_mat_hub <- traits_mat[rownames(datExpr), , drop = FALSE]
  trait_names <- colnames(traits_mat_hub)

  pick_primary_trait <- function() {
    cfg_primary <- as.character(W$export$hub_genes$primary_trait %||% "")
    if (nzchar(cfg_primary)) {
      if (cfg_primary %in% trait_names) return(cfg_primary)
      stop("配置的 primary_trait 不在 trait_matrix 中：", cfg_primary)
    }

    if (!is.na(auto_primary_trait_name) && nzchar(auto_primary_trait_name) && auto_primary_trait_name %in% trait_names) {
      return(auto_primary_trait_name)
    }

    candidates <- as_char_vec(W$export$hub_genes$primary_trait_candidates)
    for (nm in candidates) {
      if (nm %in% trait_names) return(nm)
    }

    return(trait_names[1])
  }

  primary_trait_name <- pick_primary_trait()
  log_msg(logger, "INFO", paste0("Primary trait for hub ranking: ", primary_trait_name))

  trait_vec <- as.numeric(traits_mat_hub[[primary_trait_name]])
  names(trait_vec) <- rownames(traits_mat_hub)

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
  kME_df$gene_id <- rownames(kME_df)
  rownames(kME_df) <- kME_df$gene_id

  hub_base <- gene_module_df
  hub_base$GS <- as.numeric(GS[hub_base$gene_id])
  hub_base$GS_p <- as.numeric(GS_p[hub_base$gene_id])

  hub_df <- cbind(hub_base, kME_df[hub_base$gene_id, setdiff(colnames(kME_df), "gene_id"), drop = FALSE])
  data.table::fwrite(hub_df, file.path(SUBDIRS$hub, "gene_level_GS_kME.tsv"), sep = "\t", quote = FALSE)

  top_n <- as.integer(W$export$hub_genes$top_n_per_module)
  rank_by <- as_char_vec(W$export$hub_genes$rank_by)

  map_color_to_ME <- function(color) paste0("ME", color)
  get_kME_col_for_module <- function(me_name) paste0("kME_", me_name)

  mods <- sort(unique(hub_df$module_color))
  mods <- mods[mods != "grey"]
  hub_list_all <- list()

  for (mc in mods) {
    sub <- hub_df[hub_df$module_color == mc, , drop = FALSE]
    me_col <- get_kME_col_for_module(map_color_to_ME(mc))
    if (!me_col %in% colnames(sub)) next
    if (nrow(sub) == 0) next

    score_kME <- abs(sub[[me_col]])
    score_GS  <- abs(sub$GS)

    if (length(rank_by) == 0) {
      sub <- sub[order(score_kME, score_GS, decreasing = TRUE), , drop = FALSE]
    } else if (length(rank_by) == 1) {
      if (rank_by[1] == "GS") {
        sub <- sub[order(score_GS, score_kME, decreasing = TRUE), , drop = FALSE]
      } else {
        sub <- sub[order(score_kME, score_GS, decreasing = TRUE), , drop = FALSE]
      }
    } else {
      key1 <- if (rank_by[1] == "GS") score_GS else score_kME
      key2 <- if (rank_by[2] == "GS") score_GS else score_kME
      sub <- sub[order(key1, key2, decreasing = TRUE), , drop = FALSE]
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
  }
}

used_params <- list(
  config_yaml = CONFIG_YAML,
  wgcna = W,
  matched_samples = rownames(datExpr),
  n_samples = nrow(datExpr),
  n_genes_final = ncol(datExpr),
  chosen_power = chosen_power,
  method = method,
  threads = W$threads,
  auto_primary_trait = auto_primary_trait_name
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

log_msg(logger, "INFO", "DONE: WGCNA finished successfully.")

