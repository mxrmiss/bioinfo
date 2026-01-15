#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# =============================================================================
# 08b_any_vs_other.R
#
# 在 08a_enrich_module.py 写完 base tasks.tsv 之后、调用 08_g_enrich.R 之前运行：
#   - 生成 one-vs-rest（任意 target 组织 vs 其它所有组织合并）的输入文件：
#       results/08_enrich/<target>_vs_other/{all.list,up.list,down.list,universe.list,gsea_ranks.tsv,meta.tsv,deg.tsv}
#   - 并把任务合并进 results/08_enrich/tasks.tsv（不改变原有任务，只增加/更新 *_vs_other 的 label 行）
#
# 配置读取优先级（皇上硬要求）：
#   1) 优先从 config.yaml 的 one_vs_rest 块读取；
#   2) 若 config.yaml 中不存在 one_vs_rest（或关键字段缺失），才使用本脚本顶部“用户自定义参数区”的默认值。
#
# 硬约束（皇上需求清单）：
# - 仅在 config 开启时运行（enable=true）；否则正常退出（status=0，不算错）。
# - 全程相对路径；不生成任何 .bak_* 文件（严禁）。
# - tasks.tsv 合并必须“幂等”：重复运行不会无限追加同一 label 的 all/up/down 三行。
# - GSEA 支撑：必须准备 gsea_ranks.tsv 等输入。
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

# --------------------------- 用户自定义参数区（仅在 config 缺失 one_vs_rest 时启用） ---------------------------

USER_ENABLE <- FALSE
USER_TARGETS <- c("foot")          # 目标组织名：必须能匹配 data/samples.tsv 的 group 值（建议小写）
USER_LABEL_SUFFIX <- "_vs_other"   # 输出 label 后缀
USER_REPLACE_EXISTING_LABEL <- TRUE

# 是否加入 block（paired/批次）控制：design = ~ block + group2
USER_USE_BLOCK <- FALSE
USER_BLOCK_COLUMN <- ""           # 例如 "batch" 或 "individual"；留空表示自动猜（batch/individual）

# DEG 阈值（用于 up/down list 供 ORA）
USER_DEG_LFC <- 1.0
USER_DEG_FDR <- 0.05

# GSEA 排名分数来源：stat 或 log2fc
USER_GSEA_SCORE_FROM <- "stat"

# 最少样本数（target 与 other）
USER_MIN_SAMPLES_PER_GROUP <- 3

# 低表达预过滤（加速/稳定）：rowSums(counts) >= 10
USER_MIN_ROW_SUM_COUNTS <- 10

# -----------------------------------------------------------------------------------------------------------

log_msg <- function(level = "INFO", ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(ts, "[08b_any_vs_other]", sprintf("[%s]", level), paste0(...), "\n")
}

stop_with <- function(...) {
  log_msg("ERROR", ...)
  stop(paste0(...), call. = FALSE)
}

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

cfg_get <- function(cfg, keys, default = NULL) {
  cur <- cfg
  for (k in keys) {
    if (is.null(cur) || !is.list(cur) || is.null(cur[[k]])) return(default)
    cur <- cur[[k]]
  }
  cur
}

tolower_trim <- function(x) tolower(trimws(as.character(x)))

# --------------------------- 读取 config.yaml ---------------------------

cfg_path <- "config.yaml"
cfg <- list()
if (file.exists(cfg_path)) {
  cfg <- yaml::read_yaml(cfg_path)
} else {
  log_msg("WARNING", "未找到 config.yaml，将仅使用脚本顶部用户参数区。")
  cfg <- list()
}

# dirs
dir_matrix <- cfg_get(cfg, c("dirs", "matrix"), "results/05_matrix")
dir_annot  <- cfg_get(cfg, c("dirs", "annot"),  "results/07_annot")
dir_enrich <- cfg_get(cfg, c("dirs", "enrich"), "results/08_enrich")

# data
samples_tsv <- cfg_get(cfg, c("data", "samples_tsv"), "data/samples.tsv")

# 矩阵文件（沿用项目约定）
counts_fp <- file.path(dir_matrix, "counts", "gene_counts.tsv")
tpm_fp    <- file.path(dir_matrix, "tpms", "gene_tpm.tsv")

# 注释宇宙文件（沿用项目约定）
gene2go_fp   <- file.path(dir_annot, "gene2go.tsv")
gene2path_fp <- file.path(dir_annot, "gene2pathway.tsv")

# tasks.tsv
tasks_fp <- file.path(dir_enrich, "tasks.tsv")

# inputs 目录（外包专用）：只确保存在，不覆盖（皇上要求）
inputs_dir <- cfg_get(cfg, c("background", "ora_inputs_dir"), file.path(dir_enrich, "inputs"))

# one_vs_rest 块（核心：若缺失则用 USER_*）
one_cfg <- cfg_get(cfg, c("one_vs_rest"), NULL)

enable <- if (!is.null(one_cfg) && !is.null(one_cfg$enable)) isTRUE(one_cfg$enable) else isTRUE(USER_ENABLE)
targets <- if (!is.null(one_cfg) && !is.null(one_cfg$targets)) as.character(one_cfg$targets) else as.character(USER_TARGETS)
label_suffix <- if (!is.null(one_cfg) && !is.null(one_cfg$label_suffix)) as.character(one_cfg$label_suffix) else as.character(USER_LABEL_SUFFIX)
replace_existing <- if (!is.null(one_cfg) && !is.null(one_cfg$replace_existing_label)) isTRUE(one_cfg$replace_existing_label) else isTRUE(USER_REPLACE_EXISTING_LABEL)

use_block <- if (!is.null(one_cfg) && !is.null(one_cfg$use_block)) isTRUE(one_cfg$use_block) else isTRUE(USER_USE_BLOCK)
block_col <- if (!is.null(one_cfg) && !is.null(one_cfg$block_column)) as.character(one_cfg$block_column) else as.character(USER_BLOCK_COLUMN)

deg_lfc <- cfg_get(cfg, c("deg", "lfc"), USER_DEG_LFC)
deg_fdr <- cfg_get(cfg, c("deg", "fdr"), USER_DEG_FDR)

gsea_score_from <- cfg_get(cfg, c("gsea", "score_from"), USER_GSEA_SCORE_FROM)
gsea_score_from <- tolower_trim(gsea_score_from)

min_samples <- cfg_get(cfg, c("deg", "min_samples_per_group"), USER_MIN_SAMPLES_PER_GROUP)
min_row_sum <- cfg_get(cfg, c("one_vs_rest", "min_row_sum_counts"), USER_MIN_ROW_SUM_COUNTS)
if (is.null(min_row_sum) || is.na(min_row_sum)) min_row_sum <- USER_MIN_ROW_SUM_COUNTS

targets <- tolower_trim(targets)
targets <- targets[nzchar(targets)]
targets <- unique(targets)

log_msg("INFO", "enable=", enable)
log_msg("INFO", "targets=", if (length(targets) > 0) paste(targets, collapse = ",") else "NA")
log_msg("INFO", "label_suffix=", label_suffix)
log_msg("INFO", "deg_lfc=", deg_lfc, " deg_fdr=", deg_fdr, " gsea_score_from=", gsea_score_from)
log_msg("INFO", "use_block=", use_block, " block_column=", if (nzchar(block_col)) block_col else "(auto/none)")
log_msg("INFO", "min_samples_per_group=", min_samples, " min_row_sum_counts=", min_row_sum)

# 若未启用，则直接退出（给 08a 自动调度用：不算错误）
if (!isTRUE(enable)) {
  log_msg("INFO", "one_vs_rest 未启用（enable=false），08b 直接退出。")
  quit(save = "no", status = 0)
}
if (length(targets) == 0) {
  log_msg("WARNING", "one_vs_rest 已启用但 targets 为空，08b 直接退出。")
  quit(save = "no", status = 0)
}

# --------------------------- 运行前检查 + 建目录 ---------------------------

if (!file.exists(samples_tsv)) stop_with("找不到 samples.tsv：", samples_tsv)
if (!file.exists(tasks_fp))    stop_with("找不到 tasks.tsv（必须先由 08a_enrich_module.py 生成）：", tasks_fp)
if (!file.exists(counts_fp))   stop_with("找不到 gene_counts.tsv：", counts_fp)

# inputs 文件夹：不存在就创建，存在不覆盖（皇上要求）
ensure_dir(inputs_dir)
log_msg("INFO", "已检查/创建 inputs 文件夹（存在不覆盖）：", inputs_dir)

# 注释宇宙：缺一个会导致背景覆盖率不可信，但仍可运行（按你的项目常规应当都存在）
if (!file.exists(gene2go_fp))   log_msg("WARNING", "缺少 gene2go.tsv：", gene2go_fp)
if (!file.exists(gene2path_fp)) log_msg("WARNING", "缺少 gene2pathway.tsv：", gene2path_fp)

# 需要 DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop_with("缺少 R 包 DESeq2（Bioconductor）。请先安装后再运行 08b。")
}
suppressPackageStartupMessages(library(DESeq2))

# --------------------------- 读入 samples / counts / tpm（tpm 可选） ---------------------------

smp <- data.table::fread(samples_tsv, sep = "\t", header = TRUE)
if (!all(c("sample", "group") %in% colnames(smp))) {
  stop_with("samples.tsv 必须包含列：sample, group；实际列：", paste(colnames(smp), collapse = ","))
}
smp[, sample := as.character(sample)]
smp[, group  := tolower_trim(group)]

# block 列自动猜
guess_block_col <- function(dt) {
  for (cand in c("batch", "individual", "ind", "Batch", "Individual", "Ind")) {
    if (cand %in% colnames(dt)) return(cand)
  }
  return(NA_character_)
}
if (isTRUE(use_block) && !nzchar(block_col)) {
  bc <- guess_block_col(smp)
  if (!is.na(bc)) block_col <- bc
}

cnt_dt <- data.table::fread(counts_fp, sep = "\t", header = TRUE)
gene_col <- if ("gene_id" %in% colnames(cnt_dt)) "gene_id" else colnames(cnt_dt)[1]
cnt_dt[, gene_id := as.character(get(gene_col))]
if (gene_col != "gene_id") cnt_dt[, (gene_col) := NULL]

count_samples <- setdiff(colnames(cnt_dt), "gene_id")
keep_samples <- intersect(smp$sample, count_samples)
if (length(keep_samples) < 4) {
  stop_with("counts 与 samples 匹配样本过少：matched=", length(keep_samples),
            "；请检查 samples.tsv 的 sample 值是否与 gene_counts.tsv 列名一致。")
}
smp2 <- smp[sample %in% keep_samples]
cnt2 <- cnt_dt[, c("gene_id", smp2$sample), with = FALSE]

# tpm 可选：若不存在则当作全 0（仍可跑，只是 detectable 规则退化为 counts>0）
tpm_dt <- NULL
if (file.exists(tpm_fp)) {
  tpm_dt <- data.table::fread(tpm_fp, sep = "\t", header = TRUE)
  tpm_gene_col <- if ("gene_id" %in% colnames(tpm_dt)) "gene_id" else colnames(tpm_dt)[1]
  tpm_dt[, gene_id := as.character(get(tpm_gene_col))]
  if (tpm_gene_col != "gene_id") tpm_dt[, (tpm_gene_col) := NULL]
  tpm_keep <- intersect(colnames(tpm_dt), c("gene_id", smp2$sample))
  tpm_dt <- tpm_dt[, tpm_keep, with = FALSE]
} else {
  log_msg("WARNING", "未找到 gene_tpm.tsv，将仅使用 counts>0 作为 detectable 规则：", tpm_fp)
}

# --------------------------- 注释宇宙 gene_universe = gene2go ∪ gene2pathway ---------------------------

annot_universe <- character()
if (file.exists(gene2go_fp)) {
  g2g <- data.table::fread(gene2go_fp, sep = "\t", header = TRUE)
  if ("gene_id" %in% colnames(g2g)) annot_universe <- c(annot_universe, as.character(g2g$gene_id))
}
if (file.exists(gene2path_fp)) {
  g2k <- data.table::fread(gene2path_fp, sep = "\t", header = TRUE)
  if ("gene_id" %in% colnames(g2k)) annot_universe <- c(annot_universe, as.character(g2k$gene_id))
}
annot_universe <- unique(annot_universe)
log_msg("INFO", "annot_universe size=", length(annot_universe))

# --------------------------- tasks.tsv 读入（最后写回） ---------------------------

tasks <- data.table::fread(tasks_fp, sep = "\t", header = TRUE)
need_cols <- c("label","type","test_set","test_file","background_file","universe_meta","outdir","n_deg_all")
miss <- setdiff(need_cols, colnames(tasks))
if (length(miss) > 0) {
  stop_with("tasks.tsv 缺少必要列：", paste(miss, collapse = ", "))
}

# --------------------------- 工具函数 ---------------------------

write_list <- function(fp, vec) {
  dt <- data.table(gene_id = as.character(vec))
  fwrite(dt, fp, sep = "\t", quote = FALSE, na = "NA")
}

calc_detectable_and_universe <- function(count_mat, tpm_mat, gene_ids, annot_univ) {
  # detectable：任一样本 counts>0 或 tpm>0
  if (is.null(tpm_mat)) {
    det_flag <- rowSums(count_mat > 0) > 0
  } else {
    det_flag <- (rowSums(count_mat > 0) > 0) | (rowSums(tpm_mat > 0) > 0)
  }
  detectable <- gene_ids[det_flag]
  universe <- intersect(detectable, annot_univ)
  list(detectable = detectable, universe = universe)
}

# 幂等合并：按 (label, test_set) 去重
# - replace_existing=TRUE：先删该 label 的 all/up/down，再加入新行
# - replace_existing=FALSE：若该 label 的某个 test_set 已存在，则不再追加（避免重复）
merge_tasks_idempotent <- function(tasks_dt, new_rows_dt, replace_existing) {
  stopifnot(all(c("label","test_set") %in% colnames(tasks_dt)))
  stopifnot(all(c("label","test_set") %in% colnames(new_rows_dt)))

  lab <- unique(new_rows_dt$label)
  if (length(lab) != 1) stop_with("内部错误：new_rows_dt label 不唯一：", paste(lab, collapse=","))

  target_sets <- unique(new_rows_dt$test_set)

  if (isTRUE(replace_existing)) {
    # 删除该 label 的目标 test_set 行，再追加新行（替换语义）
    tasks_dt <- tasks_dt[!(label == lab & test_set %in% target_sets)]
    tasks_dt <- rbindlist(list(tasks_dt, new_rows_dt), use.names = TRUE, fill = TRUE)
    return(tasks_dt)
  }

  # replace_existing=FALSE：只补齐缺失的 test_set 行（保证幂等，不重复追加）
  exist_key <- paste(tasks_dt$label, tasks_dt$test_set, sep = "||")
  new_key <- paste(new_rows_dt$label, new_rows_dt$test_set, sep = "||")
  keep_new <- !(new_key %in% exist_key)
  if (any(keep_new)) {
    tasks_dt <- rbindlist(list(tasks_dt, new_rows_dt[keep_new]), use.names = TRUE, fill = TRUE)
  }
  return(tasks_dt)
}

# --------------------------- 主循环：每个 target 做 one-vs-rest ---------------------------

append_rows_all <- list()

for (target in targets) {
  label <- paste0(target, label_suffix)
  outdir_label <- file.path(dir_enrich, label)
  ensure_dir(outdir_label)

  # 构建 group2：target vs other
  smp2[, group2 := ifelse(group == target, target, "other")]

  n_t <- sum(smp2$group2 == target)
  n_o <- sum(smp2$group2 == "other")

  log_msg("INFO", "Target=", target, " label=", label, " n_target=", n_t, " n_other=", n_o)

  if (n_t < min_samples || n_o < min_samples) {
    log_msg("WARNING", "样本数不足，跳过：target=", target, "（需要每组至少 ", min_samples, "）")
    next
  }

  # counts 矩阵（gene x sample）
  count_mat <- as.matrix(cnt2[, -1, with = FALSE])
  rownames(count_mat) <- cnt2$gene_id
  suppressWarnings(storage.mode(count_mat) <- "numeric")
  count_mat[is.na(count_mat)] <- 0

  # 预过滤：rowSums(counts) >= min_row_sum
  keep_gene <- rowSums(count_mat) >= as.numeric(min_row_sum)
  count_mat <- count_mat[keep_gene, , drop = FALSE]
  gene_ids <- rownames(count_mat)

  # tpm 矩阵（可选）
  tpm_mat <- NULL
  if (!is.null(tpm_dt)) {
    tpm_dt2 <- tpm_dt[gene_id %in% gene_ids, ]
    setkey(tpm_dt2, gene_id)
    tpm_dt2 <- tpm_dt2[gene_ids]
    # 用 smp2$sample 顺序对齐（与 count_mat 列顺序一致）
    tpm_mat <- as.matrix(tpm_dt2[, smp2$sample, with = FALSE])
    rownames(tpm_mat) <- tpm_dt2$gene_id
    suppressWarnings(storage.mode(tpm_mat) <- "numeric")
    tpm_mat[is.na(tpm_mat)] <- 0
  }

  # detectable + universe
  du <- calc_detectable_and_universe(
    count_mat = count_mat,
    tpm_mat = tpm_mat,
    gene_ids = gene_ids,
    annot_univ = annot_universe
  )
  detectable <- du$detectable
  universe <- du$universe

  universe_fp <- file.path(outdir_label, "universe.list")
  write_list(universe_fp, sort(unique(universe)))

  coverage <- if (length(detectable) > 0) length(universe) / length(detectable) else 0

  meta_fp <- file.path(outdir_label, "meta.tsv")
  meta_dt <- data.table(
    label = label,
    n_detectable = length(detectable),
    n_annot_mapped = length(universe),
    universe_size = length(universe),
    coverage = as.numeric(sprintf("%.4f", coverage)),
    detectable_rule = cfg_get(cfg, c("background","detectable"), "TPM>0_or_count>0_in>=1_sample"),
    samples_used = paste(smp2$sample, collapse = ",")
  )
  fwrite(meta_dt, meta_fp, sep = "\t", quote = FALSE, na = "NA")

  # DESeq2：target vs other
  coldata <- data.frame(
    row.names = smp2$sample,
    group2 = factor(smp2$group2, levels = c("other", target))
  )

  design_formula <- ~ group2
  if (isTRUE(use_block)) {
    if (!nzchar(block_col) || !(block_col %in% colnames(smp2))) {
      log_msg("WARNING", "use_block=TRUE 但找不到 block 列，回退为 ~ group2。block_col=", block_col)
    } else {
      coldata$block <- factor(as.character(smp2[[block_col]]))
      design_formula <- ~ block + group2
      log_msg("INFO", "DESeq2 design: ~ block + group2 （block_col=", block_col, "）")
    }
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat[, smp2$sample, drop = FALSE]),
    colData = coldata,
    design = design_formula
  )

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group2", target, "other"))

  res_dt <- as.data.table(res)
  res_dt[, gene_id := rownames(res)]
  setcolorder(res_dt, c("gene_id", setdiff(colnames(res_dt), "gene_id")))

  deg_fp <- file.path(outdir_label, "deg.tsv")
  fwrite(res_dt, deg_fp, sep = "\t", quote = FALSE, na = "NA")

  # up/down/all lists（供 ORA）
  up <- res_dt[!is.na(padj) & padj <= deg_fdr & !is.na(log2FoldChange) & log2FoldChange >= deg_lfc, gene_id]
  down <- res_dt[!is.na(padj) & padj <= deg_fdr & !is.na(log2FoldChange) & log2FoldChange <= -deg_lfc, gene_id]
  up <- sort(unique(as.character(up)))
  down <- sort(unique(as.character(down)))
  all_deg <- sort(unique(c(up, down)))

  write_list(file.path(outdir_label, "up.list"), up)
  write_list(file.path(outdir_label, "down.list"), down)
  write_list(file.path(outdir_label, "all.list"), all_deg)

  n_deg_all <- length(all_deg)
  log_msg("INFO", "DEG counts: label=", label, " all=", n_deg_all, " up=", length(up), " down=", length(down))

  # gsea_ranks.tsv（gene_id, score, rank_source）
  score <- NULL
  rank_source <- NULL
  if (gsea_score_from == "stat") {
    score <- res_dt$stat
    rank_source <- "DESeq2.stat"
  } else if (gsea_score_from %in% c("log2fc","log2foldchange","log2_fold_change")) {
    score <- res_dt$log2FoldChange
    rank_source <- "DESeq2.log2FoldChange"
  } else {
    score <- res_dt$stat
    rank_source <- "DESeq2.stat"
  }

  ranks_dt <- data.table(
    gene_id = res_dt$gene_id,
    score = as.numeric(score),
    rank_source = rank_source
  )
  ranks_dt <- ranks_dt[!is.na(score) & !is.na(gene_id) & nzchar(gene_id)]
  setorder(ranks_dt, -score)

  fwrite(ranks_dt, file.path(outdir_label, "gsea_ranks.tsv"),
         sep = "\t", quote = FALSE, na = "NA")

  # 生成 tasks 新行（保持原表头）
  new_rows <- data.table(
    label = label,
    type = "rna",
    test_set = c("all","up","down"),
    test_file = c(file.path(outdir_label,"all.list"),
                  file.path(outdir_label,"up.list"),
                  file.path(outdir_label,"down.list")),
    background_file = universe_fp,
    universe_meta = meta_fp,
    outdir = outdir_label,
    n_deg_all = as.character(n_deg_all)
  )

  append_rows_all[[label]] <- new_rows
}

if (length(append_rows_all) == 0) {
  log_msg("WARNING", "未生成任何 one-vs-rest 任务（可能 targets 为空或样本数不足）。")
  quit(save = "no", status = 0)
}

# --------------------------- 合并写回 tasks.tsv（幂等、无备份） ---------------------------

for (lab in names(append_rows_all)) {
  tasks <- merge_tasks_idempotent(tasks, append_rows_all[[lab]], replace_existing = replace_existing)
}

# 最终再做一次严格去重：按 (label, test_set) 去重，保留第一条（避免不期望覆盖）
# 若 replace_existing=TRUE，上面已删除旧行再追加；此处主要防御异常重复。
setkey(tasks, label, test_set)
tasks <- unique(tasks, by = c("label", "test_set"))

fwrite(tasks, tasks_fp, sep = "\t", quote = FALSE, na = "NA")

log_msg("INFO", "tasks.tsv updated: merged labels = ", paste(names(append_rows_all), collapse = ", "))
log_msg("INFO", "===== 08b_any_vs_other.R done =====")

