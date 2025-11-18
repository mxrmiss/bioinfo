#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# go_local.R —— 离线导出 GO ID → 英文名称（term name）字典（只产原材料 TSV，不画图）
# 运行方式：Rscript scripts/go_local.R
# 产物：results/07_annot/go/term2name.tsv（两列：go_id, term_name）
# 日志：logs/go_local.log（固定文件名，无时间戳，追加写入）
# 依赖：yaml、readr、dplyr、stringr、tidyr、AnnotationDbi、GO.db
# 约定：所有参数集中在脚本顶部；若 config.yaml 存在则优先读取覆盖默认值

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(AnnotationDbi)
  library(GO.db)
})

# =========================
# 1) 顶部参数（可手动修改）
# =========================
DEFAULTS <- list(
  paths = list(
    anno_dir = "results/07_annot/go",              # [改动] 统一到约定目录
    gene2go  = "ref/annotations/gene2go.tsv",      # 由 07_prepare_emapper_annotations.py 产出（未改）
    logs_dir = "logs"                               # 日志目录（固定文件名：go_local.log）
  ),
  go_local = list(
    export_all = FALSE                              # TRUE：导出 GO.db 全量；FALSE：仅导出 gene2go 中出现的 GO
  )
)

# =========================
# 2) 读取 config.yaml（若存在则覆盖默认）
# =========================
cfg <- DEFAULTS
if (file.exists("config.yaml")) {
  y <- yaml::read_yaml("config.yaml")
  if (!is.null(y$paths))         cfg$paths         <- modifyList(cfg$paths,         y$paths)
  if (!is.null(y$go_local))      cfg$go_local      <- modifyList(cfg$go_local,      y$go_local)
}

anno_dir <- cfg$paths$anno_dir
gene2go  <- cfg$paths$gene2go
logs_dir <- cfg$paths$logs_dir
export_all <- isTRUE(cfg$go_local$export_all)

dir.create(anno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 3) 简易日志
# =========================
LOG <- file.path(logs_dir, "go_local.log")
.log <- function(level, ...) {
  cat(paste0("[", level, "] ", paste0(..., collapse = "")), "\n", file = LOG, append = TRUE)
  message(paste0("[", level, "] ", paste0(..., collapse = "")))
}
log_info  <- function(...) .log("INFO",  ...)
log_warn  <- function(...) .log("WARN",  ...)
log_error <- function(...) .log("ERR",   ...)

log_info("==== go_local 启动 ====")
log_info("anno_dir: ", anno_dir)
log_info("gene2go : ", gene2go)
log_info("export_all = ", export_all)

# =========================
# 4) 工具函数
# =========================
norm_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

extract_go_from_df <- function(df) {
  # 在整张 gene2go 表中抓取形如 GO:0000000 的 token
  tokens <- unlist(strsplit(paste0(unlist(df), collapse = " "), "[,;| \t]+"))
  tokens <- tokens[grepl("^GO:[0-9]{7}$", tokens)]
  unique(tokens)
}

# =========================
# 5) 构建 GO 号候选集
# =========================
go_ids <- character()

if (export_all) {
  log_info("export_all=TRUE：从 GO.db 读取全部 GOID")
  all_terms <- as.list(GO.db::GOTERM)
  go_ids <- names(all_terms)
  go_ids <- unique(go_ids[!is.na(go_ids)])
  log_info("全部 GOID 数量：", length(go_ids))
} else {
  log_info("严格模式：仅导出 gene2go.tsv 中出现的 GO 号")
  if (!file.exists(gene2go)) {
    log_error("未找到 gene2go.tsv：", gene2go)
    quit(status = 2)
  }
  g2g <- tryCatch(readr::read_tsv(gene2go, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(g2g)) {
    log_error("读取 gene2go.tsv 失败：", gene2go)
    quit(status = 2)
  }
  cand <- NULL
  if (all(c("gene_id","go_id") %in% colnames(g2g))) {
    cand <- unique(g2g$go_id)
  } else {
    log_warn("gene2go.tsv 未检测到 'go_id' 列，尝试从整表内容抽取 GO:xxxxxxx")
    cand <- extract_go_from_df(g2g)
  }
  go_ids <- unique(cand[!is.na(cand)])
  if (!length(go_ids)) {
    log_error("未在 gene2go.tsv 中解析到任何 GO:xxxxxxx；请检查文件格式与列名")
    quit(status = 2)
  }
  log_info("在 gene2go.tsv 中解析到 GOID 数量：", length(go_ids))
}

# =========================
# 6) 查询 GO.db → TERM
# =========================
log_info("开始从 GO.db 查询 TERM（分块稳态查询）")
res <- tibble(GO = character(), name = character())

chunk_select <- function(v, k = 5000L) {
  n <- length(v); if (n == 0) return(list())
  split(v, ceiling(seq_along(v) / k))
}
chunks <- chunk_select(go_ids, 5000L)

for (chs in chunks) {
  suppressWarnings({
    tt <- suppressMessages(select(GO.db, keys = chs, keytype = "GOID", columns = c("GOID", "TERM")))
  })
  if (is.null(tt) || !nrow(tt)) next
  qi <- tt %>% transmute(GO = .data$GOID, name = .data$TERM) %>% distinct(GO, .keep_all = TRUE)
  res <- bind_rows(res, qi)
}
res <- res %>%
  distinct(GO, .keep_all = TRUE) %>%
  arrange(factor(GO, levels = go_ids))

# 若 GO.db 中缺失某些 ID，补空名占位，保证键全集
if (nrow(res) < length(go_ids)) {
  lacks <- setdiff(go_ids, res$GO)
  if (length(lacks)) {
    msg <- paste0(
      "下列 GOID 在 GO.db 中未找到，将以空名称占位（仅显示前 10 个）：",
      paste(head(lacks, 10), collapse = ", "),
      ifelse(length(lacks) > 10, " ...", "")
    )
    log_warn(msg)
    res <- bind_rows(res, tibble(GO = lacks, name = ""))
    res <- res %>% arrange(factor(GO, levels = go_ids))
  }
}

# =========================
# 7) 写出（仅改表头与落点）
# =========================
out_fp <- file.path(anno_dir, "term2name.tsv")
# [改动] 写出前统一列名为约定口径（仅改表头，不改数据）
res <- res %>% dplyr::rename(go_id = GO, term_name = name)

tryCatch({
  readr::write_tsv(res, out_fp)
  log_info("写出：", norm_path(out_fp), " （", nrow(res), " 行）")
  log_info("DONE: go_local 完成（模式：", ifelse(export_all, "export_all=TRUE", "严格模式"), "）")
}, error = function(e) {
  log_error("写出失败：", conditionMessage(e))
  quit(status = 2)
})