#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# go_local.R —— 离线导出 GO ID → 英文名称（term name）字典（只产原材料 TSV，不画图）
# 运行方式：Rscript scripts/go_local.R
# 产物：ref/annotations/go_terms.tsv（两列：GO, name）
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
    anno_dir = "ref/annotations",                # go_terms.tsv 输出目录
    gene2go  = "ref/annotations/gene2go.tsv",    # 由 07_prepare_emapper_annotations.py 产出
    logs_dir = "logs"                             # 日志目录（固定文件名：go_local.log）
  ),
  go_local = list(
    export_all = FALSE                            # TRUE：导出 GO.db 全量；FALSE：仅导出 gene2go 中出现的 GO
  )
)

# =========================
# 2) 小工具：路径/合并/切块
# =========================
`%||%` <- function(a, b) if (is.null(a)) b else a

norm_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

merge_list_deep <- function(x, y) {
  for (n in names(y)) {
    if (!is.list(y[[n]]) || is.null(x[[n]])) x[[n]] <- y[[n]] else x[[n]] <- merge_list_deep(x[[n]], y[[n]])
  }
  x
}

chunk_select <- function(ids, chunk = 5000L) {
  n <- length(ids)
  if (n == 0) return(list())
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
  lapply(idx, function(i) ids[i])
}

extract_go_from_df <- function(df) {
  all_txt <- df %>%
    mutate(across(everything(), as.character)) %>%
    unite("__tmp__", everything(), sep = "||", remove = FALSE) %>%
    pull("__tmp__")
  unique(str_extract_all(all_txt, "GO:\\d{7}") %>% unlist())
}

# =========================
# 3) 读取 config.yaml 并合并参数
# =========================
cfg_path <- "config.yaml"
if (!file.exists(cfg_path)) stop("缺少配置文件：config.yaml（请在项目根目录运行或复制配置文件）")
cfg_file <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
if (is.null(cfg_file)) stop("config.yaml 解析失败")
CFG <- merge_list_deep(DEFAULTS, cfg_file)

anno_dir <- CFG$paths$anno_dir %||% DEFAULTS$paths$anno_dir
gene2go  <- CFG$paths$gene2go  %||% DEFAULTS$paths$gene2go
logs_dir <- CFG$paths$logs_dir %||% DEFAULTS$paths$logs_dir
export_all <- isTRUE(CFG$go_local$export_all %||% DEFAULTS$go_local$export_all)

# 统一绝对路径，保证稳定落点
anno_dir <- norm_path(file.path(getwd(), anno_dir))
gene2go  <- norm_path(file.path(getwd(), gene2go))
logs_dir <- norm_path(file.path(getwd(), logs_dir))

dir.create(anno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 4) 日志系统（固定文件名 + 流式屏幕输出）
# =========================
log_path <- file.path(logs_dir, "go_local.log")
.log_con <- file(log_path, open = "at")  # 追加写入
on.exit({ try(close(.log_con), silent = TRUE) }, add = TRUE)

.now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

.log_write <- function(level, msg) {
  line <- sprintf("[%s] [%s] %s", .now(), level, msg)
  # 写日志文件
  try(writeLines(text = line, con = .log_con, sep = "\n"), silent = TRUE)
  # 屏幕实时输出
  cat(line, "\n", sep = "")
}

log_info  <- function(...)  .log_write("INFO",  paste0(...))
log_warn  <- function(...)  .log_write("WARN",  paste0(...))
log_error <- function(...)  .log_write("ERROR", paste0(...))

log_info("go_local 启动；日志文件：", norm_path(log_path))
log_info("使用配置：", cfg_path)
log_info("anno_dir = ", anno_dir)
log_info("gene2go  = ", gene2go)
log_info("logs_dir = ", logs_dir)
log_info("export_all = ", export_all)

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
    log_error("未找到 gene2go 文件：", gene2go, "（先运行 07_prepare_emapper_annotations.py 产出 gene2go.tsv）")
    quit(status = 2)
  }
  g2g <- suppressWarnings(readr::read_tsv(gene2go, col_types = readr::cols(.default = readr::col_character())))
  if (!nrow(g2g)) {
    log_error("gene2go.tsv 为空：", gene2go)
    quit(status = 2)
  }
  cand <- character()
  if ("GO" %in% names(g2g)) {
    cand <- str_extract_all(paste(g2g$GO, collapse = ";"), "GO:\\d{7}") %>% unlist()
  } else {
    log_warn("gene2go.tsv 未检测到 'GO' 列，尝试从整表内容抽取 GO:xxxxxxx")
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
chunks <- chunk_select(go_ids, chunk = 5000L)

for (i in seq_along(chunks)) {
  keys_i <- chunks[[i]]
  log_info(sprintf("查询块 %d/%d，size=%d", i, length(chunks), length(keys_i)))
  qi <- tryCatch(
    AnnotationDbi::select(GO.db, keys = keys_i, keytype = "GOID", columns = c("TERM")),
    error = function(e) {
      log_warn("AnnotationDbi::select 出错，本块跳过：", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(qi) && nrow(qi)) {
    qi <- qi %>% transmute(GO = .data$GOID, name = .data$TERM %||% "")
    res <- bind_rows(res, qi)
  }
}

# 去重并按输入顺序排列
res <- res %>%
  distinct(GO, .keep_all = TRUE) %>%
  mutate(GO = as.character(GO)) %>%
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
# 7) 写出
# =========================
out_fp <- file.path(anno_dir, "go_terms.tsv")
tryCatch({
  readr::write_tsv(res, out_fp)
  log_info("写出：", norm_path(out_fp), " （", nrow(res), " 行）")
  log_info("DONE: go_local 完成（模式：", ifelse(export_all, "export_all=TRUE", "严格模式"), "）")
}, error = function(e) {
  log_error("写出失败：", conditionMessage(e))
  quit(status = 2)
})

