#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# go_local.R —— 从 GO.db 离线导出 GO 编号对应的英文名称（仅生成名称字典）
# 设计要点：
#   1) 默认优先读取项目根目录下的 config.yaml（以配置为准），脚本内默认仅作后备。
#   2) 严格依赖 ref/annotations/gene2go.tsv 中实际出现的 GO 号；若未找到该文件则报错退出，不做兜底。
#   3) 输出两列表：GO, name 到 ref/annotations/go_terms.tsv，供 08 富集脚本拼接 term_name。
#   4) 离线工作：使用 GO.db + AnnotationDbi，不联网，不自动安装依赖。
#
# 依赖包：yaml、readr、dplyr、stringr、AnnotationDbi、GO.db
# 若缺包，请按提示手工安装（示例）：
#   install.packages("BiocManager")
#   BiocManager::install(c("GO.db","AnnotationDbi"))
#   install.packages(c("yaml","readr","dplyr","stringr"))

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(stringr)
})

# ============ 基础工具函数 ============
`%||%` <- function(a, b) if (is.null(a)) b else a

merge_list_deep <- function(a, b){
  # 递归合并列表：b 覆盖 a
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (!is.list(a) || !is.list(b)) return(b)
  keys <- union(names(a), names(b))
  out <- list()
  for (k in keys){
    out[[k]] <- merge_list_deep(a[[k]], b[[k]])
  }
  out
}

norm_path <- function(x, fallback){
  if (is.null(x) || !length(x)) return(fallback)
  x <- as.character(x[1])
  if (!nzchar(x)) return(fallback)
  x
}

ensure_dir <- function(p){
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}

timestamp_now <- function(){
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

# ============ 脚本默认参数（会被 config.yaml 覆盖） ============
DEFAULTS <- list(
  paths = list(
    anno_dir = "ref/annotations",
    logs_dir = "logs"
  ),
  go_local = list(
    export_all = FALSE # 严格模式：默认 FALSE 表示必须基于 gene2go.tsv 中出现的 GO 号导出
  )
)

# ============ 读取 config.yaml 并合并参数 ============
cfg_path <- "config.yaml"
if (!file.exists(cfg_path)) {
  stop("缺少配置文件：config.yaml（请在项目根目录运行）")
}
cfg_file <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
if (is.null(cfg_file)) stop("config.yaml 解析失败")

CFG <- merge_list_deep(DEFAULTS, cfg_file)

anno_dir <- norm_path(CFG$paths$anno_dir %||% "ref/annotations", "ref/annotations")
logs_dir <- norm_path(CFG$paths$logs_dir %||% "logs", "logs")
export_all <- isTRUE(CFG$go_local$export_all)

ensure_dir(anno_dir)
ensure_dir(logs_dir)

# 日志文件（同时镜像到屏幕）
log_fp <- file.path(logs_dir, sprintf("07b_go_local.%s.log", timestamp_now()))
zz <- file(log_fp, open = "wt", encoding = "UTF-8")
sink(zz, split = TRUE)
sink(zz, split = TRUE, type = "message")
on.exit({
  try(sink(type="message"), silent=TRUE)
  try(sink(), silent=TRUE)
  try(close(zz), silent=TRUE)
}, add = TRUE)

message("[INFO] go_local.R 启动")
message("[INFO] anno_dir = ", anno_dir)
message("[INFO] logs_dir = ", logs_dir)
message("[INFO] export_all = ", export_all)

# ============ 依赖检查 ============
pkgs_needed <- c("AnnotationDbi","GO.db")
pk_missing <- pkgs_needed[!sapply(pkgs_needed, requireNamespace, quietly = TRUE)]
if (length(pk_missing)) {
  stop(
    "缺少依赖包：", paste(pk_missing, collapse = ", "),
    "\n请手工安装，例如：\n",
    "  install.packages('BiocManager')\n",
    "  BiocManager::install(c('GO.db','AnnotationDbi'))\n"
  )
}
library(AnnotationDbi)
library(GO.db)

# ============ 读取 gene2go.tsv 并抽取 GO 号 ============
g2go_fp <- file.path(anno_dir, "gene2go.tsv")
if (!export_all) {
  if (!file.exists(g2go_fp)) {
    stop("严格模式：未发现 ", g2go_fp, "。请先运行 07_prepare_emapper_annotations.py 生成 gene2go.tsv，或将 go_local$export_all 设为 TRUE。")
  }
}

go_ids <- character()

if (!export_all) {
  message("[INFO] 读取：", g2go_fp)
  # 宽松读取（仅使用前两列，避免无关列名差异）
  g2go <- suppressMessages(readr::read_tsv(g2go_fp, col_types = cols(.default = "c")))
  if (!nrow(g2go)) stop("gene2go.tsv 为空：", g2go_fp)

  # 兼容列名：优先识别 GO 字段名；否则使用第二列
  if ("GO" %in% names(g2go)) {
    go_col <- "GO"
  } else {
    if (ncol(g2go) < 2) stop("gene2go.tsv 列数不足，至少需要两列（gene_id, GO）。")
    go_col <- names(g2go)[2]
    message("[WARN] gene2go.tsv 未发现名为 'GO' 的列，将使用第二列：", go_col)
  }

  go_ids <- unique(g2go[[go_col]])
  go_ids <- go_ids[!is.na(go_ids) & nzchar(go_ids)]
  # 仅保留形如 GO:0000001 的条目
  go_ids <- go_ids[str_detect(go_ids, "^GO:\\d{7}$")]
  go_ids <- unique(go_ids)

  if (!length(go_ids)) {
    stop("gene2go.tsv 中未解析到任何合法 GO 编号（形如 GO:0000001）。")
  }
  message("[INFO] 从 gene2go.tsv 解析到 GO 条目数：", length(go_ids))
} else {
  # 导出 GO.db 全部 GO→名称（可能体量较大）
  message("[INFO] export_all = TRUE，将导出 GO.db 中全部条目。")
  # 使用 GO.db 内置 key 列
  go_ids <- keys(GO.db, keytype = "GOID")
  go_ids <- unique(go_ids[str_detect(go_ids, "^GO:\\d{7}$")])
  message("[INFO] GO.db 中可用 GO 条目数：", length(go_ids))
  if (!length(go_ids)) stop("在 GO.db 中未找到任何 GO:xxxxxxx 条目。")
}

# ============ 查询 GO 名称 ============
# 使用 AnnotationDbi::select，从 GOID → TERM
chunk_select <- function(ids, chunk = 5000L){
  out <- list()
  n <- length(ids)
  if (n <= 0) return(tibble(GO = character(), name = character()))
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
  for (ii in idx){
    keys_part <- ids[ii]
    tbl <- AnnotationDbi::select(GO.db,
                                 keys = keys_part,
                                 columns = c("GOID","TERM"),
                                 keytype = "GOID")
    tbl <- as_tibble(tbl) %>%
      transmute(GO = as.character(GOID), name = as.character(TERM)) %>%
      distinct()
    out[[length(out)+1]] <- tbl
    message("[INFO] 查询进度：", max(ii), " / ", n)
  }
  bind_rows(out)
}

res <- chunk_select(go_ids, chunk = 5000L) %>%
  mutate(name = if_else(is.na(name), "", name)) %>%
  distinct(GO, .keep_all = TRUE) %>%
  arrange(GO)

# 若有缺名条目，填空字符串并提示
missing_n <- sum(!res$GO %in% go_ids) # 正常应为 0
if (missing_n > 0) {
  message("[WARN] 检测到 ", missing_n, " 条非目标 ID，已忽略。")
}

# 确保覆盖所有请求的 GO（即使 GO.db 中缺失 TERM，也要给出空名称）
if (!export_all) {
  lacks <- setdiff(go_ids, res$GO)
  if (length(lacks)) {
    message("[WARN] 有 ", length(lacks), " 个 GO 在 GO.db 中未找到 TERM，将以空名称补齐（可接受）：示例 ", paste(head(lacks, 5), collapse = ", "), ifelse(length(lacks) > 5, " ...", ""))
    res <- bind_rows(res, tibble(GO = lacks, name = ""))
  }
}

# ============ 写出 ============
out_fp <- file.path(anno_dir, "go_terms.tsv")
readr::write_tsv(res, out_fp)
message("[OK] 写出：", out_fp, " （", nrow(res), " 行）")

# 概要信息
message(sprintf("[DONE] go_local.R 完成。导出模式：%s；输出路径：%s",
                ifelse(export_all, "export_all=TRUE（全部 GO）", "严格模式（仅 gene2go 中出现的 GO）"),
                out_fp))

