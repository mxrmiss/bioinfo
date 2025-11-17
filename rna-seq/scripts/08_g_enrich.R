#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# 08_g_enrich.R — GO/KEGG 富集（兼容 DEG 与外部基因集；仅产 TSV 原材料）
# -------------------------------
# 主要特性：
#   * 背景 = 表达基因（RNA-seq）或用户提供的 background.list（PSG/CAFE）
#   * 结果列：set_id, term_name, overlap, set_size, pvalue, padj, bg_size, sig_size
#              + hit_in_set, hit_in_sig, gene_ratio, bg_ratio
#   * 词典：ref/annotations/go_terms.tsv（GO,name）；kegg_pathway.tsv（pathway_id,name）
#   * 输入模式：
#       - enrich$input = "deg"   （默认）读取 results/deg/tables/*.deseq2.tsv
#       - enrich$input = "lists" 读取用户 sets_long.tsv 或若干 *.list + background.list
#   * 不画图，不出 xlsx（xlsx 在 09 脚本生成）
#
suppressPackageStartupMessages({
  library(jsonlite); library(yaml); library(readr); library(dplyr)
  library(tidyr); library(stringr); library(purrr)
})

# ---------- 工具 ----------
getp <- function(lst, ..., default=NULL){
  ref <- lst
  for (k in list(...)){
    if (is.null(ref)) return(default)
    if (!(is.list(ref) || is.environment(ref))) return(default)
    v <- ref[[as.character(k)]]
    if (is.null(v)) return(default)
    ref <- v
  }
  ref
}
norm_path <- function(x, fallback){
  if (is.null(x) || length(x)==0) return(as.character(fallback))
  x <- as.character(x)[1]; if (!nzchar(x)) return(as.character(fallback)); x
}
norm_num <- function(x, fallback){
  if (is.null(x) || length(x)==0) return(as.numeric(fallback))
  y <- suppressWarnings(as.numeric(x[1])); if (is.na(y)) return(as.numeric(fallback)); y
}
get_int <- function(x, default){
  v <- suppressWarnings(as.integer(x))
  if (is.null(v) || length(v)==0 || is.na(v) || v < 1) return(as.integer(default)) else v
}
ensure_dir <- function(p){ dir.create(p, recursive=TRUE, showWarnings=FALSE); p }

core_id <- function(x){
  x <- as.character(x)
  x <- sub('^[^|]+\\|', '', x)
  x <- sub('\\.exon\\d+$', '', x)
  x <- sub('\\.cds\\d+$',  '', x)
  x <- sub('\\.\\d+$',     '', x)
  x
}

make_empty <- function(){
  tibble(
    set_id = character(), term_name = character(),
    overlap = integer(), set_size = integer(),
    pvalue = double(),   padj = double(),
    bg_size = integer(), sig_size = integer(),
    hit_in_set = double(), hit_in_sig = double(),
    gene_ratio = character(), bg_ratio = character()
  )
}

# ---------- 读取配置 ----------
cfg_path <- "config.yaml"
if (!file.exists(cfg_path)) stop("缺少 config.yaml（请在项目根目录运行）")
cfg <- tryCatch(yaml::read_yaml(cfg_path), error=function(e) NULL)
if (is.null(cfg)) stop("config.yaml 解析失败")

# 允许通过环境变量覆盖核心参数（由 08_enrich_module.py 注入）
env_json <- Sys.getenv("RNASEQ_CONFIG_JSON", unset = "")
if (nzchar(env_json)) {
  env_cfg <- tryCatch(jsonlite::fromJSON(env_json, simplifyVector = TRUE), error = function(e) NULL)
  if (!is.null(env_cfg)) {
    # 仅覆盖我们使用的层级
    cfg$paths   <- modifyList(cfg$paths %||% list(),   env_cfg$paths   %||% list())
    cfg$deg     <- modifyList(cfg$deg %||% list(),     env_cfg$deg     %||% list())
    cfg$enrich  <- modifyList(cfg$enrich %||% list(),  env_cfg$enrich  %||% list())
  }
}

`%||%` <- function(a,b) if (is.null(a)) b else a

paths <- cfg$paths %||% list()
anno_dir   <- norm_path(paths$anno_dir,   "ref/annotations")
deg_dir    <- norm_path(paths$deg_dir,    "results/deg")
matrix_dir <- norm_path(paths$matrix_dir, "results/matrix")
enrich_dir <- norm_path(paths$enrich_dir, "results/enrich")

ensure_dir(file.path(enrich_dir, "tables"))
ensure_dir(file.path(enrich_dir, "genes"))

alpha_cfg  <- norm_num(getp(cfg,"deg","alpha", 0.05), 0.05)
lfc_cfg    <- norm_num(getp(cfg,"deg","lfc",   1.0),  1.0)

min_set <- get_int(getp(cfg,"enrich","min_set_size", 10),   10)
max_set <- get_int(getp(cfg,"enrich","max_set_size", 500), 500)
p_cut   <- norm_num(getp(cfg,"enrich","pval_cutoff", 0.05), 0.05)

padj_m_raw <- getp(cfg,"enrich","p_adjust","BH")
if (is.null(padj_m_raw) || length(padj_m_raw)==0) padj_m_raw <- "BH"
padj_m <- as.character(padj_m_raw)[1]
valid_m <- c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")
if (!(padj_m %in% valid_m)) padj_m <- "BH"

kegg_ex_glob <- isTRUE(getp(cfg,"enrich","kegg","exclude_global_pathways", TRUE))
kegg_req_ko  <- isTRUE(getp(cfg,"enrich","kegg","require_has_ko", FALSE))
kegg_dedup   <- isTRUE(getp(cfg,"enrich","kegg","dedup_gene_pathway", TRUE))
global_pw <- c("ko01100","ko00000","ko00001","ko00002","ko00003","ko00004")

input_mode <- tolower(as.character(getp(cfg,"enrich","input","deg")))
if (!(input_mode %in% c("deg","lists"))) input_mode <- "deg"

message(sprintf(
  "[INFO] Enrich params: mode=%s | alpha=%.3g lfc=%.3g p_cut=%.3g padj=%s min_set=%d max_set=%d | KEGG exclude_global=%s require_has_ko=%s dedup=%s",
  input_mode, alpha_cfg, lfc_cfg, p_cut, padj_m, min_set, max_set,
  kegg_ex_glob, kegg_req_ko, kegg_dedup
))

# ---------- 注释与“名称映射” ----------
g2go_fp <- file.path(anno_dir, "gene2go.tsv")
g2pw_fp <- file.path(anno_dir, "gene2pathway.tsv")
g2ko_fp <- file.path(anno_dir, "gene2ko.tsv")
go_name_fp <- file.path(anno_dir, "go_terms.tsv")
pw_name_fp <- file.path(anno_dir, "kegg_pathway.tsv")

if (!file.exists(g2go_fp)) stop("缺少注释文件：", g2go_fp)
if (!file.exists(g2pw_fp)) stop("缺少注释文件：", g2pw_fp)

g2go <- read_tsv(g2go_fp, show_col_types = FALSE) %>%
  transmute(gene_id = core_id(gene_id), GO = GO) %>%
  filter(nzchar(gene_id), nzchar(GO)) %>% distinct()

g2pw <- read_tsv(g2pw_fp, show_col_types = FALSE) %>%
  transmute(gene_id = core_id(gene_id), pathway_id = pathway_id) %>%
  filter(nzchar(gene_id), nzchar(pathway_id)) %>% distinct()

if (kegg_ex_glob){
  g2pw <- g2pw %>% filter(!(pathway_id %in% global_pw))
}
if (kegg_req_ko){
  if (!file.exists(g2ko_fp)) stop("需要 gene2ko.tsv，但未找到：", g2ko_fp)
  g2ko <- read_tsv(g2ko_fp, show_col_types = FALSE) %>%
    transmute(gene_id = core_id(gene_id)) %>%
    filter(nzchar(gene_id)) %>% distinct()
  g2pw <- g2pw %>% semi_join(g2ko, by="gene_id")
}
if (kegg_dedup){
  g2pw <- g2pw %>% distinct(gene_id, pathway_id, .keep_all = FALSE)
}

go_names <- NULL
if (file.exists(go_name_fp)){
  go_names <- read_tsv(go_name_fp, show_col_types = FALSE) %>%
    transmute(GO = GO, term_name = as.character(name)) %>% distinct()
} else {
  message("[INFO] 未提供 go_terms.tsv（GO→name），将仅输出 GO ID。")
}
pw_names <- NULL
if (file.exists(pw_name_fp)){
  pw_names <- read_tsv(pw_name_fp, show_col_types = FALSE) %>%
    transmute(pathway_id = pathway_id, term_name = as.character(name)) %>% distinct()
} else {
  message("[INFO] 未提供 kegg_pathway.tsv（pathway_id→name），将仅输出 KEGG 编号。")
}

# ---------- 背景集 ----------
bg_from_matrix <- function(){
  fp <- file.path(matrix_dir, "meta", "expressed_genes.tsv")
  if (!file.exists(fp)) stop("缺少背景基因清单：", fp)
  tb <- read_tsv(fp, show_col_types = FALSE)
  if ("gene_id" %in% names(tb)) core_id(tb$gene_id) else core_id(tb[[1]])
}

bg_set <- character()
labels <- character()
sig_picker <- NULL  # function(label, subset) -> gene vector

if (input_mode == "deg"){
  # DEG 模式
  bg_set <- unique(bg_from_matrix())

  deg_tbl_dir <- file.path(deg_dir, "tables")
  if (!dir.exists(deg_tbl_dir)) stop("缺少目录：", deg_tbl_dir)
  deg_files <- list.files(deg_tbl_dir, pattern="\\.deseq2\\.tsv$", full.names=TRUE)
  if (!length(deg_files)) stop("未发现 *.deseq2.tsv")

  deg_map <- list()
  for (fp in deg_files){
    lab <- sub("\\.deseq2\\.tsv$", "", basename(fp))
    tb <- read_tsv(fp, show_col_types = FALSE)
    need <- c("gene_id","log2FoldChange","padj","sig")
    miss <- setdiff(need, names(tb))
    if (length(miss)) stop("DEG 文件缺列：", paste(miss, collapse=", "), " @", fp)
    tb <- tb %>%
      mutate(gene_id = core_id(as.character(gene_id)),
             padj = ifelse(is.na(padj), 1, padj),
             sig  = ifelse(is.na(sig), "NS", sig),
             keep = (padj <= alpha_cfg) & (abs(log2FoldChange) >= lfc_cfg))
    deg_map[[lab]] <- tb
  }
  labels <- names(deg_map)
  sig_picker <- function(label, subset){
    d <- deg_map[[label]]
    if (subset=="All")  return(unique(d$gene_id[d$keep]))
    if (subset=="Up")   return(unique(d$gene_id[d$keep & d$sig=="Up"]))
    if (subset=="Down") return(unique(d$gene_id[d$keep & d$sig=="Down"]))
    character()
  }

} else {
  # lists 模式（PSG/CAFE）
  lists <- getp(cfg,"enrich","lists", list())
  bg_fp <- norm_path(lists$background, "results/custom_sets/background.list")
  if (!file.exists(bg_fp)) stop("lists 模式缺少背景：", bg_fp)
  bg_set <- unique(core_id(read_lines(bg_fp)))

  # sets_long 优先；否则扫描目录中的 *.list
  sets_long_fp <- lists$sets_long
  sets_dir     <- lists$dir %||% "results/custom_sets"
  sets_tbl <- NULL
  if (!is.null(sets_long_fp) && nzchar(sets_long_fp) && file.exists(sets_long_fp)){
    sets_tbl <- read_tsv(sets_long_fp, show_col_types = FALSE)
    need <- c("set","subset","gene_id")
    miss <- setdiff(need, names(sets_tbl))
    if (length(miss)) stop("sets_long 缺列：", paste(miss, collapse=", "))
    sets_tbl <- sets_tbl %>%
      mutate(set = as.character(set),
             subset = ifelse(is.na(subset)|!nzchar(subset),"All",as.character(subset)),
             gene_id = core_id(as.character(gene_id))) %>%
      filter(nzchar(set), nzchar(subset), nzchar(gene_id))
  } else {
    # 读取 *.list（文件名：<set>_<subset>.list 或 <set>.list 默认为 All）
    lf <- list.files(sets_dir, pattern="\\.list$", full.names=TRUE)
    if (!length(lf)) stop("lists 模式未发现任何 .list：", sets_dir)
    items <- lapply(lf, function(f){
      nm <- sub("\\.list$","", basename(f))
      m <- regexec("^(.+?)_(All|Up|Down)$", nm)
      s <- regmatches(nm, m)[[1]]
      if (length(s)==3){
        set <- s[2]; subset <- s[3]
      } else {
        set <- nm; subset <- "All"
      }
      tibble(set=set, subset=subset, gene_id=core_id(read_lines(f)))
    })
    sets_tbl <- bind_rows(items) %>% filter(nzchar(gene_id))
  }
  labels <- unique(sets_tbl$set)
  sig_picker <- function(label, subset){
    sets_tbl %>% filter(set==label, subset==subset) %>% pull(gene_id) %>% unique()
  }
}

bg_set <- unique(bg_set)
bg_size <- length(bg_set)

# ---------- 富集主函数 ----------
enrich_one <- function(sig_genes, g2set, set_col, name_lut=NULL){
  sig <- unique(core_id(sig_genes))
  universe <- unique(core_id(bg_set))

  ann <- g2set %>% filter(gene_id %in% universe)
  sets <- split(ann$gene_id, ann[[set_col]])
  if (!length(sets)) return(make_empty())

  lens <- vapply(sets, length, 1L)
  keep_ids <- names(sets)[lens >= min_set & lens <= max_set]
  sets <- sets[keep_ids]

  sig <- intersect(sig, universe)
  k <- length(sig); M <- length(universe)
  if (!length(sets) || k==0) return(make_empty())

  res <- lapply(names(sets), function(sid){
    S <- unique(sets[[sid]]); K <- length(S)
    x <- sum(sig %in% S)
    p <- phyper(q = x-1, m = K, n = M-K, k = k, lower.tail = FALSE)
    tibble(set_id = sid, overlap = x, set_size = K, pvalue = as.numeric(p))
  }) %>% bind_rows()

  if (!nrow(res)) return(make_empty())

  res$pvalue[!is.finite(res$pvalue) | is.na(res$pvalue)] <- 1
  res <- res %>%
    mutate(
      padj = p.adjust(pvalue, method = padj_m),
      bg_size = M, sig_size = k,
      hit_in_set = ifelse(set_size>0, overlap / set_size, 0),
      hit_in_sig = ifelse(sig_size>0, overlap / sig_size, 0),
      gene_ratio = paste0(overlap, ":", sig_size),
      bg_ratio   = paste0(set_size, ":", bg_size)
    ) %>%
    arrange(padj, pvalue)

  if (!is.null(name_lut) && nrow(name_lut)){
    key <- names(name_lut)[1]
    res <- res %>% left_join(name_lut, by = c("set_id" = key)) %>%
      mutate(term_name = ifelse(is.na(term_name), "", term_name)) %>%
      select(set_id, term_name, everything())
  } else {
    res <- res %>% mutate(term_name = "") %>% select(set_id, term_name, everything())
  }
  res
}

write_enrich <- function(tb, out_fp){
  if (is.null(tb) || !nrow(tb)) tb <- make_empty()
  write_tsv(tb, out_fp)
  message("[OK] 写出：", out_fp)
}

emit_gene_tables <- function(type_tag, subset_tag, sig_genes, g2set, res_tbl, name_lut, out_dir){
  # 仅输出“发生过命中”的条目集合
  used <- res_tbl %>% filter(overlap > 0) %>% pull(set_id) %>% unique()
  if (!length(used)) return(invisible(NULL))
  sig <- unique(core_id(sig_genes))
  ann <- g2set %>% filter((!!sym(names(g2set)[2])) %in% used, gene_id %in% sig)

  if (!is.null(name_lut) && nrow(name_lut)){
    key <- names(name_lut)[1]
    ann <- ann %>% left_join(name_lut, by = setNames(list(key), names(g2set)[2])) %>%
      mutate(term_name = ifelse(is.na(term_name), "", term_name))
  } else {
    ann <- ann %>% mutate(term_name = "")
  }

  # 长表
  long_tb <- ann %>%
    transmute(type = type_tag, subset = subset_tag,
              set_id = (!!sym(names(g2set)[2])), term_name = term_name,
              gene_id = gene_id) %>%
    distinct()
  write_tsv(long_tb, file.path(out_dir, sprintf("%s_genes_%s_%s_long.tsv", tolower(type_tag), subset_tag, "all")))

  # 宽表（每个 term 一行，命中基因以逗号拼接）
  wide_tb <- long_tb %>%
    group_by(type, subset, set_id, term_name) %>%
    summarise(genes = paste(sort(unique(gene_id)), collapse=","), .groups="drop")
  write_tsv(wide_tb, file.path(out_dir, sprintf("%s_genes_%s_%s_byterm.tsv", tolower(type_tag), subset_tag, "all")))
}

# ---------- 主循环：逐 label × {All,Up,Down} ----------
for (lab in labels){
  for (subset in c("All","Up","Down")){
    sig <- sig_picker(lab, subset)
    if (!length(sig)) next

    cov_go <- length(intersect(bg_set, unique(g2go$gene_id)))
    cov_pw <- length(intersect(bg_set, unique(g2pw$gene_id)))
    message(sprintf("[CHECK] %s/%s: bg=%d, GO_covered=%d, PW_covered=%d, sig=%d",
                    lab, subset, length(bg_set), cov_go, cov_pw, length(sig)))

    # GO
    go_res <- enrich_one(sig, g2go, "GO", name_lut = go_names)
    write_enrich(go_res, file.path(enrich_dir, "tables", sprintf("go_enrich_%s_%s.tsv", lab, subset)))

    # KEGG
    pw_res <- enrich_one(sig, g2pw, "pathway_id", name_lut = pw_names)
    write_enrich(pw_res, file.path(enrich_dir, "tables", sprintf("kegg_enrich_%s_%s.tsv", lab, subset)))

    # 附：命中基因长/宽表（仅为辅助查看；不做下游输入）
    emit_gene_tables("GO", subset, sig, g2go, go_res, go_names, file.path(enrich_dir, "genes"))
    emit_gene_tables("KEGG", subset, sig, g2pw, pw_res, pw_names, file.path(enrich_dir, "genes"))
  }
}

message("富集分析完成 ✅")

