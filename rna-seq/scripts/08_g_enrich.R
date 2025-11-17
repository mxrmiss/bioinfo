#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# 08_g_enrich.R —— GO/KEGG 富集（稳态完整版；兼容 DEG 与外部基因集，仅产 TSV 原材料）
suppressPackageStartupMessages({
  library(yaml); library(readr); library(dplyr); library(stringr); library(tidyr); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
getp <- function(lst, ..., default = NULL){
  keys <- list(...); x <- lst
  for(k in keys){ if(is.null(x)) break; if(!is.list(x)){x<-NULL;break}; x <- x[[k]] }
  if (is.null(x) || (is.atomic(x) && length(x)==0)) return(default); x
}
as_scalar_chr <- function(x){ if (is.null(x)||length(x)==0) return(NA_character_); tolower(trimws(as.character(x)[1])) }
as_scalar_num <- function(x){ if (is.null(x)||length(x)==0) return(NA_real_); suppressWarnings(as.numeric(x[1])) }
ensure_dir <- function(p){ dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
core_id <- function(x){ x <- as.character(x); x <- sub('^[^|]+\\|','',x); x <- sub('\\.exon\\d+$','',x); x <- sub('\\.cds\\d+$','',x); x <- sub('\\.\\d+$','',x); x }
norm_num <- function(x, d){ y <- as_scalar_num(x); if(is.na(y)) d else y }
get_int <- function(x,d){ y <- as_scalar_num(x); y <- if(is.na(y)) d else floor(y); if(is.na(y)) d else y }

# ---- 读取配置 ----
cfg_path <- "config.yaml"; if(!file.exists(cfg_path)) stop("[ERR] 未找到配置文件 config.yaml（请在项目根目录运行）")
cfg <- tryCatch(yaml::read_yaml(cfg_path), error=function(e) NULL); if(is.null(cfg)) stop("[ERR] config.yaml 解析失败")

paths <- list(
  anno_dir   = getp(cfg,"paths","anno_dir",    default="ref/annotations"),
  deg_dir    = getp(cfg,"paths","deg_dir",     default="results/deg"),
  enrich_dir = getp(cfg,"paths","enrich_dir",  default="results/enrich"),
  matrix_dir = getp(cfg,"paths","matrix_dir",  default="results/matrix"),
  logs_dir   = getp(cfg,"paths","logs_dir",    default="logs")
)
paths <- lapply(paths, function(p) normalizePath(file.path(getwd(), p), winslash="/", mustWork=FALSE))

input_mode <- as_scalar_chr(getp(cfg,"enrich","input",       default="deg"))
# NOTE: p_adjust 特殊处理，不使用 as_scalar_chr（避免被强制小写）
p_adjust_raw <- getp(cfg,"enrich","p_adjust", default="BH")
p_cutoff   <- as_scalar_num(getp(cfg,"enrich","pval_cutoff", default=0.05))
min_set    <- get_int(getp(cfg,"enrich","min_set_size",      default=10), 10)
max_set    <- get_int(getp(cfg,"enrich","max_set_size",      default=500), 500)
source_opt <- as_scalar_chr(getp(cfg,"enrich","source",      default="emapper_only"))

deg_alpha  <- norm_num(getp(cfg,"deg","alpha", default=0.05), 0.05)
deg_lfc    <- norm_num(getp(cfg,"deg","lfc",   default=1.0),  1.0)

# ---- 规范化 p_adjust：BH/BY 用大写，其余小写 ----
normalize_padj <- function(x){
  if (is.null(x) || length(x)==0) return("BH")
  s <- trimws(as.character(x[1]))
  if (toupper(s) %in% c("BH","BY")) return(toupper(s))
  s_low <- tolower(s)
  valid_low <- c("holm","hochberg","hommel","bonferroni","fdr","none")
  if (s_low %in% valid_low) return(s_low)
  "BH"
}
p_adjust <- normalize_padj(p_adjust_raw)

# ---- 参数校验与回显 ----
if (is.na(input_mode) || !(input_mode %in% c("deg","lists"))) { message("[WARN] enrich.input 非法/空，回退 'deg'"); input_mode <- "deg" }
valid_mixed <- c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")
if (!(p_adjust %in% valid_mixed)) { message("[WARN] enrich.p_adjust 非法/空，回退 'BH'"); p_adjust <- "BH" }
if (is.na(p_cutoff) || p_cutoff<=0 || p_cutoff>1) { message("[WARN] enrich.pval_cutoff 非法/空，回退 0.05"); p_cutoff <- 0.05 }
if (is.na(min_set) || min_set<1){ message("[WARN] enrich.min_set_size 非法/空，回退 10"); min_set <- 10 }
if (is.na(max_set) || max_set<min_set){ message("[WARN] enrich.max_set_size 非法/空，回退 500"); max_set <- max(500, min_set) }

message("[init] Using config: ", normalizePath(cfg_path, winslash="/", mustWork=FALSE))
message("[init] paths.anno_dir   = ", paths$anno_dir)
message("[init] paths.deg_dir    = ", paths$deg_dir)
message("[init] paths.enrich_dir = ", paths$enrich_dir)
message("[init] enrich.input     = ", input_mode)
message("[init] enrich.p_adjust  = ", p_adjust)
message("[init] enrich.pval_cut  = ", p_cutoff)
message("[init] enrich.min/max   = ", paste0(min_set,"/",max_set))
message("[init] deg.alpha/lfc    = ", paste0(deg_alpha,"/",deg_lfc))
message("[init] source           = ", source_opt)

# ---- 名称与映射 ----
go_name_fp <- file.path(paths$anno_dir, "go_terms.tsv"); if(!file.exists(go_name_fp)) stop("[ERR] 缺少 ", go_name_fp)
go_names <- suppressWarnings(readr::read_tsv(go_name_fp, col_types="cc")) %>% distinct(GO, .keep_all=TRUE)

pw_name_fp <- file.path(paths$anno_dir, "kegg_pathway.tsv"); if(!file.exists(pw_name_fp)) stop("[ERR] 缺少 ", pw_name_fp)
pw_names <- suppressWarnings(readr::read_tsv(pw_name_fp, col_types="cc"))
if (!("pathway_id" %in% names(pw_names))) {
  if ("pathway" %in% names(pw_names)) names(pw_names)[names(pw_names)=="pathway"] <- "pathway_id"
  if ("id" %in% names(pw_names)) names(pw_names)[names(pw_names)=="id"] <- "pathway_id"
}
pw_names <- pw_names %>% distinct(pathway_id, .keep_all=TRUE)

pw_t2g_fp <- file.path(paths$anno_dir, "term2gene_kegg_pathway.tsv"); if(!file.exists(pw_t2g_fp)) stop("[ERR] 缺少 ", pw_t2g_fp)
g2pw <- suppressWarnings(readr::read_tsv(pw_t2g_fp, col_types="cc"))
if (!("pathway_id" %in% names(g2pw))) {
  if ("pathway" %in% names(g2pw)) names(g2pw)[names(g2pw)=="pathway"] <- "pathway_id"
}
g2pw <- g2pw %>% transmute(term_id = .data$pathway_id, gene = core_id(.data$gene)) %>% distinct()

# ---- GO term2gene：优先 term2gene_go.tsv；否则由 gene2go.tsv 现算（使用 gene_id 列）----
go_t2g_fp <- file.path(paths$anno_dir, "term2gene_go.tsv")
if (file.exists(go_t2g_fp)) {
  g2go <- suppressWarnings(readr::read_tsv(go_t2g_fp, col_types="cc")) %>%
    transmute(term_id = .data$GO, gene = core_id(.data$gene)) %>% distinct()
} else {
  g2g_fp <- file.path(paths$anno_dir, "gene2go.tsv"); if(!file.exists(g2g_fp)) stop("[ERR] 缺少 ", g2g_fp, "（或提供 term2gene_go.tsv）")
  df_go <- suppressWarnings(readr::read_tsv(g2g_fp, col_types = cols(.default=col_character())))
  if (!all(c("gene_id","GO") %in% names(df_go))) {
    stop("[ERR] gene2go.tsv 需包含表头：gene_id\tGO；当前列为：", paste(names(df_go), collapse=","))
  }
  g2go <- df_go %>%
    transmute(term_id = .data$GO, gene = core_id(.data$gene_id)) %>%
    filter(!is.na(term_id) & !is.na(gene) & nzchar(term_id) & nzchar(gene)) %>%
    distinct()
}

# ---- 背景 ----
bg_file_candidates <- c(file.path(paths$enrich_dir,"background.list"), file.path(paths$deg_dir,"background.list"))
bg_vec <- NULL; for(fp in bg_file_candidates){ if(file.exists(fp)){ bg_vec <- readr::read_lines(fp); break } }
if (is.null(bg_vec)) { message("[info] 未提供 background.list，使用 GO/KEGG 映射联合基因全集作为背景"); bg_vec <- union(g2go$gene, g2pw$gene) }
bg_vec <- core_id(unique(bg_vec)); bg_vec <- bg_vec[!is.na(bg_vec) & nzchar(bg_vec)]

# ---- 输入读取 ----
read_deg_tables <- function(deg_dir, alpha, lfc){
  tab_dir <- file.path(deg_dir, "tables"); if(!dir.exists(tab_dir)) stop("[ERR] DEG 表目录不存在：", tab_dir)
  fs <- list.files(tab_dir, pattern="\\.deseq2\\.tsv$", full.names=TRUE); if(!length(fs)) stop("[ERR] 未找到 deseq2 结果：", tab_dir)
  res <- list()
  for (f in fs) {
    df <- suppressWarnings(readr::read_tsv(f, col_types = cols(.default=col_character())))
    # 统一列名：gene_id→gene；Log2FoldChange→log2FoldChange；兼容 Gene/id/FDR/log2FC
    if (!("gene" %in% names(df)) && "gene_id" %in% names(df)) names(df)[names(df)=="gene_id"] <- "gene"
    if (!("log2FoldChange" %in% names(df)) && "Log2FoldChange" %in% names(df)) names(df)[names(df)=="Log2FoldChange"] <- "log2FoldChange"
    if (!("gene" %in% names(df))) {
      if ("Gene" %in% names(df)) names(df)[names(df)=="Gene"] <- "gene"
      if ("id"   %in% names(df)) names(df)[names(df)=="id"]   <- "gene"
    }
    if (!("padj" %in% names(df)) && "FDR" %in% names(df)) names(df)[names(df)=="FDR"] <- "padj"
    if (!("log2FoldChange" %in% names(df)) && "log2FC" %in% names(df)) names(df)[names(df)=="log2FC"] <- "log2FoldChange"
    if (!("gene" %in% names(df))) stop("[ERR] DEG 表缺少基因列。请提供列名为 'gene'；当前列：", paste(names(df), collapse=","))
    if (!("log2FoldChange" %in% names(df))) stop("[ERR] DEG 表缺少 log2FoldChange 列；当前列：", paste(names(df), collapse=","))

    df <- df %>% mutate(padj = suppressWarnings(as.numeric(padj)),
                        log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
                        gene = core_id(gene)) %>% filter(!is.na(gene) & nzchar(gene))
    lab <- sub("\\.deseq2\\.tsv$", "", basename(f))
    all_genes <- unique(df$gene)
    sig <- df %>% filter(!is.na(padj) & padj <= alpha & !is.na(log2FoldChange))
    sig_all <- unique(sig$gene); sig_up <- unique(sig$gene[sig$log2FoldChange >=  lfc]); sig_dn <- unique(sig$gene[sig$log2FoldChange <= -lfc])
    res[[lab]] <- list(all = intersect(all_genes, bg_vec), up = intersect(sig_up, bg_vec), down = intersect(sig_dn, bg_vec), sig_all = intersect(sig_all, bg_vec))
  }
  res
}

read_long_sets <- function(enrich_dir){
  long_fp <- file.path(enrich_dir,"sets_long.tsv"); out <- list()
  if (file.exists(long_fp)) {
    df <- suppressWarnings(readr::read_tsv(long_fp, col_types="ccc"))
    if (!all(c("label","subset","gene") %in% names(df))) stop("[ERR] sets_long.tsv 需包含列：label, subset, gene")
    df <- df %>% mutate(gene = core_id(gene)) %>% filter(nzchar(gene))
    for (lb in unique(df$label)) out[[lb]] <- list(all = unique(core_id(df$gene[df$label==lb])))
    return(out)
  }
  lst_dir <- file.path(enrich_dir,"sets"); if(!dir.exists(lst_dir)) stop("[ERR] 未找到 sets_long.tsv，且缺少 sets/ 目录：", enrich_dir)
  fs <- list.files(lst_dir, pattern="\\.list$", full.names=TRUE); if(!length(fs)) stop("[ERR] sets/ 目录下未找到 *.list")
  for (f in fs){ lb <- sub("\\.list$","",basename(f)); g <- readr::read_lines(f); out[[lb]] <- list(all = core_id(unique(g))) }
  out
}

# ---- 富集主程 ----
enrich_one <- function(sig_genes, t2g_df, id_col="term_id", name_lut=NULL){
  if (length(sig_genes)==0) return(tibble(set_id=character(), term_name=character(), overlap=integer(), set_size=integer(),
                                          pvalue=double(), padj=double(), bg_size=integer(), sig_size=integer(),
                                          hit_in_set=character(), hit_in_sig=character(), gene_ratio=double(), bg_ratio=double()))
  sig_genes <- unique(sig_genes); bg <- unique(bg_vec); bg <- unique(bg[bg %in% union(t2g_df$gene, sig_genes)])
  N <- length(bg); n <- length(intersect(sig_genes, bg)); if(N==0 || n==0){ message("[warn] 背景为空或无交集，跳过富集"); return(tibble()) }
  term_sets <- t2g_df %>% filter(.data$gene %in% bg) %>% group_by(.data[[id_col]]) %>%
    summarise(genes=list(unique(.data$gene)), .groups="drop") %>% rename(set_id = !!id_col)
  if (!nrow(term_sets)) return(tibble())
  term_sets <- term_sets %>% mutate(set_size = lengths(genes)) %>% filter(set_size >= min_set, set_size <= max_set)
  if (!nrow(term_sets)) return(tibble())
  res <- term_sets %>% mutate(hit_in_set_list = purrr::map(genes, ~ intersect(.x, sig_genes)),
                              overlap = lengths(hit_in_set_list), bg_size = N, sig_size = n) %>%
    filter(overlap > 0) %>%
    mutate(pvalue = phyper(q = overlap-1, m = set_size, n = bg_size - set_size, k = sig_size, lower.tail = FALSE),
           hit_in_set = vapply(hit_in_set_list, function(x) paste(x, collapse=";"), character(1)),
           hit_in_sig = hit_in_set,
           gene_ratio = ifelse(set_size>0, overlap/set_size, 0),
           bg_ratio   = ifelse(bg_size>0, set_size/bg_size, 0)) %>%
    select(set_id, overlap, set_size, pvalue, bg_size, sig_size, hit_in_set, hit_in_sig, gene_ratio, bg_ratio)
  if (!nrow(res)) return(tibble())
  if (!is.null(name_lut)) { names(name_lut) <- c("set_id","term_name"); res <- res %>% left_join(name_lut %>% distinct(set_id, term_name), by="set_id") }
  else { res <- res %>% mutate(term_name = NA_character_) }
  res <- res %>% mutate(padj = p.adjust(pvalue, method = p_adjust)) %>% relocate(term_name, .after=set_id) %>% arrange(pvalue, padj)
  res
}

write_enrich <- function(df, fp){
  if (is.null(df) || !nrow(df)) {
    readr::write_tsv(tibble(set_id=character(), term_name=character(), overlap=integer(), set_size=integer(), pvalue=double(), padj=double(),
                            bg_size=integer(), sig_size=integer(), hit_in_set=character(), hit_in_sig=character(), gene_ratio=double(), bg_ratio=double()), fp)
  } else readr::write_tsv(df, fp)
}

emit_gene_tables <- function(type, subset, sig, t2g, enrich_df, name_lut, out_dir){
  if (missing(enrich_df) || is.null(enrich_df) || !nrow(enrich_df)) return(invisible(NULL))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  hit_in_set <- enrich_df %>% select(set_id, term_name, hit_in_set) %>%
    mutate(genes = str_split(hit_in_set, ";")) %>% select(-hit_in_set) %>% tidyr::unnest(genes) %>% filter(nzchar(genes)) %>% distinct()
  hit_in_sig <- hit_in_set
  readr::write_tsv(hit_in_set, file.path(out_dir, sprintf("%s_%s__hit_in_set__%s.tsv", type, subset, label_cur)))
  readr::write_tsv(hit_in_sig, file.path(out_dir, sprintf("%s_%s__hit_in_sig__%s.tsv", type, subset, label_cur)))
}

# ---- 运行 ----
dir.create(paths$enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$enrich_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$enrich_dir, "genes"),  recursive = TRUE, showWarnings = FALSE)

if (input_mode == "deg") {
  message("[mode] 使用 DEG 模式")
  sets <- read_deg_tables(paths$deg_dir, deg_alpha, deg_lfc)
  for (label_cur in names(sets)) {
    S <- sets[[label_cur]]
    for (subset in c("all","up","down")) {
      sig <- S[[subset]]
      go_res <- enrich_one(sig, g2go, "term_id", name_lut = rename(go_names, term_id = GO, name = name))
      write_enrich(go_res, file.path(paths$enrich_dir,"tables", sprintf("go_enrich_%s_%s.tsv", label_cur, subset)))
      pw_res <- enrich_one(sig, g2pw, "term_id", name_lut = rename(pw_names, term_id = pathway_id, name = name))
      write_enrich(pw_res, file.path(paths$enrich_dir,"tables", sprintf("kegg_enrich_%s_%s.tsv", label_cur, subset)))
      emit_gene_tables("GO", subset, sig, g2go, go_res, go_names, file.path(paths$enrich_dir,"genes"))
      emit_gene_tables("KEGG", subset, sig, g2pw, pw_res, pw_names, file.path(paths$enrich_dir,"genes"))
    }
  }
} else {
  message("[mode] 使用 lists 模式")
  sets <- read_long_sets(paths$enrich_dir)
  for (label_cur in names(sets)) {
    sig <- sets[[label_cur]]$all; subset <- "all"
    go_res <- enrich_one(sig, g2go, "term_id", name_lut = rename(go_names, term_id = GO, name = name))
    write_enrich(go_res, file.path(paths$enrich_dir,"tables", sprintf("go_enrich_%s_%s.tsv", label_cur, subset)))
    pw_res <- enrich_one(sig, g2pw, "term_id", name_lut = rename(pw_names, term_id = pathway_id, name = name))
    write_enrich(pw_res, file.path(paths$enrich_dir,"tables", sprintf("kegg_enrich_%s_%s.tsv", label_cur, subset)))
    emit_gene_tables("GO", subset, sig, g2go, go_res, go_names, file.path(paths$enrich_dir,"genes"))
    emit_gene_tables("KEGG", subset, sig, g2pw, pw_res, pw_names, file.path(paths$enrich_dir,"genes"))
  }
}

message("[done] 富集分析完成 ✅")

