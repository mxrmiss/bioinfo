#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# 08_g_enrich.R —— GO/KEGG 富集（背景口径稳态 + BP/CC/MF & by-term 宽表）
suppressPackageStartupMessages({
  library(yaml); library(readr); library(dplyr); library(stringr); library(tidyr); library(purrr)
})

`%||%` <- function(a,b) if (is.null(a)) b else a
getp <- function(lst, ..., default=NULL){
  ks <- list(...); x <- lst
  for(k in ks){ if(is.null(x)) break; if(!is.list(x)){x<-NULL;break}; x <- x[[k]] }
  if (is.null(x) || (is.atomic(x) && length(x)==0)) return(default); x
}
as_scalar_chr <- function(x){ if (is.null(x)||length(x)==0) return(NA_character_); tolower(trimws(as.character(x)[1])) }
as_scalar_num <- function(x){ if (is.null(x)||length(x)==0) return(NA_real_); suppressWarnings(as.numeric(x[1])) }
ensure_dir <- function(p){ dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
core_id <- function(x){ x <- as.character(x); x <- sub('^[^|]+\\|','',x); x <- sub('\\.exon\\d+$','',x); x <- sub('\\.cds\\d+$','',x); x <- sub('\\.\\d+$','',x); x }
norm_num <- function(x, d){ y <- as_scalar_num(x); if(is.na(y)) d else y }
get_int <- function(x,d){ y <- as_scalar_num(x); y <- if(is.na(y)) d else floor(y); if(is.na(y)) d else y }

# ---------- 读取配置 ----------
cfg_path <- "config.yaml"; if(!file.exists(cfg_path)) stop("[ERR] 未找到 config.yaml（请在项目根目录运行）")
cfg <- tryCatch(yaml::read_yaml(cfg_path), error=function(e) NULL); if(is.null(cfg)) stop("[ERR] config.yaml 解析失败")

paths <- list(
  anno_dir   = getp(cfg,"paths","anno_dir",    default="ref/annotations"),
  deg_dir    = getp(cfg,"paths","deg_dir",     default="results/deg"),
  enrich_dir = getp(cfg,"paths","enrich_dir",  default="results/enrich"),
  logs_dir   = getp(cfg,"paths","logs_dir",    default="logs")
)
paths <- lapply(paths, function(p) normalizePath(file.path(getwd(), p), winslash="/", mustWork=FALSE))

input_mode <- as_scalar_chr(getp(cfg,"enrich","input",       default="deg"))

# p_adjust：BH/BY 大写，其它小写
normalize_padj <- function(x){
  if (is.null(x) || length(x)==0) return("BH")
  s <- trimws(as.character(x[1]))
  if (toupper(s) %in% c("BH","BY")) return(toupper(s))
  s_low <- tolower(s)
  valid_low <- c("holm","hochberg","hommel","bonferroni","fdr","none")
  if (s_low %in% valid_low) return(s_low)
  "BH"
}
p_adjust <- normalize_padj(getp(cfg,"enrich","p_adjust", default="BH"))

p_cutoff   <- as_scalar_num(getp(cfg,"enrich","pval_cutoff", default=0.05))
min_set    <- get_int(getp(cfg,"enrich","min_set_size",      default=10), 10)
max_set    <- get_int(getp(cfg,"enrich","max_set_size",      default=500), 500)
source_opt <- as_scalar_chr(getp(cfg,"enrich","source",      default="emapper_only"))
deg_alpha  <- norm_num(getp(cfg,"deg","alpha", default=0.05), 0.05)
deg_lfc    <- norm_num(getp(cfg,"deg","lfc",   default=1.0),  1.0)

# ---------- 参数回显 ----------
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

# ---------- 读取名称/映射 ----------
go_name_fp <- file.path(paths$anno_dir, "go_terms.tsv"); if(!file.exists(go_name_fp)) stop("[ERR] 缺少 ", go_name_fp)
go_names <- suppressWarnings(readr::read_tsv(go_name_fp, col_types=cols(.default=col_character()))) %>% distinct(GO, .keep_all=TRUE)
# 可能存在 ontology / namespace 列
go_names$namespace <- go_names$namespace %||% go_names$ontology %||% NA_character_

pw_name_fp <- file.path(paths$anno_dir, "kegg_pathway.tsv"); if(!file.exists(pw_name_fp)) stop("[ERR] 缺少 ", pw_name_fp)
pw_names <- suppressWarnings(readr::read_tsv(pw_name_fp, col_types="cc"))
if (!("pathway_id" %in% names(pw_names))) {
  if ("pathway" %in% names(pw_names)) names(pw_names)[names(pw_names)=="pathway"] <- "pathway_id"
  if ("id"      %in% names(pw_names)) names(pw_names)[names(pw_names)=="id"]      <- "pathway_id"
}
pw_names <- pw_names %>% distinct(pathway_id, .keep_all=TRUE)

pw_t2g_fp <- file.path(paths$anno_dir, "term2gene_kegg_pathway.tsv"); if(!file.exists(pw_t2g_fp)) stop("[ERR] 缺少 ", pw_t2g_fp)
g2pw <- suppressWarnings(readr::read_tsv(pw_t2g_fp, col_types="cc"))
if (!("pathway_id" %in% names(g2pw))) {
  if ("pathway" %in% names(g2pw)) names(g2pw)[names(g2pw)=="pathway"] <- "pathway_id"
}
g2pw <- g2pw %>% transmute(term_id = .data$pathway_id, gene = core_id(.data$gene)) %>% distinct()

go_t2g_fp <- file.path(paths$anno_dir, "term2gene_go.tsv")
if (file.exists(go_t2g_fp)) {
  g2go <- suppressWarnings(readr::read_tsv(go_t2g_fp, col_types="cc")) %>%
    transmute(term_id = .data$GO, gene = core_id(.data$gene)) %>% distinct()
} else {
  g2g_fp <- file.path(paths$anno_dir, "gene2go.tsv"); if(!file.exists(g2g_fp)) stop("[ERR] 缺少 ", g2g_fp, "（或提供 term2gene_go.tsv）")
  df_go <- suppressWarnings(readr::read_tsv(g2g_fp, col_types = cols(.default=col_character())))
  if (!all(c("gene_id","GO") %in% names(df_go))) {
    stop("[ERR] gene2go.tsv 需包含列：gene_id, GO；当前列：", paste(names(df_go), collapse=","))
  }
  g2go <- df_go %>%
    transmute(term_id = .data$GO, gene = core_id(.data$gene_id)) %>%
    filter(!is.na(term_id) & !is.na(gene) & nzchar(term_id) & nzchar(gene)) %>% distinct()
}

# ---------- 背景集合（C 优先，其次 B） ----------
bg_file <- file.path(paths$enrich_dir, "background.list")
bg_all <- if (file.exists(bg_file)) {
  message("[bg] 使用统一背景：", bg_file)
  core_id(unique(readr::read_lines(bg_file)))
} else {
  message("[bg] 未提供 background.list：按体系统一背景（GO/KEGG 各自映射全集）")
  NULL
}
bg_all <- bg_all[!is.na(bg_all) & nzchar(bg_all)]

bg_go   <- if (!is.null(bg_all)) bg_all else unique(g2go$gene)
bg_pw   <- if (!is.null(bg_all)) bg_all else unique(g2pw$gene)

# ---------- 输入读取 ----------
read_deg_tables <- function(deg_dir, alpha, lfc){
  tab_dir <- file.path(deg_dir, "tables"); if(!dir.exists(tab_dir)) stop("[ERR] DEG 表目录不存在：", tab_dir)
  fs <- list.files(tab_dir, pattern="\\.deseq2\\.tsv$", full.names=TRUE); if(!length(fs)) stop("[ERR] 未找到 deseq2 结果：", tab_dir)
  res <- list()
  for (f in fs) {
    df <- suppressWarnings(readr::read_tsv(f, col_types = cols(.default=col_character())))
    if (!("gene" %in% names(df)) && "gene_id" %in% names(df)) names(df)[names(df)=="gene_id"] <- "gene"
    if (!("log2FoldChange" %in% names(df)) && "Log2FoldChange" %in% names(df)) names(df)[names(df)=="Log2FoldChange"] <- "log2FoldChange"
    if (!("padj" %in% names(df)) && "FDR" %in% names(df)) names(df)[names(df)=="FDR"] <- "padj"
    if (!("log2FoldChange" %in% names(df)) && "log2FC" %in% names(df)) names(df)[names(df)=="log2FC"] <- "log2FoldChange"
    if (!("gene" %in% names(df))) stop("[ERR] DEG 表缺少基因列（gene / gene_id）")
    if (!("log2FoldChange" %in% names(df))) stop("[ERR] DEG 表缺少 log2FoldChange 列")

    df <- df %>% mutate(padj = suppressWarnings(as.numeric(padj)),
                        log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
                        gene = core_id(gene)) %>% filter(!is.na(gene) & nzchar(gene))
    lab <- sub("\\.deseq2\\.tsv$", "", basename(f))
    all_genes <- unique(df$gene)
    sig <- df %>% filter(!is.na(padj) & padj <= alpha & !is.na(log2FoldChange))
    sig_all <- unique(sig$gene); sig_up <- unique(sig$gene[sig$log2FoldChange >=  lfc]); sig_dn <- unique(sig$gene[sig$log2FoldChange <= -lfc])
    res[[lab]] <- list(all = all_genes, up = sig_up, down = sig_dn, sig_all = sig_all)
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

# ---------- 富集主程（使用固定背景） ----------
enrich_one <- function(sig_genes, t2g_df, bg_vec, id_col="term_id", name_lut=NULL){
  # 背景固定（按体系或统一背景）；不并入子集
  bg <- unique(core_id(bg_vec))
  sig_genes <- unique(core_id(sig_genes))
  # term 集合限定在背景内
  term_sets <- t2g_df %>%
    filter(.data$gene %in% bg) %>%
    group_by(.data[[id_col]]) %>%
    summarise(genes=list(unique(.data$gene)), .groups="drop") %>%
    rename(set_id = !!id_col)
  if (!nrow(term_sets)) return(tibble())

  term_sets <- term_sets %>% mutate(set_size = lengths(genes)) %>%
    filter(set_size >= min_set, set_size <= max_set)

  hit_in_set_list <- term_sets$genes
  hit_in_sig_list <- purrr::map(hit_in_set_list, ~ intersect(.x, sig_genes))

  overlap <- lengths(hit_in_sig_list)
  N <- length(bg); n <- length(sig_genes)
  if (N==0 || n==0) return(tibble())

  res <- term_sets %>%
    transmute(
      set_id = set_id,
      overlap = overlap,
      set_size = set_size,
      pvalue = phyper(q = overlap-1, m = set_size, n = N - set_size, k = n, lower.tail = FALSE),
      bg_size = N,
      sig_size = n,
      hit_in_set = vapply(hit_in_set_list, function(x) paste(x, collapse=";"), character(1)),
      hit_in_sig = vapply(hit_in_sig_list, function(x) paste(x, collapse=";"), character(1)),
      gene_ratio = ifelse(set_size>0, overlap/set_size, 0),
      bg_ratio   = ifelse(N>0, set_size/N, 0)
    )

  if (!is.null(name_lut)) {
    names(name_lut) <- c("set_id","term_name")
    res <- res %>% left_join(name_lut %>% distinct(set_id, term_name), by="set_id")
  } else res <- res %>% mutate(term_name = NA_character_)

  res %>% relocate(term_name, .after=set_id) %>%
    mutate(padj = p.adjust(pvalue, method = p_adjust)) %>%
    arrange(pvalue, padj)
}

write_enrich <- function(df, fp){
  if (is.null(df) || !nrow(df)) {
    readr::write_tsv(tibble(set_id=character(), term_name=character(), overlap=integer(), set_size=integer(), pvalue=double(), padj=double(),
                            bg_size=integer(), sig_size=integer(), hit_in_set=character(), hit_in_sig=character(), gene_ratio=double(), bg_ratio=double()), fp)
  } else readr::write_tsv(df, fp)
}

# by-term 宽表（把 genes 列写一行一 term）
write_by_term <- function(df_hits, fp){
  if (is.null(df_hits) || !nrow(df_hits)) {
    readr::write_tsv(tibble(set_id=character(), term_name=character(), genes=character()), fp)
  } else {
    df2 <- df_hits %>% select(set_id, term_name, genes = hit_in_sig)
    readr::write_tsv(df2, fp)
  }
}

# ---------- 目录 ----------
dir.create(paths$enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$enrich_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$enrich_dir, "genes"),  recursive = TRUE, showWarnings = FALSE)

# ---------- 跑起来 ----------
if (input_mode == "deg") {
  message("[mode] 使用 DEG 模式")
  sets <- read_deg_tables(paths$deg_dir, deg_alpha, deg_lfc)

  for (label_cur in names(sets)) {
    S <- sets[[label_cur]]
    for (subset in c("all","up","down")) {
      sig <- S[[subset]]

      # GO 富集（体系统一背景）
      go_res <- enrich_one(sig, g2go, bg_go, "term_id",
                           name_lut = rename(go_names, set_id = GO, term_name = name))
      write_enrich(go_res, file.path(paths$enrich_dir,"tables", sprintf("go_enrich_%s_%s.tsv", label_cur, subset)))

      # KEGG 富集（体系统一背景）
      pw_res <- enrich_one(sig, g2pw, bg_pw, "term_id",
                           name_lut = rename(pw_names, set_id = pathway_id, term_name = name))
      write_enrich(pw_res, file.path(paths$enrich_dir,"tables", sprintf("kegg_enrich_%s_%s.tsv", label_cur, subset)))

      # ---- by-term 宽表：GO 的 BP/CC/MF 拆分 + KEGG pathway ----
      # 先把 GO 结果带上 namespace
      if (nrow(go_res)) {
        go_res_ns <- go_res %>% left_join(rename(go_names, set_id=GO, ns = namespace), by="set_id")
        bp <- go_res_ns %>% filter(tolower(ns) %in% c("bp","biological_process"))
        cc <- go_res_ns %>% filter(tolower(ns) %in% c("cc","cellular_component"))
        mf <- go_res_ns %>% filter(tolower(ns) %in% c("mf","molecular_function"))

        write_by_term(bp, file.path(paths$enrich_dir,"genes", sprintf("GO_BP_%s_genes_by_term.tsv", subset)))
        write_by_term(cc, file.path(paths$enrich_dir,"genes", sprintf("GO_CC_%s_genes_by_term.tsv", subset)))
        write_by_term(mf, file.path(paths$enrich_dir,"genes", sprintf("GO_MF_%s_genes_by_term.tsv", subset)))
      } else {
        # 仍输出空表，便于 09 装订
        write_by_term(NULL, file.path(paths$enrich_dir,"genes", sprintf("GO_BP_%s_genes_by_term.tsv", subset)))
        write_by_term(NULL, file.path(paths$enrich_dir,"genes", sprintf("GO_CC_%s_genes_by_term.tsv", subset)))
        write_by_term(NULL, file.path(paths$enrich_dir,"genes", sprintf("GO_MF_%s_genes_by_term.tsv", subset)))
      }

      # KEGG pathway by-term（使用交集）
      write_by_term(pw_res, file.path(paths$enrich_dir,"genes", sprintf("KEGG_pathway_%s_genes_by_term.tsv", subset)))

      # 兼容保留：原有命中长表（若您已在上游生成，会被 09 读取；此处不再额外写）
    }
  }

} else {
  message("[mode] 使用 lists 模式")
  sets <- read_long_sets(paths$enrich_dir)
  for (label_cur in names(sets)) {
    sig <- sets[[label_cur]]$all; subset <- "all"
    go_res <- enrich_one(sig, g2go, bg_go, "term_id",
                         name_lut = rename(go_names, set_id = GO, term_name = name))
    write_enrich(go_res, file.path(paths$enrich_dir,"tables", sprintf("go_enrich_%s_%s.tsv", label_cur, subset)))
    pw_res <- enrich_one(sig, g2pw, bg_pw, "term_id",
                         name_lut = rename(pw_names, set_id = pathway_id, term_name = name))
    write_enrich(pw_res, file.path(paths$enrich_dir,"tables", sprintf("kegg_enrich_%s_%s.tsv", label_cur, subset)))

    # lists 模式也写 by-term
    if (nrow(go_res)) {
      go_res_ns <- go_res %>% left_join(rename(go_names, set_id=GO, ns = namespace), by="set_id")
      bp <- go_res_ns %>% filter(tolower(ns) %in% c("bp","biological_process"))
      cc <- go_res_ns %>% filter(tolower(ns) %in% c("cc","cellular_component"))
      mf <- go_res_ns %>% filter(tolower(ns) %in% c("mf","molecular_function"))
      write_by_term(bp, file.path(paths$enrich_dir,"genes", "GO_BP_all_genes_by_term.tsv"))
      write_by_term(cc, file.path(paths$enrich_dir,"genes", "GO_CC_all_genes_by_term.tsv"))
      write_by_term(mf, file.path(paths$enrich_dir,"genes", "GO_MF_all_genes_by_term.tsv"))
    } else {
      write_by_term(NULL, file.path(paths$enrich_dir,"genes", "GO_BP_all_genes_by_term.tsv"))
      write_by_term(NULL, file.path(paths$enrich_dir,"genes", "GO_CC_all_genes_by_term.tsv"))
      write_by_term(NULL, file.path(paths$enrich_dir,"genes", "GO_MF_all_genes_by_term.tsv"))
    }
    write_by_term(pw_res, file.path(paths$enrich_dir,"genes", "KEGG_pathway_all_genes_by_term.tsv"))
  }
}

message("[done] 富集分析完成 ✅")

