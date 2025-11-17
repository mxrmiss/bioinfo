#!/usr/bin/env Rscript
# =============================================================================
# g_enrich.R  v9.9-royalGOclassic-xlsx
# - 保持原有逻辑/输出不变；只做“增量”：基因清单表 + 目录梳理 + 每对比导出Excel
# - 仍然：所有参数在顶部CFG修改；无命令行参数
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(clusterProfiler); library(ggplot2); library(purrr); library(scales)
  library(grid)
  has_GOdb        <<- requireNamespace("GO.db", quietly = TRUE) &&
                      requireNamespace("AnnotationDbi", quietly = TRUE)
  has_systemfonts <<- requireNamespace("systemfonts", quietly = TRUE)
  has_ragg        <<- requireNamespace("ragg", quietly = TRUE)
  has_writexl     <<- requireNamespace("writexl", quietly = TRUE)  # 用于导出 Excel
})

options(warn = 1)
SAFE_MIN <- .Machine$double.xmin
`%||%` <- function(a,b) if (!is.null(a)) a else b
pt_to_mm <- function(pt) pt/2.84527559

# ╔═══════════════ 皇上御用配置 ═══════════════╗
CFG <- list(
  # 基础路径
  deg_dir   = "results/deg",
  counts    = "results/matrix/gene_counts.tsv",

  # 注释资源
  go_tsv    = "ref/go_gene2go.tsv",
  kegg_t2g  = "ref/term2gene_kegg_pathway.tsv",
  kegg_t2n  = "ref/term2name_kegg_pathway.tsv",
  kegg_g2k  = "ref/kegg_gene2ko.tsv",

  # 模块功能（已默认关闭，保持现状）
  mod_t2g   = "ref/term2gene_kegg_module.tsv",
  mod_t2n   = "ref/term2name_kegg_module.tsv",
  enable_module = FALSE,

  # 富集与作图
  padj = 0.05, lfc = 1,
  top = 20, per_cat = 10,
  do_bp = TRUE, do_cc = TRUE, do_mf = TRUE,
  merge_all = TRUE,

  # 画图参数
  fig_single_w=6.5, fig_single_h=6.0,
  fig_double_w=8.0, fig_double_h=7.0,
  png_dpi=600, base_font_size=11, note_font_size_pt=9.0,
  font_family = "Arial", line_width=0.45, tick_len_mm=1.6,
  label_trunc_width = 45,
  color_up="#ECA8A9", color_down="#74AED4", color_ns="#D3E2B7",
  color_purple="#CFAFD4", color_yellow="#F7C97E",
  add_titles=FALSE, y_axis_title=FALSE,

  # ★新增：输出增强开关
  write_term_genes = TRUE,   # 为每个富集结果写出基因清单（长/宽两份）
  export_xlsx      = TRUE    # 为每个对比导出一个Excel工作簿（需 writexl 包）
)
# ╚════════════════════════════════════════════╝


# -------------------------------- 工具：字体/ID/清洗 -------------------------------
preflight_arial <- function(){
  if (!has_systemfonts) { cat("[G][提示] 未安装 systemfonts；用系统 sans。\n"); return(invisible(FALSE)) }
  mf <- suppressWarnings(systemfonts::match_fonts(CFG$font_family))
  ok <- is.data.frame(mf) && nrow(mf)>0 && any(!is.na(mf$path)&nzchar(mf$path))
  if (!ok) cat(sprintf("[G][提示] 未找到 '%s'，将用系统 sans。\n", CFG$font_family))
  invisible(ok)
}

read_prefixes <- function(){
  if (file.exists("ref/id_prefixes.txt")) {
    x <- readLines("ref/id_prefixes.txt", warn = FALSE); x <- gsub("\\s+","",x); x[nzchar(x)]
  } else character()
}
ID_PREFIXES <- read_prefixes()
core_id <- function(x){
  x <- as.character(x); x <- sub("\\|.*$","",x); x <- sub("\\.[0-9]+$","",x)
  if (length(ID_PREFIXES)>0) for (p in ID_PREFIXES) if (nzchar(p)) {
    pe <- gsub("([][^$.|+(){}])","\\\\\\1",p); x <- sub(paste0("^",pe),"",x)
  }
  x
}
clean_term_name <- function(s){
  s <- as.character(s)
  s <- gsub("\\s*\\[EC[: ]?[0-9\\.-]+\\]\\s*","",s,perl=TRUE)
  s <- gsub("\\s*\\(EC[: ]?[0-9\\.-]+\\)\\s*","",s,perl=TRUE)
  s <- gsub("\\s+"," ", s)
  s <- trimws(s)
  ifelse(is.na(s)|s=="","<NA>",s)
}
is_go_or_ko_id <- function(x){
  x <- as.character(x)
  grepl("^GO:\\d{7}$", x) | grepl("^ko\\d{5}$", x, ignore.case = TRUE)
}

save_png_pdf <- function(p, filebase, w, h, dpi=CFG$png_dpi){
  if (has_ragg) { ragg::agg_png(paste0(filebase,".png"), width=w, height=h, units="in", res=dpi, background="white"); print(p); dev.off()
  } else { png(paste0(filebase,".png"), width=w, height=h, units="in", res=dpi, bg="white", type="cairo"); print(p); dev.off() }
  tryCatch({ grDevices::cairo_pdf(paste0(filebase,".pdf"), width=w, height=h, family=CFG$font_family); print(p); dev.off() },
           error=function(e){ grDevices::pdf(paste0(filebase,".pdf"), width=w, height=h, family=CFG$font_family); print(p); dev.off() })
}


# -------------------------------- 读 GO -----------------------------------------
.read_go_meta <- function(ids){
  ids <- unique(ids)
  if (!has_GOdb || !length(ids)) return(NULL)
  tryCatch({
    m <- suppressMessages(
      AnnotationDbi::select(
        GO.db::GO.db,
        keys     = ids,
        keytype  = "GOID",
        columns  = c("TERM","ONTOLOGY")
      )
    )
    if (!is.null(m) && nrow(m)) {
      m <- as_tibble(m) %>% transmute(go_id = GOID, term = as.character(TERM), ont = toupper(ONTOLOGY)) %>% distinct()
      return(m)
    } else NULL
  }, error=function(e) NULL)
}

read_go_term2 <- function(path){
  x <- suppressMessages(read_tsv(path, col_types="cccccc", show_col_types=FALSE, progress=FALSE))
  colnames(x) <- tolower(colnames(x))
  names(x)[1:3] <- c("gene_id","go_id","go_name_en")[seq_len(min(3, ncol(x)))]
  x$gene_id <- core_id(x$gene_id)
  if ("go_name_en" %in% names(x)) x$go_name_en <- clean_term_name(x$go_name_en)

  dbm <- .read_go_meta(unique(x$go_id))
  t2g <- x %>% select(go_id, gene_id) %>% distinct()

  t2n <- NULL
  if (!is.null(dbm) || ("go_name_en" %in% names(x))) {
    file_t2n <- if ("go_name_en" %in% names(x)) x %>% select(go_id, go_name_en) %>% distinct() else NULL
    db_t2n   <- if (!is.null(dbm)) dbm %>% select(go_id, go_name_en = term) %>% distinct() else NULL
    t2n <- full_join(file_t2n, db_t2n, by = "go_id") %>%
      transmute(go_id, go_name_en = clean_term_name(coalesce(.data$go_name_en.y, .data$go_name_en.x))) %>%
      filter(!is.na(go_name_en) & go_name_en != "") %>% distinct()
  }

  t2o <- NULL
  if (!is.null(dbm)) {
    t2o <- dbm %>% select(go_id, ont) %>% filter(ont %in% c("BP","CC","MF")) %>% distinct()
  } else {
    ont_col <- dplyr::case_when(
      "ontology" %in% names(x) ~ "ontology",
      "aspect"   %in% names(x) ~ "aspect",
      "namespace"%in% names(x) ~ "namespace",
      TRUE ~ NA_character_
    )
    if (!is.na(ont_col)) {
      t2o <- x %>% transmute(go_id, ont = toupper(.data[[ont_col]])) %>% filter(ont %in% c("BP","CC","MF")) %>% distinct()
    }
  }
  list(t2g = t2g, t2n = t2n, t2o = t2o)
}


# -------------------------------- KEGG 清洗 -------------------------------------
norm_kegg_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^ko:",  "", x, ignore.case = TRUE)
  x <- gsub("^map",  "ko", x, ignore.case = TRUE)
  x <- gsub("^md:",  "", x, ignore.case = TRUE)
  toupper(x)
}

prepare_kegg_t2g <- function(t2g_raw, gene2ko_path) {
  nm <- names(t2g_raw)
  if (length(nm) < 2) stop("[G][错误] term2gene 文件格式错误，应为两列。")
  df <- t2g_raw %>% transmute(ID = norm_kegg_id(.data[[nm[1]]]),
                              gene = core_id(.data[[nm[2]]])) %>%
    filter(nzchar(ID) & nzchar(gene)) %>% distinct()
  prop_ko <- mean(grepl("^K\\d{5}$", toupper(df$gene)), na.rm = TRUE)
  if (prop_ko > 0.5 && file.exists(gene2ko_path)) {
    g2k <- read_tsv(gene2ko_path, show_col_types = FALSE, col_types = "cc") %>%
      transmute(gene = core_id(..1), ko = toupper(..2)) %>%
      filter(nzchar(gene) & grepl("^K\\d{5}$", ko)) %>% distinct()
    df <- df %>% mutate(ko = toupper(gene)) %>% inner_join(g2k, by = "ko") %>%
      transmute(ID = ID, gene = gene.y) %>% distinct()
  }
  df
}


# -------------------------------- 作图（现代风格） --------------------------------
bubble_plot <- function(df, topn=CFG$top, is_kegg=FALSE){
  df$p.adjust <- as.numeric(df$p.adjust); df <- df %>% filter(is.finite(p.adjust))
  parts <- strsplit(as.character(df$GeneRatio), "/")
  # —— 修正：明确返回类型为字符，避免 vapply 类型报错
  num <- suppressWarnings(as.numeric(vapply(parts, `[`, FUN.VALUE = character(1), 1L)))
  den <- suppressWarnings(as.numeric(vapply(parts, `[`, FUN.VALUE = character(1), 2L)))
  df$GeneRatioNum <- ifelse(is.finite(num/den), num/den, NA_real_)
  df <- df %>% arrange(p.adjust, desc(GeneRatioNum)) %>% head(topn)

  df$Description <- clean_term_name(df$Description)
  df$ID <- as.character(df$ID)
  desc <- df$Description; desc[desc == "<NA>"] <- NA
  lk <- if (isTRUE(is_kegg) && any(grepl("^M\\d{5}$", df$ID)) && exists("MOD_NAME_LK", inherits=TRUE)) {
          MOD_NAME_LK
        } else if (isTRUE(is_kegg) && exists("KO_NAME_LK", inherits=TRUE)) {
          KO_NAME_LK
        } else if (exists("GO_NAME_LK", inherits=TRUE)) {
          GO_NAME_LK
        } else NULL
  mapv <- if (!is.null(lk)) as.character(lk[df$ID]) else NA_character_
  mapv[is.na(mapv) | mapv=="" | mapv=="<NA>"] <- NA
  df$Label <- dplyr::coalesce(mapv, desc, df$ID)
  df$Label[is.na(df$Label) | df$Label=="<NA>"] <- df$ID[is.na(df$Label) | df$Label=="<NA>"]
  df$Label <- stringr::str_trunc(df$Label, CFG$label_trunc_width)

  df$Key <- paste0(df$ID,"::",df$Description)
  df$Key <- factor(df$Key, levels = rev(unique(df$Key)))

  ylab <- if (CFG$y_axis_title) if (is_kegg) "Pathway" else "GO term" else NULL
  xmax <- max(df$GeneRatioNum, na.rm=TRUE); if (!is.finite(xmax)||xmax<=0) xmax <- 1

  ggplot(df, aes(x=GeneRatioNum, y=Key, size=Count, color=p.adjust)) +
    geom_point(shape=16, alpha=.85) +
    scale_y_discrete(labels = setNames(df$Label, df$Key), drop = TRUE) +
    scale_x_continuous(name="GeneRatio", expand=expansion(mult=c(0.02,0.12))) +
    coord_cartesian(xlim=c(0, xmax*1.20), clip="off") +
    scale_size_continuous(name = "Count", range = c(1.8, 6.0),
                          breaks = scales::pretty_breaks(n = 4),
                          guide  = guide_legend(order = 1)) +
    scale_color_gradientn(
      name  = "p.adjust",
      colors = c("#FF0000","#FF7F00","#FFFF00","#00FF00","#00FFFF","#0000FF","#8B00FF"),
      guide  = guide_colorbar(reverse = TRUE, order = 2, barwidth = .7, barheight = 6)
    ) +
    labs(x="GeneRatio", y=ylab) +
    theme_classic(base_size=CFG$base_font_size, base_family=CFG$font_family) +
    theme(axis.text=element_text(color="black", family=CFG$font_family),
          axis.title=element_text(family=CFG$font_family),
          panel.border=element_rect(color="black", fill=NA, linewidth=CFG$line_width),
          legend.position="right",
          plot.margin=margin(10,60,10,40))
}

bar_plot <- function(df, topn=CFG$top, fill_hex=CFG$color_down, is_kegg=FALSE){
  df$p.adjust <- as.numeric(df$p.adjust)
  df <- df %>% arrange(p.adjust, desc(Count)) %>% slice_head(n=topn)

  df$Description <- clean_term_name(df$Description)
  df$ID <- as.character(df$ID)
  desc <- df$Description; desc[desc == "<NA>"] <- NA
  lk <- if (isTRUE(is_kegg) && any(grepl("^M\\d{5}$", df$ID)) && exists("MOD_NAME_LK", inherits=TRUE)) {
          MOD_NAME_LK
        } else if (isTRUE(is_kegg) && exists("KO_NAME_LK", inherits=TRUE)) {
          KO_NAME_LK
        } else if (exists("GO_NAME_LK", inherits=TRUE)) {
          GO_NAME_LK
        } else NULL
  mapv <- if (!is.null(lk)) as.character(lk[df$ID]) else NA_character_
  mapv[is.na(mapv) | mapv=="" | mapv=="<NA>"] <- NA
  df$Label <- dplyr::coalesce(mapv, desc, df$ID)
  df$Label[is.na(df$Label) | df$Label=="<NA>"] <- df$ID[is.na(df$Label) | df$Label=="<NA>"]
  df$Label <- stringr::str_trunc(df$Label, CFG$label_trunc_width)

  df$Key   <- paste0(df$ID,"::",df$Description)
  df$neglog10FDR <- -log10(pmax(df$p.adjust, SAFE_MIN))
  df$Key <- factor(df$Key, levels = rev(unique(df$Key)))

  ylab <- if (CFG$y_axis_title) if (is_kegg) "Pathway" else "GO term" else NULL

  xmax  <- max(df$neglog10FDR, na.rm=TRUE); if (!is.finite(xmax) || xmax <= 0) xmax <- 1
  nudge <- max(0.04 * xmax, 0.08)

  ggplot(df, aes(x=neglog10FDR, y=Key)) +
    geom_col(fill=fill_hex, width=.8) +
    geom_text(aes(label=Count),
              hjust=0, nudge_x=nudge,
              size=pt_to_mm(CFG$note_font_size_pt),
              color="black", family=CFG$font_family) +
    scale_y_discrete(labels = setNames(df$Label, df$Key), drop = TRUE) +
    labs(x="-log10(FDR)", y=ylab) +
    coord_cartesian(xlim=c(0, xmax * 1.32), clip="off") +
    theme_classic(base_size=CFG$base_font_size, base_family=CFG$font_family) +
    theme(axis.text=element_text(color="black", family=CFG$font_family),
          axis.line=element_line(color="black"),
          panel.border=element_rect(color="black", fill=NA, linewidth=CFG$line_width),
          plot.margin=margin(10,80,10,40))
}


# ------------------------------- 结果修复与基因表 ---------------------------------
repair_description <- function(df, name2){
  # 空就原样返回
  if (is.null(df) || !nrow(df)) return(df)

  # 若 name2 不是标准 ID/Description，也尝试从 go 格式适配
  if (!is.null(name2) && !all(c("ID","Description") %in% names(name2))) {
    if (all(c("go_id","go_name_en") %in% names(name2))) {
      name2 <- name2 %>% transmute(ID = go_id, Description = go_name_en)
    } else {
      # name2 无法适配就仅对 df 做最小清洗
      df$Description <- clean_term_name(df$Description %||% df$ID)
      return(df)
    }
  }

  # 现在的 name2：ID, Description；做一份去重并清洗
  n2 <- name2 %>%
    transmute(ID = as.character(ID),
              Description = clean_term_name(as.character(Description))) %>%
    filter(!is.na(ID) & nzchar(ID)) %>%
    distinct(ID, .keep_all = TRUE)

  # df 中必要列的最小保障
  if (!"ID" %in% names(df)) df$ID <- as.character(df$ID %||% NA_character_)
  if (!"Description" %in% names(df)) df$Description <- NA_character_

  df <- df %>%
    mutate(ID = as.character(ID),
           Description = as.character(Description),
           # 如果 Description 看起来像 ID（GO:xxxxxxx / koXXXXX），先置空待会儿用字典回填
           Description = ifelse(is_go_or_ko_id(Description), NA_character_, Description)) %>%
    left_join(n2, by = "ID", suffix = c("", "..n2")) %>%
    mutate(
      Description = clean_term_name(coalesce(Description, `Description..n2`, ID))
    ) %>%
    select(-`Description..n2`)

  # 兜底：仍为空的用 ID 顶上
  bad <- is.na(df$Description) | df$Description == "<NA>"
  if (any(bad)) df$Description[bad] <- df$ID[bad]

  df
}

# ★新增：将富集结果拆出“长表”基因清单
.extract_gene_long <- function(df){
  if (is.null(df) || !nrow(df) || !"geneID" %in% names(df)) return(tibble())
  df %>%
    transmute(
      ID,
      Description = clean_term_name(Description),
      gene = strsplit(as.character(geneID), "/")
    ) %>% tidyr::unnest(gene) %>%
    transmute(ID = as.character(ID),
              Description = as.character(Description),
              gene = core_id(gene)) %>%
    distinct()
}

# ★新增：写出两张基因表（长表 + 宽表）
.write_term_gene_tables <- function(df, outfile_prefix, genes_dir){
  if (is.null(genes_dir) || !dir.exists(genes_dir)) return(invisible(NULL))
  g_long <- .extract_gene_long(df)
  fp_long <- file.path(genes_dir, paste0(outfile_prefix, "_genes.tsv"))
  fp_wide <- file.path(genes_dir, paste0(outfile_prefix, "_genes_by_term.tsv"))
  if (!nrow(g_long)) {
    write_tsv(tibble(), fp_long)
    write_tsv(tibble(), fp_wide)
    return(invisible(NULL))
  }
  write_tsv(g_long, fp_long)
  g_wide <- g_long %>%
    group_by(ID, Description) %>%
    summarise(GeneCount = n(),
              Genes = paste(gene, collapse = ";"),
              .groups = "drop")
  write_tsv(g_wide, fp_wide)
  invisible(NULL)
}


# -------------------------------- 富集执行 ----------------------------------------
do_enrich <- function(gset, term2, name2, color_hex, outfile_prefix,
                      plot_dir, out_dir, universe=NULL, is_kegg=FALSE,
                      deg=NULL, gcol=NULL, lfc_col=NULL, padj_col=NULL,
                      tables_dir=NULL, genes_dir=NULL){

  if (length(gset) < 10) {
    # 样本偏少 -> GSEA 兜底/空表
    if (!is.null(deg) && !is.null(gcol) && !is.null(lfc_col) && !is.null(padj_col)) {
      ranked <- deg %>% filter(!is.na(.data[[padj_col]])) %>%
        mutate(score = sign(.data[[lfc_col]]) * -log10(pmax(.data[[padj_col]], SAFE_MIN))) %>%
        arrange(desc(score)) %>% distinct(!!sym(gcol), .keep_all=TRUE)
      rnks <- ranked$score; names(rnks) <- ranked[[gcol]]
      gs <- tryCatch({ GSEA(geneList=rnks, TERM2GENE=term2, TERM2NAME=name2, pAdjustMethod="BH") }, error=function(e) NULL)
      df <- if (!is.null(gs) && nrow(as.data.frame(gs))>0) as.data.frame(gs) else tibble()
      write_tsv(df, file.path(out_dir, paste0(outfile_prefix, "_GSEA.tsv")))
      if (nrow(df)) {
        p <- clusterProfiler::dotplot(gs, showCategory=min(CFG$top, nrow(df))) + theme_classic(base_size=CFG$base_font_size)
        save_png_pdf(p, file.path(plot_dir, paste0(outfile_prefix,"_GSEA_dotplot")), CFG$fig_single_w, CFG$fig_single_h)
      }
      # 同步 tables/ 与 genes/
      if (!is.null(tables_dir)) write_tsv(df, file.path(tables_dir, paste0(outfile_prefix, ".tsv")))
      if (isTRUE(CFG$write_term_genes) && !is.null(genes_dir)) {
        if (nrow(df)) .write_term_gene_tables(df, outfile_prefix, genes_dir)
        else { .write_term_gene_tables(tibble(), outfile_prefix, genes_dir) }
      }
    } else {
      df <- tibble()
      write_tsv(df, file.path(out_dir, paste0(outfile_prefix, ".tsv")))
      if (!is.null(tables_dir)) write_tsv(df, file.path(tables_dir, paste0(outfile_prefix, ".tsv")))
      if (isTRUE(CFG$write_term_genes) && !is.null(genes_dir)) .write_term_gene_tables(df, outfile_prefix, genes_dir)
    }
    return(invisible(NULL))
  }

  er <- tryCatch({
    enricher(gene=gset, TERM2GENE=term2, TERM2NAME=name2,
             universe=universe, pvalueCutoff=1, qvalueCutoff=1,
             pAdjustMethod="BH", minGSSize=3, maxGSSize=10000)
  }, error=function(e) NULL)
  df <- as.data.frame(er)
  if (!nrow(df)) {
    write_tsv(tibble(), file.path(out_dir, paste0(outfile_prefix, ".tsv")))
    if (!is.null(tables_dir)) write_tsv(tibble(), file.path(tables_dir, paste0(outfile_prefix, ".tsv")))
    if (isTRUE(CFG$write_term_genes) && !is.null(genes_dir)) .write_term_gene_tables(tibble(), outfile_prefix, genes_dir)
    return(invisible(NULL))
  }

  name2_std <- if (all(c("ID","Description") %in% names(name2))) name2 else {
    if (all(c("go_id","go_name_en") %in% names(name2))) name2 %>% transmute(ID=go_id, Description=go_name_en)
    else NULL
  }
  df <- repair_description(df, name2_std)

  if (isTRUE(is_kegg)) {
    df <- df %>% filter(!grepl("(?i)infection|cancer|cardiomyopathy|insulin|neurodegenerative|virus", Description))
  }

  write_tsv(df, file.path(out_dir, paste0(outfile_prefix, ".tsv")))
  if (!is.null(tables_dir)) write_tsv(df, file.path(tables_dir, paste0(outfile_prefix, ".tsv")))
  if (isTRUE(CFG$write_term_genes) && !is.null(genes_dir)) .write_term_gene_tables(df, outfile_prefix, genes_dir)

  save_png_pdf(bubble_plot(df, CFG$top, is_kegg), file.path(plot_dir, paste0(outfile_prefix,"_bubble")), CFG$fig_single_w, CFG$fig_single_h)
  save_png_pdf(bar_plot(df, CFG$top, color_hex, is_kegg),  file.path(plot_dir, paste0(outfile_prefix,"_bar")),    CFG$fig_single_w, CFG$fig_single_h)
}


# -------------------------------- 主流程 ------------------------------------------
preflight_arial()

cat("[G] 选用 counts：", CFG$counts, "\n")
counts <- read_tsv(CFG$counts, show_col_types=FALSE)
gene_col <- if ("gene" %in% names(counts)) "gene" else names(counts)[1]
universe_all <- unique(core_id(counts[[gene_col]]))
cat("[G] universe size = ", length(universe_all), "\n", sep="")

# 读 GO（固定返回 list: t2g/t2n/t2o）
go <- read_go_term2(CFG$go_tsv)
cat("[G] GO 名称：", if (has_GOdb) "GO.db 优先 + 本地回补\n" else "仅文件第三列/ID 兜底\n")

# 建 GO NAME 字典（优先 GO.db，再用文件第三列补空位）
GO_NAME_LK <<- NULL
{
  ids_for_name <- if (!is.null(go$t2g) && nrow(go$t2g)) unique(go$t2g$go_id) else character()
  if (has_GOdb && length(ids_for_name)) {
    dbm <- .read_go_meta(ids_for_name)
    if (!is.null(dbm) && nrow(dbm)) GO_NAME_LK <<- setNames(clean_term_name(dbm$term), dbm$go_id)
  }
  if (is.null(GO_NAME_LK) && !is.null(go$t2n)) {
    GO_NAME_LK <<- setNames(clean_term_name(go$t2n$go_name_en), go$t2n$go_id)
  } else if (!is.null(GO_NAME_LK) && !is.null(go$t2n)) {
    fill <- setNames(clean_term_name(go$t2n$go_name_en), go$t2n$go_id)
    miss <- is.na(GO_NAME_LK) | GO_NAME_LK=="" | GO_NAME_LK=="<NA>"
    if (any(miss)) GO_NAME_LK[miss] <- fill[names(GO_NAME_LK)[miss]]
  }
  if (!is.null(GO_NAME_LK)) GO_NAME_LK[is.na(GO_NAME_LK) | GO_NAME_LK=="" | GO_NAME_LK=="<NA>"] <- NA_character_
}

# 依据本体拆分
go_list <- list(
  t2g = if (!is.null(go$t2g)) go$t2g %>% transmute(ID=go_id, gene=gene_id) %>% distinct() else tibble(ID=character(), gene=character()),
  t2n = if (!is.null(go$t2n)) go$t2n %>% transmute(ID=go_id, Description=go_name_en) %>% distinct() else NULL,
  t2o = if (!is.null(go$t2o)) go$t2o %>% transmute(ID=go_id, ONT=ont) %>% distinct() else NULL
)
go_bp <- if (!is.null(go_list$t2o) && nrow(go_list$t2o)) {
  keep <- go_list$t2o %>% filter(ONT=="BP") %>% select(ID) %>% distinct()
  list(t2g=inner_join(go_list$t2g, keep, by="ID"),
       t2n=if (!is.null(go_list$t2n)) inner_join(go_list$t2n, keep, by="ID") else NULL)
} else list(t2g=go_list$t2g, t2n=go_list$t2n)
go_cc <- if (!is.null(go_list$t2o) && nrow(go_list$t2o)) {
  keep <- go_list$t2o %>% filter(ONT=="CC") %>% select(ID) %>% distinct()
  list(t2g=inner_join(go_list$t2g, keep, by="ID"),
       t2n=if (!is.null(go_list$t2n)) inner_join(go_list$t2n, keep, by="ID") else NULL)
} else list(t2g=go_list$t2g, t2n=go_list$t2n)
go_mf <- if (!is.null(go_list$t2o) && nrow(go_list$t2o)) {
  keep <- go_list$t2o %>% filter(ONT=="MF") %>% select(ID) %>% distinct()
  list(t2g=inner_join(go_list$t2g, keep, by="ID"),
       t2n=if (!is.null(go_list$t2n)) inner_join(go_list$t2n, keep, by="ID") else NULL)
} else list(t2g=go_list$t2g, t2n=go_list$t2n)

# KEGG Pathway 准备
kegg_path_ok <- file.exists(CFG$kegg_t2g) && file.exists(CFG$kegg_t2n) && file.info(CFG$kegg_t2g)$size>0
if (kegg_path_ok){
  kegg_t2g_raw <- read_tsv(CFG$kegg_t2g, col_types="cc", col_names=FALSE, show_col_types=FALSE)
  kegg_t2n     <- read_tsv(CFG$kegg_t2n, col_types="cc", col_names=FALSE, show_col_types=FALSE)
  colnames(kegg_t2n) <- c("ID","Description")
  kegg_t2n <- kegg_t2n %>% mutate(ID = norm_kegg_id(ID), Description = clean_term_name(Description)) %>% distinct()
  KO_NAME_LK <<- setNames(kegg_t2n$Description, kegg_t2n$ID)
  KO_NAME_LK[is.na(KO_NAME_LK) | KO_NAME_LK=="" | KO_NAME_LK=="<NA>"] <- NA_character_
  cat("[G] KEGG Pathway 映射：已启用（离线，ID 已统一）\n")
} else { kegg_t2g_raw <- NULL; kegg_t2n <- NULL; cat("[G][提示] KEGG Pathway 映射缺失或为空。\n") }

# 覆盖率提示
cov_go <- length(intersect(universe_all, unique(go_list$t2g$gene))) / max(1, length(universe_all))
cat(sprintf("[G] Annotation coverage: GO=%.2f%%\n", 100*cov_go))
if (!is.null(kegg_t2g_raw)){
  kegg_t2g <- prepare_kegg_t2g(kegg_t2g_raw, CFG$kegg_g2k)
  cov_kegg_path <- length(intersect(universe_all, unique(kegg_t2g$gene))) / max(1, length(universe_all))
  cat(sprintf("[G] Annotation coverage (KEGG Pathway) = %.2f%%\n", 100*cov_kegg_path))
}

# 遍历对比目录
all_dirs <- list.dirs(CFG$deg_dir, recursive=FALSE)
contrast_dirs <- all_dirs[!startsWith(basename(all_dirs), "_")]
if (!length(contrast_dirs)) stop("[G][致命] 未在 ", CFG$deg_dir, " 下发现对比目录。")

for (cdir in contrast_dirs){
  cname <- basename(cdir)
  deg_path <- if (file.exists(file.path(cdir,"DEG.tsv"))) file.path(cdir,"DEG.tsv") else file.path(cdir,"deg.tsv")
  if (!file.exists(deg_path)) { cat("[G][跳过] ", cname, " 缺少 DEG.tsv/deg.tsv\n", sep=""); next }

  out_dir  <- file.path(cdir, "enrich"); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
  # ★新增：统一副本与基因清单目录
  tables_dir <- file.path(out_dir, "tables"); dir.create(tables_dir, recursive=TRUE, showWarnings=FALSE)
  genes_dir  <- file.path(out_dir, "genes");  dir.create(genes_dir,  recursive=TRUE, showWarnings=FALSE)

  cat("[G] 处理对比：", cname, "  使用文件：", basename(deg_path), "\n", sep="")
  deg <- read_tsv(deg_path, show_col_types=FALSE)
  gcol <- if ("gene" %in% names(deg)) "gene" else names(deg)[1]
  lfc_col <- grep("log2|log2fc|log2foldchange|lfc", names(deg), value=TRUE, ignore.case=TRUE)[1]
  padj_col<- grep("padj|p\\.adjust|fdr|qvalue|qval", names(deg), value=TRUE, ignore.case=TRUE)[1]
  if (is.na(lfc_col) || is.na(padj_col)) { cat("[G][跳过] 无法识别 LFC/Padj 列：", cname, "\n"); next }

  deg <- deg %>% mutate(gene=core_id(.data[[gcol]]),
                        lfc=suppressWarnings(as.numeric(.data[[lfc_col]])),
                        padj=suppressWarnings(as.numeric(.data[[padj_col]]))) %>% filter(!is.na(padj))
  universe <- intersect(universe_all, unique(deg$gene))
  deg_f <- deg %>% filter(padj <= CFG$padj & abs(lfc) >= CFG$lfc)
  up    <- deg_f %>% filter(lfc>0) %>% pull(gene) %>% unique()
  down  <- deg_f %>% filter(lfc<0) %>% pull(gene) %>% unique()
  allg  <- unique(c(up,down))
  cat(sprintf("[G] Up=%d  Down=%d  (universe=%d)\n", length(up), length(down), length(universe)))

  # GO 富集
  if (CFG$do_bp) { do_enrich(up,   go_bp$t2g,   if (!is.null(go_bp$t2n)) go_bp$t2n else NULL, CFG$color_up,     "GO_BP_up",
                             plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir)
                   do_enrich(down, go_bp$t2g,   if (!is.null(go_bp$t2n)) go_bp$t2n else NULL, CFG$color_down,  "GO_BP_down",
                             plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir) }
  if (CFG$do_cc) { do_enrich(up,   go_cc$t2g,   if (!is.null(go_cc$t2n)) go_cc$t2n else NULL, CFG$color_purple, "GO_CC_up",
                             plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir)
                   do_enrich(down, go_cc$t2g,   if (!is.null(go_cc$t2n)) go_cc$t2n else NULL, CFG$color_purple, "GO_CC_down",
                             plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir) }
  if (CFG$do_mf) { do_enrich(up,   go_mf$t2g,   if (!is.null(go_mf$t2n)) go_mf$t2n else NULL, CFG$color_yellow, "GO_MF_up",
                             plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir)
                   do_enrich(down, go_mf$t2g,   if (!is.null(go_mf$t2n)) go_mf$t2n else NULL, CFG$color_yellow, "GO_MF_down",
                             plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir) }

  # KEGG up/down（保持原有 all_updown 不变）
  if (!is.null(kegg_t2g_raw)) {
    if (!exists("kegg_t2g") || is.null(kegg_t2g)) kegg_t2g <- prepare_kegg_t2g(kegg_t2g_raw, CFG$kegg_g2k)
    do_enrich(up,   kegg_t2g, kegg_t2n, CFG$color_up,   "KEGG_pathway_up",
              plot_dir, out_dir, universe, TRUE, deg, "gene","lfc","padj", tables_dir, genes_dir)
    do_enrich(down, kegg_t2g, kegg_t2n, CFG$color_down, "KEGG_pathway_down",
              plot_dir, out_dir, universe, TRUE, deg, "gene","lfc","padj", tables_dir, genes_dir)
  } else {
    cat("[G][提示] 跳过 KEGG_pathway_up/down（Pathway 映射为空）。\n")
  }

  # 合并（上+下）
  if (isTRUE(CFG$merge_all) && length(allg) >= 10) {
    do_enrich(allg, go_list$t2g, if (!is.null(go_list$t2n)) go_list$t2n else NULL, CFG$color_ns, "GO_all_updown",
              plot_dir, out_dir, universe, FALSE, deg, "gene","lfc","padj", tables_dir, genes_dir)
    if (!is.null(kegg_t2g_raw)) {
      do_enrich(allg, prepare_kegg_t2g(kegg_t2g_raw, CFG$kegg_g2k), kegg_t2n, CFG$color_ns, "KEGG_pathway_all_updown",
                plot_dir, out_dir, universe, TRUE, deg, "gene","lfc","padj", tables_dir, genes_dir)
    } else cat("[G][提示] 跳过 KEGG_pathway_all_updown（Pathway 映射为空）。\n")
  }

  # —— GO 三分面合并柱（保持原逻辑，仅使用 out_dir 内文件）
  try({
    pick <- c("ID","Description","p.adjust","Count")
    rd <- function(f){ df <- suppressMessages(read_tsv(f, show_col_types=FALSE)); if (!nrow(df)) return(df); df %>% select(any_of(pick)) }
    files_bp <- list.files(out_dir, pattern="(?i)^GO_BP_(up|down)\\.tsv$", full.names=TRUE)
    files_cc <- list.files(out_dir, pattern="(?i)^GO_CC_(up|down)\\.tsv$", full.names=TRUE)
    files_mf <- list.files(out_dir, pattern="(?i)^GO_MF_(up|down)\\.tsv$", full.names=TRUE)
    df <- bind_rows(
      if (length(files_bp)) bind_rows(lapply(files_bp, rd)) %>% mutate(Category = "BP") else NULL,
      if (length(files_cc)) bind_rows(lapply(files_cc, rd)) %>% mutate(Category = "CC") else NULL,
      if (length(files_mf)) bind_rows(lapply(files_mf, rd)) %>% mutate(Category = "MF") else NULL
    )
    if (!is.null(df) && nrow(df)) {
      df <- df %>% mutate(Description=clean_term_name(Description),
                          padj_safe=pmax(ifelse(is.na(p.adjust),1,p.adjust), SAFE_MIN))
      d <- df %>% group_by(Category, ID) %>%
        summarise(Description=first(Description[order(padj_safe, -as.numeric(Count))]),
                  Count=max(suppressWarnings(as.numeric(Count)), na.rm=TRUE),
                  padj=min(padj_safe, na.rm=TRUE), .groups="drop") %>%
        mutate(score=-log10(padj)) %>% filter(is.finite(score))

      per_cat <- max(1L, as.integer(CFG$per_cat))
      d_top <- d %>% group_by(Category) %>% slice_max(order_by = score, n = per_cat, with_ties = FALSE) %>%
        ungroup() %>% mutate(
          ID   = as.character(ID),
          desc = Description,
          desc = ifelse(desc == "<NA>", NA_character_, desc),
          dict_raw = if (!is.null(GO_NAME_LK)) unname(GO_NAME_LK[ID]) else NA_character_,
          dict = ifelse(is.na(dict_raw) | dict_raw == "" | dict_raw == "<NA>", NA_character_, dict_raw),
          Label = stringr::str_trunc(dplyr::coalesce(desc, dict, ID), CFG$label_trunc_width),
          Key   = paste(Category, ID, sep = "::")
        ) %>%
        arrange(Category, desc(score)) %>%
        group_by(Category) %>% mutate(Key = factor(Key, levels = rev(Key))) %>% ungroup()

      if (nrow(d_top)) {
        fill_map <- c(BP=CFG$color_up, CC=CFG$color_down, MF=CFG$color_purple)
        xmax <- max(d_top$score, na.rm=TRUE); if (!is.finite(xmax)||xmax<=0) xmax <- 1
        nudge <- max(0.04 * xmax, 0.08)
        ylab <- if (CFG$y_axis_title) "GO term" else NULL
        p <- ggplot(d_top, aes(x=score, y=Key, fill=Category)) +
          scale_x_continuous(expand=expansion(mult=c(0.02,0.05))) +
          geom_col(width=0.8) +
          geom_text(aes(label=Count), hjust=0, nudge_x=nudge,
                    family=CFG$font_family, size=pt_to_mm(CFG$note_font_size_pt)) +
          scale_y_discrete(labels=setNames(d_top$Label, d_top$Key), drop = TRUE) +
          facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
          scale_fill_manual(values=fill_map, guide="none") +
          labs(x="-log10(FDR)", y=ylab) +
          coord_cartesian(xlim=c(0, xmax*1.32), clip="off") +
          theme_classic(base_size=CFG$base_font_size) +
          theme(strip.background=element_blank(),
                strip.text=element_text(face="bold", family=CFG$font_family, size=11),
                axis.line.x=element_line(color="black", linewidth=CFG$line_width, lineend="square"),
                axis.line.y=element_line(color="black", linewidth=CFG$line_width, lineend="square"),
                axis.ticks=element_line(color="black", linewidth=CFG$line_width),
                axis.ticks.length=unit(CFG$tick_len_mm,"mm"),
                axis.text=element_text(color="black", family=CFG$font_family),
                axis.title=element_text(family=CFG$font_family),
                panel.border=element_rect(color="black", fill=NA, linewidth=CFG$line_width),
                panel.spacing=unit(0.3,"lines"),
                plot.margin=margin(10,40,10,40))
        save_png_pdf(p, file.path(plot_dir, "GO_combined_bar"), CFG$fig_double_w, CFG$fig_double_h)
      }
    }
  }, silent=TRUE)
}

# ------------------------------- 汇总表与导出 -------------------------------------
collect_tables <- function(pattern){
  files <- list.files(CFG$deg_dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  if (!length(files)) return(tibble())
  map_dfr(files, function(f){
    contrast <- basename(dirname(dirname(f)))
    df <- suppressMessages(read_tsv(f, show_col_types=FALSE))
    if (!nrow(df)) return(tibble())
    df %>% transmute(contrast=contrast, ID=.data$ID,
                     Description=clean_term_name(.data$Description),
                     FDR=suppressWarnings(as.numeric(.data$p.adjust)),
                     score=-log10(pmax(suppressWarnings(as.numeric(p.adjust)), SAFE_MIN)))
  })
}

summary_dir <- file.path(CFG$deg_dir, "_summary"); dir.create(summary_dir, recursive=TRUE, showWarnings=FALSE)
go_bp_sum    <- collect_tables("enrich/GO_BP_(up|down)\\.tsv$")
go_cc_sum    <- collect_tables("enrich/GO_CC_(up|down)\\.tsv$")
go_mf_sum    <- collect_tables("enrich/GO_MF_(up|down)\\.tsv$")
kegg_sum     <- collect_tables("enrich/KEGG_(pathway_)?(up|down)\\.tsv$")
go_all_sum   <- collect_tables("enrich/GO_all_updown\\.tsv$")
kegg_all_sum <- collect_tables("enrich/KEGG_(pathway|module)_all_updown\\.tsv$")
if (nrow(go_bp_sum))     write_tsv(go_bp_sum,     file.path(summary_dir, "enrich_top_terms_GO_BP.tsv"))
if (nrow(go_cc_sum))     write_tsv(go_cc_sum,     file.path(summary_dir, "enrich_top_terms_GO_CC.tsv"))
if (nrow(go_mf_sum))     write_tsv(go_mf_sum,     file.path(summary_dir, "enrich_top_terms_GO_MF.tsv"))
if (nrow(kegg_sum))      write_tsv(kegg_sum,      file.path(summary_dir, "enrich_top_terms_KEGG.tsv"))
if (nrow(go_all_sum))    write_tsv(go_all_sum,    file.path(summary_dir, "enrich_top_terms_GO_ALL.tsv"))
if (nrow(kegg_all_sum))  write_tsv(kegg_all_sum,  file.path(summary_dir, "enrich_top_terms_KEGG_ALL.tsv"))

# ★新增：收集所有对比的“长表”基因清单（GO + KEGG）
.collect_genes_long <- function(pattern){
  files <- list.files(CFG$deg_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(tibble())
  map_dfr(files, function(f){
    # f: .../<contrast>/enrich/genes/<prefix>_genes.tsv
    contrast <- basename(dirname(dirname(dirname(f))))
    df <- suppressMessages(read_tsv(f, show_col_types = FALSE))
    if (!nrow(df)) return(tibble())
    if (!all(c("ID","Description","gene") %in% names(df))) return(tibble())
    df %>% transmute(contrast = contrast,
                     ID = as.character(ID),
                     Description = clean_term_name(Description),
                     gene = core_id(gene))
  })
}
go_long   <- .collect_genes_long("enrich/genes/GO_.*_genes\\.tsv$")
kegg_long <- .collect_genes_long("enrich/genes/KEGG_pathway_.*_genes\\.tsv$")
if (nrow(go_long))   write_tsv(go_long,   file.path(summary_dir, "term2gene_GO_long.tsv"))
if (nrow(kegg_long)) write_tsv(kegg_long, file.path(summary_dir, "term2gene_KEGG_long.tsv"))

cat("[G] 汇总完成。\n")

# ★新增：每对比导出 Excel（多Sheet），完全可选
.export_contrast_xlsx <- function(contrast_dir){
  if (!isTRUE(CFG$export_xlsx)) return(invisible(NULL))
  if (!has_writexl) { cat("[G][提示] 未安装 writexl，跳过 Excel 导出。\n"); return(invisible(NULL)) }

  enrich_dir <- file.path(contrast_dir, "enrich")
  tables_dir <- file.path(enrich_dir,  "tables")
  genes_dir  <- file.path(enrich_dir,  "genes")
  if (!dir.exists(enrich_dir)) return(invisible(NULL))

  # 收集该对比的表
  tsvs <- list.files(enrich_dir, pattern="\\.tsv$", full.names=TRUE)
  tsvs_tables <- if (dir.exists(tables_dir)) list.files(tables_dir, pattern="\\.tsv$", full.names=TRUE) else character()
  genes_long  <- if (dir.exists(genes_dir))  list.files(genes_dir,  pattern="_genes\\.tsv$", full.names=TRUE) else character()
  genes_wide  <- if (dir.exists(genes_dir))  list.files(genes_dir,  pattern="_genes_by_term\\.tsv$", full.names=TRUE) else character()

  # 组装 Sheet：Summary + GO/KEGG 结果 + 基因清单
  read_safe <- function(f){
    df <- tryCatch(suppressMessages(read_tsv(f, show_col_types=FALSE)), error=function(e) tibble())
    if (!nrow(df)) return(tibble())
    # 防止 writexl 因 list/unknown 列报错：全部转为可写类型
    df <- mutate(df, across(everything(), as.character))
    df
  }

  wb <- list()

  # Summary（该对比 enrich/*.tsv 合并）
  if (length(tsvs)) {
    sm <- map_dfr(tsvs, function(f){
      tag <- tools::file_path_sans_ext(basename(f))
      df <- read_safe(f)
      if (!nrow(df)) return(tibble())
      df %>% mutate(file = tag)
    })
    if (nrow(sm)) wb[["Summary"]] <- sm
  }

  # GO & KEGG（tables/ 中的“标准化副本”作为主结果）
  add_sheet_from <- function(files, sheet){
    if (!length(files)) return()
    dd <- map_dfr(files, function(f){
      tag <- tools::file_path_sans_ext(basename(f))
      df <- read_safe(f)
      if (!nrow(df)) return(tibble())
      df %>% mutate(file = tag)
    })
    if (nrow(dd)) wb[[sheet]] <<- dd
  }
  add_sheet_from(tsvs_tables[grepl("^.*/GO_",  tsvs_tables)], "GO_results")
  add_sheet_from(tsvs_tables[grepl("KEGG_pathway", tsvs_tables)], "KEGG_results")

  # 基因清单
  if (length(genes_long)) {
    gl <- map_dfr(genes_long, function(f){
      tag <- tools::file_path_sans_ext(basename(f))
      df <- read_safe(f); if (!nrow(df)) return(tibble())
      df %>% mutate(file = tag)
    })
    if (nrow(gl)) wb[["Genes_long"]] <- gl
  }
  if (length(genes_wide)) {
    gw <- map_dfr(genes_wide, function(f){
      tag <- tools::file_path_sans_ext(basename(f))
      df <- read_safe(f); if (!nrow(df)) return(tibble())
      df %>% mutate(file = tag)
    })
    if (nrow(gw)) wb[["Genes_by_term"]] <- gw
  }

  if (length(wb)) {
    out_xlsx <- file.path(enrich_dir, paste0(basename(contrast_dir), ".xlsx"))
    writexl::write_xlsx(wb, out_xlsx)
    cat("[G] Excel 导出：", out_xlsx, "\n", sep="")
  }
}

# 对每个对比输出 Excel
if (isTRUE(CFG$export_xlsx) && has_writexl) {
  for (cdir in contrast_dirs) .export_contrast_xlsx(cdir)
} else if (isTRUE(CFG$export_xlsx) && !has_writexl) {
  cat("[G][提示] 开启了 export_xlsx 但系统未安装 writexl；请先安装。\n")
}

# 收工
invisible(NULL)
