#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# 08_g_enrich.R —— GO/KEGG ORA + GSEA 执行模块
#
# 功能：
#   1) 读取 config.yaml 与 results/08_enrich/tasks.tsv；
#   2) 读取注释词典（gene2go/gene2pathway/GO term2name/obsolete_map/KEGG term2gene/term2name）；
#   3) 对每个任务（label + test_set）执行 ORA（GO + KEGG），输出统一表头的 by-term 表；
#   4) 生成 GO_sig_up / GO_sig_down（GO 三本体纵合 + FDR 过滤）；
#   5) 若存在 gene_id → gene_name 词典，则在 ORA 结果中增加 gene_names 列；
#   6) 若 config.gsea.enable = TRUE，则对每个 RNA label 运行 GSEA（GO + KEGG），并输出绘图原材料；
#   7) 汇总所有 label 的 ORA + GSEA 显著条目数，写出 results/08_enrich/summary.tsv。
#
# 说明：
# - GO obsolete_map.tsv：若缺失或为空，则跳过替换，直接使用 gene2go。
# - 本脚本不输出 test_set=all 的 ORA 文件（你流水线 tasks.tsv 仍可包含 all 任务，但不会写出对应表格）。

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
  library(GO.db)
  library(AnnotationDbi)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

log_msg <- function(level = "INFO", ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(...)
  cat(ts, "[08_enrich]", sprintf("[%s]", level), txt, "\n")
}

stop_with <- function(...) {
  log_msg("ERROR", ...)
  stop(paste0(...), call. = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

read_gene_list <- function(path) {
  # 兼容有 header (gene_id) 或无 header 的单列文件
  if (!file.exists(path)) {
    return(character())
  }
  dt <- tryCatch(
    data.table::fread(path, sep = "\t", header = TRUE),
    error = function(e) NULL
  )
  if (!is.null(dt) && "gene_id" %in% colnames(dt)) {
    g <- dt$gene_id
    g <- g[!is.na(g)]
    return(unique(g))
  }
  dt2 <- data.table::fread(path, sep = "\t", header = FALSE)
  if (ncol(dt2) < 1L) return(character())
  g <- dt2[[1]]
  g <- g[!is.na(g)]
  unique(g)
}

gene_ids_to_names <- function(gene_ids, gene_name_map) {
  if (is.null(gene_name_map)) return(rep(NA_character_, length(gene_ids)))
  out <- unname(gene_name_map[gene_ids])
  out[is.na(out)] <- NA_character_
  out
}

parse_ratio_vec <- function(x) {
  sapply(
    strsplit(as.character(x), "/"),
    function(v) {
      if (length(v) != 2) return(NA_real_)
      num <- suppressWarnings(as.numeric(v[1]))
      den <- suppressWarnings(as.numeric(v[2]))
      if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
      num / den
    }
  )
}

format_go_sub <- function(dt_sub, ontology, test_set, universe_size,
                          minGS, maxGS, gene_name_map) {
  if (is.null(dt_sub) || nrow(dt_sub) == 0L) {
    return(data.table(
      term_id        = character(),
      term_name      = character(),
      test_set       = character(),
      ontology       = character(),
      gene_ids       = character(),
      gene_names     = character(),
      gene_count     = integer(),
      bg_size        = integer(),
      gene_ratio     = numeric(),
      bg_ratio       = numeric(),
      p_value        = numeric(),
      p_adjust       = numeric(),
      universe_size  = integer(),
      min_gs         = integer(),
      max_gs         = integer()
    ))
  }

  gene_list <- strsplit(dt_sub$geneID, "/")
  gene_ids_vec <- vapply(
    gene_list,
    function(x) paste(sort(unique(x)), collapse = ";"),
    character(1)
  )
  gene_names_vec <- vapply(
    gene_list,
    function(x) {
      nms <- gene_ids_to_names(x, gene_name_map)
      if (all(is.na(nms))) return(NA_character_)
      valid <- !is.na(nms)
      if (!any(valid)) return(NA_character_)
      paste(unique(nms[valid]), collapse = ";")
    },
    character(1)
  )

  gene_ratio_num <- parse_ratio_vec(dt_sub$GeneRatio)
  bg_ratio_num   <- parse_ratio_vec(dt_sub$BgRatio)

  data.table(
    term_id       = dt_sub$ID,
    term_name     = dt_sub$Description,
    test_set      = test_set,
    ontology      = ontology,
    gene_ids      = gene_ids_vec,
    gene_names    = gene_names_vec,
    gene_count    = dt_sub$Count,
    bg_size       = universe_size,
    gene_ratio    = gene_ratio_num,
    bg_ratio      = bg_ratio_num,
    p_value       = dt_sub$pvalue,
    p_adjust      = dt_sub$p_adjust,
    universe_size = universe_size,
    min_gs        = minGS,
    max_gs        = maxGS
  )
}

# KEGG 输出表：固定列顺序 + gene_names + ontology=KEGG
format_kegg_res <- function(dt, test_set, universe_size,
                            minGS, maxGS, gene_name_map, count_mode) {
  # 列顺序约定：
  # pathway_id, term_name, test_set, count_mode,
  # gene_ids, gene_names, gene_count, bg_size, gene_ratio, bg_ratio,
  # p_value, p_adjust, universe_size, min_gs, max_gs, ontology
  if (is.null(dt) || nrow(dt) == 0L) {
    return(data.table(
      pathway_id    = character(),
      term_name     = character(),
      test_set      = character(),
      count_mode    = character(),
      gene_ids      = character(),
      gene_names    = character(),
      gene_count    = integer(),
      bg_size       = integer(),
      gene_ratio    = numeric(),
      bg_ratio      = numeric(),
      p_value       = numeric(),
      p_adjust      = numeric(),
      universe_size = integer(),
      min_gs        = integer(),
      max_gs        = integer(),
      ontology      = character()
    ))
  }

  gene_list <- strsplit(dt$geneID, "/")
  gene_ids_vec <- vapply(
    gene_list,
    function(x) paste(sort(unique(x)), collapse = ";"),
    character(1)
  )
  gene_names_vec <- vapply(
    gene_list,
    function(x) {
      nms <- gene_ids_to_names(x, gene_name_map)
      if (all(is.na(nms))) return(NA_character_)
      valid <- !is.na(nms)
      if (!any(valid)) return(NA_character_)
      paste(unique(nms[valid]), collapse = ";")
    },
    character(1)
  )
  gene_ratio_num <- parse_ratio_vec(dt$GeneRatio)
  bg_ratio_num   <- parse_ratio_vec(dt$BgRatio)

  data.table(
    pathway_id    = dt$ID,
    term_name     = dt$Description,
    test_set      = test_set,
    count_mode    = count_mode,
    gene_ids      = gene_ids_vec,
    gene_names    = gene_names_vec,
    gene_count    = dt$Count,
    bg_size       = universe_size,
    gene_ratio    = gene_ratio_num,
    bg_ratio      = bg_ratio_num,
    p_value       = dt$pvalue,
    p_adjust      = dt$p.adjust,
    universe_size = universe_size,
    min_gs        = minGS,
    max_gs        = maxGS,
    ontology      = "KEGG"
  )
}


filter_ora_sig <- function(dt, fdr_cutoff) {
  if (is.null(dt) || nrow(dt) == 0L) return(dt)
  dt[!is.na(p_adjust) & p_adjust <= fdr_cutoff]
}

#------------------------- 读取配置与任务 -------------------------#

config_file <- "config.yaml"
if (!file.exists(config_file)) {
  stop_with("找不到配置文件：", config_file)
}
cfg <- yaml::read_yaml(config_file)

dirs_cfg <- cfg$dirs %||% list()
dir_enrich <- dirs_cfg$enrich %||% "results/08_enrich"
dir_annot  <- dirs_cfg$annot  %||% "results/07_annot"
dir_maps   <- dirs_cfg$maps   %||% "results/03_maps"

task_path <- file.path(dir_enrich, "tasks.tsv")
if (!file.exists(task_path)) {
  stop_with("找不到任务清单：", task_path)
}
tasks <- data.table::fread(task_path, sep = "\t", header = TRUE)
if (nrow(tasks) == 0L) {
  stop_with("tasks.tsv 为空，没有需要执行的 ORA 任务。")
}
required_task_cols <- c("label", "type", "test_set", "test_file",
                        "background_file", "universe_meta", "outdir", "n_deg_all")
miss_task_cols <- setdiff(required_task_cols, colnames(tasks))
if (length(miss_task_cols) > 0L) {
  stop_with("tasks.tsv 缺少列：", paste(miss_task_cols, collapse = ", "))
}

#------------------------- 读取注释资源 -------------------------#

annot_dir <- dir_annot

gene2go_fp        <- file.path(annot_dir, "gene2go.tsv")
gene2pathway_fp   <- file.path(annot_dir, "gene2pathway.tsv")
go_term2name_fp   <- file.path(annot_dir, "go", "term2name.tsv")
go_obsolete_fp    <- file.path(annot_dir, "go", "obsolete_map.tsv")
kegg_term2gene_fp <- file.path(annot_dir, "kegg", "term2gene.tsv")
kegg_term2name_fp <- file.path(annot_dir, "kegg", "term2name.tsv")

if (!file.exists(gene2go_fp))        stop_with("缺少 gene2go.tsv：", gene2go_fp)
if (!file.exists(gene2pathway_fp))   stop_with("缺少 gene2pathway.tsv：", gene2pathway_fp)
if (!file.exists(go_term2name_fp))   stop_with("缺少 go/term2name.tsv：", go_term2name_fp)
if (!file.exists(kegg_term2gene_fp)) stop_with("缺少 kegg/term2gene.tsv：", kegg_term2gene_fp)
if (!file.exists(kegg_term2name_fp)) stop_with("缺少 kegg/term2name.tsv：", kegg_term2name_fp)

gene2go        <- data.table::fread(gene2go_fp, sep = "\t", header = TRUE)
gene2pathway   <- data.table::fread(gene2pathway_fp, sep = "\t", header = TRUE)
go_term2name   <- data.table::fread(go_term2name_fp, sep = "\t", header = TRUE)
kegg_term2gene <- data.table::fread(kegg_term2gene_fp, sep = "\t", header = TRUE)
kegg_term2name <- data.table::fread(kegg_term2name_fp, sep = "\t", header = TRUE)

#------------------------- gene_name 词典（可选） -------------------------#

gene_name_map <- NULL
gene_info_fp <- file.path(dir_maps, "gene_info.tsv")
if (file.exists(gene_info_fp)) {
  gi <- data.table::fread(gene_info_fp, sep = "\t", header = TRUE)
  if (all(c("gene_id", "gene_name") %in% colnames(gi))) {
    gene_name_map <- gi$gene_name
    names(gene_name_map) <- gi$gene_id
    log_msg("INFO", "已加载 gene_id→gene_name 词典，条目数=", length(gene_name_map))
  } else {
    log_msg("WARNING", "gene_info.tsv 不含 gene_id/gene_name 列，忽略 gene_name 映射。")
  }
} else {
  log_msg("INFO", "未检测到 gene_info.tsv，将不输出 gene_names 列（值为 NA）。")
}

#------------------------- 参数：ORA / GSEA -------------------------#

go_minGS <- cfg$go$minGS %||% 10L
go_maxGS <- cfg$go$maxGS %||% 500L
go_padj_method <- cfg$go$p_adjust_method %||% "BH"
go_obsolete_policy <- cfg$go$obsolete_policy %||% "replace_or_consider"

enrich_fdr <- cfg$enrich$fdr %||% 0.05
enrich_output_sig_sorted <- cfg$enrich$output_sig_sorted %||% TRUE

# ORA 正式输出不在 enricher() 内按显著性提前截断；后续统一用 p_adjust 判断显著。
ora_pvalue_cutoff <- 1.0
ora_qvalue_cutoff <- 1.0

kegg_count_mode <- cfg$kegg$count_mode %||% "by_gene"
kegg_padj_method <- cfg$kegg$p_adjust_method %||% "BH"

gsea_cfg <- cfg$gsea %||% list()
gsea_enable <- isTRUE(gsea_cfg$enable)
gsea_score_from <- tolower(gsea_cfg$score_from %||% "stat")
gsea_minGS <- gsea_cfg$minGS %||% 10L
gsea_maxGS <- gsea_cfg$maxGS %||% 500L
gsea_padj_method <- gsea_cfg$p_adjust_method %||% "BH"
gsea_fdr <- gsea_cfg$fdr %||% 0.25
gsea_ontologies <- gsea_cfg$ontologies
if (is.null(gsea_ontologies)) {
  gsea_ontologies <- c("BP")
}

#------------------------- GO obsolete 处理 -------------------------#

if (!file.exists(go_obsolete_fp)) {
  log_msg("WARNING", "未找到 go/obsolete_map.tsv，将跳过 obsolete GO term 处理，直接使用 gene2go。")
  gene2go_clean <- gene2go
} else {
  go_obsolete <- data.table::fread(go_obsolete_fp, sep = "\t", header = TRUE)
  if (nrow(go_obsolete) == 0L) {
    log_msg("WARNING", "go/obsolete_map.tsv 为空，将跳过 obsolete GO term 处理，直接使用 gene2go。")
    gene2go_clean <- gene2go
  } else if (!all(c("old_go_id", "action", "new_go_id") %in% colnames(go_obsolete))) {
    log_msg("WARNING", "go/obsolete_map.tsv 列名不完整（需 old_go_id/action/new_go_id），将跳过 obsolete 处理，直接使用 gene2go。")
    gene2go_clean <- gene2go
  } else {
    obsolete_map <- copy(go_obsolete)
    setnames(obsolete_map, c("old_go_id", "action", "new_go_id"),
             c("old", "action", "new"))
    obsolete_map[, old := as.character(old)]
    obsolete_map[, action := as.character(action)]
    obsolete_map[, new := as.character(new)]

    gene2go_clean <- gene2go[, .(gene_id, go_id)]
    gene2go_clean[, go_id := as.character(go_id)]

    gene2go_clean <- merge(
      gene2go_clean,
      obsolete_map,
      by.x = "go_id",
      by.y = "old",
      all.x = TRUE
    )
    gene2go_clean[, final_go := go_id]

    idx_replace <- !is.na(action) & action == "replace"
    gene2go_clean[idx_replace & nzchar(new), final_go := new]

    idx_consider <- !is.na(action) & action == "consider"
    if (go_obsolete_policy == "replace_or_consider") {
      gene2go_clean[idx_consider & nzchar(new), final_go := new]
    }

    idx_remove <- !is.na(action) & action == "remove"
    gene2go_clean <- gene2go_clean[!idx_remove]
    gene2go_clean <- gene2go_clean[!is.na(final_go) & nzchar(final_go)]
    gene2go_clean <- unique(gene2go_clean[, .(gene_id, go_id = final_go)])

    log_msg("INFO", "GO obsolete 处理完成：原始行数=", nrow(gene2go),
            "，处理后行数=", nrow(gene2go_clean))
  }
}

#------------------------- 组装 TERM2GENE / TERM2NAME -------------------------#

go_term2gene <- gene2go_clean[, .(go_id, gene_id)]
colnames(go_term2gene) <- c("term_id", "gene_id")
colnames(go_term2name) <- c("term_id", "term_name")

kegg_term2gene2 <- kegg_term2gene[, .(pathway_id, gene_id)]
colnames(kegg_term2gene2) <- c("term_id", "gene_id")
colnames(kegg_term2name) <- c("term_id", "term_name")

# GO ontology 从 GO.db 查
go_ids_all <- unique(go_term2name$term_id)
ont_df <- AnnotationDbi::select(
  GO.db,
  keys = go_ids_all,
  keytype = "GOID",
  columns = c("ONTOLOGY")
)
ont_df <- unique(as.data.table(ont_df))
setnames(ont_df, c("GOID", "ONTOLOGY"), c("term_id", "ontology"))
go_ontology_map <- ont_df$ontology
names(go_ontology_map) <- ont_df$term_id

# GO ORA 使用 BP/CC/MF 独立 TERM2GENE，避免三本体混合富集后再拆分。
go_term2gene_annot <- merge(
  go_term2gene,
  unique(ont_df[ontology %in% c("BP", "CC", "MF"), .(term_id, ontology)]),
  by = "term_id",
  all.x = FALSE
)
go_term2name_annot <- merge(
  go_term2name,
  unique(ont_df[ontology %in% c("BP", "CC", "MF"), .(term_id, ontology)]),
  by = "term_id",
  all.x = FALSE
)

go_term2gene_by_ont <- list()
go_term2name_by_ont <- list()
for (ont in c("BP", "CC", "MF")) {
  go_term2gene_by_ont[[ont]] <- go_term2gene_annot[ontology == ont, .(term_id, gene_id)]
  go_term2name_by_ont[[ont]] <- go_term2name_annot[ontology == ont, .(term_id, term_name)]
}

# GO / KEGG 的 gene 宇宙
go_gene_universe   <- unique(go_term2gene$gene_id)
kegg_gene_universe <- unique(kegg_term2gene2$gene_id)

#------------------------- 主循环：按 label 处理 ORA + GSEA -------------------------#

summary_rows <- list()
labels_unique <- unique(tasks$label)

for (lb in labels_unique) {
  tasks_lb <- tasks[label == lb]
  log_msg("INFO", "开始处理 label=", lb, "，任务数=", nrow(tasks_lb))

  meta_path <- unique(tasks_lb$universe_meta)
  if (length(meta_path) != 1L) {
    log_msg("WARNING", "label=", lb, " 的 universe_meta 不唯一，将使用第一项：", meta_path[1])
  }
  meta_path <- meta_path[1]
  if (!file.exists(meta_path)) {
    log_msg("WARNING", "label=", lb, " 缺少 meta.tsv：", meta_path,
            "，summary 中 universe_size/coverage 记为 NA。")
    meta_dt <- NULL
    universe_size_meta <- NA_integer_
    coverage_meta <- NA_real_
  } else {
    meta_dt <- data.table::fread(meta_path, sep = "\t", header = TRUE)
    universe_size_meta <- meta_dt$universe_size[1]
    coverage_meta <- meta_dt$coverage[1]
  }

  # -------- ORA：对每个 test_set 运行 GO + KEGG，并写 by-term 表 --------
  for (i in seq_len(nrow(tasks_lb))) {
    row <- tasks_lb[i]
    label  <- row$label
    type   <- row$type
    test_set <- row$test_set
    test_file <- row$test_file
    bg_file   <- row$background_file
    outdir    <- row$outdir

    ensure_dir(outdir)
    log_msg("INFO", "ORA: label=", label, ", test_set=", test_set,
            "，type=", type, "，test_file=", test_file)

    test_genes <- read_gene_list(test_file)
    bg_genes   <- read_gene_list(bg_file)

    if (length(bg_genes) == 0L) {
      log_msg("WARNING", "label=", label, ", test_set=", test_set,
              " 的背景集合为空，GO/KEGG 将产生空表。")
    }

    universe_go    <- intersect(bg_genes, go_gene_universe)
    universe_kegg  <- intersect(bg_genes, kegg_gene_universe)
    test_go        <- intersect(test_genes, universe_go)
    test_kegg      <- intersect(test_genes, universe_kegg)

    universe_size_go   <- if (!is.na(universe_size_meta)) universe_size_meta else length(universe_go)
    universe_size_kegg <- if (!is.na(universe_size_meta)) universe_size_meta else length(universe_kegg)

    # ---- GO ORA ----
    go_res_list <- list()
    for (ont in c("BP", "CC", "MF")) {
      go_res_list[[ont]] <- format_go_sub(
        dt_sub        = NULL,
        ontology      = ont,
        test_set      = test_set,
        universe_size = universe_size_go,
        minGS         = go_minGS,
        maxGS         = go_maxGS,
        gene_name_map = gene_name_map
      )
    }

    if (length(test_go) >= 1L && length(universe_go) >= go_minGS) {
      for (ont in c("BP", "CC", "MF")) {
        term2gene_ont <- go_term2gene_by_ont[[ont]]
        term2name_ont <- go_term2name_by_ont[[ont]]
        if (is.null(term2gene_ont) || nrow(term2gene_ont) == 0L) {
          next
        }

        ont_gene_universe <- unique(term2gene_ont$gene_id)
        universe_go_ont <- intersect(universe_go, ont_gene_universe)
        test_go_ont <- intersect(test_go, universe_go_ont)

        if (length(test_go_ont) < 1L || length(universe_go_ont) < go_minGS) {
          next
        }

        eg <- tryCatch(
          clusterProfiler::enricher(
            gene          = test_go_ont,
            universe      = universe_go_ont,
            TERM2GENE     = term2gene_ont,
            TERM2NAME     = term2name_ont,
            pAdjustMethod = go_padj_method,
            minGS         = go_minGS,
            maxGS         = go_maxGS,
            pvalueCutoff  = ora_pvalue_cutoff,
            qvalueCutoff  = ora_qvalue_cutoff
          ),
          error = function(e) {
            log_msg("WARNING", "GO ", ont, " enricher 失败：label=", label,
                    ", test_set=", test_set, "；原因：", e$message)
            NULL
          }
        )

        if (!is.null(eg) && nrow(as.data.frame(eg)) > 0L) {
          df <- as.data.table(as.data.frame(eg))
          df[, p_adjust := stats::p.adjust(pvalue, method = go_padj_method)]
          go_res_list[[ont]] <- format_go_sub(
            dt_sub        = df,
            ontology      = ont,
            test_set      = test_set,
            universe_size = universe_size_go,
            minGS         = go_minGS,
            maxGS         = go_maxGS,
            gene_name_map = gene_name_map
          )
        }
      }
    }

    for (ont in c("BP", "CC", "MF")) {
      out_fp <- file.path(outdir, sprintf("GO_%s_by_term_%s.tsv", ont, test_set))
      if (!identical(test_set, "all")) {
        go_res_out <- filter_ora_sig(go_res_list[[ont]], enrich_fdr)
        data.table::fwrite(go_res_out, file = out_fp, sep = "\t",
                           quote = FALSE, na = "NA")
      }
    }

    # ---- KEGG ORA ----
    if (length(test_kegg) >= 1L && length(universe_kegg) >= go_minGS) {
      ek <- tryCatch(
        clusterProfiler::enricher(
          gene          = test_kegg,
          universe      = universe_kegg,
          TERM2GENE     = kegg_term2gene2,
          TERM2NAME     = kegg_term2name,
          pAdjustMethod = kegg_padj_method,
          minGS         = go_minGS,
          maxGS         = go_maxGS,
          pvalueCutoff  = ora_pvalue_cutoff,
          qvalueCutoff  = ora_qvalue_cutoff
        ),
        error = function(e) {
          log_msg("WARNING", "KEGG enricher 失败：label=", label,
                  ", test_set=", test_set, "；原因：", e$message)
          NULL
        }
      )
      if (!is.null(ek) && nrow(as.data.frame(ek)) > 0L) {
        df_k <- as.data.table(as.data.frame(ek))
        kegg_res <- format_kegg_res(
          dt           = df_k,
          test_set     = test_set,
          universe_size= universe_size_kegg,
          minGS        = go_minGS,
          maxGS        = go_maxGS,
          gene_name_map= gene_name_map,
          count_mode   = kegg_count_mode
        )
      } else {
        kegg_res <- format_kegg_res(
          dt           = NULL,
          test_set     = test_set,
          universe_size= universe_size_kegg,
          minGS        = go_minGS,
          maxGS        = go_maxGS,
          gene_name_map= gene_name_map,
          count_mode   = kegg_count_mode
        )
      }
    } else {
      kegg_res <- format_kegg_res(
        dt           = NULL,
        test_set     = test_set,
        universe_size= universe_size_kegg,
        minGS        = go_minGS,
        maxGS        = go_maxGS,
        gene_name_map= gene_name_map,
        count_mode   = kegg_count_mode
      )
    }

    kegg_out_fp <- file.path(outdir, sprintf("KEGG_by_term_%s.tsv", test_set))
    if (!identical(test_set, "all")) {
      kegg_res_out <- filter_ora_sig(kegg_res, enrich_fdr)
      data.table::fwrite(kegg_res_out, file = kegg_out_fp, sep = "\t",
                         quote = FALSE, na = "NA")
    }
  }

  # -------- 生成 GO_sig_up / GO_sig_down（按 label） --------
  outdir_label <- tasks_lb$outdir[1]
  ensure_dir(outdir_label)

  for (set_name in c("up", "down")) {
    go_files <- file.path(
      outdir_label,
      sprintf("GO_%s_by_term_%s.tsv", c("BP", "CC", "MF"), set_name)
    )
    go_dt_list <- lapply(go_files, function(fp) {
      if (!file.exists(fp)) return(NULL)
      dt <- data.table::fread(fp, sep = "\t", header = TRUE)
      if (nrow(dt) == 0L) return(NULL)
      dt
    })
    go_dt_list <- Filter(Negate(is.null), go_dt_list)
    if (length(go_dt_list) == 0L) {
      empty_dt <- data.table(
        term_id = character(),
        term_name = character(),
        test_set = character(),
        ontology = character(),
        gene_ids = character(),
        gene_names = character(),
        gene_count = integer(),
        bg_size = integer(),
        gene_ratio = numeric(),
        bg_ratio = numeric(),
        p_value = numeric(),
        p_adjust = numeric(),
        universe_size = integer(),
        min_gs = integer(),
        max_gs = integer()
      )
      out_sig_fp <- file.path(outdir_label, sprintf("GO_sig_%s.tsv", set_name))
      data.table::fwrite(empty_dt, file = out_sig_fp, sep = "\t",
                         quote = FALSE, na = "NA")
      next
    }
    merged <- data.table::rbindlist(go_dt_list, use.names = TRUE, fill = TRUE)
    sig <- merged[p_adjust <= enrich_fdr]
    if (nrow(sig) > 0L && isTRUE(enrich_output_sig_sorted)) {
      sig <- sig[order(p_adjust, -gene_count)]
    }
    out_sig_fp <- file.path(outdir_label, sprintf("GO_sig_%s.tsv", set_name))
    data.table::fwrite(sig, file = out_sig_fp, sep = "\t",
                       quote = FALSE, na = "NA")
  }

  # -------- ORA 显著数统计（summary 用） --------
  go_sig_all_n <- 0L
  go_sig_up_n <- 0L
  go_sig_down_n <- 0L
  for (set_name in c("all", "up", "down")) {
    fp <- file.path(outdir_label, sprintf("GO_sig_%s.tsv", set_name))
    if (file.exists(fp)) {
      dt <- data.table::fread(fp, sep = "\t", header = TRUE)
      n_sig <- nrow(dt)
      if (set_name == "all") go_sig_all_n <- n_sig
      if (set_name == "up")  go_sig_up_n  <- n_sig
      if (set_name == "down")go_sig_down_n<- n_sig
    }
  }

  kegg_sig_all_n <- 0L
  kegg_sig_up_n  <- 0L
  kegg_sig_down_n<- 0L
  for (set_name in c("all", "up", "down")) {
    fp <- file.path(outdir_label, sprintf("KEGG_by_term_%s.tsv", set_name))
    if (!file.exists(fp)) next
    dt <- data.table::fread(fp, sep = "\t", header = TRUE)
    if (nrow(dt) == 0L) next
    n_sig <- nrow(dt[p_adjust <= enrich_fdr])
    if (set_name == "all") kegg_sig_all_n <- n_sig
    if (set_name == "up")  kegg_sig_up_n  <- n_sig
    if (set_name == "down")kegg_sig_down_n<- n_sig
  }

  # -------- GSEA：若开启且为 RNA label --------
  go_gsea_sig_n <- NA_integer_
  kegg_gsea_sig_n <- NA_integer_

  type_lb <- as.character(tasks_lb$type[1])
  if (gsea_enable && identical(type_lb, "rna")) {
    ranks_fp <- file.path(outdir_label, "gsea_ranks.tsv")
    if (!file.exists(ranks_fp)) {
      log_msg("INFO", "label=", lb, " 未找到 gsea_ranks.tsv，跳过 GSEA。")
    } else {
      log_msg("INFO", "GSEA: label=", lb, "，使用排名文件：", ranks_fp)
      ranks_dt <- data.table::fread(ranks_fp, sep = "\t", header = TRUE)
      if (!all(c("gene_id", "score") %in% colnames(ranks_dt))) {
        log_msg("WARNING", "gsea_ranks.tsv 格式不正确（需要 gene_id, score），跳过 label=", lb)
      } else {
        ranks_dt <- ranks_dt[!is.na(score)]
        ranks_dt <- ranks_dt[!is.na(gene_id) & nzchar(gene_id)]
        if (nrow(ranks_dt) >= 2L) {
          ranks_dt <- ranks_dt[order(-score)]
          geneList <- ranks_dt$score
          names(geneList) <- ranks_dt$gene_id

          gsea_dir <- file.path(outdir_label, "gsea")
          ensure_dir(gsea_dir)
          curves_dir <- file.path(gsea_dir, "curves")
          ensure_dir(curves_dir)

          # geneList：用于复现 GSEA 曲线与命中位置
          geneList_dt <- data.table(
            gene_id = names(geneList),
            score   = as.numeric(geneList)
          )
          data.table::fwrite(geneList_dt, file = file.path(gsea_dir, "geneList.tsv"),
                             sep = "\t", quote = FALSE, na = "NA")
          saveRDS(geneList, file = file.path(gsea_dir, "geneList.rds"))

          # 记录关键参数与运行环境（方便复现）
          params_dt <- data.table(
            param = c("label", "score_from", "minGS", "maxGS", "p_adjust_method", "gsea_fdr", "ontologies"),
            value = c(
              as.character(lb),
              as.character(gsea_score_from),
              as.character(gsea_minGS),
              as.character(gsea_maxGS),
              as.character(gsea_padj_method),
              as.character(gsea_fdr),
              paste(as.character(gsea_ontologies), collapse = ",")
            )
          )
          data.table::fwrite(params_dt, file = file.path(gsea_dir, "gsea_params.tsv"),
                             sep = "\t", quote = FALSE, na = "NA")
          suppressWarnings({
            capture.output(sessionInfo(), file = file.path(gsea_dir, "sessionInfo.txt"))
          })

          # 文件名清洗：避免不同系统/工具链对特殊字符不兼容
          safe_name <- function(x) {
            x2 <- gsub("[:/\\\\\\s]+", "_", as.character(x))
            x2 <- gsub("[^A-Za-z0-9._-]+", "_", x2)
            x2
          }

          # 生成 running ES 曲线数据（与经典 GSEA running enrichment plot 对应）
          build_curve_dt <- function(geneList_vec, geneset_vec, exponent = 1) {
            N <- length(geneList_vec)
            all_genes <- names(geneList_vec)
            hits_flag <- all_genes %in% geneset_vec
            Nh <- sum(hits_flag)
            Nm <- N - Nh
            if (Nh == 0 || Nm == 0) {
              return(data.table(
                rank       = seq_len(N),
                gene_id    = all_genes,
                score      = as.numeric(geneList_vec),
                is_hit     = as.integer(hits_flag),
                running_ES = rep(0, N),
                hit_index  = as.integer(ifelse(hits_flag, cumsum(hits_flag), NA))
              ))
            }
            weights <- abs(geneList_vec)^exponent
            weights_hit <- weights * hits_flag
            sum_hit <- sum(weights_hit)
            if (sum_hit == 0) {
              return(data.table(
                rank       = seq_len(N),
                gene_id    = all_genes,
                score      = as.numeric(geneList_vec),
                is_hit     = as.integer(hits_flag),
                running_ES = rep(0, N),
                hit_index  = as.integer(ifelse(hits_flag, cumsum(hits_flag), NA))
              ))
            }
            step <- ifelse(
              hits_flag,
              weights_hit / sum_hit,
              -1 / Nm
            )
            running_ES <- cumsum(step)
            data.table(
              rank       = seq_len(N),
              gene_id    = all_genes,
              score      = as.numeric(geneList_vec),
              is_hit     = as.integer(hits_flag),
              running_ES = as.numeric(running_ES),
              hit_index  = as.integer(ifelse(hits_flag, cumsum(hits_flag), NA))
            )
          }

          # ---- GSEA-GO ----
          go_gsea_sig_n <- NA_integer_
          if (length(intersect(names(geneList), go_gene_universe)) >= gsea_minGS) {
            go_gsea <- tryCatch(
              clusterProfiler::GSEA(
                geneList      = geneList,
                TERM2GENE     = go_term2gene,
                TERM2NAME     = go_term2name,
                pAdjustMethod = gsea_padj_method,
                minGS         = gsea_minGS,
                maxGS         = gsea_maxGS,
                pvalueCutoff  = 1.0,
                verbose       = FALSE
              ),
              error = function(e) {
                log_msg("WARNING", "GSEA-GO 失败：label=", lb, "；原因：", e$message)
                NULL
              }
            )
            if (!is.null(go_gsea) && nrow(as.data.frame(go_gsea)) > 0L) {
              dt <- as.data.table(as.data.frame(go_gsea))
              dt[, ontology := go_ontology_map[as.character(ID)]]
              dt <- dt[ontology %in% gsea_ontologies]
              if (nrow(dt) > 0L) {
                if (!"qvalues" %in% colnames(dt)) {
                  dt[, qvalues := stats::p.adjust(`p.adjust`, method = "BH")]
                }
                go_gsea_sig_n <- nrow(dt[qvalues <= gsea_fdr])
                go_gsea_out <- data.table(
                  term_id          = dt$ID,
                  term_name        = dt$Description,
                  ontology         = dt$ontology,
                  nes              = dt$NES,
                  enrichment_score = dt$enrichmentScore,
                  p_value          = dt$pvalue,
                  p_adjust         = dt$p.adjust,
                  q_value          = dt$qvalues,
                  core_enriched    = dt$core_enrichment,
                  size             = dt$setSize
                )
                data.table::fwrite(go_gsea_out, file = file.path(gsea_dir, "GO_gsea.tsv"),
                                   sep = "\t", quote = FALSE, na = "NA")
                saveRDS(go_gsea, file = file.path(gsea_dir, "GO_gsea.rds"))

                gs_list_go <- tryCatch(go_gsea@geneSets, error = function(e) NULL)
                if (!is.null(gs_list_go) && length(gs_list_go) > 0L) {
                  gs_dt_go <- data.table(
                    term_id  = names(gs_list_go),
                    gene_ids = vapply(gs_list_go, function(x) paste(unique(x), collapse = ";"), character(1))
                  )
                  data.table::fwrite(gs_dt_go, file = file.path(gsea_dir, "GO_genesets_used.tsv"),
                                     sep = "\t", quote = FALSE, na = "NA")

                  dt_top_go <- copy(dt)
                  dt_top_go[, absNES := abs(NES)]
                  setorder(dt_top_go, qvalues, -absNES)
                  dt_top_go[, absNES := NULL]
                  if (nrow(dt_top_go) > 10L) dt_top_go <- dt_top_go[1:10]

                  for (tid in dt_top_go$ID) {
                    gs_vec <- gs_list_go[[as.character(tid)]]
                    if (is.null(gs_vec) || length(gs_vec) == 0L) next
                    curve_dt <- build_curve_dt(geneList, gs_vec, exponent = 1)
                    out_curve_fp <- file.path(curves_dir, sprintf("GO_%s.tsv", safe_name(tid)))
                    data.table::fwrite(curve_dt, file = out_curve_fp, sep = "\t",
                                       quote = FALSE, na = "NA")
                  }
                } else {
                  empty_gs <- data.table(term_id = character(), gene_ids = character())
                  data.table::fwrite(empty_gs, file = file.path(gsea_dir, "GO_genesets_used.tsv"),
                                     sep = "\t", quote = FALSE, na = "NA")
                }
              } else {
                go_gsea_sig_n <- 0L
                empty_dt <- data.table(
                  term_id          = character(),
                  term_name        = character(),
                  ontology         = character(),
                  nes              = numeric(),
                  enrichment_score = numeric(),
                  p_value          = numeric(),
                  p_adjust         = numeric(),
                  q_value          = numeric(),
                  core_enriched    = character(),
                  size             = integer()
                )
                data.table::fwrite(empty_dt, file = file.path(gsea_dir, "GO_gsea.tsv"),
                                   sep = "\t", quote = FALSE, na = "NA")
                saveRDS(go_gsea, file = file.path(gsea_dir, "GO_gsea.rds"))
                empty_gs <- data.table(term_id = character(), gene_ids = character())
                data.table::fwrite(empty_gs, file = file.path(gsea_dir, "GO_genesets_used.tsv"),
                                   sep = "\t", quote = FALSE, na = "NA")
              }
            }
          }

          # ---- GSEA-KEGG ----
          kegg_gsea_sig_n <- NA_integer_
          if (length(intersect(names(geneList), kegg_gene_universe)) >= gsea_minGS) {
            kegg_gsea <- tryCatch(
              clusterProfiler::GSEA(
                geneList      = geneList,
                TERM2GENE     = kegg_term2gene2,
                TERM2NAME     = kegg_term2name,
                pAdjustMethod = gsea_padj_method,
                minGS         = gsea_minGS,
                maxGS         = gsea_maxGS,
                pvalueCutoff  = 1.0,
                verbose       = FALSE
              ),
              error = function(e) {
                log_msg("WARNING", "GSEA-KEGG 失败：label=", lb, "；原因：", e$message)
                NULL
              }
            )
            if (!is.null(kegg_gsea) && nrow(as.data.frame(kegg_gsea)) > 0L) {
              dtk <- as.data.table(as.data.frame(kegg_gsea))
              if (!"qvalues" %in% colnames(dtk)) {
                dtk[, qvalues := stats::p.adjust(`p.adjust`, method = "BH")]
              }
              kegg_gsea_sig_n <- nrow(dtk[qvalues <= gsea_fdr])
              kegg_gsea_out <- data.table(
                term_id          = dtk$ID,
                term_name        = dtk$Description,
                ontology         = "KEGG",
                nes              = dtk$NES,
                enrichment_score = dtk$enrichmentScore,
                p_value          = dtk$pvalue,
                p_adjust         = dtk$p.adjust,
                q_value          = dtk$qvalues,
                core_enriched    = dtk$core_enrichment,
                size             = dtk$setSize,
                count_mode       = kegg_count_mode
              )
              data.table::fwrite(kegg_gsea_out, file = file.path(gsea_dir, "KEGG_gsea.tsv"),
                                 sep = "\t", quote = FALSE, na = "NA")
              saveRDS(kegg_gsea, file = file.path(gsea_dir, "KEGG_gsea.rds"))

              gs_list_kegg <- tryCatch(kegg_gsea@geneSets, error = function(e) NULL)
              if (!is.null(gs_list_kegg) && length(gs_list_kegg) > 0L) {
                gs_dt_kegg <- data.table(
                  term_id  = names(gs_list_kegg),
                  gene_ids = vapply(gs_list_kegg, function(x) paste(unique(x), collapse = ";"), character(1))
                )
                data.table::fwrite(gs_dt_kegg, file = file.path(gsea_dir, "KEGG_genesets_used.tsv"),
                                   sep = "\t", quote = FALSE, na = "NA")

                dtk_top <- copy(dtk)
                dtk_top[, absNES := abs(NES)]
                setorder(dtk_top, qvalues, -absNES)
                dtk_top[, absNES := NULL]
                if (nrow(dtk_top) > 10L) dtk_top <- dtk_top[1:10]

                for (tid in dtk_top$ID) {
                  gs_vec <- gs_list_kegg[[as.character(tid)]]
                  if (is.null(gs_vec) || length(gs_vec) == 0L) next
                  curve_dt <- build_curve_dt(geneList, gs_vec, exponent = 1)
                  out_curve_fp <- file.path(curves_dir, sprintf("KEGG_%s.tsv", safe_name(tid)))
                  data.table::fwrite(curve_dt, file = out_curve_fp, sep = "\t",
                                     quote = FALSE, na = "NA")
                }
              } else {
                empty_gs <- data.table(term_id = character(), gene_ids = character())
                data.table::fwrite(empty_gs, file = file.path(gsea_dir, "KEGG_genesets_used.tsv"),
                                   sep = "\t", quote = FALSE, na = "NA")
              }
            } else {
              kegg_gsea_sig_n <- 0L
              empty_dt <- data.table(
                term_id          = character(),
                term_name        = character(),
                ontology         = character(),
                nes              = numeric(),
                enrichment_score = numeric(),
                p_value          = numeric(),
                p_adjust         = numeric(),
                q_value          = numeric(),
                core_enriched    = character(),
                size             = integer(),
                count_mode       = character()
              )
              data.table::fwrite(empty_dt, file = file.path(gsea_dir, "KEGG_gsea.tsv"),
                                 sep = "\t", quote = FALSE, na = "NA")
              saveRDS(kegg_gsea, file = file.path(gsea_dir, "KEGG_gsea.rds"))
              empty_gs <- data.table(term_id = character(), gene_ids = character())
              data.table::fwrite(empty_gs, file = file.path(gsea_dir, "KEGG_genesets_used.tsv"),
                                 sep = "\t", quote = FALSE, na = "NA")
            }
          }

        } else {
          log_msg("WARNING", "label=", lb, " 的 GSEA 排名有效行数不足，跳过 GSEA。")
        }
      }
    }
  }

  n_deg_all_lb <- as.character(tasks_lb$n_deg_all[1])

  summary_rows[[lb]] <- data.table(
    label          = lb,
    type           = type_lb,
    n_deg_all      = n_deg_all_lb,
    go_sig_all     = go_sig_all_n,
    go_sig_up      = go_sig_up_n,
    go_sig_down    = go_sig_down_n,
    kegg_sig_all   = kegg_sig_all_n,
    kegg_sig_up    = kegg_sig_up_n,
    kegg_sig_down  = kegg_sig_down_n,
    go_gsea_sig    = go_gsea_sig_n,
    kegg_gsea_sig  = kegg_gsea_sig_n,
    universe_size  = if (!is.na(universe_size_meta)) universe_size_meta else NA_integer_,
    coverage       = if (!is.na(coverage_meta)) coverage_meta else NA_real_
  )
}

#------------------------- 写 summary.tsv -------------------------#

summary_dt <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
summary_fp <- file.path(dir_enrich, "summary.tsv")
data.table::fwrite(summary_dt, file = summary_fp, sep = "\t",
                   quote = FALSE, na = "NA")

log_msg("INFO", "已写出富集汇总表：", summary_fp)
log_msg("INFO", "===== 08_g_enrich.R 执行完成 =====")

