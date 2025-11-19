#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# 08_g_enrich.R —— GO/KEGG ORA + GSEA 执行模块（vNext + GSEA 最终版）
#
# 职责：
#   1) 读取 config.yaml 与 results/08_enrich/tasks.tsv；
#   2) 读取注释词典（gene2go/gene2pathway/GO term2name/obsolete_map/KEGG term2gene/term2name）；
#   3) 对每个任务（label + test_set）执行 ORA（GO + KEGG），输出统一表头的 by-term 表；
#   4) 为每个 label 生成 GO_sig_all/up/down（GO 三本体纵合 + FDR 过滤）；
#   5) 若存在 gene_id → gene_name 词典，则在 ORA 结果中增加 gene_names 列；
#   6) 若 config.gsea.enable = TRUE，则对每个 RNA label 运行 GSEA：
#        - 读 results/08_enrich/<label>/gsea_ranks.tsv
#        - 对 GO（可配置本体）和 KEGG 执行 GSEA
#        - 输出 gsea/GO_gsea.tsv / gsea/KEGG_gsea.tsv
#   7) 汇总所有 label 的 ORA + GSEA 显著条目数，写出 results/08_enrich/summary.tsv。

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

format_kegg_res <- function(dt, test_set, universe_size,
                            minGS, maxGS, gene_name_map, count_mode) {
  if (is.null(dt) || nrow(dt) == 0L) {
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
      max_gs         = integer(),
      count_mode     = character()
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
    term_id       = dt$ID,
    term_name     = dt$Description,
    test_set      = test_set,
    ontology      = "KEGG",
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
    count_mode    = count_mode
  )
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
if (!file.exists(go_obsolete_fp))    stop_with("缺少 go/obsolete_map.tsv：", go_obsolete_fp)
if (!file.exists(kegg_term2gene_fp)) stop_with("缺少 kegg/term2gene.tsv：", kegg_term2gene_fp)
if (!file.exists(kegg_term2name_fp)) stop_with("缺少 kegg/term2name.tsv：", kegg_term2name_fp)

gene2go        <- data.table::fread(gene2go_fp, sep = "\t", header = TRUE)
gene2pathway   <- data.table::fread(gene2pathway_fp, sep = "\t", header = TRUE)
go_term2name   <- data.table::fread(go_term2name_fp, sep = "\t", header = TRUE)
go_obsolete    <- data.table::fread(go_obsolete_fp, sep = "\t", header = TRUE)
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

#------------------------- GO obsolete 处理 + ontology 映射 -------------------------#

go_minGS <- cfg$go$minGS %||% 10L
go_maxGS <- cfg$go$maxGS %||% 500L
go_padj_method <- cfg$go$p_adjust_method %||% "BH"
go_obsolete_policy <- cfg$go$obsolete_policy %||% "replace_or_consider"

enrich_fdr <- cfg$enrich$fdr %||% 0.05
enrich_output_sig_sorted <- cfg$enrich$output_sig_sorted %||% TRUE

kegg_count_mode <- cfg$kegg$count_mode %||% "by_gene"
kegg_padj_method <- cfg$kegg$p_adjust_method %||% "BH"

# GSEA 配置
gsea_cfg <- cfg$gsea %||% list()
gsea_enable <- isTRUE(gsea_cfg$enable)
gsea_score_from <- tolower(gsea_cfg$score_from %||% "stat")
gsea_minGS <- gsea_cfg$minGS %||% 10L
gsea_maxGS <- gsea_cfg$maxGS %||% 500L
gsea_padj_method <- gsea_cfg$p_adjust_method %||% "BH"
gsea_fdr <- gsea_cfg$fdr %||% 0.25
gsea_ontologies <- gsea_cfg$ontologies
if (is.null(gsea_ontologies)) {
  gsea_ontologies <- c("BP")  # 默认只看 BP
}

# GO obsolete 处理
if (!all(c("old_go_id", "action", "new_go_id") %in% colnames(go_obsolete))) {
  log_msg("WARNING", "go/obsolete_map.tsv 列不完整，将不进行 obsolete 替换。")
  gene2go_clean <- gene2go
} else {
  obsolete_map <- go_obsolete
  colnames(obsolete_map) <- c("old", "action", "new")

  gene2go_clean <- gene2go[, .(gene_id, go_id)]
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

go_term2gene <- gene2go_clean[, .(go_id, gene_id)]
colnames(go_term2gene) <- c("term_id", "gene_id")
colnames(go_term2name) <- c("term_id", "term_name")

kegg_term2gene2 <- kegg_term2gene[, .(pathway_id, gene_id)]
colnames(kegg_term2gene2) <- c("term_id", "gene_id")
colnames(kegg_term2name) <- c("term_id", "term_name")

# GO ontology 从 GO.db 查（term2name 本身没有 ontology 列）
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

# GO / KEGG 的 gene 宇宙
go_gene_universe <- unique(go_term2gene$gene_id)
kegg_gene_universe <- unique(kegg_term2gene2$gene_id)

#------------------------- 主循环：按 label 处理 ORA + GSEA -------------------------#

summary_rows <- list()
labels_unique <- unique(tasks$label)

for (lb in labels_unique) {
  tasks_lb <- tasks[label == lb]
  log_msg("INFO", "开始处理 label=", lb, "，任务数=", nrow(tasks_lb))

  # meta 信息（背景）
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

    # universe_size：优先 meta，其次实际 universe
    universe_size_go <- if (!is.na(universe_size_meta)) universe_size_meta else length(universe_go)
    universe_size_kegg <- if (!is.na(universe_size_meta)) universe_size_meta else length(universe_kegg)

    # ---- GO ORA ----
    go_res_list <- list()
    if (length(test_go) >= 1L && length(universe_go) >= go_minGS) {
      eg <- tryCatch(
        clusterProfiler::enricher(
          gene          = test_go,
          universe      = universe_go,
          TERM2GENE     = go_term2gene,
          TERM2NAME     = go_term2name,
          pAdjustMethod = go_padj_method,
          minGS         = go_minGS,
          maxGS         = go_maxGS,
          qvalueCutoff  = 1.0
        ),
        error = function(e) {
          log_msg("WARNING", "GO enricher 失败：label=", label,
                  ", test_set=", test_set, "；原因：", e$message)
          NULL
        }
      )
      if (!is.null(eg) && nrow(as.data.frame(eg)) > 0L) {
        df <- as.data.table(as.data.frame(eg))
        # 添加 ontology 信息
        df[, ontology := go_ontology_map[as.character(ID)]]
        # 按本体拆分并重新做 BH
        for (ont in c("BP", "CC", "MF")) {
          sub <- df[ontology == ont]
          if (nrow(sub) == 0L) {
            go_res_list[[ont]] <- format_go_sub(
              dt_sub        = NULL,
              ontology      = ont,
              test_set      = test_set,
              universe_size = universe_size_go,
              minGS         = go_minGS,
              maxGS         = go_maxGS,
              gene_name_map = gene_name_map
            )
          } else {
            sub[, p_adjust := stats::p.adjust(pvalue, method = go_padj_method)]
            go_res_list[[ont]] <- format_go_sub(
              dt_sub        = sub,
              ontology      = ont,
              test_set      = test_set,
              universe_size = universe_size_go,
              minGS         = go_minGS,
              maxGS         = go_maxGS,
              gene_name_map = gene_name_map
            )
          }
        }
      } else {
        # 无结果
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
      }
    } else {
      # 样本太少或 universe 太小，直接空表
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
    }

    # 写 GO by-term
    for (ont in c("BP", "CC", "MF")) {
      out_fp <- file.path(outdir, sprintf("GO_%s_by_term_%s.tsv", ont, test_set))
      data.table::fwrite(go_res_list[[ont]], file = out_fp, sep = "\t",
                         quote = FALSE, na = "NA")
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
          qvalueCutoff  = 1.0
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
    data.table::fwrite(kegg_res, file = kegg_out_fp, sep = "\t",
                       quote = FALSE, na = "NA")
  } # end for tasks_lb (各 test_set)

  # -------- 生成 GO_sig_all/up/down（按 label） --------
  outdir_label <- tasks_lb$outdir[1]
  ensure_dir(outdir_label)

  for (set_name in c("all", "up", "down")) {
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
        min_gs