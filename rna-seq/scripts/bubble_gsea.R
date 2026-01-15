#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# =============================================================================
# bubble_gsea.R
# GSEA bubble plot (GO / KEGG) for results/08_enrich/<label>/gsea/*.tsv
#
# 输入（由 08_g_enrich.R 产生）：
#   - results/08_enrich/<label>/gsea/GO_gsea.tsv
#   - results/08_enrich/<label>/gsea/KEGG_gsea.tsv
#
# 输出（本脚本产生）：
#   - results/figure/gsea/<label>/bubble/
#       <label>.GO.gsea.bubble.*.pdf/png
#       <label>.KEGG.gsea.bubble.*.pdf/png
#
# 设计目标（顶刊常见风格）：
#   - 白底、无网格
#   - 四边封闭 panel border（不是半开放轴线）
#   - 只在“正负 NES 同时展示”时才画 x=0 参考线
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(scales)
  library(viridis)
})

# =============================================================================
# 1) 用户参数区（皇上只需要改这里）
# =============================================================================

# -------------------------
# [A] 输入输出路径
# -------------------------

# GSEA 结果根目录（相对路径，换机器可直接跑）
ENRICH_DIR <- "results/08_enrich"

# 图片输出根目录（皇上指定：results/figure/gsea）
OUT_BASE_DIR <- "results/figure/gsea"

# 想画哪些 label？
# - 推荐：只画一个，例如 c("foot_vs_other")
# - 也可以写 "ALL"：自动扫描 ENRICH_DIR 下所有含 gsea 结果的 label
LABELS_TO_PLOT <- c("foot_vs_other")

# 要画哪些类型？可选：c("GO","KEGG") 或只画其中一个
KINDS_TO_PLOT <- c("GO", "KEGG")

# -------------------------
# [B] 展示策略（Top N / 正负 NES）
# -------------------------

# 默认展示 Top 20（皇上已同意）
TOP_N <- 20L

# 只研究斧足：默认只画正 NES（富集在 foot 一侧）
# 可选：
#   "pos"  只画 NES > 0
#   "neg"  只画 NES < 0
#   "both" 正负都画（此时会自动画 x=0 参考线）
NES_SIDE <- "pos"

# GSEA 过滤阈值：
# - GSEA 传统上常用 FDR<=0.25 做“探索性发现”（比 ORA 更宽松）
# - 你也可以改成 0.05，更严格但更容易“没结果”
GSEA_FDR_CUTOFF <- 0.25

# -------------------------
# [C] 标注与外观（顶刊风格：封闭边框、少冗余）
# -------------------------

# GO 是否在 term 前加 [BP]/[CC]/[MF] 前缀（你的图里就是这样）
GO_PREFIX_ONTOLOGY <- TRUE

# 通路/术语换行宽度（字符数）
TERM_WRAP_WIDTH <- 55

# 是否显示标题（投稿图一般可以开，也可以关）
SHOW_TITLE <- TRUE

# 字体（尽量用系统安全字体，避免 cairo/pdf 字体坑）
FONT_FAMILY <- "sans"
BASE_SIZE <- 14

# 颜色：viridis 家族（色盲友好）
# 可选： "viridis" "magma" "plasma" "inferno" "cividis"
COLOR_MAP_OPTION <- "viridis"
VIRIDIS_N <- 256

# 点大小映射：setSize 的范围
POINT_SIZE_RANGE <- c(3, 9)

# 输出分辨率
PNG_DPI <- 300

# =============================================================================
# 2) 工具函数（无需修改）
# =============================================================================

log_msg <- function(level = "INFO", ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(ts, "[bubble_gsea]", sprintf("[%s]", level), paste0(...), "\n")
}

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

safe_name <- function(x) {
  x2 <- gsub("[:/\\\\\\s]+", "_", as.character(x))
  x2 <- gsub("[^A-Za-z0-9._-]+", "_", x2)
  x2
}

sanitize_viridis_option <- function(opt) {
  valid_opts <- c("viridis", "magma", "plasma", "inferno", "cividis")
  if (!(opt %in% valid_opts)) {
    log_msg("WARNING", "COLOR_MAP_OPTION 非法：", opt, "；回退为 viridis")
    return("viridis")
  }
  opt
}

# 生成对数色标 breaks（给 FDR 用）
generate_fixed_breaks <- function(min_val, max_val, n = 5) {
  if (!is.finite(min_val) || min_val <= 0) min_val <- 1e-300
  if (!is.finite(max_val) || max_val <= 0) max_val <- 0.25
  if (min_val == max_val) max_val <- max_val * 10
  log_seq <- seq(log10(max_val), log10(min_val), length.out = n)
  10^log_seq
}

# =============================================================================
# 3) 读入 & 整理 GSEA 数据
# =============================================================================

read_gsea_tsv <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  dt <- tryCatch(read_tsv(fp, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt
}

prep_gsea_dt <- function(df, kind = c("GO", "KEGG")) {
  kind <- match.arg(kind)

  # 兼容列名（以 08_g_enrich.R 输出为准）
  # GO/KEGG_gsea.tsv 应包含：term_name / ontology(仅GO) / nes / p_adjust / size(setSize)
  out <- df %>%
    mutate(
      nes = as.numeric(.data$nes),
      fdr = as.numeric(.data$p_adjust),
      setSize = as.numeric(.data$size),
      term_name = as.character(.data$term_name)
    ) %>%
    filter(
      !is.na(nes),
      !is.na(fdr),
      !is.na(setSize),
      is.finite(nes),
      is.finite(fdr),
      is.finite(setSize)
    )

  if (nrow(out) == 0) return(NULL)

  # 只保留显著（按你设定的阈值）
  out <- out %>% filter(fdr <= GSEA_FDR_CUTOFF)
  if (nrow(out) == 0) return(NULL)

  # NES 方向筛选
  if (NES_SIDE == "pos") out <- out %>% filter(nes > 0)
  if (NES_SIDE == "neg") out <- out %>% filter(nes < 0)
  if (nrow(out) == 0) return(NULL)

  # GO 前缀：[BP]/[CC]/[MF]
  if (kind == "GO" && GO_PREFIX_ONTOLOGY && ("ontology" %in% colnames(out))) {
    out <- out %>%
      mutate(
        ontology = as.character(.data$ontology),
        term_show = ifelse(!is.na(ontology) & nzchar(ontology),
                           paste0("[", ontology, "] ", term_name),
                           term_name)
      )
  } else {
    out <- out %>% mutate(term_show = term_name)
  }

  # Top N 选择：按 FDR 升序；同 FDR 时优先 abs(NES) 大的（更“有方向性”）
  out <- out %>%
    arrange(fdr, desc(abs(nes))) %>%
    slice_head(n = TOP_N)

  # y 轴顺序：按 NES 从小到大排列，让 NES 最大的出现在最上面（更直观）
  out <- out %>%
    mutate(
      term_wrapped = str_wrap(term_show, width = TERM_WRAP_WIDTH)
    ) %>%
    arrange(nes) %>%
    mutate(term_wrapped = factor(term_wrapped, levels = unique(term_wrapped)))

  out
}

# =============================================================================
# 4) 作图函数（封闭边框、少冗余）
# =============================================================================

plot_gsea_bubble <- function(df_ready, title_txt) {
  if (is.null(df_ready) || nrow(df_ready) == 0) return(NULL)

  # 颜色：viridis 连续色板（对数变换显示 FDR）
  COLOR_MAP_OPTION2 <- sanitize_viridis_option(COLOR_MAP_OPTION)
  color_vec <- viridis(n = VIRIDIS_N, option = COLOR_MAP_OPTION2, direction = -1)

  fdr_min <- min(df_ready$fdr, na.rm = TRUE)
  fdr_max <- max(df_ready$fdr, na.rm = TRUE)
  plot_breaks <- generate_fixed_breaks(fdr_min, fdr_max, n = 5)
  plot_limits <- c(fdr_min, fdr_max)

  # 是否画 x=0 参考线：
  # - 只有在正负 NES 同时展示时才有必要
  # - 只画正 NES 时，x=0 线会和边框/轴叠一起，显得“多余竖线”
  show_zero_line <- (NES_SIDE == "both") && any(df_ready$nes < 0) && any(df_ready$nes > 0)

  p <- ggplot(df_ready, aes(x = nes, y = term_wrapped)) +
    geom_point(aes(size = setSize, color = fdr), alpha = 0.95) +

    scale_color_gradientn(
      colours = color_vec,
      trans   = "log10",
      breaks  = plot_breaks,
      limits  = plot_limits,
      labels  = label_scientific(digits = 2),
      oob     = squish,
      name    = "FDR"
    ) +

    scale_size_continuous(
      range = POINT_SIZE_RANGE,
      name  = "setSize"
    ) +

    labs(
      x = "NES",
      y = NULL,
      title = if (isTRUE(SHOW_TITLE)) title_txt else NULL
    ) +

    theme_bw(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      # 顶刊常见：无网格
      panel.grid = element_blank(),

      # 关键：封闭边框（四边都有）
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),

      # 不用“半开放轴线”
      axis.line = element_blank(),

      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),

      # 右侧图例更紧凑一点
      legend.key.height = unit(0.8, "cm")
    )

  if (show_zero_line) {
    p <- p + geom_vline(xintercept = 0, linewidth = 0.5)
  }

  p
}

# =============================================================================
# 5) 主流程
# =============================================================================

# label 自动发现
labels_final <- LABELS_TO_PLOT
if (length(labels_final) == 1 && identical(labels_final, "ALL")) {
  # 只扫描 ENRICH_DIR 下一级子目录
  subdirs <- list.dirs(ENRICH_DIR, full.names = FALSE, recursive = FALSE)
  # 只保留存在 gsea/GO_gsea.tsv 或 gsea/KEGG_gsea.tsv 的 label
  keep <- c()
  for (lb in subdirs) {
    go_fp <- file.path(ENRICH_DIR, lb, "gsea", "GO_gsea.tsv")
    kk_fp <- file.path(ENRICH_DIR, lb, "gsea", "KEGG_gsea.tsv")
    if (file.exists(go_fp) || file.exists(kk_fp)) keep <- c(keep, lb)
  }
  labels_final <- unique(keep)
}

log_msg("INFO", "labels to plot: ", paste(labels_final, collapse = ", "))

for (lb in labels_final) {

  out_dir <- file.path(OUT_BASE_DIR, lb, "bubble")
  ensure_dir(out_dir)

  # 逐类处理
  if ("GO" %in% KINDS_TO_PLOT) {
    go_fp <- file.path(ENRICH_DIR, lb, "gsea", "GO_gsea.tsv")
    log_msg("INFO", "Reading GO: ", go_fp)
    go_dt <- read_gsea_tsv(go_fp)
    go_ready <- if (!is.null(go_dt)) prep_gsea_dt(go_dt, kind = "GO") else NULL

    if (!is.null(go_ready)) {
      ttl <- paste0(lb, " | GO GSEA (Top ", TOP_N, ")")
      p <- plot_gsea_bubble(go_ready, ttl)
      if (!is.null(p)) {
        base <- paste0(safe_name(lb), ".GO.gsea.bubble.top", TOP_N, ".", NES_SIDE)
        ggsave(file.path(out_dir, paste0(base, ".pdf")), p,
               width = 11, height = 8, device = cairo_pdf)
        ggsave(file.path(out_dir, paste0(base, ".png")), p,
               width = 11, height = 8, dpi = PNG_DPI, bg = "white")
        log_msg("INFO", "Saved GO bubble: ", file.path(out_dir, base))
      }
    } else {
      log_msg("WARNING", "GO: no data to plot for label=", lb, " (after filtering).")
    }
  }

  if ("KEGG" %in% KINDS_TO_PLOT) {
    kk_fp <- file.path(ENRICH_DIR, lb, "gsea", "KEGG_gsea.tsv")
    log_msg("INFO", "Reading KEGG: ", kk_fp)
    kk_dt <- read_gsea_tsv(kk_fp)
    kk_ready <- if (!is.null(kk_dt)) prep_gsea_dt(kk_dt, kind = "KEGG") else NULL

    if (!is.null(kk_ready)) {
      ttl <- paste0(lb, " | KEGG GSEA (Top ", TOP_N, ")")
      p <- plot_gsea_bubble(kk_ready, ttl)
      if (!is.null(p)) {
        base <- paste0(safe_name(lb), ".KEGG.gsea.bubble.top", TOP_N, ".", NES_SIDE)
        ggsave(file.path(out_dir, paste0(base, ".pdf")), p,
               width = 11, height = 8, device = cairo_pdf)
        ggsave(file.path(out_dir, paste0(base, ".png")), p,
               width = 11, height = 8, dpi = PNG_DPI, bg = "white")
        log_msg("INFO", "Saved KEGG bubble: ", file.path(out_dir, base))
      }
    } else {
      log_msg("WARNING", "KEGG: no data to plot for label=", lb, " (after filtering).")
    }
  }
}

log_msg("INFO", "All done.")

