#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

CONFIG_YAML <- "config.yaml"

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}

need_pkg(c("data.table", "ggplot2", "yaml"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_mkdir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_yaml_safe <- function(path) {
  if (!file.exists(path)) stop("找不到配置文件：", path)
  yaml::read_yaml(path)
}

require_cfg <- function(x, name) {
  if (is.null(x)) stop("配置缺失：", name)
  x
}

log_info <- function(...) {
  cat(sprintf("[%s] [INFO] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

open_png <- function(path, width_in, height_in, res_dpi) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width_in, height = height_in, units = "in", res = res_dpi)
  } else {
    png(filename = path, width = width_in * res_dpi, height = height_in * res_dpi, res = res_dpi, type = "cairo")
  }
}

open_pdf <- function(path, width_in, height_in) {
  if (capabilities("cairo")) {
    cairo_pdf(file = path, width = width_in, height = height_in)
  } else {
    pdf(file = path, width = width_in, height = height_in, useDingbats = FALSE)
  }
}

read_mt_tsv <- function(path) {
  if (!file.exists(path)) stop("文件不存在：", path)
  dt <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  if (ncol(dt) < 2) stop("表太小：", path)

  if ("module" %in% colnames(dt)) {
    rn <- as.character(dt$module)
    mat_dt <- dt[, setdiff(colnames(dt), "module"), drop = FALSE]
  } else {
    rn <- as.character(dt[[1]])
    mat_dt <- dt[, -1, drop = FALSE]
  }

  if (anyDuplicated(rn) > 0) stop("模块行名重复，请检查：", path)

  mat <- as.matrix(mat_dt)
  rownames(mat) <- rn
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat
}

read_trait_metadata <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  dt <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  if (!nrow(dt)) return(NULL)
  cn <- colnames(dt)
  if (!"trait" %in% cn) {
    alt <- intersect(c("trait_name", "name", "column", "trait_id"), cn)
    if (length(alt) > 0) colnames(dt)[match(alt[1], cn)] <- "trait"
  }
  if (!"trait" %in% colnames(dt)) return(NULL)
  dt$trait <- as.character(dt$trait)
  dt
}

fmt_r_vec <- function(x) {
  x <- as.numeric(x)
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- sprintf("%.2f", x[ok])
  out
}

fmt_fdr_vec <- function(x) {
  x <- as.numeric(x)
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  if (!any(ok)) return(out)

  small <- ok & x < 1e-3
  mid   <- ok & !small

  if (any(small)) out[small] <- format(x[small], scientific = TRUE, digits = 2)
  if (any(mid))   out[mid]   <- sprintf("%.3f", x[mid])

  out
}

title_case_words <- function(x) {
  w <- strsplit(tolower(x), "\\s+")[[1]]
  w <- w[nzchar(w)]
  if (length(w) == 0) return("")
  paste(paste0(toupper(substring(w, 1, 1)), substring(w, 2)), collapse = " ")
}

fmt_trait_label_one <- function(x, label_map = NULL, onehot_prefix = "group_") {
  x_raw <- as.character(x)

  if (!is.null(label_map) && x_raw %in% names(label_map)) {
    return(as.character(label_map[[x_raw]]))
  }

  x2 <- x_raw
  if (!is.null(onehot_prefix) && nzchar(onehot_prefix)) {
    x2 <- gsub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\]|\\\\)", "\\\\\\1", onehot_prefix)), "", x2)
  }
  x2 <- gsub("_", " ", x2)
  x2 <- trimws(x2)
  title_case_words(x2)
}

fmt_trait_label_vec <- function(v, label_map = NULL, onehot_prefix = "group_") {
  vapply(v, fmt_trait_label_one, FUN.VALUE = character(1), label_map = label_map, onehot_prefix = onehot_prefix)
}

infer_trait_type <- function(trait_name, meta_row = NULL, onehot_prefix = "group_") {
  if (!is.null(meta_row) && nrow(meta_row) > 0) {
    for (col_nm in c("trait_type", "type", "kind", "category")) {
      if (col_nm %in% colnames(meta_row)) {
        val <- tolower(trimws(as.character(meta_row[[col_nm]][1])))
        if (nzchar(val)) return(val)
      }
    }
  }

  nm <- as.character(trait_name)
  if (!is.null(onehot_prefix) && nzchar(onehot_prefix) && startsWith(nm, onehot_prefix)) {
    return("group_onehot")
  }
  if (grepl("(^|_)binary$", nm, ignore.case = TRUE)) return("binary")
  if (grepl("(^|_)(score|index|axis)$", nm, ignore.case = TRUE)) return("score")
  if (grepl("(^|_)(order|rank)$", nm, ignore.case = TRUE)) return("ordered")
  "unknown"
}

is_discrete_vector <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  length(ux) <= 2
}

pick_sort_trait <- function(sort_trait, cor_mat, fdr_mat, trait_meta = NULL, onehot_prefix = "group_") {
  cols <- colnames(cor_mat)

  if (!is.null(sort_trait) && nzchar(as.character(sort_trait))) {
    st <- as.character(sort_trait)
    if (st %in% cols) return(st)

    hit <- grep(st, cols, fixed = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])

    stop("sort_trait 不在 module_trait_cor.tsv 列名中：", st,
         "\n可用列名：", paste(cols, collapse = ", "))
  }

  score_df <- data.frame(
    trait = cols,
    best_abs_cor = apply(abs(cor_mat), 2, function(v) suppressWarnings(max(v, na.rm = TRUE))),
    best_fdr = apply(fdr_mat, 2, function(v) suppressWarnings(min(v, na.rm = TRUE))),
    stringsAsFactors = FALSE
  )

  score_df$best_abs_cor[!is.finite(score_df$best_abs_cor)] <- -Inf
  score_df$best_fdr[!is.finite(score_df$best_fdr)] <- Inf
  score_df$trait_type <- vapply(score_df$trait, function(tt) {
    meta_row <- NULL
    if (!is.null(trait_meta)) meta_row <- trait_meta[trait_meta$trait == tt, , drop = FALSE]
    infer_trait_type(tt, meta_row = meta_row, onehot_prefix = onehot_prefix)
  }, FUN.VALUE = character(1))

  score_df$type_priority <- ifelse(score_df$trait_type %in% c("score", "ordered", "continuous", "numeric"), 1,
                                   ifelse(score_df$trait_type %in% c("binary", "group_onehot", "categorical"), 2, 3))

  score_df$is_onehot <- !is.null(onehot_prefix) & nzchar(onehot_prefix) & startsWith(score_df$trait, onehot_prefix)
  score_df$is_binary_name <- grepl("(^|_)binary$", score_df$trait, ignore.case = TRUE)
  score_df$name_penalty <- ifelse(score_df$is_onehot, 2L, ifelse(score_df$is_binary_name, 1L, 0L))
  score_df$signif_flag <- ifelse(is.finite(score_df$best_fdr) & score_df$best_fdr <= 0.05, 1L, 0L)

  ord <- order(
    -score_df$signif_flag,
    score_df$type_priority,
    score_df$name_penalty,
    score_df$best_fdr,
    -score_df$best_abs_cor,
    score_df$trait
  )
  score_df <- score_df[ord, , drop = FALSE]
  chosen <- score_df$trait[1]

  log_info("Auto-selected sort_trait: ", chosen,
           " (type=", score_df$trait_type[1],
           ", best_fdr=", ifelse(is.finite(score_df$best_fdr[1]), format(score_df$best_fdr[1], scientific = TRUE, digits = 3), "NA"),
           ", best_abs_cor=", ifelse(is.finite(score_df$best_abs_cor[1]), sprintf("%.3f", score_df$best_abs_cor[1]), "NA"), ")")

  chosen
}

cfg_all <- read_yaml_safe(CONFIG_YAML)

cfg_plot_all <- require_cfg(cfg_all$wgcna_plot, "wgcna_plot")
cfg_plot <- require_cfg(cfg_plot_all$module_trait_heatmap, "wgcna_plot.module_trait_heatmap")
cfg_wgcna <- require_cfg(cfg_all$wgcna, "wgcna")

INPUT_COR_TSV <- require_cfg(cfg_plot$input_cor_tsv, "wgcna_plot.module_trait_heatmap.input_cor_tsv")
INPUT_FDR_TSV <- require_cfg(cfg_plot$input_fdr_tsv, "wgcna_plot.module_trait_heatmap.input_fdr_tsv")
INPUT_ASSIGN_TSV <- cfg_plot$input_assign_tsv %||% NULL
INPUT_TRAIT_META_TSV <- cfg_plot$input_trait_meta_tsv %||% NULL
OUTDIR <- require_cfg(cfg_plot$outdir, "wgcna_plot.module_trait_heatmap.outdir")

if (is.null(INPUT_TRAIT_META_TSV) || !nzchar(as.character(INPUT_TRAIT_META_TSV))) {
  inferred_meta <- file.path(dirname(dirname(INPUT_COR_TSV)), "traits", "trait_metadata.tsv")
  if (file.exists(inferred_meta)) INPUT_TRAIT_META_TSV <- inferred_meta
}

APPEND_MODULE_SIZE <- isTRUE(cfg_plot$append_module_size %||% TRUE)
SHOW_TITLE <- isTRUE(cfg_plot$show_title %||% FALSE)
MAIN_TITLE <- as.character(cfg_plot$main_title %||% "Module-trait relationship")

PNG_WIDTH_IN  <- as.numeric(require_cfg(cfg_plot$png_width_in,  "wgcna_plot.module_trait_heatmap.png_width_in"))
PNG_HEIGHT_IN <- as.numeric(require_cfg(cfg_plot$png_height_in, "wgcna_plot.module_trait_heatmap.png_height_in"))
PNG_RES_DPI   <- as.numeric(require_cfg(cfg_plot$png_res_dpi,   "wgcna_plot.module_trait_heatmap.png_res_dpi"))
PDF_WIDTH_IN  <- as.numeric(require_cfg(cfg_plot$pdf_width_in,  "wgcna_plot.module_trait_heatmap.pdf_width_in"))
PDF_HEIGHT_IN <- as.numeric(require_cfg(cfg_plot$pdf_height_in, "wgcna_plot.module_trait_heatmap.pdf_height_in"))

FILL_LIMIT <- as.numeric(require_cfg(cfg_plot$fill_limit, "wgcna_plot.module_trait_heatmap.fill_limit"))
FILL_LOW   <- as.character(require_cfg(cfg_plot$fill_low,  "wgcna_plot.module_trait_heatmap.fill_low"))
FILL_MID   <- as.character(require_cfg(cfg_plot$fill_mid,  "wgcna_plot.module_trait_heatmap.fill_mid"))
FILL_HIGH  <- as.character(require_cfg(cfg_plot$fill_high, "wgcna_plot.module_trait_heatmap.fill_high"))

SORT_TRAIT <- cfg_plot$sort_trait %||% NULL
HIDE_TRAITS <- unlist(cfg_plot$hide_traits %||% list(), use.names = FALSE)
HIDE_TRAITS <- as.character(HIDE_TRAITS)
TRAIT_LABEL_MAP <- cfg_plot$trait_label_map %||% NULL

FONT_FAMILY <- as.character(cfg_plot$font_family %||% "Arial")
BASE_FONT_SIZE <- as.numeric(cfg_plot$base_font_size %||% 12)
CELL_TEXT_SIZE <- as.numeric(cfg_plot$cell_text_size %||% 3.0)
AXIS_TEXT_SIZE <- as.numeric(cfg_plot$axis_text_size %||% 11)
LEGEND_TITLE_SIZE <- as.numeric(cfg_plot$legend_title_size %||% 11)
LEGEND_TEXT_SIZE <- as.numeric(cfg_plot$legend_text_size %||% 10)
TITLE_SIZE <- as.numeric(cfg_plot$title_size %||% 13)

X_TEXT_ANGLE <- as.numeric(cfg_plot$x_text_angle %||% 45)
X_TEXT_HJUST <- as.numeric(cfg_plot$x_text_hjust %||% 1)
X_TEXT_VJUST <- as.numeric(cfg_plot$x_text_vjust %||% 1)
Y_TEXT_SIZE <- as.numeric(cfg_plot$y_text_size %||% AXIS_TEXT_SIZE)
X_TEXT_SIZE <- as.numeric(cfg_plot$x_text_size %||% AXIS_TEXT_SIZE)

ONEHOT_PREFIX <- cfg_wgcna$traits$group_onehot_prefix %||% cfg_wgcna$traits$onehot_prefix %||% "group_"

safe_mkdir(OUTDIR)

cor_mat <- read_mt_tsv(INPUT_COR_TSV)
fdr_mat <- read_mt_tsv(INPUT_FDR_TSV)
trait_meta <- read_trait_metadata(INPUT_TRAIT_META_TSV)

common_rows <- intersect(rownames(cor_mat), rownames(fdr_mat))
common_cols <- intersect(colnames(cor_mat), colnames(fdr_mat))
if (length(common_rows) < 2 || length(common_cols) < 1) {
  stop("cor 与 fdr 矩阵行列对不上，请检查输入文件。")
}

cor_mat <- cor_mat[common_rows, common_cols, drop = FALSE]
fdr_mat <- fdr_mat[common_rows, common_cols, drop = FALSE]

module_size <- NULL
if (isTRUE(APPEND_MODULE_SIZE) && !is.null(INPUT_ASSIGN_TSV) && file.exists(INPUT_ASSIGN_TSV)) {
  assign_dt <- data.table::fread(INPUT_ASSIGN_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  if ("module_color" %in% colnames(assign_dt)) {
    tab <- table(assign_dt$module_color)
    module_size <- as.integer(tab)
    names(module_size) <- names(tab)
  }
}

sort_trait <- pick_sort_trait(SORT_TRAIT, cor_mat, fdr_mat, trait_meta = trait_meta, onehot_prefix = ONEHOT_PREFIX)
ord <- order(abs(cor_mat[, sort_trait]), decreasing = TRUE, na.last = TRUE)
cor_mat <- cor_mat[ord, , drop = FALSE]
fdr_mat <- fdr_mat[rownames(cor_mat), , drop = FALSE]

row_lab_raw <- rownames(cor_mat)
row_lab <- gsub("^ME", "", row_lab_raw)

if (!is.null(module_size)) {
  nn <- module_size[row_lab]
  nn[is.na(nn)] <- NA_integer_
  row_lab <- ifelse(is.na(nn), row_lab, paste0(row_lab, " (n=", nn, ")"))
}

col_lab_raw <- colnames(cor_mat)
col_lab <- fmt_trait_label_vec(col_lab_raw, label_map = TRAIT_LABEL_MAP, onehot_prefix = ONEHOT_PREFIX)

keep_cols <- rep(TRUE, length(col_lab_raw))
if (length(HIDE_TRAITS) > 0) {
  keep_cols <- !(col_lab_raw %in% HIDE_TRAITS | col_lab %in% HIDE_TRAITS)
}

cor_mat_plot <- cor_mat[, keep_cols, drop = FALSE]
fdr_mat_plot <- fdr_mat[, keep_cols, drop = FALSE]
col_lab_plot <- col_lab[keep_cols]

if (ncol(cor_mat_plot) < 1) {
  stop("隐藏 trait 后没有剩余列可绘图，请检查 hide_traits 设置。")
}

df_long <- expand.grid(
  module = row_lab,
  trait  = col_lab_plot,
  stringsAsFactors = FALSE
)
df_long$r   <- as.vector(cor_mat_plot)
df_long$fdr <- as.vector(fdr_mat_plot)
df_long$label <- paste0(fmt_r_vec(df_long$r), "\n(", fmt_fdr_vec(df_long$fdr), ")")

df_long$module <- factor(df_long$module, levels = row_lab)
df_long$trait  <- factor(df_long$trait, levels = col_lab_plot)

theme_set(theme_classic(base_size = BASE_FONT_SIZE, base_family = FONT_FAMILY))

p <- ggplot(df_long, aes(x = trait, y = module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = CELL_TEXT_SIZE, lineheight = 0.95, family = FONT_FAMILY) +
  scale_fill_gradient2(
    limits = c(-FILL_LIMIT, FILL_LIMIT),
    low = FILL_LOW,
    mid = FILL_MID,
    high = FILL_HIGH
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  theme(
    axis.text.x = element_text(size = X_TEXT_SIZE, angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST, vjust = X_TEXT_VJUST),
    axis.text.y = element_text(size = Y_TEXT_SIZE),
    axis.ticks  = element_blank(),
    axis.line   = element_blank(),
    plot.margin = margin(8, 8, 8, 8),
    legend.title = element_text(size = LEGEND_TITLE_SIZE),
    legend.text  = element_text(size = LEGEND_TEXT_SIZE)
  )

if (isTRUE(SHOW_TITLE)) {
  p <- p + ggtitle(MAIN_TITLE) +
    theme(plot.title = element_text(size = TITLE_SIZE, hjust = 0.5))
}

png_out <- file.path(OUTDIR, "11a_module_trait_heatmap.png")
pdf_out <- file.path(OUTDIR, "11a_module_trait_heatmap.pdf")

open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI)
print(p)
invisible(dev.off())

open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN)
print(p)
invisible(dev.off())

log_info("sort_trait used for ordering modules: ", sort_trait)
cat("[OK]", png_out, "\n")
cat("[OK]", pdf_out, "\n")

