#!/usr/bin/env Rscript
# =============================================================================
# 11a_module_trait_heatmap.R
# WGCNA Module–Trait heatmap（每格显示 r + (FDR)）
#
# 本版更新（按皇上旨意）：
# - 去掉横纵轴外框线（axis.line / ticks 全关）
# - 配色更清新：Mint–White–Peach（固定 [-1,1]）
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 顶部参数区（手动改这里） ---------------------------

INPUT_COR_TSV     <- "results/10_wgcna/assoc/module_trait_cor.tsv"
INPUT_FDR_TSV     <- "results/10_wgcna/assoc/module_trait_padj.tsv"
INPUT_ASSIGN_TSV  <- "results/10_wgcna/modules/gene_module_assignments.tsv" # 可不存在

OUTDIR <- "results/10_wgcna/figures"

PRIMARY_TRAIT_FOR_SORT <- "target_vs_rest"
APPEND_MODULE_SIZE <- TRUE

SHOW_TITLE <- FALSE
MAIN_TITLE <- "Module–trait relationship"

PNG_WIDTH_IN  <- 10
PNG_HEIGHT_IN <- 7
PNG_RES_DPI   <- 300
PDF_WIDTH_IN  <- 10
PDF_HEIGHT_IN <- 7

FILL_LIMIT <- 1

# 清新配色（负相关→正相关）
# Mint(负) - White(0) - Peach(正)
FILL_LOW  <- "#7FD3C8"
FILL_MID  <- "white"
FILL_HIGH <- "#F79D93"

# ------------------------------------------------------------------------------

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}
need_pkg(c("data.table", "ggplot2"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

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

fmt_trait_label_one <- function(x) {
  x <- as.character(x)
  x <- gsub("^group_", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)

  if (tolower(x) %in% c("target vs rest", "target_vs_rest")) return("Foot vs Rest")

  w <- strsplit(tolower(x), "\\s+")[[1]]
  w <- paste0(toupper(substring(w, 1, 1)), substring(w, 2))
  paste(w, collapse = " ")
}
fmt_trait_label_vec <- function(v) vapply(v, fmt_trait_label_one, FUN.VALUE = character(1))

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

  if (anyDuplicated(rn) > 0) stop("行名重复（module 重复），请检查：", path)

  mat <- as.matrix(mat_dt)
  rownames(mat) <- rn
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (anyNA(mat)) warning("矩阵中出现 NA（可能是非数值字符导致），请检查输入：", path)
  mat
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

# =============================================================================
# 主流程
# =============================================================================

safe_mkdir(OUTDIR)

cor_mat <- read_mt_tsv(INPUT_COR_TSV)
fdr_mat <- read_mt_tsv(INPUT_FDR_TSV)

common_rows <- intersect(rownames(cor_mat), rownames(fdr_mat))
common_cols <- intersect(colnames(cor_mat), colnames(fdr_mat))
if (length(common_rows) < 2 || length(common_cols) < 2) {
  stop("cor 与 fdr 矩阵行/列对不上，请检查输入文件：\n",
       " - ", INPUT_COR_TSV, "\n",
       " - ", INPUT_FDR_TSV)
}
cor_mat <- cor_mat[common_rows, common_cols, drop = FALSE]
fdr_mat <- fdr_mat[common_rows, common_cols, drop = FALSE]

module_size <- NULL
if (isTRUE(APPEND_MODULE_SIZE) && file.exists(INPUT_ASSIGN_TSV)) {
  assign_dt <- data.table::fread(INPUT_ASSIGN_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  if ("module_color" %in% colnames(assign_dt)) {
    tab <- table(assign_dt$module_color)
    module_size <- as.integer(tab)
    names(module_size) <- names(tab)
  }
}

sort_trait <- PRIMARY_TRAIT_FOR_SORT
if (!sort_trait %in% colnames(cor_mat)) {
  hint <- grep(sort_trait, colnames(cor_mat), fixed = TRUE, value = TRUE)
  if (length(hint) > 0) sort_trait <- hint[1] else {
    stop("PRIMARY_TRAIT_FOR_SORT 不在 cor 矩阵列名中：", PRIMARY_TRAIT_FOR_SORT,
         "\n可用列名包括：", paste(colnames(cor_mat), collapse = ", "))
  }
}
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
col_lab <- fmt_trait_label_vec(col_lab_raw)

df_long <- expand.grid(
  module = row_lab,
  trait  = col_lab,
  stringsAsFactors = FALSE
)
df_long$r   <- as.vector(cor_mat)
df_long$fdr <- as.vector(fdr_mat)
df_long$label <- paste0(fmt_r_vec(df_long$r), "\n(", fmt_fdr_vec(df_long$fdr), ")")

df_long$module <- factor(df_long$module, levels = row_lab)
df_long$trait  <- factor(df_long$trait,  levels = col_lab)

base_family <- "sans"
theme_set(theme_classic(base_size = 12, base_family = base_family))

p <- ggplot(df_long, aes(x = trait, y = module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3.0, lineheight = 0.95) +
  scale_fill_gradient2(
    limits = c(-FILL_LIMIT, FILL_LIMIT),
    low = FILL_LOW, mid = FILL_MID, high = FILL_HIGH
  ) +
  labs(x = NULL, y = NULL, fill = "Correlation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks  = element_blank(),

    # ✅ 去掉横纵轴外框线（关键）
    axis.line   = element_blank(),

    # 让整体更“清爽”
    plot.margin = margin(8, 8, 8, 8),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

if (isTRUE(SHOW_TITLE)) p <- p + ggtitle(MAIN_TITLE)

png_out <- file.path(OUTDIR, "11a_module_trait_heatmap.png")
pdf_out <- file.path(OUTDIR, "11a_module_trait_heatmap.pdf")

open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI); print(p); invisible(dev.off())
open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN); print(p); invisible(dev.off())

cat("[OK]", png_out, "\n")
cat("[OK]", pdf_out, "\n")

