#!/usr/bin/env Rscript
# =============================================================================
# 11e_dendrogram_module_colors.R
# Gene dendrogram + module colors（Dynamic tree cut vs Merged dynamic）
# ✅ 自动只输出“包含指定模块颜色”的那个 block（默认：lightgreen）
# ✅ 新增：输出 pickSoftThreshold 两联图（Scale independence + Mean connectivity）
#    - 只读 10 脚本产物 TSV：results/10_wgcna/power/pickSoftThreshold.fitIndices.tsv
#    - 缺 TSV 直接报错（皇上要求：不做兜底）
# ✅ 自动加线：
#    - 左图红色水平线：scale_free_R2_target（从 wgcna_objects.rds 读取）
#    - 两图竖虚线+点：chosen_power（从 wgcna_objects.rds 读取）
#
# 输出：
# - results/10_wgcna/figures/11e_pickSoftThreshold.png/pdf
# - results/10_wgcna/figures/11e_dendrogram_module_colors.<target>.png/pdf
# - 同目录输出 block 诊断表：11e_block_summary.<target>.tsv
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 顶部参数区（手动改这里） ---------------------------

WGCNA_RDS <- "results/10_wgcna/rds/wgcna_objects.rds"
OUTDIR    <- "results/10_wgcna/figures"

# 10 脚本输出的 pickSoftThreshold 原材料（只用 TSV；缺了就 stop）
SFT_FITINDICES_TSV <- "results/10_wgcna/power/pickSoftThreshold.fitIndices.tsv"

# 你要锁定的模块颜色（支持 "lightgreen" 或 "MElightgreen"）
TARGET_MODULE_COLOR <- "lightgreen"

# 如果多个 block 都含该颜色：选“该颜色基因数最多”的 block
PICK_STRATEGY <- "max_count"  # 可选："max_count" / "first"

SHOW_TITLE <- FALSE
MAIN_TITLE <- "Cluster dendrogram"

# dendrogram 图尺寸
PNG_WIDTH_IN  <- 8
PNG_HEIGHT_IN <- 4
PNG_RES_DPI   <- 300
PDF_WIDTH_IN  <- 8
PDF_HEIGHT_IN <- 4

# pickSoftThreshold 两联图尺寸（皇上指定：每个面板接近方形）
SFT_PNG_WIDTH_IN  <- 8
SFT_PNG_HEIGHT_IN <- 4
SFT_PNG_RES_DPI   <- 300
SFT_PDF_WIDTH_IN  <- 8
SFT_PDF_HEIGHT_IN <- 4

# 数字标注大小（红色 power 数字）
SFT_LABEL_CEX <- 1.10

FONT_FAMILY <- "Arial"

# ------------------------------------------------------------------------------

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}
need_pkg(c("WGCNA", "data.table"))

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

open_png <- function(path, w, h, res) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = w, height = h, units = "in", res = res)
  } else {
    png(path, width = w * res, height = h * res, res = res, type = "cairo")
  }
}

open_pdf <- function(path, w, h) {
  if (capabilities("cairo")) cairo_pdf(path, width = w, height = h)
  else pdf(path, width = w, height = h, useDingbats = FALSE)
}

# =============================================================================
# 0) 基本检查：目录 + 必要输入（皇上要求：缺就直接报错）
# =============================================================================

safe_mkdir(OUTDIR)

if (!file.exists(WGCNA_RDS)) {
  stop("找不到：", WGCNA_RDS, "（请先运行 10_wgcna_modules.R 生成 wgcna_objects.rds）")
}
if (!file.exists(SFT_FITINDICES_TSV)) {
  stop(
    "找不到：", SFT_FITINDICES_TSV, "\n",
    "皇上说过：11e 之前一定会运行 10；请先运行 10_wgcna_modules.R 生成该 TSV。"
  )
}

# =============================================================================
# 1) 读取 wgcna_objects.rds（用于拿到 chosen_power 和 R2 阈值）
# =============================================================================

obj <- readRDS(WGCNA_RDS)

# chosen_power：10 脚本保存的最终 power
chosen_power <- obj$chosen_power %||% NULL
if (is.null(chosen_power) || length(chosen_power) != 1 || is.na(as.numeric(chosen_power))) {
  stop("wgcna_objects.rds 内缺少合法的 chosen_power（请确认 10 脚本是否正常完成并保存对象）。")
}
chosen_power <- as.numeric(chosen_power)

# R2 阈值：10 脚本用于选 power 的阈值（scale_free_R2_target）
r2_target <- NULL
if (!is.null(obj$W) && !is.null(obj$W$power) && !is.null(obj$W$power$scale_free_R2_target)) {
  r2_target <- as.numeric(obj$W$power$scale_free_R2_target)
}
if (is.null(r2_target) || length(r2_target) != 1 || is.na(r2_target)) {
  stop("wgcna_objects.rds 内缺少 W$power$scale_free_R2_target（请确认 10 脚本默认参数/配置是否包含该项）。")
}

# =============================================================================
# 2) 画 pickSoftThreshold 两联图（只读 TSV + 自动线）
# =============================================================================

plot_soft_threshold_from_tsv <- function(tsv_path, r2_target, chosen_power, label_cex = 0.9) {
  fit <- data.table::fread(tsv_path, sep = "\t", header = TRUE, data.table = FALSE)

  need_cols <- c("Power", "SFT.R.sq", "mean.k.")
  miss <- setdiff(need_cols, colnames(fit))
  if (length(miss) > 0) {
    stop(
      "pickSoftThreshold.fitIndices.tsv 缺少列：", paste(miss, collapse = ", "),
      "。当前列名：", paste(colnames(fit), collapse = ", ")
    )
  }

  fit$Power    <- as.numeric(fit$Power)
  fit$SFT.R.sq <- as.numeric(fit$SFT.R.sq)
  fit$mean.k.  <- as.numeric(fit$mean.k.)

  if (anyNA(fit$Power) || anyNA(fit$SFT.R.sq) || anyNA(fit$mean.k.)) {
    stop("fitIndices.tsv 中存在无法转为数值的内容（NA），请检查 TSV。")
  }

  # 为了让面板更“方”，边距不要太肥
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(
    mfrow = c(1, 2),
    mar = c(4.2, 4.6, 2.6, 1.2),
    cex.axis = 1.40,
    cex.lab = 1.53,
    cex.main = 1.66
  )

  par(family = FONT_FAMILY)

  # ---------------- 左图：Scale independence ----------------
  x <- fit$Power
  y <- fit$SFT.R.sq

  # 让顶部不贴边，避免红字挤到框线上
  ylim_top <- min(1, max(y, r2_target) + 0.06)
  ylim_bot <- min(y, 0, r2_target) - 0.06

  plot(
    x, y,
    xlab = "Soft threshold (power)",
    ylab = "Scale-free topology fit (signed R²)",
    type = "n",
    main = "Scale independence",
    ylim = c(ylim_bot, ylim_top)
  )

  # 红色水平线：R2 阈值
  abline(h = r2_target, col = "red", lwd = 1)

  # 竖虚线：chosen power（两幅图都画）
  abline(v = chosen_power, col = "grey40", lty = 2, lwd = 1)

  # 红色数字标注（power）
  text(x, y, labels = fit$Power, col = "red", cex = label_cex)


  # ---------------- 右图：Mean connectivity ----------------
  x2 <- fit$Power
  y2 <- fit$mean.k.

  plot(
    x2, y2,
    xlab = "Soft threshold (power)",
    ylab = "Mean connectivity",
    type = "n",
    main = "Mean connectivity",
    ylim = c(min(y2) - 0.03 * diff(range(y2)), max(y2) + 0.03 * diff(range(y2)))
  )

  abline(v = chosen_power, col = "grey40", lty = 2, lwd = 1)
  text(x2, y2, labels = fit$Power, col = "red", cex = label_cex)

}

sft_png_out <- file.path(OUTDIR, "11e_pickSoftThreshold.png")
sft_pdf_out <- file.path(OUTDIR, "11e_pickSoftThreshold.pdf")

open_png(sft_png_out, SFT_PNG_WIDTH_IN, SFT_PNG_HEIGHT_IN, SFT_PNG_RES_DPI)
plot_soft_threshold_from_tsv(SFT_FITINDICES_TSV, r2_target = r2_target, chosen_power = chosen_power, label_cex = SFT_LABEL_CEX)
dev.off()

open_pdf(sft_pdf_out, SFT_PDF_WIDTH_IN, SFT_PDF_HEIGHT_IN)
plot_soft_threshold_from_tsv(SFT_FITINDICES_TSV, r2_target = r2_target, chosen_power = chosen_power, label_cex = SFT_LABEL_CEX)
dev.off()

# =============================================================================
# 3) dendrogram + module colors（你原来的 11e 主逻辑：基本不动）
# =============================================================================

net <- obj$net
if (is.null(net)) stop("wgcna_objects.rds 内缺少 net。")
if (is.null(net$dendrograms)) stop("net 内缺少 dendrograms。")
if (is.null(net$blockGenes)) stop("net 内缺少 blockGenes（无法稳妥对齐 block）。")

datExpr <- obj$datExpr
if (is.null(datExpr)) stop("wgcna_objects.rds 内缺少 datExpr（用于取基因名）。")

all_genes <- colnames(datExpr)
if (is.null(all_genes) || length(all_genes) < 10) stop("datExpr 基因列名异常。")

# 颜色对象：可能是 vector（全基因）或 list（按 block）
merged_raw   <- net$colors
unmerged_raw <- net$unmergedColors
if (is.null(merged_raw) || is.null(unmerged_raw)) stop("net 内缺少 colors/unmergedColors。")

# 目标模块颜色规范化
target <- as.character(TARGET_MODULE_COLOR)
target <- gsub("^ME", "", target, ignore.case = TRUE)
target <- tolower(target)

# 为每个 block 构建“对齐后的颜色向量”
dends   <- net$dendrograms
blocks  <- net$blockGenes
n_blocks <- length(dends)

get_block_vec <- function(x_raw, idx, block_i) {
  if (is.list(x_raw)) {
    if (length(x_raw) < block_i) stop("colors list 长度小于 block 数：", block_i)
    v <- x_raw[[block_i]]
  } else {
    v <- x_raw[idx]
  }
  v
}

block_info <- data.frame(
  block = seq_len(n_blocks),
  n_genes = NA_integer_,
  has_target = FALSE,
  target_count = 0L,
  stringsAsFactors = FALSE
)

aligned_colors_list <- vector("list", n_blocks)

for (i in seq_len(n_blocks)) {
  dend <- dends[[i]]
  idx  <- blocks[[i]]

  if (is.null(idx) || length(idx) < 2) next
  idx <- as.integer(idx)

  genes_block <- all_genes[idx]
  block_info$n_genes[i] <- length(genes_block)

  mer <- get_block_vec(merged_raw, idx, i)
  unm <- get_block_vec(unmerged_raw, idx, i)

  mer_col <- WGCNA::labels2colors(mer)
  unm_col <- WGCNA::labels2colors(unm)

  names(mer_col) <- genes_block
  names(unm_col) <- genes_block

  lab <- dend$labels
  if (is.null(lab)) {
    lab <- genes_block
    dend$labels <- lab
  }

  mer2 <- mer_col[lab]
  unm2 <- unm_col[lab]

  if (any(is.na(mer2)) || any(is.na(unm2))) {
    if (length(lab) == length(genes_block)) {
      mer2 <- unname(mer_col)
      unm2 <- unname(unm_col)
      if (length(mer2) != length(lab)) stop("block ", i, "：颜色与 dend labels 仍不匹配。")
    } else {
      stop("block ", i, "：dendrogram labels 与 blockGenes 无法对齐，请检查 WGCNA 输出结构。")
    }
  }

  aligned_colors_list[[i]] <- list(dend = dend, unm = unm2, mer = mer2)

  has_t <- any(tolower(mer2) == target, na.rm = TRUE)
  block_info$has_target[i] <- has_t
  if (has_t) block_info$target_count[i] <- sum(tolower(mer2) == target, na.rm = TRUE)
}

diag_out <- file.path(OUTDIR, paste0("11e_block_summary.", target, ".tsv"))
data.table::fwrite(block_info, diag_out, sep = "\t", quote = FALSE)

# 选定要画的 block
cand <- which(block_info$has_target)

pick_block <- function() {
  if (length(cand) == 0) return(NA_integer_)
  if (length(cand) == 1) return(cand[1])
  if (tolower(PICK_STRATEGY) == "first") return(cand[1])
  tc <- block_info$target_count[cand]
  cand[which.max(tc)]
}

best_i <- pick_block()

if (is.na(best_i)) {
  stop(
    "在任何 block 的 Merged dynamic 中都找不到模块颜色：", target, "\n",
    "已输出诊断表：", diag_out, "\n",
    "请把 TARGET_MODULE_COLOR 改成实际存在的颜色名（例如 turquoise/blue/brown/lightgreen 等）。"
  )
}

x <- aligned_colors_list[[best_i]]
if (is.null(x)) stop("选中的 block 无有效对齐数据：block=", best_i)

dend_best <- x$dend
unm_best  <- x$unm
mer_best  <- x$mer

# 作图 & 输出（Dynamic tree cut vs Merged dynamic）
png_out <- file.path(OUTDIR, paste0("11e_dendrogram_module_colors.", target, ".png"))
pdf_out <- file.path(OUTDIR, paste0("11e_dendrogram_module_colors.", target, ".pdf"))

plot_one <- function() {
  par(family = FONT_FAMILY)
  WGCNA::plotDendroAndColors(
    dendro = dend_best,
    colors = cbind(unm_best, mer_best),
    groupLabels = c("Dynamic tree cut", "Merged dynamic"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = if (isTRUE(SHOW_TITLE)) MAIN_TITLE else ""
  )
}

open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI); plot_one(); dev.off()
open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN); plot_one(); dev.off()

cat("[OK] pickSoftThreshold plot:\n")
cat("[OK] ", sft_png_out, " (R2_target=", r2_target, ", chosen_power=", chosen_power, ")\n", sep = "")
cat("[OK] ", sft_pdf_out, " (R2_target=", r2_target, ", chosen_power=", chosen_power, ")\n", sep = "")
cat("[OK] Selected block:", best_i, " (target=", target, ", count=", block_info$target_count[best_i], ")\n", sep = "")
cat("[OK] ", png_out, "\n", sep = "")
cat("[OK] ", pdf_out, "\n", sep = "")
cat("[OK] ", diag_out, "\n", sep = "")

