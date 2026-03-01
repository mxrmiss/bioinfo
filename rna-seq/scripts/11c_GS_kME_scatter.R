#!/usr/bin/env Rscript
# =============================================================================
# 11c_GS_kME_scatter.R
# GS vs kME scatter（论文级：保留方向 + 回归线 + 角落 r/p/n + 可选基因标注）
#
# 设计目标（皇上痛点对症下药）：
# - 不依赖 config.yaml，顶部参数一把梭
# - 换模块/换颜色不容易炸：kME 列名自动兜底匹配
# - 不再输出烦人的 geom_smooth() message
# - ggrepel 可选（没装也能跑）
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 顶部参数区（手动改这里） ---------------------------

INPUT_GS_KME_TSV <- "results/10_wgcna/hub/gene_level_GS_kME.tsv"
OUTDIR <- "results/10_wgcna/figures"

# 要画哪些模块（module_color；不含 "ME" 前缀）
MODULE_COLORS <- c("lightgreen")

# 主 trait 标签（留空则自动从文件 primary_trait 推断并美化；如 target_vs_rest -> Foot vs Rest）
PRIMARY_TRAIT_LABEL <- ""

# 是否使用绝对值（默认 FALSE：保留方向，符合你要“讲故事”的口径）
USE_ABS <- FALSE

# 相关性与 p 值：默认 Pearson（稳且无需额外依赖）
# 若你特别想 bicor，可把 COR_METHOD 改成 "bicor"（需要 WGCNA 包）
COR_METHOD <- "pearson"  # "pearson" / "bicor"

# 每个模块标注多少个基因（按 |kME|+|GS| 排序）
N_LABEL <- 12

# 是否显示标题（默认关）
SHOW_TITLE <- FALSE

# 图片尺寸
PNG_WIDTH_IN  <- 6
PNG_HEIGHT_IN <- 6
PNG_RES_DPI   <- 300
PDF_WIDTH_IN  <- 6
PDF_HEIGHT_IN <- 6

# 输出文件名前缀
OUT_PREFIX <- "11c_GS_kME"

# ------------------------------------------------------------------------------

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}

need_pkg(c("data.table", "ggplot2"))

has_repel <- requireNamespace("ggrepel", quietly = TRUE)
has_wgcna <- requireNamespace("WGCNA", quietly = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
if (has_repel) suppressPackageStartupMessages(library(ggrepel))

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

open_png <- function(path, w, h, res) {
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(path, w, h, units = "in", res = res)
  else png(path, width = w * res, height = h * res, res = res, type = "cairo")
}
open_pdf <- function(path, w, h) {
  if (capabilities("cairo")) cairo_pdf(path, width = w, height = h)
  else pdf(path, width = w, height = h, useDingbats = FALSE)
}

# 把 trait 名美化成更像人话的标签
fmt_trait_label <- function(x) {
  x <- as.character(x)
  x <- gsub("^binary__", "", x)
  x <- gsub("^group_", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)

  if (tolower(x) %in% c("target vs rest", "target vs. rest", "target vs rest ")) return("Foot vs Rest")

  words <- strsplit(tolower(x), "\\s+")[[1]]
  words <- paste0(toupper(substring(words, 1, 1)), substring(words, 2))
  paste(words, collapse = " ")
}

# 自动寻找某个 module_color 对应的 kME 列名
# 首选：kME_MElightgreen（来自 10 脚本：kME_ + ME列名）
guess_kme_col <- function(dt_cols, module_color) {
  # 常规：kME_MElightgreen
  cand1 <- paste0("kME_ME", module_color)
  if (cand1 %in% dt_cols) return(cand1)

  # 兜底：有些人可能输出成 kME_MElightgreen（同上）或 kME_lightgreen（少了 ME）
  cand2 <- paste0("kME_", module_color)
  if (cand2 %in% dt_cols) return(cand2)

  # 再兜底：用正则找包含 module_color 的 kME 列
  pat <- paste0("^kME_.*", module_color, "$")
  hit <- grep(pat, dt_cols, value = TRUE)
  if (length(hit) >= 1) return(hit[1])

  return(NA_character_)
}

# =============================================================================
# 主流程
# =============================================================================

safe_mkdir(OUTDIR)

if (!file.exists(INPUT_GS_KME_TSV)) stop("找不到：", INPUT_GS_KME_TSV)

dt <- data.table::fread(INPUT_GS_KME_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

need_cols <- c("gene_id", "module_color", "GS")
if (!all(need_cols %in% colnames(dt))) {
  stop("gene_level_GS_kME.tsv 缺少必要列：", paste(setdiff(need_cols, colnames(dt)), collapse = ", "))
}

# 自动推断 trait label
trait_label <- PRIMARY_TRAIT_LABEL
if (!nzchar(trait_label)) {
  if ("primary_trait" %in% colnames(dt)) {
    trait_label <- fmt_trait_label(dt$primary_trait[which(!is.na(dt$primary_trait))[1]])
    if (!nzchar(trait_label) || is.na(trait_label)) trait_label <- "Trait"
  } else {
    trait_label <- "Trait"
  }
}

# cor method 检查
cor_method <- tolower(COR_METHOD)
if (cor_method == "bicor" && !has_wgcna) {
  message("[WARN] 你选择了 COR_METHOD=bicor，但当前环境未安装 WGCNA，将自动改用 Pearson。")
  cor_method <- "pearson"
}

base_family <- "Arial"
theme_set(theme_classic(base_size = 14, base_family = base_family))

for (mc in MODULE_COLORS) {
  sub <- dt[dt$module_color == mc, , drop = FALSE]
  if (nrow(sub) < 5) {
    message("[WARN] 模块 ", mc, " 基因数太少（<5），跳过。")
    next
  }

  kme_col <- guess_kme_col(colnames(dt), mc)
  if (is.na(kme_col)) stop("找不到模块 ", mc, " 对应的 kME 列名（期望类似 kME_MElightgreen）。")

  sub$kME <- as.numeric(sub[[kme_col]])
  sub$GS2 <- as.numeric(sub$GS)

  ok <- is.finite(sub$kME) & is.finite(sub$GS2)
  sub <- sub[ok, , drop = FALSE]
  if (nrow(sub) < 5) {
    message("[WARN] 模块 ", mc, " 有效点太少（<5），跳过。")
    next
  }

  if (isTRUE(USE_ABS)) {
    sub$kME_plot <- abs(sub$kME)
    sub$GS_plot  <- abs(sub$GS2)
  } else {
    sub$kME_plot <- sub$kME
    sub$GS_plot  <- sub$GS2
  }

  # 选择要标注的点：综合分数
  sub$score <- abs(sub$kME) + abs(sub$GS2)
  lab <- sub[order(sub$score, decreasing = TRUE), , drop = FALSE]
  lab <- head(lab, n = min(N_LABEL, nrow(lab)))

  # 相关性 + p 值
  if (cor_method == "bicor") {
    r <- as.numeric(WGCNA::bicor(sub$GS_plot, sub$kME_plot, use = "pairwise.complete.obs"))
    # 用 cor.test 给 p 值（bicor 的精确 p 值没必要硬抠；这里用 Pearson p 值做展示即可）
    ct <- suppressWarnings(stats::cor.test(sub$GS_plot, sub$kME_plot, method = "pearson"))
    pval <- ct$p.value
  } else {
    ct <- suppressWarnings(stats::cor.test(sub$GS_plot, sub$kME_plot, method = "pearson"))
    r <- unname(ct$estimate)
    pval <- ct$p.value
  }

  ann_p <- ifelse(pval < 1e-3, format(pval, scientific = TRUE, digits = 2), sprintf("%.3f", pval))
  ann <- sprintf("n = %d\nr = %.2f\np = %s", nrow(sub), r, ann_p)

  # 作图
  p <- ggplot(sub, aes(x = GS_plot, y = kME_plot)) +
    geom_point(size = 2.2, alpha = 0.75) +
    suppressMessages(geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9)) +
    labs(
      x = paste0("GS (", trait_label, ")"),
      y = paste0("kME (to ME", mc, ")")
    ) +
    theme(
      axis.ticks = element_blank()
    )

  # 角落注释
  p <- p + annotate(
    "text",
    x = min(sub$GS_plot, na.rm = TRUE),
    y = max(sub$kME_plot, na.rm = TRUE),
    label = ann,
    hjust = 0, vjust = 1,
    size = 4
  )

  # 标注基因
  if (nrow(lab) > 0) {
    if (has_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = lab,
        aes(x = if (isTRUE(USE_ABS)) abs(GS2) else GS2,
            y = if (isTRUE(USE_ABS)) abs(kME) else kME,
            label = gene_id),
        size = 3.8,
        max.overlaps = Inf,
        min.segment.length = 0
      )
    } else {
      p <- p + geom_text(
        data = lab,
        aes(x = if (isTRUE(USE_ABS)) abs(GS2) else GS2,
            y = if (isTRUE(USE_ABS)) abs(kME) else kME,
            label = gene_id),
        size = 3.0,
        vjust = -0.6
      )
    }
  }

  if (isTRUE(SHOW_TITLE)) p <- p + ggtitle(paste0("ME", mc))

  png_out <- file.path(OUTDIR, paste0(OUT_PREFIX, "_ME", mc, ".png"))
  pdf_out <- file.path(OUTDIR, paste0(OUT_PREFIX, "_ME", mc, ".pdf"))

  open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI); print(p); invisible(dev.off())
  open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN); print(p); invisible(dev.off())

  cat("[OK]", png_out, "\n")
  cat("[OK]", pdf_out, "\n")
}

