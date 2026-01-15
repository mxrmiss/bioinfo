#!/usr/bin/env Rscript
# =============================================================================
# 11f_module_network_force.R
# WGCNA 模块网络图（力导向布局，证据型/可复现/论文展示版）
#
# 修复点（2026-01-09+）：
# - adjacency / TOM 在少数环境可能丢 dimnames（rownames/colnames）
# - 本版强制补齐 dimnames，规避筛边 idx 映射失败
#
# 图形规范（2026-01-09+）：
# - 删除横纵轴边框（axis line）
# - 去除圆圈边框（节点无描边）
# - 默认标题关闭
# - 图上不显示 module / edge rule 文本
# - 去除标签旁“黑点/短杠”（关闭 ggrepel 的 segment）
# - 修正类似 Sco04g08230- 末尾多余字符（仅清洗 label，不改 gene_id）
# - hub 节点更突出（无描边，靠 size boost）
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# 0) 顶部参数区（皇上只需要改这里）
# =============================================================================

WGCNA_RDS <- "results/10_wgcna/rds/wgcna_objects.rds"
HUB_TSV   <- "results/10_wgcna/hub/hub_genes_by_module.tsv"
OUTDIR    <- "results/10_wgcna/figures"

MODULE_COLOR <- "lightgreen"

# 显示基因数量
MAX_NODES <- 80
# 视觉上突出核心基因
TOP_HUB_N <- 15

EDGE_METHOD          <- "top_n"   # "quantile" / "top_n"
EDGE_KEEP_QUANTILE   <- 0.85
EDGE_MIN_EDGES       <- 60
EDGE_LOWEST_QUANTILE <- 0.80
# 显示边数量
EDGE_TOP_N           <- 200

RANDOM_SEED <- 20260109

# ---------------- 标签数量开关（保留） ----------------
LABEL_MODE  <- "dense"        # "paper" / "balanced" / "dense"
# 显示基因id
LABEL_TOP_N <- 50     # 若填数字（如 18），覆盖 LABEL_MODE

USE_LABEL_REPEL      <- TRUE
MAX_LABEL_OVERLAPS   <- 50

# ✅ 默认标题关闭
SHOW_TITLE <- FALSE

# ✅ hub 视觉强调（不加描边，只加一点点 size）
HUB_SIZE_BOOST <- 1.35   # 建议 1.2~1.6，越大 hub 越显眼

PNG_WIDTH_IN  <- 7
PNG_HEIGHT_IN <- 7
PNG_RES_DPI   <- 600
PDF_WIDTH_IN  <- 7
PDF_HEIGHT_IN <- 7

EXPORT_NODE_EDGE_TSV <- TRUE

# =============================================================================
# 1) 依赖与工具函数
# =============================================================================

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss) > 0) stop("缺少 R 包：", paste(miss, collapse = ", "), "。请先安装。")
}
need_pkg(c("data.table", "WGCNA", "igraph", "ggplot2"))

has_repel <- requireNamespace("ggrepel", quietly = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(igraph)
  library(ggplot2)
})
if (has_repel) suppressPackageStartupMessages(library(ggrepel))

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

open_png <- function(path, w, h, res) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = w, height = h, units = "in", res = res)
  } else {
    png(filename = path, width = w * res, height = h * res, res = res, type = "cairo")
  }
}
open_pdf <- function(path, w, h) {
  if (capabilities("cairo")) {
    cairo_pdf(file = path, width = w, height = h)
  } else {
    pdf(file = path, width = w, height = h, useDingbats = FALSE)
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

resolve_label_n <- function(mode, override_n) {
  if (!is.na(override_n) && is.numeric(override_n) && override_n >= 1) return(as.integer(override_n))
  mode <- tolower(as.character(mode %||% "paper"))
  if (mode == "paper") return(12L)
  if (mode == "balanced") return(15L)
  if (mode == "dense") return(20L)
  12L
}

# ✅ 关键稳妥：强制补齐矩阵 dimnames
ensure_square_dimnames <- function(mat, names_vec) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) != ncol(mat)) stop("Expected a square matrix, got: ", nrow(mat), "x", ncol(mat))
  if (length(names_vec) != nrow(mat)) {
    stop("Dimnames length mismatch: length(names_vec)=", length(names_vec),
         " but matrix is ", nrow(mat), "x", ncol(mat))
  }
  rownames(mat) <- names_vec
  colnames(mat) <- names_vec
  mat
}

# ✅ 只清洗“显示用标签”，不改真实 gene_id（避免影响匹配/导出）
clean_label_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  # 去掉常见的小点/项目符号（有些字体会渲染成黑点）
  x <- gsub("[•·●◦∙]", "", x, perl = TRUE)
  # 去掉首尾非字母数字点下划线的字符（比如末尾 '-'）
  x <- gsub("^[^A-Za-z0-9_.]+", "", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9_.]+$", "", x, perl = TRUE)
  # 再保险：去掉中间的非常规字符（gene id 正常不该出现）
  x <- gsub("[^A-Za-z0-9_.]", "", x, perl = TRUE)
  x
}

pick_edges_by_quantile <- function(tom,
                                  q_start = 0.95,
                                  min_edges = 60,
                                  q_lowest = 0.80) {
  if (is.null(rownames(tom)) || is.null(colnames(tom))) {
    stop("TOM 必须有 rownames/colnames。")
  }
  n <- nrow(tom)
  if (n < 3) stop("模块基因太少，无法构网。")

  ut <- upper.tri(tom, diag = FALSE)
  w <- tom[ut]
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) stop("TOM 中没有可用的正权重（全 0？）。")

  qs <- seq(q_start, q_lowest, by = -0.01)
  best_idx <- NULL
  best_thr <- NA_real_

  for (q in qs) {
    thr <- as.numeric(quantile(w, probs = q, na.rm = TRUE))
    idx <- which(ut & (tom >= thr), arr.ind = TRUE)

    if (nrow(idx) >= min_edges) {
      best_idx <- idx
      best_thr <- thr
      break
    }
    if (is.null(best_idx) && nrow(idx) > 0) {
      best_idx <- idx
      best_thr <- thr
    }
  }

  if (is.null(best_idx) || nrow(best_idx) == 0) stop("按分位数筛边失败：得到 0 条边。")

  edges <- data.frame(
    from = rownames(tom)[best_idx[, 1]],
    to   = colnames(tom)[best_idx[, 2]],
    weight = tom[best_idx],
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$from < edges$to, , drop = FALSE]
  edges <- edges[is.finite(edges$weight) & edges$weight > 0, , drop = FALSE]

  list(edges = edges, threshold = best_thr)
}

pick_edges_by_topn <- function(tom, top_n = 200) {
  if (is.null(rownames(tom)) || is.null(colnames(tom))) {
    stop("TOM 必须有 rownames/colnames。")
  }
  ut <- upper.tri(tom, diag = FALSE)
  w <- tom[ut]
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) stop("TOM 中没有可用的正权重（全 0？）。")

  ord <- order(w, decreasing = TRUE)
  k <- min(length(ord), as.integer(top_n))
  keep_w <- w[ord[seq_len(k)]]
  thr <- min(keep_w)

  idx <- which(ut & (tom >= thr), arr.ind = TRUE)

  edges <- data.frame(
    from = rownames(tom)[idx[, 1]],
    to   = colnames(tom)[idx[, 2]],
    weight = tom[idx],
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$from < edges$to, , drop = FALSE]
  edges <- edges[order(edges$weight, decreasing = TRUE), , drop = FALSE]
  edges <- head(edges, n = k)

  list(edges = edges, threshold = thr)
}

# =============================================================================
# 2) 读取对象，准备模块基因
# =============================================================================

safe_mkdir(OUTDIR)
set.seed(RANDOM_SEED)

if (!file.exists(WGCNA_RDS)) stop("找不到：", WGCNA_RDS)
obj <- readRDS(WGCNA_RDS)

datExpr <- obj$datExpr
W <- obj$W
chosen_power <- obj$chosen_power
moduleColors <- obj$moduleColors

if (is.null(datExpr) || is.null(W) || is.null(chosen_power) || is.null(moduleColors)) {
  stop("wgcna_objects.rds 缺少 datExpr/W/chosen_power/moduleColors，无法构建网络。")
}

all_genes <- colnames(datExpr)
names(moduleColors) <- all_genes

module_genes <- names(moduleColors)[moduleColors == MODULE_COLOR]
module_genes <- intersect(module_genes, all_genes)
if (length(module_genes) < 10) stop("模块基因数太少：", length(module_genes), "（<10 不建议画网络）")

hub <- NULL
hub_genes <- character()
if (file.exists(HUB_TSV)) {
  hub <- data.table::fread(HUB_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  if (all(c("gene_id", "module_color", "kME_primary_module", "GS") %in% colnames(hub))) {
    hub_sub <- hub[hub$module_color == MODULE_COLOR, , drop = FALSE]
    if (nrow(hub_sub) > 0) {
      hub_sub$abs_score <- abs(hub_sub$kME_primary_module) + abs(hub_sub$GS)
      hub_sub <- hub_sub[order(hub_sub$abs_score, decreasing = TRUE), , drop = FALSE]
      hub_top <- head(hub_sub, n = min(TOP_HUB_N, nrow(hub_sub)))
      hub_genes <- intersect(unique(hub_top$gene_id), module_genes)
    }
  }
}

# =============================================================================
# 3) adjacency + TOM（模块内）并强制补齐 dimnames
# =============================================================================

expr_m <- datExpr[, module_genes, drop = FALSE]
gene_names <- colnames(expr_m)
if (is.null(gene_names) || length(gene_names) != ncol(expr_m)) {
  stop("expr_m 缺少列名（基因名）。请检查 datExpr / module_genes 是否正常。")
}

netType <- W$network$networkType %||% "signed"
corType <- tolower(W$network$corType %||% "bicor")
useObs  <- as.character(W$network$useObs %||% "pairwise.complete.obs")
maxPOut <- as.numeric(W$network$maxPOutliers %||% 0.05)

corFnc <- if (corType == "bicor") "bicor" else "cor"
corOpt <- if (corType == "bicor") list(use = useObs, maxPOutliers = maxPOut) else list(use = useObs)

adj <- WGCNA::adjacency(expr_m,
                        power = as.numeric(chosen_power),
                        type = netType,
                        corFnc = corFnc,
                        corOptions = corOpt)

adj <- ensure_square_dimnames(adj, gene_names)

TOMType <- W$network$TOMType %||% "signed"
tom <- WGCNA::TOMsimilarity(adj, TOMType = TOMType)

tom <- ensure_square_dimnames(tom, gene_names)
diag(tom) <- 0

# =============================================================================
# 4) 选边
# =============================================================================

edge_info <- NULL
if (tolower(EDGE_METHOD) == "top_n") {
  edge_info <- pick_edges_by_topn(tom, top_n = EDGE_TOP_N)
} else {
  edge_info <- pick_edges_by_quantile(tom,
                                      q_start = EDGE_KEEP_QUANTILE,
                                      min_edges = EDGE_MIN_EDGES,
                                      q_lowest = EDGE_LOWEST_QUANTILE)
}

edges <- edge_info$edges
edge_thr <- edge_info$threshold

if (nrow(edges) < 10) stop("可用边太少（<10）。建议降低 EDGE_KEEP_QUANTILE 或减小 EDGE_MIN_EDGES。")

# =============================================================================
# 5) 构图 + 限制规模 + 选择子图
# =============================================================================

g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = gene_names))

comp <- igraph::components(g)
memb <- comp$membership
comp_sizes <- comp$csize

if (length(hub_genes) > 0) {
  hub_in_comp <- tapply(V(g)$name %in% hub_genes, memb, sum)
  best_comp <- as.integer(names(hub_in_comp)[which.max(hub_in_comp)])
} else {
  best_comp <- which.max(comp_sizes)
}

keep_v <- V(g)[memb == best_comp]
g2 <- igraph::induced_subgraph(g, keep_v)

if (igraph::vcount(g2) > MAX_NODES) {
  deg <- igraph::degree(g2)
  vname <- names(deg)
  is_hub <- vname %in% hub_genes
  score <- as.numeric(deg) + 1000 * as.numeric(is_hub)
  keep2 <- names(sort(score, decreasing = TRUE))[seq_len(MAX_NODES)]
  g2 <- igraph::induced_subgraph(g2, keep2)
}

kme_map <- setNames(rep(0, length(gene_names)), gene_names)
if (!is.null(hub) && all(c("gene_id", "module_color", "kME_primary_module") %in% colnames(hub))) {
  hub_sub2 <- hub[hub$module_color == MODULE_COLOR, , drop = FALSE]
  if (nrow(hub_sub2) > 0) {
    tmp <- setNames(hub_sub2$kME_primary_module, hub_sub2$gene_id)
    kme_map[names(tmp)] <- tmp
  }
}
V(g2)$kME <- kme_map[V(g2)$name]
V(g2)$kME[is.na(V(g2)$kME)] <- 0
V(g2)$is_hub <- V(g2)$name %in% hub_genes
V(g2)$degree <- igraph::degree(g2)

# =============================================================================
# 6) 标签
# =============================================================================

label_n <- resolve_label_n(LABEL_MODE, LABEL_TOP_N)
label_n <- min(label_n, igraph::vcount(g2))

nodes_df <- data.frame(
  name = V(g2)$name,
  is_hub = V(g2)$is_hub,
  kME = V(g2)$kME,
  degree = V(g2)$degree,
  stringsAsFactors = FALSE
)
nodes_df$score <- 1000 * as.numeric(nodes_df$is_hub) + abs(nodes_df$kME) + 0.05 * nodes_df$degree
nodes_df <- nodes_df[order(nodes_df$score, decreasing = TRUE), , drop = FALSE]
label_genes <- head(nodes_df$name, n = label_n)

# =============================================================================
# 7) 布局与作图（无轴线/无描边/默认无标题/不显示 module 与 edge rule）
# =============================================================================

layout_xy <- igraph::layout_with_fr(g2)
colnames(layout_xy) <- c("x", "y")

png_out <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".png"))
pdf_out <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".pdf"))
node_out <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".nodes.tsv"))
edge_out <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".edges.tsv"))

base_family <- "sans"
theme_set(theme_classic(base_size = 14, base_family = base_family))

ed <- igraph::as_data_frame(g2, what = "edges")
ed$weight <- E(g2)$weight

nd <- data.frame(
  name = V(g2)$name,
  x = layout_xy[, 1],
  y = layout_xy[, 2],
  kME = V(g2)$kME,
  degree = V(g2)$degree,
  is_hub = V(g2)$is_hub,
  stringsAsFactors = FALSE
)

# 节点基础大小：优先用 |kME|，否则用 degree
nd$base_size <- abs(nd$kME)
if (all(nd$base_size == 0)) nd$base_size <- nd$degree

# ✅ hub size boost（不加描边）
nd$plot_size <- nd$base_size * ifelse(nd$is_hub, as.numeric(HUB_SIZE_BOOST), 1.0)

# ✅ label 文本清洗（只影响显示，不改 name）
nd$label_raw <- ifelse(nd$name %in% label_genes, nd$name, "")
nd$label <- ifelse(nd$label_raw != "", clean_label_text(nd$label_raw), "")

pos <- nd[, c("name", "x", "y")]
ed <- merge(ed, pos, by.x = "from", by.y = "name", all.x = TRUE)
ed <- merge(ed, pos, by.x = "to", by.y = "name", all.x = TRUE, suffixes = c(".from", ".to"))

p <- ggplot() +
  geom_segment(
    data = ed,
    aes(x = x.from, y = y.from, xend = x.to, yend = y.to, linewidth = weight),
    alpha = 0.35,
    color = "grey40"
  ) +
  geom_point(
    data = nd,
    aes(x = x, y = y, size = plot_size),
    shape = 16,
    color = MODULE_COLOR,
    alpha = 0.9
  ) +
  scale_linewidth(range = c(0.2, 1.4)) +
  scale_size(range = c(2.5, 9)) +
  labs(
    # ✅ 即使打开标题，也不出现 module/edge rule（只给一个中性标题）
    title = if (isTRUE(SHOW_TITLE)) "Co-expression network" else NULL,
    subtitle = NULL
  ) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

# ✅ 关键：去除 ggrepel 的 segment（黑点/短杠的主要来源）
if (isTRUE(USE_LABEL_REPEL) && has_repel) {
  p <- p + ggrepel::geom_text_repel(
    data = nd[nd$label != "", , drop = FALSE],
    aes(x = x, y = y, label = label),
    size = 4,
    min.segment.length = Inf,  # 强制不画 segment
    segment.color = NA,        # 双保险：彻底关闭 segment
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = MAX_LABEL_OVERLAPS
  )
} else {
  p <- p + geom_text(
    data = nd[nd$label != "", , drop = FALSE],
    aes(x = x, y = y, label = label),
    size = 3.6,
    vjust = -0.7
  )
}

open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI); print(p); dev.off()
open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN); print(p); dev.off()

cat("[OK]", png_out, "\n")
cat("[OK]", pdf_out, "\n")

# =============================================================================
# 8) 导出 node/edge TSV
# =============================================================================

if (isTRUE(EXPORT_NODE_EDGE_TSV)) {
  nd_out <- data.frame(
    node = V(g2)$name,
    module_color = MODULE_COLOR,
    kME = V(g2)$kME,
    degree = V(g2)$degree,
    is_hub = V(g2)$is_hub,
    is_labeled = V(g2)$name %in% label_genes,
    stringsAsFactors = FALSE
  )

  ed_out <- igraph::as_data_frame(g2, what = "edges")
  ed_out$weight <- E(g2)$weight
  colnames(ed_out)[1:2] <- c("from", "to")

  data.table::fwrite(nd_out, node_out, sep = "\t", quote = FALSE)
  data.table::fwrite(ed_out, edge_out, sep = "\t", quote = FALSE)

  cat("[OK]", node_out, "\n")
  cat("[OK]", edge_out, "\n")
}

cat("[DONE] 11f network finished.\n")

