#!/usr/bin/env Rscript
# =============================================================================
# 11f_module_network_force.R
# WGCNA 模块网络图（力导向布局，证据型/可复现/论文展示版）
#
# 2026-02-04 最终增量改进（皇上指定口径）：
# - ✅ 强制读取 data/samples.tsv（缺失直接停止；不提供任何兜底变量）
# - ✅ 顶部参数可指定目标组织（TARGET_GROUP，例如 foot）
# - ✅ 节点颜色：mean(log2(TPM+1))（在目标组织的重复样本上先 log2(TPM+1)，再取平均）
# - ✅ 同时读取 results/10_wgcna/lists/module_<color>.list 第2列 gene_name 作为标签
#   * 若 gene_name 缺失/为空：不显示标签（但节点仍保留，保证网络结构不被悄悄改变）
# - ✅ 保持所有路径为相对路径
#
# 2026-02-04 版面改进（论文观感增强）：
# - ✅ 布局支持 KK / DrL（更均匀，减少节点拥挤）
# - ✅ PCA 旋转布局 + coord_equal，避免竖长条
# - ✅ 图例标题缩短：只显示 "log2(TPM+1)"（mean 与 TARGET_GROUP 在图注说明）
#
# 2026-02-04 终端输出清爽化（本次按皇上要求新增）：
# - ✅ 只加一次坐标系（避免 Coordinate system already present...）
# - ✅ dev.off() 不回显（避免 null device 1）
# - ✅ 小清新连续渐变配色（可在顶部一键切换）
#
# 2026-02-05 皇上追加改进：
# - ✅ 字体参数移动到顶部参数区（FONT_FAMILY / BASE_FONT_SIZE）
# - ✅ 若标签英文名重复：自动追加 -1/-2/...（只作用于显示 label，不改变 gene_id 和网络结构）
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# 0) 顶部参数区（皇上只需要改这里）
# =============================================================================

WGCNA_RDS <- "results/10_wgcna/rds/wgcna_objects.rds"
HUB_TSV   <- "results/10_wgcna/hub/hub_genes_by_module.tsv"
OUTDIR    <- "results/10_wgcna/figures"

MODULE_COLOR <- "lightgreen"

# ✅ 模块 list（第2列 gene_name 用于显示标签）
MODULE_LIST <- paste0("results/10_wgcna/lists/module_", MODULE_COLOR, ".list")

# ✅ TPM 与样本表（强制存在）
TPM_TSV     <- "results/05_matrix/tpms/gene_tpm.tsv"
SAMPLES_TSV <- "data/samples.tsv"

# ✅ 指定目标组织（例如 foot/gill/mantle...）
TARGET_GROUP <- "foot"

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
LABEL_TOP_N <- 50             # 若填数字（如 18），覆盖 LABEL_MODE

USE_LABEL_REPEL      <- TRUE
MAX_LABEL_OVERLAPS   <- 50

# ✅ 默认标题关闭
SHOW_TITLE <- FALSE

# ✅ hub 视觉强调（不加描边，只加一点点 size）
HUB_SIZE_BOOST <- 1.35   # 建议 1.2~1.6

# ✅ 版面与布局（新增）
LAYOUT_METHOD <- "drl"   # "kk" / "drl"
LAYOUT_KK_MAXITER <- 2000
LAYOUT_EXPAND <- 1.35    # 节点间距整体拉开（1.0=不变，建议 1.1~1.4）
ROTATE_LAYOUT <- TRUE    # PCA 旋转，避免竖长条
COORD_EQUAL   <- TRUE    # 等比例坐标，避免拉伸
PANEL_PADDING_FRAC <- 0.08  # 四周留白比例

# ✅ 小清新配色方案（本次新增）
# - mint: 蓝绿薄荷（默认推荐，清爽耐看）
# - lavender: 青紫雾面（更高级）
# - sunrise: 低值浅暖、高值清青（更活泼但仍克制）
COLOR_SCHEME <- "mint"   # "mint" / "lavender" / "sunrise"

# ✅ 字体放到顶部，方便皇上一键修改
FONT_FAMILY   <- "Arial"
BASE_FONT_SIZE <- 14

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

cap_first <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
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

# ✅ 只清洗“显示用标签”，不改真实 gene_id（允许中间 '-'）
clean_label_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[•·●◦∙]", "", x, perl = TRUE)
  x <- gsub("^[^A-Za-z0-9_.-]+", "", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9_.-]+$", "", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9_.-]", "", x, perl = TRUE)
  x
}

# ✅ 重复英文名自动追加 -1/-2/...（只作用于显示 label，不改变 gene_id）
# 修复点：不要用 ave(df$key, ...)（字符会把结果也变成字符），改为对纯数字向量分组计数
make_unique_labels <- function(labels, keys) {
  labels <- as.character(labels)
  keys <- as.character(keys)
  out <- labels

  idx <- which(nzchar(labels))
  if (length(idx) <= 1) return(out)

  df <- data.frame(
    idx = idx,
    label = labels[idx],
    key = keys[idx],
    stringsAsFactors = FALSE
  )

  # 为了可复现：同名时按 key（gene_id）排序决定谁保留原名
  df <- df[order(df$label, df$key), , drop = FALSE]

  # 关键：用纯数字向量做 ave 的输入，保证 rank_in_label 是整数
  df$rank_in_label <- ave(seq_len(nrow(df)), df$label, FUN = seq_along)
  df$rank_in_label <- as.integer(df$rank_in_label)

  df$suffix <- ifelse(df$rank_in_label == 1L, "", paste0("-", df$rank_in_label - 1L))
  df$label_u <- paste0(df$label, df$suffix)

  out[df$idx] <- df$label_u
  out
}

# ✅ PCA 旋转布局，让长轴更接近水平
rotate_layout_pca <- function(xy) {
  xy <- as.matrix(xy)
  if (nrow(xy) < 3) return(xy)
  xyc <- scale(xy, center = TRUE, scale = FALSE)
  s <- svd(xyc)
  rot <- s$v
  xyr <- xyc %*% rot
  xyr
}

# ✅ 计算布局：kk / drl
compute_layout <- function(g, method = "drl", kk_maxiter = 2000) {
  m <- tolower(as.character(method))
  if (m == "kk") {
    return(igraph::layout_with_kk(g, maxiter = as.integer(kk_maxiter)))
  }
  if (m == "drl") {
    return(igraph::layout_with_drl(g))
  }
  igraph::layout_with_fr(g)
}

# ✅ 小清新连续渐变色标（本次新增）
get_color_scale <- function(scheme, legend_title) {
  sc <- tolower(as.character(scheme %||% "mint"))

  if (sc == "lavender") {
    return(scale_color_gradientn(
      colours = c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#6A51A3"),
      name = legend_title
    ))
  }

  if (sc == "sunrise") {
    return(scale_color_gradientn(
      colours = c("#FFF2CC", "#F7D89C", "#BFE6D0", "#66C2A4", "#2CA25F"),
      name = legend_title
    ))
  }

  scale_color_gradientn(
    colours = c("#E9F7F2", "#CDEEE6", "#A7E1D8", "#7CCBD3", "#4DA9D8"),
    name = legend_title
  )
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

# =============================================================================
# 2.1) 读取模块 list（gene_id -> gene_name，用于标签显示）
# =============================================================================

if (!file.exists(MODULE_LIST)) stop("找不到模块 list：", MODULE_LIST)

ml <- data.table::fread(MODULE_LIST, sep = "\t", header = FALSE, data.table = FALSE,
                        fill = TRUE, quote = "", colClasses = "character")
if (ncol(ml) < 1) stop("模块 list 为空或格式异常：", MODULE_LIST)

colnames(ml)[1] <- "gene_id"
if (ncol(ml) >= 2) {
  colnames(ml)[2] <- "gene_name"
} else {
  ml$gene_name <- NA_character_
}

ml$gene_id <- trimws(as.character(ml$gene_id))
ml$gene_name <- trimws(as.character(ml$gene_name))
ml$gene_name[is.na(ml$gene_name)] <- ""

gene_name_map <- setNames(ml$gene_name, ml$gene_id)

# =============================================================================
# 2.2) hub（可选）
# =============================================================================

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
# 5.1) 强制读取 samples.tsv 并从 TPM 中提取目标组织表达（mean(log2(TPM+1))）
# =============================================================================

if (!file.exists(SAMPLES_TSV)) stop("找不到 samples.tsv：", SAMPLES_TSV)
if (!file.exists(TPM_TSV)) stop("找不到 TPM 矩阵：", TPM_TSV)

smp <- data.table::fread(SAMPLES_TSV, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
need_cols <- c("sample", "group")
if (!all(need_cols %in% colnames(smp))) {
  stop("samples.tsv 必须包含列：sample, group。当前列为：", paste(colnames(smp), collapse = ", "))
}

target_samples <- smp$sample[as.character(smp$group) == as.character(TARGET_GROUP)]
target_samples <- unique(as.character(target_samples))
target_samples <- target_samples[nzchar(target_samples)]

if (length(target_samples) == 0) {
  stop("samples.tsv 中找不到 group == '", TARGET_GROUP, "' 的样本。请检查 TARGET_GROUP 或 samples.tsv 内容。")
}

tpm_header <- names(data.table::fread(TPM_TSV, sep = "\t", header = TRUE, nrows = 0, data.table = FALSE, check.names = FALSE))
if (length(tpm_header) < 2) stop("TPM 表头异常（列数 < 2）：", TPM_TSV)
if (!("gene_id" %in% tpm_header)) {
  stop("TPM 文件必须包含列名 gene_id（第一列）。当前表头为：", paste(tpm_header, collapse = ", "))
}

target_cols <- intersect(target_samples, tpm_header)
if (length(target_cols) == 0) {
  stop("TPM 文件中找不到 TARGET_GROUP 对应样本列。\n",
       "TARGET_GROUP='", TARGET_GROUP, "' samples: ", paste(target_samples, collapse = ", "), "\n",
       "TPM header (first 30): ", paste(head(tpm_header, 30), collapse = ", "))
}

tpm_sub <- data.table::fread(TPM_TSV, sep = "\t", header = TRUE, data.table = TRUE,
                            select = c("gene_id", target_cols), check.names = FALSE)

for (cc in target_cols) {
  suppressWarnings(tpm_sub[, (cc) := log2(as.numeric(get(cc)) + 1)])
}
tpm_sub[, expr_color := rowMeans(.SD, na.rm = TRUE), .SDcols = target_cols]

genes_in_graph <- V(g2)$name
tpm_sub <- tpm_sub[gene_id %in% genes_in_graph]

missing_expr <- setdiff(genes_in_graph, tpm_sub$gene_id)
if (length(missing_expr) > 0) {
  stop("TPM 文件中缺少网络子图中的部分基因（无法严格上色）。缺失数量=",
       length(missing_expr), "。示例：", paste(head(missing_expr, 10), collapse = ", "))
}

expr_map <- setNames(as.numeric(tpm_sub$expr_color), as.character(tpm_sub$gene_id))

# =============================================================================
# 6) 标签（仅使用 gene_name；缺失 gene_name 则不显示标签）
# =============================================================================

label_n <- resolve_label_n(LABEL_MODE, LABEL_TOP_N)
label_n <- min(label_n, igraph::vcount(g2))

nodes_df <- data.frame(
  name = V(g2)$name,
  is_hub = V(g2)$is_hub,
  kME = V(g2)$kME,
  degree = V(g2)$degree,
  gene_name = unname(gene_name_map[V(g2)$name] %||% ""),
  stringsAsFactors = FALSE
)
nodes_df$gene_name[is.na(nodes_df$gene_name)] <- ""
nodes_df$score <- 1000 * as.numeric(nodes_df$is_hub) + abs(nodes_df$kME) + 0.05 * nodes_df$degree
nodes_df <- nodes_df[order(nodes_df$score, decreasing = TRUE), , drop = FALSE]

cand <- nodes_df[nzchar(nodes_df$gene_name), , drop = FALSE]
if (nrow(cand) == 0) {
  label_genes <- character()
} else {
  label_genes <- head(cand$name, n = min(label_n, nrow(cand)))
}

# =============================================================================
# 7) 布局与作图（无轴线/无描边/默认无标题）
# =============================================================================

layout_xy <- compute_layout(g2, method = LAYOUT_METHOD, kk_maxiter = LAYOUT_KK_MAXITER)
colnames(layout_xy) <- c("x", "y")

if (isTRUE(ROTATE_LAYOUT)) {
  layout_xy <- rotate_layout_pca(layout_xy)
}

layout_xy <- layout_xy * as.numeric(LAYOUT_EXPAND)

png_out  <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".png"))
pdf_out  <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".pdf"))
node_out <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".nodes.tsv"))
edge_out <- file.path(OUTDIR, paste0("11f_network_", MODULE_COLOR, ".edges.tsv"))

theme_set(theme_classic(base_size = as.numeric(BASE_FONT_SIZE), base_family = as.character(FONT_FAMILY)))

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

nd$base_size <- abs(nd$kME)
if (all(nd$base_size == 0)) nd$base_size <- nd$degree
nd$plot_size <- nd$base_size * ifelse(nd$is_hub, as.numeric(HUB_SIZE_BOOST), 1.0)

nd$expr_color <- as.numeric(unname(expr_map[nd$name]))

nd$gene_name <- unname(gene_name_map[nd$name] %||% "")
nd$gene_name[is.na(nd$gene_name)] <- ""
nd$label_raw <- ifelse(nd$name %in% label_genes, nd$gene_name, "")
nd$label <- ifelse(nzchar(nd$label_raw), clean_label_text(nd$label_raw), "")

# ✅ 重复英文名唯一化（只改显示 label）
nd$label_display <- make_unique_labels(nd$label, nd$name)

pos <- nd[, c("name", "x", "y")]
ed <- merge(ed, pos, by.x = "from", by.y = "name", all.x = TRUE)
ed <- merge(ed, pos, by.x = "to", by.y = "name", all.x = TRUE, suffixes = c(".from", ".to"))

# ✅ 图例标题只显示 log2(TPM+1)（mean 与 TARGET_GROUP 放图注说明）
legend_title <- "log2(TPM+1)"

xrange <- range(nd$x, finite = TRUE)
yrange <- range(nd$y, finite = TRUE)
xpad <- diff(xrange) * as.numeric(PANEL_PADDING_FRAC)
ypad <- diff(yrange) * as.numeric(PANEL_PADDING_FRAC)
if (!is.finite(xpad) || xpad == 0) xpad <- 1
if (!is.finite(ypad) || ypad == 0) ypad <- 1

p <- ggplot() +
  geom_segment(
    data = ed,
    aes(x = x.from, y = y.from, xend = x.to, yend = y.to, linewidth = weight),
    alpha = 0.35,
    color = "grey40"
  ) +
  geom_point(
    data = nd,
    aes(x = x, y = y, size = plot_size, color = expr_color),
    shape = 16,
    alpha = 0.95
  ) +
  scale_linewidth(range = c(0.2, 1.4), guide = "none") +
  scale_size(range = c(2.5, 9), guide = "none") +
  get_color_scale(COLOR_SCHEME, legend_title) +
  labs(
    title = if (isTRUE(SHOW_TITLE)) "Co-expression network" else NULL,
    subtitle = NULL
  ) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

# ✅ 关键：只加一次坐标系（本次修复 Coordinate system 提示）
if (isTRUE(COORD_EQUAL)) {
  p <- p + coord_equal(
    xlim = c(xrange[1] - xpad, xrange[2] + xpad),
    ylim = c(yrange[1] - ypad, yrange[2] + ypad)
  )
} else {
  p <- p + coord_cartesian(
    xlim = c(xrange[1] - xpad, xrange[2] + xpad),
    ylim = c(yrange[1] - ypad, yrange[2] + ypad)
  )
}

if (isTRUE(USE_LABEL_REPEL) && has_repel) {
  p <- p + ggrepel::geom_text_repel(
    data = nd[nzchar(nd$label_display), , drop = FALSE],
    aes(x = x, y = y, label = label_display),
    size = 4,
    min.segment.length = Inf,
    segment.color = NA,
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = MAX_LABEL_OVERLAPS
  )
} else {
  p <- p + geom_text(
    data = nd[nzchar(nd$label_display), , drop = FALSE],
    aes(x = x, y = y, label = label_display),
    size = 3.6,
    vjust = -0.7
  )
}

# ✅ dev.off 不回显（本次消灭 null device 1）
open_png(png_out, PNG_WIDTH_IN, PNG_HEIGHT_IN, PNG_RES_DPI); print(p); invisible(dev.off())
open_pdf(pdf_out, PDF_WIDTH_IN, PDF_HEIGHT_IN); print(p); invisible(dev.off())

cat("[OK]", png_out, "\n")
cat("[OK]", pdf_out, "\n")

# =============================================================================
# 8) 导出 node/edge TSV
# =============================================================================

if (isTRUE(EXPORT_NODE_EDGE_TSV)) {
  nd_out <- data.frame(
    node = V(g2)$name,
    module_color = MODULE_COLOR,
    gene_name = unname(gene_name_map[V(g2)$name] %||% ""),
    label_display = nd$label_display[match(V(g2)$name, nd$name)],
    target_group = as.character(TARGET_GROUP),
    expr_mean_log2_tpm1 = unname(expr_map[V(g2)$name]),
    kME = V(g2)$kME,
    degree = V(g2)$degree,
    is_hub = V(g2)$is_hub,
    is_labeled = V(g2)$name %in% label_genes,
    stringsAsFactors = FALSE
  )
  nd_out$gene_name[is.na(nd_out$gene_name)] <- ""
  nd_out$label_display[is.na(nd_out$label_display)] <- ""
  nd_out$expr_mean_log2_tpm1 <- as.numeric(nd_out$expr_mean_log2_tpm1)

  ed_out <- igraph::as_data_frame(g2, what = "edges")
  ed_out$weight <- E(g2)$weight
  colnames(ed_out)[1:2] <- c("from", "to")

  data.table::fwrite(nd_out, node_out, sep = "\t", quote = FALSE)
  data.table::fwrite(ed_out, edge_out, sep = "\t", quote = FALSE)

  cat("[OK]", node_out, "\n")
  cat("[OK]", edge_out, "\n")
}

cat("[DONE] 11f network finished.\n")

