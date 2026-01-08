#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(igraph)
  library(ggraph)
})

# =========================
# 皇上在这里手填参数
# =========================

INPUT_TSV <- "input/KEGG_by_term_test.tsv"

# 固定唯一输出目录（短名、无时间戳；每次运行先删再建）
OUT_PARENT_DIR <- "output"
OUT_DIR_NAME   <- "kegg_net"

# 选前多少个 term（按 p_adjust 从小到大）
TOP_N_TERMS <- 12

# gene 节点过滤：只保留连接 term 数 >= 1（建议别动）
MIN_GENE_DEGREE <- 1

# 只标共享 gene（degree>=2）
LABEL_GENE_DEGREE_AT_LEAST <- 2

# 不再额外强行标 top genes（避免拥挤）
LABEL_TOP_GENES <- 0

# term 是否都标
LABEL_TERMS <- TRUE

# 标题
TITLE_ON   <- TRUE
TITLE_TEXT <- "KEGG gene–concept network (unique expanded OGs)"

# 字体（默认 Arial；若系统不可用会自动 fallback）
FONT_FAMILY <- "Arial"

# 输出
DPI <- 300
FIG_WIDTH_IN  <- 12
FIG_HEIGHT_IN <- 7

# term 节点：大小与显著性映射
TERM_NODE_SIZE_MIN <- 5
TERM_NODE_SIZE_MAX <- 12
TERM_FILL_LOW  <- "#95C8F2"
TERM_FILL_HIGH <- "#9E90E6"

# gene 节点样式
GENE_NODE_SIZE  <- 1.6
GENE_NODE_COLOR <- "#9A9A9A"

# 边样式
EDGE_ALPHA <- 0.45
EDGE_WIDTH <- 0.35

# FR 布局迭代次数（更大更紧）
FR_NITER <- 2500

# 标签字号
TERM_LABEL_SIZE <- 3.0
GENE_LABEL_SIZE <- 3.0

# 输入表列名（按你的 pipeline）
COL_PATHWAY_ID <- "pathway_id"
COL_TERM_NAME  <- "term_name"
COL_GENE_IDS   <- "gene_ids"
COL_GENE_COUNT <- "gene_count"
COL_P_ADJUST   <- "p_adjust"

GENE_ID_SEP <- ";"


# =========================
# 工具函数
# =========================

mkdir_p <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

reset_dir <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop(paste0("无法创建输出目录：", path), call. = FALSE)
}

stop_if_missing_cols <- function(df, cols) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(
      paste0(
        "输入表缺少必要列：", paste(miss, collapse = ", "), "\n",
        "当前列名：", paste(colnames(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

parse_gene_ids <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || toupper(x) == "NA") return(character(0))
  parts <- unlist(strsplit(x, GENE_ID_SEP, fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  parts[!duplicated(parts)]
}

resolve_font <- function(requested) {
  cand <- unique(c(requested, "Arial", "Liberation Sans", "DejaVu Sans", "sans"))
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    fam <- unique(systemfonts::system_fonts()$family)
    for (f in cand) if (f %in% fam) return(f)
  }
  "sans"
}

enable_showtext <- function(dpi) {
  if (requireNamespace("showtext", quietly = TRUE)) {
    suppressPackageStartupMessages(library(showtext))
    showtext_opts(dpi = dpi)
    showtext_auto(enable = TRUE)
    return(TRUE)
  }
  return(FALSE)
}

pick_png_device <- function() {
  if (requireNamespace("ragg", quietly = TRUE)) return(ragg::agg_png)
  if (capabilities("cairo")) return("cairo-png")
  return("png")
}

pick_pdf_device <- function() {
  if (capabilities("cairo")) return(cairo_pdf)
  return("pdf")
}

save_plot_safe <- function(p, png_file, pdf_file, w, h, dpi) {
  ggsave(png_file, plot = p, width = w, height = h, dpi = dpi, bg = "white", device = pick_png_device())
  ggsave(pdf_file, plot = p, width = w, height = h, bg = "white", device = pick_pdf_device())
}

# 20 色马卡龙（你的皇家配色）
ROYAL_MACARON <- c(
  "#F79D93","#95C8F2","#F6CD96","#9FD5CB","#F6C6E7",
  "#A99BEF","#F5B07E","#A7D9F7","#F3A6C9","#B6E0B6",
  "#FBE3A8","#8FE3F2","#F4B1A8","#8FB1F2","#F0D07A",
  "#7FD3C8","#F7B6D2","#C3A3F5","#D6E7B5","#9E90E6"
)

# =========================
# 主流程
# =========================

FONT_USED <- resolve_font(FONT_FAMILY)

OUT_DIR <- file.path(OUT_PARENT_DIR, OUT_DIR_NAME)
mkdir_p(OUT_PARENT_DIR)
reset_dir(OUT_DIR)

message("========== KEGG gene–concept network (star-cluster) ==========")
message("[IN ] ", INPUT_TSV)
message("[OUT] ", OUT_DIR)

if (!file.exists(INPUT_TSV)) {
  stop(paste0("找不到输入文件：", INPUT_TSV, "（请在 magic/ 下运行）"), call. = FALSE)
}

df <- read_tsv(INPUT_TSV, show_col_types = FALSE, progress = FALSE)
stop_if_missing_cols(df, c(COL_PATHWAY_ID, COL_TERM_NAME, COL_GENE_IDS, COL_GENE_COUNT, COL_P_ADJUST))

df2 <- df %>%
  mutate(
    p_adjust_num = suppressWarnings(as.numeric(.data[[COL_P_ADJUST]])),
    gene_count_num = suppressWarnings(as.numeric(.data[[COL_GENE_COUNT]]))
  ) %>%
  arrange(p_adjust_num) %>%
  slice_head(n = TOP_N_TERMS) %>%
  mutate(
    term_id    = paste0("TERM::", .data[[COL_PATHWAY_ID]]),
    term_label = as.character(.data[[COL_TERM_NAME]]),  # ✅ 只用 term_name，不带 KO
    score      = -log10(pmax(p_adjust_num, 1e-300))
  ) %>%
  select(term_id, term_label, score, p_adjust_num, gene_count_num, all_of(COL_GENE_IDS))

term_tbl <- df2 %>%
  transmute(
    term_id = term_id,
    term_label = term_label,        # ✅ legend/label 都用 term_label
    category = term_label,          # ✅ 边颜色/legend 分类也用 term_label（不带 KO）
    score = score,
    p_adjust = p_adjust_num,
    gene_count = ifelse(is.na(gene_count_num), 0, gene_count_num),
    gene_ids_raw = .data[[COL_GENE_IDS]]
  )

edges <- term_tbl %>%
  mutate(gene_vec = lapply(gene_ids_raw, parse_gene_ids)) %>%
  select(term_id, category, gene_vec) %>%
  tidyr::unnest(gene_vec) %>%
  rename(gene_id = gene_vec) %>%
  mutate(gene_node = paste0("GENE::", gene_id)) %>%
  transmute(from = term_id, to = gene_node)

if (nrow(edges) == 0) {
  stop("gene_ids 解析为空：请检查 gene_ids 列是否为 ';' 分隔，且不是 NA。", call. = FALSE)
}

# 建图（无向二部图）
g0 <- graph_from_data_frame(edges, directed = FALSE)

V(g0)$kind <- ifelse(startsWith(V(g0)$name, "TERM::"), "term", "gene")
V(g0)$label <- V(g0)$name
V(g0)$category <- NA_character_
V(g0)$score <- NA_real_
V(g0)$gene_count <- NA_real_

# term 属性回填
term_attr <- term_tbl %>% select(term_id, term_label, category, score, gene_count) %>% distinct()
term_idx <- which(V(g0)$kind == "term")
term_names <- V(g0)$name[term_idx]
m <- match(term_names, term_attr$term_id)

V(g0)$label[term_idx] <- term_attr$term_label[m]     # ✅ term 节点标签 = term_name
V(g0)$category[term_idx] <- term_attr$category[m]
V(g0)$score[term_idx] <- term_attr$score[m]
V(g0)$gene_count[term_idx] <- term_attr$gene_count[m]

# gene 标签：显示 gene id
gene_idx <- which(V(g0)$kind == "gene")
V(g0)$label[gene_idx] <- sub("^GENE::", "", V(g0)$name[gene_idx])

# gene degree 过滤
deg0 <- degree(g0)
V(g0)$degree <- as.integer(deg0[V(g0)$name])

keep_genes <- V(g0)$name[V(g0)$kind == "gene" & V(g0)$degree >= MIN_GENE_DEGREE]
keep_terms <- V(g0)$name[V(g0)$kind == "term"]
g <- induced_subgraph(g0, vids = intersect(V(g0)$name, c(keep_terms, keep_genes)))

if (vcount(g) == 0 || ecount(g) == 0) {
  stop("过滤后网络为空：请把 MIN_GENE_DEGREE 调小（比如 1）或把 TOP_N_TERMS 调大。", call. = FALSE)
}

# 把 category 写进边属性（每条边必有一个 term 端点）
ee <- ends(g, E(g), names = TRUE)
edge_cat <- character(nrow(ee))
for (i in seq_len(nrow(ee))) {
  a <- ee[i, 1]
  b <- ee[i, 2]
  if (startsWith(a, "TERM::")) {
    edge_cat[i] <- V(g)$category[match(a, V(g)$name)]
  } else {
    edge_cat[i] <- V(g)$category[match(b, V(g)$name)]
  }
}
E(g)$category <- edge_cat

# category 调色板（按 term_label 分配颜色）
cats <- sort(unique(na.omit(term_attr$category)))
pal <- rep(ROYAL_MACARON, length.out = length(cats))
cat_color <- setNames(pal, cats)

# 标签策略：term 全标；gene 只标共享（degree>=2）
deg_g <- degree(g)
gene_nodes <- names(deg_g)[startsWith(names(deg_g), "GENE::")]
gene_deg <- deg_g[gene_nodes]

shared_gene_nodes <- names(gene_deg)[gene_deg >= LABEL_GENE_DEGREE_AT_LEAST]
top_gene_nodes <- character(0)
if (is.numeric(LABEL_TOP_GENES) && LABEL_TOP_GENES > 0) {
  top_gene_nodes <- names(sort(gene_deg, decreasing = TRUE))[seq_len(min(LABEL_TOP_GENES, length(gene_nodes)))]
}

V(g)$show_label <- FALSE
if (LABEL_TERMS) V(g)$show_label[V(g)$kind == "term"] <- TRUE
V(g)$show_label[V(g)$kind == "gene" & V(g)$name %in% union(shared_gene_nodes, top_gene_nodes)] <- TRUE

# 字体防炸（不再打印 [1] TRUE）
invisible(enable_showtext(DPI))

# FR 力导向布局：星团风格关键
set.seed(1)
lay <- create_layout(g, layout = "fr", niter = FR_NITER)

p <- ggraph(lay) +
  geom_edge_link0(
    aes(color = category),
    alpha = EDGE_ALPHA,
    linewidth = EDGE_WIDTH,
    show.legend = TRUE
  ) +
  geom_node_point(
    data = function(d) d %>% filter(kind == "gene"),
    color = GENE_NODE_COLOR,
    size = GENE_NODE_SIZE,
    alpha = 0.95
  ) +
  geom_node_point(
    data = function(d) d %>% filter(kind == "term"),
    aes(size = gene_count, fill = score),
    shape = 21,
    color = "#333333",
    stroke = 0.35,
    alpha = 0.98
  ) +
  scale_fill_gradient(low = TERM_FILL_LOW, high = TERM_FILL_HIGH, name = "-log10(p_adjust)") +
  scale_size_continuous(range = c(TERM_NODE_SIZE_MIN, TERM_NODE_SIZE_MAX), name = "term gene_count") +
  scale_color_manual(values = cat_color, name = "category") +
  theme_classic(base_family = FONT_USED) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

use_repel <- requireNamespace("ggrepel", quietly = TRUE)

# term 标签（不带 KO）
p <- p +
  geom_node_text(
    data = function(d) d %>% filter(show_label, kind == "term"),
    aes(label = label),
    family = FONT_USED,
    size = TERM_LABEL_SIZE,
    repel = use_repel,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.15, "lines"),
    segment.alpha = 0.25,
    max.overlaps = Inf
  )

# gene 标签：只标共享 gene（degree>=2）
p <- p +
  geom_node_text(
    data = function(d) d %>% filter(show_label, kind == "gene"),
    aes(label = label),
    family = FONT_USED,
    size = GENE_LABEL_SIZE,
    repel = use_repel,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.15, "lines"),
    segment.alpha = 0.25,
    max.overlaps = Inf
  )

if (TITLE_ON && !is.na(TITLE_TEXT) && nchar(TITLE_TEXT) > 0) {
  p <- p + ggtitle(TITLE_TEXT)
}

png_file <- file.path(OUT_DIR, "kegg_cnet_star.png")
pdf_file <- file.path(OUT_DIR, "kegg_cnet_star.pdf")

# 导出表（便于复现与审稿人追溯）
edge_out <- as_tibble(as_data_frame(g, what = "edges")) %>%
  mutate(category = E(g)$category)
write_tsv(edge_out, file.path(OUT_DIR, "edges.tsv"))

node_out <- tibble(
  node = V(g)$name,
  kind = V(g)$kind,
  label = V(g)$label,
  category = V(g)$category,
  score = V(g)$score,
  gene_count = V(g)$gene_count,
  degree = as.integer(degree(g)),
  show_label = V(g)$show_label
)
write_tsv(node_out, file.path(OUT_DIR, "nodes.tsv"))

sink(file.path(OUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

save_plot_safe(p, png_file, pdf_file, FIG_WIDTH_IN, FIG_HEIGHT_IN, DPI)

message("[DONE] star-cluster cnet saved.")
message("[FONT] ", FONT_USED)
message("[OUT] ", png_file)
message("[OUT] ", pdf_file)

