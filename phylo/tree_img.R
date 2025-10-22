#!/usr/bin/env Rscript
# =============================================================================
# root_and_plot_windows.R —— 顶刊版（不拥挤·对齐美化·可加门类括注）
# =============================================================================

## ========== 1) 参数（按需改） ==============================================
lib_dir            <- file.path(Sys.getenv("USERPROFILE"), "Rlibs")
tree_path          <- "C:/Users/herol/Downloads/supermatrix.contree"
outgroups          <- c("Nematostella_vectensis")       # 可多个外群
show_branch_length <- TRUE
bootstrap_cutoff   <- 70                                # 仅显示 ≥ 阈值的自举值
font_family        <- "Arial"                           # 期刊常见无衬线
line_size_tree     <- 0.7                               # 树枝线条粗细（刊物安全值）
dpi                <- 600
png_width_px_fixed <- 3200                              # 宽度固定；高度自适应
also_pdf           <- FALSE

# 重点物种（可留空 c()；若需要高亮，bold_targets <- TRUE）
targets_bold <- c("Sinonovacula_constricta","Sinonovacula_rivularis")
bold_targets <- TRUE

# 右侧“门类括注”（按您的树 tip.label 填写；不需要可设为空 list()）
clade_strips <- list(
  list(taxa1="Nematostella_vectensis",  taxa2="Nematostella_vectensis",  label="Cnidaria"),
  list(taxa1="Homo_sapiens",            taxa2="Branchiostoma_floridae",  label="Chordata"),
  list(taxa1="Nautilus_pompilius",      taxa2="Cyclina_sinensis",        label="Mollusca"),
  list(taxa1="Capitella_teleta",        taxa2="Capitella_teleta",        label="Annelida"),
  list(taxa1="Drosophila_melanogaster", taxa2="Drosophila_melanogaster", label="Arthropoda")
)

## ========== 2) 库路径与依赖 ================================================
if (nzchar(lib_dir)) { dir.create(lib_dir, recursive=TRUE, showWarnings=FALSE); .libPaths(c(lib_dir, .libPaths())) }
ensure_pkg <- function(pkgs, bioc=FALSE){
  for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)){
    message(sprintf("正在安装 %s ...", p))
    if(bioc){
      if(!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager", repos="https://cloud.r-project.org", lib=.libPaths()[1])
      BiocManager::install(p, ask=FALSE, update=FALSE, lib=.libPaths()[1])
    } else install.packages(p, repos="https://cloud.r-project.org", lib=.libPaths()[1])
  }
}
ensure_pkg(c("ape","ragg","ggplot2"))
ensure_pkg(c("treeio","tidytree","ggtree","ggfun","yulab.utils","ggplotify"), bioc=TRUE)

suppressPackageStartupMessages({ library(ape); library(ggplot2); library(ggtree); library(ragg) })

## ========== 3) 读树 → 多外群加根 → 梯形化 → Bootstrap 清洗 ================
stopifnot(file.exists(tree_path))
tr <- read.tree(tree_path)

missing <- setdiff(outgroups, tr$tip.label)
if(length(missing)>0){
  stop(sprintf("外群未在树中：%s\n示例 tip：%s",
               paste(missing, collapse=", "), paste(head(tr$tip.label, 20), collapse=", ")))
}
rooted <- root(tr, outgroup=outgroups, resolve.root=TRUE)
rooted <- ladderize(rooted, right=TRUE)

# 自举值阈值与清洗
if(!is.null(rooted$node.label)){
  bl <- suppressWarnings(as.numeric(rooted$node.label))
  bl[is.na(bl)] <- NA_real_
  if(bootstrap_cutoff>0) bl[bl<bootstrap_cutoff] <- NA_real_
  rooted$node.label <- ifelse(is.na(bl), "", as.character(bl))
}

## ========== 4) 自适应版式参数（不拥挤的关键） =============================
n_tips <- length(rooted$tip.label)

# 根据物种数自适应：高度、字号、偏移与留白
clamp <- function(x, lo, hi) max(lo, min(hi, x))
png_height_px <- max(2400, round(n_tips * 90))         # 每个 tip ~90px，低于 2400 也不至于挤
tip_cex       <- clamp(3.2 - 0.03*(n_tips-12), 2.2, 3.2) # 物种多→字号略降（安全范围：2.2–3.2）

# 先粗算 xmax 以便确定 offset 与右侧留白
p_probe <- ggtree(rooted, layout="rectangular",
                  branch.length = if (show_branch_length) "branch.length" else "none",
                  size=line_size_tree)
xmax <- max(ggplot_build(p_probe)$data[[1]]$x, na.rm=TRUE)

label_offset <- xmax * 0.05                              # 与枝端的距离（避免撞线）
right_pad    <- xmax * 0.50                              # 右侧留白，保证对齐+门类文字

## ========== 5) 输出路径 ====================================================
out_dir <- dirname(tree_path)
fname_base <- "rooted_tree_rectangular_topjournal"
out_png <- file.path(out_dir, paste0(fname_base, ".png"))
out_pdf <- file.path(out_dir, paste0(fname_base, ".pdf"))

## ========== 6) 绘图（对齐不拥挤） ==========================================
label_map <- function(x) gsub("_"," ", x)

p <- ggtree(rooted, layout="rectangular",
            branch.length = if (show_branch_length) "branch.length" else "none",
            size=line_size_tree) +
  # 学名：右侧对齐 + 细斜体 + 统一字体
  geom_tiplab(aes(label = label_map(label)),
              size = tip_cex, fontface = "italic",
              offset = label_offset, align = TRUE,
              linesize = 0.3, linetype = "solid") +
  # Bootstrap
  geom_nodelab(aes(label = label),
               size = clamp(tip_cex*0.9, 2.0, 3.0),
               hjust = -0.10, vjust = -0.20, na.rm = TRUE) +
  theme_tree2() +
  theme(plot.margin = margin(16, 36, 14, 16),
        text = element_text(family = font_family),
        axis.title.x = element_text(margin = margin(t = 8))) +
  xlim(NA, xmax + right_pad) +
  xlab(if (show_branch_length) "Substitutions / site" else NULL)

# 重点物种覆盖加粗（可关闭）
if (bold_targets && length(intersect(targets_bold, rooted$tip.label)) > 0) {
  p <- p + geom_tiplab(data = function(d) subset(d, label %in% targets_bold),
                       aes(label = label_map(label)),
                       size = tip_cex, fontface = "bold.italic",
                       offset = label_offset, align = TRUE, linesize = 0.5)
}

# 比例尺：放左下，避免与标签冲突；同时抑制“很小的小数”露头
if (show_branch_length) {
  p <- p + ggtree::geom_treescale(x = xmax*0.10, y = -0.8, width = xmax*0.10, fontsize = 3)
}

## ========== 7) 门类括注（可选，刊物风格） ================================
add_strip <- function(p, tip1, tip2, lab){
  p + ggtree::geom_strip(
    taxa1 = tip1, taxa2 = tip2, label = lab,
    align = TRUE, color = "black", barsize = 0.6,
    fontsize = clamp(tip_cex*1.1, 2.6, 3.4), fontface = "plain",
    offset = label_offset * 1.6,   # 竖线位置
    offset.text = label_offset * 2.3) # 文字更靠右，避免压学名
}
if (length(clade_strips)) {
  for (st in clade_strips) {
    if (st$taxa1 %in% rooted$tip.label && st$taxa2 %in% rooted$tip.label) {
      p <- add_strip(p, st$taxa1, st$taxa2, st$label)
    } else {
      message(sprintf("⚠️ 跳过括注 %s（%s 或 %s 不在树中）", st$label, st$taxa1, st$taxa2))
    }
  }
}

## ========== 8) 导出（ragg 抗锯齿） ========================================
ragg::agg_png(out_png, width = png_width_px_fixed, height = png_height_px, res = dpi, background = "white")
print(p); dev.off()
if (also_pdf) {
  ggsave(out_pdf, plot = p, width = png_width_px_fixed/dpi, height = png_height_px/dpi,
         units = "in", device = cairo_pdf, dpi = dpi, bg = "white")
}
message("✅ 已保存：", out_png)
if (also_pdf) message("✅ 另存 PDF：", out_pdf)
