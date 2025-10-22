# root_and_plot_windows.R —— Windows RStudio 一键出图（支持多外群）
# 功能：支持多个外群 → 加根 → 美化树图 → 输出 PNG

## ========== 1) 参数 ==========
lib_dir <- file.path(Sys.getenv("USERPROFILE"), "Rlibs")  # 包目录，可改 "" 用默认库

tree_path <- "//wsl.localhost/Ubuntu-22.04/home/mxrmiss/project/phylo/supermatrix.contree"

# 外群可多个
outgroups <- c("A_californica_protein")   # 加州海兔 + 中华石鳖

layout_choice      <- "rectangular"  # "rectangular" 或 "circular"
show_branch_length <- TRUE

## ========== 2) 库路径 ==========
if (nzchar(lib_dir)) {
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib_dir, .libPaths()))
}

## ========== 3) 自动安装依赖 ==========
ensure_pkg <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("正在安装 %s ...", p))
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = .libPaths()[1])
        }
        BiocManager::install(p, ask = FALSE, update = FALSE, lib = .libPaths()[1])
      } else {
        install.packages(p, repos = "https://cloud.r-project.org", lib = .libPaths()[1])
      }
    }
  }
}

ensure_pkg(c("ape", "ggplot2", "ragg"))
ensure_pkg(c("ggplotify", "aplot", "ggfun", "yulab.utils"), bioc = FALSE)
ensure_pkg(c("treeio", "tidytree", "ggtree"), bioc = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(ragg)
})
have_ggtree <- requireNamespace("ggtree", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)
if (have_ggtree) {
  suppressPackageStartupMessages({
    library(ggtree)
    library(ggplot2)
  })
}

## ========== 4) 读树并加根（支持多个外群） ==========
if (!file.exists(tree_path)) stop("找不到树文件：", tree_path)

tr <- read.tree(tree_path)

missing <- setdiff(outgroups, tr$tip.label)
if (length(missing) > 0) {
  stop(sprintf("这些外群不在树里：%s\n可选 tip：%s",
               paste(missing, collapse = ", "),
               paste(head(tr$tip.label, 20), collapse = ", ")))
}

rooted <- root(tr, outgroup = outgroups, resolve.root = TRUE)

## ========== 5) 输出路径 ==========
out_dir <- dirname(tree_path)
out_png <- file.path(out_dir, if (layout_choice == "circular") "rooted_tree_circular.png" else "rooted_tree.png")

## ========== 6) 绘图 ==========
if (have_ggtree) {
  p <- ggtree(rooted,
              layout = layout_choice,
              branch.length = if (show_branch_length) "branch.length" else "none") +
    geom_tiplab(size = 3, fontface = "italic") +
    geom_nodelab(aes(label = label), size = 2.6, hjust = -0.2) +
    theme(plot.margin = margin(12,12,12,12),
          text = element_text(family = "sans"))
  
  if (show_branch_length) {
    p <- p + theme_tree2() + xlab("Substitutions / site")
  }
  
  if (layout_choice == "rectangular") {
    pdata <- p$data
    p <- p + xlim(NA, max(pdata$x, na.rm = TRUE) * 1.15)
  }
  
  ragg::agg_png(out_png, width = 2400, height = 1800, res = 300, background = "white")
  print(p)
  dev.off()
  
} else {
  ragg::agg_png(out_png, width = 2400, height = 1800, res = 300, background = "white")
  plot(rooted, show.tip.label = TRUE, cex = 0.7)
  nodelabels(rooted$node.label, frame = "n", cex = 0.6, adj = c(-0.1, -0.5))
  dev.off()
  message("⚠️ 已用 ape 生成树图（ggtree 不可用）")
}

message("✅ 已保存：", out_png)
