#!/usr/bin/env Rscript

# ==============================================================================
# MAGIC PLOTTING KIT (v4.0) - 脚本 10: evolution_tree.R
# ==============================================================================
# 功能：生成 "一树三层" 全景基因组进化图 (Tree + Orthology Bar + CAFE Heatmap)
# 级别：Nature Ecology & Evolution / MBE 标准
# 依赖：ggtreeExtra, ggnewscale
# ==============================================================================

# ------------------------------------------------------------------------------
# [用户配置区域] USER CONFIGURATION
# ------------------------------------------------------------------------------

# 1. 输入文件
TREE_FILE   <- "FigTree.tre"                 # 树文件
CAFE_FILE   <- "Base_clade_results.txt"      # CAFE 结果
COUNTS_FILE <- "Orthogroups.GeneCount.tsv"   # OrthoFinder 统计

# 2. 输出设置
OUTPUT_DIR_NAME <- "evolution_landscape"     # 输出目录名

# 3. 绘图微调
PLOT_TITLE      <- ""                        # 标题
USE_SHORT_NAME  <- TRUE                      # 是否缩写物种名 (G. species)
OFFSET_GAP      <- 0.05                      # 层级之间的间隙 (相对宽度)
BAR_WIDTH       <- 0.3                       # 柱状图层的宽度
HEATMAP_WIDTH   <- 0.15                      # 热图层的宽度

# 4. [顶刊配色库]
# Orthology Stacked Bar 配色 (Nature 风格)
COLOR_ORTHO <- c(
  "Single-copy"      = "#4DBBD5CC", # 蓝 (保守)
  "Multi-copy"       = "#E64B35CC", # 红 (扩增)
  "Species-specific" = "#00A087CC", # 绿 (特有)
  "Other"            = "#B09C85CC"  # 灰
)

# 地质背景色 (莫兰迪色系)
COLOR_GEO <- list(
  "Neogene"    = "#FFE5B4", "Paleogene"  = "#FFDAB9",
  "Cretaceous" = "#C1FFC1", "Jurassic"   = "#B0E0E6", "Triassic"   = "#D8BFD8"
)

# ==============================================================================
# [核心逻辑]
# ==============================================================================

# --- 0. 依赖检查与加载 ---
required_pkgs <- c("ape", "ggplot2", "tidyverse", "ggtree", "treeio", "ggtreeExtra", "ggnewscale", "deeptime")
for (pkg in required_pkgs) {
  if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE))) {
    stop(paste0("❌ 缺少依赖包: ", pkg, "。请先安装它：BiocManager::install('", pkg, "')"))
  }
}

# 路径准备
input_path  <- file.path("input")
output_root <- file.path("output")
if (!dir.exists("scripts")) stop("❌ 请在 magic/ 根目录下运行！")
final_output_dir <- file.path(output_root, OUTPUT_DIR_NAME)
if (!dir.exists(final_output_dir)) dir.create(final_output_dir, recursive = TRUE)

check_file <- function(f) {
  p <- file.path(input_path, f)
  if (!file.exists(p)) stop(paste("❌ 找不到文件:", p))
  return(p)
}

message("✅ 环境与依赖检查通过，开始构建进化全景图...")

# --- 1. 数据读取与清洗 ---

# [1.1] 树文件读取与智能校正
tree <- read.mcmctree(check_file(TREE_FILE))
# 智能判断：如果树高极小(<500)，假设单位为100Ma，放大100倍
max_h <- max(node.depth.edgelength(tree@phylo))
if (max_h < 500) {
  message(paste0("ℹ️  检测到树高 (", round(max_h,2), ") 较小，执行 x100 时间校正..."))
  tree@phylo$edge.length <- tree@phylo$edge.length * 100
}

# [1.2] Orthology 数据清洗 (堆叠柱状图数据)
raw_counts <- read_tsv(check_file(COUNTS_FILE), show_col_types = FALSE)
# 计算分类
classify_orthogroups <- function(df) {
  sp_cols <- setdiff(names(df), c("Orthogroup", "Total"))
  res_list <- list()
  for (sp in sp_cols) {
    counts <- df[[sp]]
    others <- df %>% select(all_of(setdiff(sp_cols, sp)))
    others_exist <- rowSums(others) > 0
    n_single <- sum(counts == 1 & others_exist)
    n_multi  <- sum(counts > 1 & others_exist)
    n_spec   <- sum(counts > 0 & !others_exist)
    n_other  <- sum(counts) - (n_single + n_multi + n_spec)
    
    # 计算百分比 (Normalized)
    total <- n_single + n_multi + n_spec + n_other
    res_list[[sp]] <- data.frame(
      Species = sp,
      Category = c("Single-copy", "Multi-copy", "Species-specific", "Other"),
      Count = c(n_single, n_multi, n_spec, n_other),
      Percent = c(n_single, n_multi, n_spec, n_other) / total * 100
    )
  }
  bind_rows(res_list)
}
bar_data <- classify_orthogroups(raw_counts)
bar_data$Category <- factor(bar_data$Category, levels = c("Single-copy", "Multi-copy", "Species-specific", "Other"))

# [1.3] CAFE 动态数据清洗 (热图数据)
cafe_raw <- readLines(check_file(CAFE_FILE))
cafe_df <- read.table(text = grep("^#", cafe_raw, value = TRUE, invert = TRUE), 
                      col.names = c("Taxon_ID", "Increase", "Decrease"))
# 提取物种名
cafe_df$Species <- gsub("<[0-9]+>", "", cafe_df$Taxon_ID)
cafe_clean <- cafe_df %>% 
  filter(Species != "") %>%
  select(Species, Increase, Decrease) %>%
  pivot_longer(cols = c("Increase", "Decrease"), names_to = "Type", values_to = "Value")

# [1.4] 名称标准化
format_name <- function(x) {
  if (USE_SHORT_NAME) gsub("^([A-Z])[a-z]+_([a-z0-9]+)$", "\\1. \\2", x) else gsub("_", " ", x)
}
tree@phylo$tip.label <- format_name(tree@phylo$tip.label)
bar_data$Species <- format_name(bar_data$Species)
cafe_clean$Species <- format_name(cafe_clean$Species)

# --- 2. 核心绘图 (One Tree, Three Layers) ---

message("🎨 正在构建图层 (ggtreeExtra)...")

# 计算地质时间
max_age <- max(node.depth.edgelength(tree@phylo))
geo_rects <- data.frame(
  period = c("Neogene", "Paleogene", "Cretaceous", "Jurassic", "Triassic"),
  start  = c(0,         23.03,       66.0,         145.0,      201.3),
  end    = c(23.03,     66.0,        145.0,        201.3,      251.9),
  color  = c(COLOR_GEO$Neogene, COLOR_GEO$Paleogene, COLOR_GEO$Cretaceous, COLOR_GEO$Jurassic, COLOR_GEO$Triassic)
) %>% filter(start < max_age)
geo_rects$end[geo_rects$end > max_age] <- max_age + 50

# --- Layer 1: 主树 (Main Tree) ---
p <- ggtree(tree, linewidth = 0.6) + 
  geom_rect(data = geo_rects, 
            aes(xmin = -end, xmax = -start, ymin = -Inf, ymax = Inf, fill = period), 
            inherit.aes = FALSE, alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = setNames(geo_rects$color, geo_rects$period)) +
  theme_tree2() + 
  scale_x_continuous(labels = abs, expand = expansion(mult = c(0, 0.05))) +
  # 物种名 (Tip Labels)
  geom_tiplab(align = TRUE, linesize = 0.3, size = 3.5, fontface = "italic", offset = max_age * 0.02) +
  xlab("Geological Time (Ma)") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(t=10))

# --- Layer 2: 基因组结构 (Stacked Bar) ---
message(" -> 添加 Layer 2: Orthology Composition...")
# 重点：p + geom_fruit
p <- p + 
  geom_fruit(
    data = bar_data,
    geom = geom_bar,
    mapping = aes(y = Species, x = Percent, fill = Category),
    orientation = "y",
    stat = "identity",
    offset = 0.25,            # 距离树末端的距离
    pwidth = BAR_WIDTH,       # 柱状图的总宽度 (相对值)
    axis.params = list(axis = "x", text.size = 2, line.size = 0, vjust = 1),
    grid.params = list()
  ) + 
  scale_fill_manual(
    name = "Orthology (%)",
    values = COLOR_ORTHO,
    guide = guide_legend(order = 1) # 图例顺序
  )

# --- Layer 3: 进化动态 (Expansion/Contraction Heatmap) ---
message(" -> 添加 Layer 3: Evolutionary Dynamics...")
# 使用 new_scale_fill 为热图开启新的颜色通道
p <- p + new_scale_fill() 

# 由于 Expansion 和 Contraction 需要不同的色阶，我们这里做一个 Trick
# 使用 log1p 处理数值，使颜色更有层次
cafe_clean$LogValue <- log1p(cafe_clean$Value)

# 绘制 Expansion (红色系)
p <- p + 
  geom_fruit(
    data = cafe_clean %>% filter(Type == "Increase"),
    geom = geom_tile,
    mapping = aes(y = Species, x = Type, fill = LogValue),
    offset = OFFSET_GAP,
    pwidth = HEATMAP_WIDTH / 2,
    axis.params = list(axis = "x", text.size = 2.5, text.angle = 90, vjust = 0.5),
  ) +
  scale_fill_gradient(
    name = "Expansion (log)",
    low = "#FFF5F0", high = "#CB181D",
    guide = guide_colorbar(order = 2, barwidth = 0.5, barheight = 3)
  )

# 绘制 Contraction (蓝色系)
p <- p + new_scale_fill() +
  geom_fruit(
    data = cafe_clean %>% filter(Type == "Decrease"),
    geom = geom_tile,
    mapping = aes(y = Species, x = Type, fill = LogValue),
    offset = 0.01, # 紧挨着 Expansion
    pwidth = HEATMAP_WIDTH / 2,
    axis.params = list(axis = "x", text.size = 2.5, text.angle = 90, vjust = 0.5),
  ) +
  scale_fill_gradient(
    name = "Contraction (log)",
    low = "#F7FBFF", high = "#084594",
    guide = guide_colorbar(order = 3, barwidth = 0.5, barheight = 3)
  )

# --- 3. 输出保存 ---
file_pdf <- file.path(final_output_dir, "Figure1_Landscape.pdf")
file_png <- file.path(final_output_dir, "Figure1_Landscape.png")

message(paste("💾 保存全景图:", file_pdf))
ggsave(file_pdf, plot = p, width = 14, height = 8, device = cairo_pdf)
ggsave(file_png, plot = p, width = 14, height = 8, dpi = 300, bg = "white")

message("✨ Figure 1 生成完毕！\n   结构：左[时间树] -> 中[同源结构 100%] -> 右[扩张/收缩热图]。")
