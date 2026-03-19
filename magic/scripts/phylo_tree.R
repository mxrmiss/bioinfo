#!/usr/bin/env Rscript

# ==============================================================================
# 脚本名称: phylo_tree.R
# 功能描述:
#   1. 读取 FigTree.tre 绘制系统发育树，并叠加 CAFE 扩张/收缩气泡。
#   2. 支持一个开关：是否将输入树中的分歧时间（branch lengths）统一乘以100。
#   3. 保留原有“树宽物理压缩”逻辑，用于排版，不改变生物学解释。
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备：加载必要的 R 包
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(ggplot2)  # 绘图核心包
    library(ggtree)   # 进化树绘制包
    library(treeio)   # 树文件读取包
    library(dplyr)    # 数据清洗工具
    library(fs)       # 文件系统工具
    library(tibble)   # 数据框工具
    library(ape)      # 基础系统发育工具
    library(stringr)  # 字符串处理工具
})

# ==============================================================================
# --- [用户参数锁定区] ---
# ==============================================================================

# [开关] 是否将输入 FigTree.tre 中的分歧时间统一乘以100
# FALSE = 不乘100，直接使用原始时间
# TRUE  = 先乘100，再继续后续计算与绘图
multiply_divtime_by_100 <- TRUE

# [缩放] 树的宽度缩放比例
# 作用：物理压缩树的长度，仅用于排版。
# 0.4 表示压缩到原长的 40%，给右侧气泡腾出空间。
# 注意：这不是时间单位变换，只是图形显示缩放。
tree_width_scale <- 0.4

# [距离] 气泡离物种名的距离 (绝对像素值)
# 作用：控制红色扩张气泡列距离树枝末端的空白大小。
species_name_space <- 100

# [距离] 两个气泡列之间的间距
# 作用：控制红色气泡和绿色气泡之间的缝隙。
bubble_col_gap <- 30

# [位置] 树上时间文字的垂直位置
# 作用：控制时间数字相对于树枝的上下偏移。1.2 表示向下沉，避免压线。
time_text_vjust <- 1.2

# [样式] 字体与颜色配置
opt_font_family <- "Arial"          # 无衬线字体
color_expansion <- "#E41A1C"        # 红色 (扩张)
color_contraction <- "#4DAF4A"      # 绿色 (收缩)
color_time_text <- "#1F78B4"        # 蓝色 (时间文字)
abbreviate_species_names <- FALSE   # 开启属名缩写 (如 Sinonovacula -> S.)

# ==============================================================================
# 2. 输入输出路径
# ==============================================================================

input_tree_path <- file.path("input", "FigTree.tre")
input_cafe_path <- file.path("input", "cafe_significant_clade_summary.tsv")
output_full_path <- file.path("output", "jindd_phylo_tree")

# 自动创建输出目录
if (!dir_exists(output_full_path)) {
    dir_create(output_full_path, recurse = TRUE)
}

# ==============================================================================
# 3. 数据读取与时间/物理缩放处理
# ==============================================================================

message("正在读取数据...")

# 读取 MCMCTree 的 NEXUS 格式树文件
tree <- read.mcmctree(input_tree_path)

# ------------------------------------------------------------------------------
# [开关逻辑] 是否将输入树中的分歧时间乘以100
# 说明：
#   这里直接修改树对象的 edge.length（枝长）。
#   这是最稳妥的做法，因为后续所有节点年龄、刻度、标签、树梢位置
#   都会自动基于新的时间体系重新计算，避免“只改标签不改树”的错误。
# ------------------------------------------------------------------------------
if (isTRUE(multiply_divtime_by_100)) {
    tree@phylo$edge.length <- tree@phylo$edge.length * 100
    message("已开启开关：输入 FigTree.tre 的分歧时间已乘以100。")
} else {
    message("开关关闭：使用 FigTree.tre 原始分歧时间。")
}

# ------------------------------------------------------------------------------
# [核心逻辑：物理缩放]
# 说明：
#   这一步只用于图形排版，不改变时间单位本身。
#   如果上面已经乘100，那么这里是在“乘100后的时间树”基础上再做物理压缩。
# ------------------------------------------------------------------------------
tree@phylo$edge.length <- tree@phylo$edge.length * tree_width_scale

# 转换为 phylo 对象进行计算
phy <- as.phylo(tree)
node_dist <- node.depth.edgelength(phy)

# scaled_tree_height: 缩放后的物理高度（用于绘图坐标定位）
scaled_tree_height <- max(node_dist)

# original_tree_height: 真实生物学时间高度（用于坐标轴和标签显示）
# 注意：这里要除以 tree_width_scale，把“物理压缩”还原掉
original_tree_height <- scaled_tree_height / tree_width_scale

# 计算每个节点用于显示的真实时间 (Age)
# 公式：(物理高度 - 节点距离) / 缩放比例
original_node_ages <- (scaled_tree_height - node_dist) / tree_width_scale

# ==============================================================================
# 4. 树表格整理与时间标签生成
# ==============================================================================

tree_tbl <- as_tibble(tree)
tree_tbl$node <- as.integer(tree_tbl$node)
tree_tbl$age_display <- original_node_ages[tree_tbl$node]
tree_tbl$is_tip_calculated <- tree_tbl$node <= length(phy$tip.label)

# --- 物种名处理函数 ---
process_label <- function(label) {
    clean <- gsub("_", " ", label)
    if (abbreviate_species_names) {
        parts <- str_split(clean, " ")[[1]]
        if (length(parts) >= 2) {
            return(paste0(toupper(substr(parts[1], 1, 1)), ". ", paste(parts[-1], collapse = " ")))
        }
    }
    return(clean)
}

# 应用名字处理
tree_tbl$display_label <- sapply(
    tree_tbl$label,
    function(x) if (x %in% phy$tip.label) process_label(x) else NA
)

# --- 时间标签生成 ---
tree_tbl_processed <- tree_tbl %>%
    mutate(
        # 提取置信区间 (HPD)
        hpd_lower = as.numeric(sapply(.data[["0.95HPD"]], function(x) if (is.null(x)) NA else x[1])),
        hpd_upper = as.numeric(sapply(.data[["0.95HPD"]], function(x) if (is.null(x)) NA else x[2])),

        # 生成标签字符串，使用真实时间 age_display
        # 说明：
        #   read.mcmctree() 读入的 0.95HPD 通常与枝长时间体系一致。
        #   因为我们前面已经直接修改了 tree 的 edge.length，
        #   所以 age_display 会随之变化。
        #   这里继续沿用原脚本写法，仅显示 age_display + 原有 HPD 字段。
        time_label = ifelse(
            !is_tip_calculated & !is.na(age_display) & !is.na(hpd_lower),
            sprintf("%.1f (%.1f-%.1f)", age_display, hpd_lower, hpd_upper),
            NA
        )
    )

tree_final <- as.treedata(tree_tbl_processed)

# ==============================================================================
# 5. 准备气泡图坐标（绝对对齐）
# ==============================================================================

# 先生成一个临时树图，提取 ggtree 计算出的叶节点 Y 坐标
p_temp <- ggtree(tree_final)
d_tips <- p_temp$data %>% filter(isTip == TRUE)
y_map <- setNames(d_tips$y, d_tips$label)

# 读取并清洗 CAFE 数据
cafe_raw <- read.table(input_cafe_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(cafe_raw) <- tolower(colnames(cafe_raw))
cafe_raw$species <- trimws(cafe_raw$species)

cafe_data <- cafe_raw %>% filter(species %in% names(y_map))

# 将树的 Y 坐标赋给气泡数据
cafe_data$y_pos <- y_map[cafe_data$species]

# [位置计算] 计算气泡的 X 坐标
# 气泡1 X坐标 = 树梢位置 + 空白距离
x_expansion <- scaled_tree_height + species_name_space

# 气泡2 X坐标 = 气泡1位置 + 间距
x_contraction <- x_expansion + bubble_col_gap

# 整理数据用于单图层绘制
exp_data <- cafe_data %>%
    select(y_pos, count = n_expanded_genes) %>%
    mutate(x_pos = x_expansion, type = "Expansion")

con_data <- cafe_data %>%
    select(y_pos, count = n_contracted_genes) %>%
    mutate(x_pos = x_contraction, type = "Contraction")

bubble_data <- bind_rows(exp_data, con_data)
size_rng <- range(bubble_data$count, na.rm = TRUE)

# ==============================================================================
# 6. 准备手绘坐标轴数据
# ==============================================================================

# 生成整百刻度 (0, 100, 200...)
original_breaks <- seq(0, floor(original_tree_height / 100) * 100, 100)

# 计算刻度在图上的物理 X 坐标
# 物理X = 缩放后树梢 - (真实时间 * 缩放比例)
breaks_x <- scaled_tree_height - (original_breaks * tree_width_scale)

# 定义轴线起点和终点
line_x_start <- 0
line_x_end <- scaled_tree_height

# ==============================================================================
# 7. 核心绘图
# ==============================================================================

message("开始绘图...")

p_final <- ggtree(tree_final, layout = "rectangular", linewidth = 0.5) +

    # [图层1: 物种名]
    geom_tiplab(aes(label = display_label), fontface = "bold.italic", offset = 5, size = 5) +

    # [图层2: 树上时间]
    geom_text(
        aes(label = time_label),
        fontface = "plain",
        size = 3.5,
        color = color_time_text,
        hjust = -0.05,
        vjust = time_text_vjust,
        na.rm = TRUE
    ) +

    # [图层3: 气泡点]
    geom_point(data = bubble_data, aes(x = x_pos, y = y_pos, size = count, color = type), alpha = 0.8) +
    scale_size_continuous(range = c(2, 8), limits = size_rng, guide = "none") +
    scale_color_manual(
        values = c("Expansion" = color_expansion, "Contraction" = color_contraction),
        guide = "none"
    ) +

    # [图层4: 气泡数字]
    geom_text(data = bubble_data, aes(x = x_pos, y = y_pos, label = count), size = 3.5, fontface = "plain", color = "black") +

    # [图层5: 气泡列标题]
    annotate(
        "text",
        x = x_expansion,
        y = length(phy$tip.label) + 1.5,
        label = "Expansion",
        color = color_expansion,
        fontface = "bold.italic",
        size = 4,
        angle = 0
    ) +
    annotate(
        "text",
        x = x_contraction,
        y = length(phy$tip.label) + 1.5,
        label = "Contraction",
        color = color_contraction,
        fontface = "bold.italic",
        size = 4,
        angle = 0
    ) +

    # [图层6: 手绘坐标轴]
    annotate("segment", x = line_x_start, xend = line_x_end, y = 0.5, yend = 0.5, linewidth = 0.5) +
    annotate("segment", x = breaks_x, xend = breaks_x, y = 0.5, yend = 0.3, linewidth = 0.5) +
    annotate("text", x = breaks_x, y = 0, label = sprintf("%.0f", original_breaks), size = 3.5, fontface = "plain") +
    annotate("text", x = line_x_end, y = 0, label = "Million years ago", hjust = -0.2, size = 4, fontface = "plain") +

    # [图层7: 画布范围]
    scale_x_continuous(limits = c(0, x_contraction + 80)) +
    scale_y_continuous(limits = c(-2, length(phy$tip.label) + 2)) +

    # [图层8: 主题清理]
    theme_void(base_family = opt_font_family) +
    theme(plot.margin = margin(t = 20, b = 20, l = 10, r = 10))

# ==============================================================================
# 8. 保存结果
# ==============================================================================

ggsave(
    file.path(output_full_path, "combined_phylogenomic_tree.pdf"),
    p_final,
    width = 18,
    height = 10,
    device = cairo_pdf
)

ggsave(
    file.path(output_full_path, "combined_phylogenomic_tree.png"),
    p_final,
    width = 18,
    height = 10,
    dpi = 300,
    bg = "white"
)

message("绘图完成。")
