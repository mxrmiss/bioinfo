#!/usr/bin/env Rscript

# ==============================================================================
# 脚本名称: phylo_tree.R (V33 - 详尽注释·最终修正版)
# 功能描述:
#   1. [核心逻辑] 采用 "单画布叠加" + "手绘坐标轴" 策略，彻底解决对齐问题。
#   2. [坐标轴] 轴线从树根(0)一直画到树梢，左端对齐；标题紧跟在 0 刻度右侧。
#   3. [样式] 严格执行去粗存精，仅物种名加粗，其余全部常规细体。
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
# --- [用户参数锁定区] (这些是你调试好的最佳参数) ---
# ==============================================================================

# [缩放] 树的宽度缩放比例
# 作用：物理压缩树的长度。0.4 表示压缩到原长的 40%，给右侧气泡腾出空间。
tree_width_scale <- 0.4

# [距离] 气泡离物种名的距离 (绝对像素值)
# 作用：控制红色扩张气泡列距离树枝末端的空白大小。
species_name_space <- 80

# [距离] 两个气泡列之间的间距
# 作用：控制红色气泡和绿色气泡之间的缝隙。
bubble_col_gap <- 30

# [位置] 树上时间文字的垂直位置
# 作用：控制时间数字相对于树枝的上下偏移。1.2 表示向下沉，避免压线。
time_text_vjust <- 1.2

# [样式] 字体与颜色配置
opt_font_family <- "Arial"          # 无衬线字体
color_expansion <- "#E41A1C"       # 红色 (扩张)
color_contraction <- "#4DAF4A"     # 绿色 (收缩)
color_time_text <- "#1F78B4"       # 蓝色 (时间文字)
abbreviate_species_names <- FALSE   # 开启属名缩写 (如 Sinonovacula -> S.)

# ==============================================================================

# 定义输入输出路径
input_tree_path <- file.path("input", "FigTree.tre")
input_cafe_path <- file.path("input", "cafe_significant_clade_summary.tsv")
output_full_path <- file.path("output", "phylo_tree")
# 自动创建输出目录
if (!dir_exists(output_full_path)) dir_create(output_full_path, recurse = TRUE)

# ------------------------------------------------------------------------------
# 2. 数据读取与物理缩放处理
# ------------------------------------------------------------------------------
message("正在读取数据...")
# 读取 MCMCTree 的 NEXUS 格式树文件
tree <- read.mcmctree(input_tree_path)

# [核心逻辑：物理缩放]
# 直接修改树对象的 edge.length (枝长) 数据。
# 这样 ggtree 在绘图时，会认为树真的变短了。
# 所有的坐标计算都将基于这个“变短后”的长度，确保气泡紧跟树梢。
tree@phylo$edge.length <- tree@phylo$edge.length * tree_width_scale

# 转换为 phylo 对象进行计算
phy <- as.phylo(tree)
node_dist <- node.depth.edgelength(phy)

# scaled_tree_height: 缩放后的物理高度 (用于绘图坐标定位)
scaled_tree_height <- max(node_dist) 
# original_tree_height: 原始的生物学高度 (用于计算刻度显示的数值)
original_tree_height <- scaled_tree_height / tree_width_scale

# 计算每个节点用于显示的原始时间 (Age)
# 公式：(物理高度 - 节点距离) / 缩放比例
original_node_ages <- (scaled_tree_height - node_dist) / tree_width_scale

# 转换为 tibble 数据框
tree_tbl <- as_tibble(tree)
tree_tbl$node <- as.integer(tree_tbl$node)
tree_tbl$age_display <- original_node_ages[tree_tbl$node] # 存入真实时间
tree_tbl$is_tip_calculated <- tree_tbl$node <= length(phy$tip.label)

# --- 物种名处理函数 ---
process_label <- function(label) {
    clean <- gsub("_", " ", label) # 去掉下划线
    if (abbreviate_species_names) {
        parts <- str_split(clean, " ")[[1]]
        # 取属名首字母大写 + 种名
        if (length(parts) >= 2) return(paste0(toupper(substr(parts[1],1,1)), ". ", paste(parts[-1], collapse=" ")))
    }
    return(clean)
}
# 应用名字处理
tree_tbl$display_label <- sapply(tree_tbl$label, function(x) if(x %in% phy$tip.label) process_label(x) else NA)

# --- 时间标签生成 ---
tree_tbl_processed <- tree_tbl %>%
  mutate(
    # 提取置信区间 (HPD)
    hpd_lower = as.numeric(sapply(.data[["0.95HPD"]], function(x) if(is.null(x)) NA else x[1])),
    hpd_upper = as.numeric(sapply(.data[["0.95HPD"]], function(x) if(is.null(x)) NA else x[2])),
    # 生成标签字符串，使用真实时间 age_display
    time_label = ifelse(!is_tip_calculated & !is.na(age_display) & !is.na(hpd_lower),
                        sprintf("%.1f (%.1f-%.1f)", age_display, hpd_lower, hpd_upper), NA)
  )
tree_final <- as.treedata(tree_tbl_processed)

# ------------------------------------------------------------------------------
# 3. 准备气泡图坐标 (绝对对齐)
# ------------------------------------------------------------------------------
# 先生成一个临时树图，提取 ggtree 计算出的叶节点 Y 坐标
# 这是实现“绝对水平对齐”的关键步骤
p_temp <- ggtree(tree_final)
d_tips <- p_temp$data %>% filter(isTip == TRUE)
y_map <- setNames(d_tips$y, d_tips$label)

# 读取并清洗 CAFE 数据
cafe_raw <- read.table(input_cafe_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
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
exp_data <- cafe_data %>% select(y_pos, count=n_expanded_genes) %>% mutate(x_pos = x_expansion, type="Expansion")
con_data <- cafe_data %>% select(y_pos, count=n_contracted_genes) %>% mutate(x_pos = x_contraction, type="Contraction")
bubble_data <- bind_rows(exp_data, con_data)
size_rng <- range(bubble_data$count, na.rm=TRUE)

# ------------------------------------------------------------------------------
# 4. 准备手绘坐标轴数据 (核心修正逻辑)
# ------------------------------------------------------------------------------
# [步骤1] 生成整百刻度 (0, 100, 200...)
# floor(...) 确保不超过最大时间
original_breaks <- seq(0, floor(original_tree_height / 100) * 100, 100)

# [步骤2] 计算刻度在图上的物理 X 坐标
# 物理X = 缩放后树梢 - (真实时间 * 缩放比例)
breaks_x <- scaled_tree_height - (original_breaks * tree_width_scale)

# [步骤3] 定义轴线的起点和终点
# 起点(左): 强制设为 0 (即树的根部位置)，确保轴线与树左对齐
line_x_start <- 0
# 终点(右): 设为树梢位置 (即 0 Ma 的位置)
line_x_end <- scaled_tree_height

# ------------------------------------------------------------------------------
# 5. 核心绘图 (叠加图层)
# ------------------------------------------------------------------------------
message("开始绘图...")

p_final <- ggtree(tree_final, layout="rectangular", linewidth=0.5) +
  
  # [图层1: 物种名] 粗斜体，offset=5 稍微离开树枝
  geom_tiplab(aes(label=display_label), fontface="bold.italic", offset=5, size=5) +
  
  # [图层2: 树上时间] 蓝色，常规体，位置下沉 (time_text_vjust)
  geom_text(aes(label=time_label), fontface="plain", size=3.5, color=color_time_text, 
            hjust=-0.05, vjust=time_text_vjust, na.rm=TRUE) +
  
  # [图层3: 气泡点]
  geom_point(data=bubble_data, aes(x=x_pos, y=y_pos, size=count, color=type), alpha=0.8) +
  scale_size_continuous(range=c(2, 8), limits=size_rng, guide="none") +
  scale_color_manual(values=c("Expansion"=color_expansion, "Contraction"=color_contraction), guide="none") +
  
  # [图层4: 气泡数字] 黑色，常规体，不加粗
  geom_text(data=bubble_data, aes(x=x_pos, y=y_pos, label=count), size=3.5, fontface="plain", color="black") +
  
  # [图层5: 气泡列标题] 粗斜体，位置悬空
  annotate("text", x=x_expansion, y=length(phy$tip.label)+1.5, label="Expansion", 
           color=color_expansion, fontface="bold.italic", size=4, angle=0) +
  annotate("text", x=x_contraction, y=length(phy$tip.label)+1.5, label="Contraction", 
           color=color_contraction, fontface="bold.italic", size=4, angle=0) +

  # [图层6: 手绘坐标轴]
  # 6.1 轴线: 黑色实线，从树根(0)画到树梢(line_x_end)
  annotate("segment", x=line_x_start, xend=line_x_end, y=0.5, yend=0.5, linewidth=0.5) +
  
  # 6.2 刻度线: 只在整百位置画
  annotate("segment", x=breaks_x, xend=breaks_x, y=0.5, yend=0.3, linewidth=0.5) +
  
  # 6.3 刻度数字: 显示整百数值
  annotate("text", x=breaks_x, y=0, label=sprintf("%.0f", original_breaks), size=3.5, fontface="plain") +
  
  # 6.4 坐标轴标题 "Million years ago"
  # x = line_x_end (树梢，即0刻度位置)
  # y = 0 (与刻度数字同一高度)
  # hjust = -0.2 (放在0的右侧，留一点点空隙，看起来紧挨着)
  annotate("text", x=line_x_end, y=0, label="Million years ago", hjust=-0.2, size=4, fontface="plain") +

  # [图层7: 画布范围]
  # 右边界延伸：包含气泡区 + 额外80单位给右边的文字
  scale_x_continuous(limits=c(0, x_contraction + 80)) +
  # 上下边界：留出空间给坐标轴和标题
  scale_y_continuous(limits=c(-2, length(phy$tip.label) + 2)) +
  
  # [图层8: 主题清理] 去掉所有默认背景和坐标轴
  theme_void(base_family = opt_font_family) +
  theme(plot.margin = margin(t=20, b=20, l=10, r=10))

# 6. 保存结果
ggsave(file.path(output_full_path, "combined_phylogenomic_tree.pdf"), p_final, width=18, height=10, device=cairo_pdf)
ggsave(file.path(output_full_path, "combined_phylogenomic_tree.png"), p_final, width=18, height=10, dpi=300, bg="white")

message("绘图完成。")
