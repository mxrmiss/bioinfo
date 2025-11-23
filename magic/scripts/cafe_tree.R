#!/usr/bin/env Rscript
# ==============================================================================
# MAGIC PLOTTING KIT - cafe_tree.R
# [功能] 绘制 Nature/Science 风格的 TimeTree
# [关键点]
#   - CAFE 数字映射到叶节点和内部节点（双列隔离：Label_Tip / Label_Int）
#   - 坐标轴使用绝对年代（Ma），从最大值到 0，无负数
#   - 不使用裁剪（coord_cartesian），避免任何图形被截断
# ==============================================================================

# --- A. 用户配置区 (User Configuration) ---
CONFIG <- list(
  INPUT_TREE_FILE = "FigTree.tre",
  INPUT_CAFE_FILE = "Base_clade_results.txt",

  OUTPUT_DIR_NAME = "cafe_tree",
  OUTPUT_FILE_NAME = "cafe_tree",

  SHOW_TITLE = TRUE,
  TITLE_TEXT = "Gene Family Evolution",

  SHOW_HPD_BARS = FALSE,

  ABBREVIATE_GENUS = TRUE,

  # [绘图参数]
  FONT_FAMILY = "sans",
  FONT_SIZE_TIP = 3.5,     # 物种名大小
  FONT_SIZE_NUM = 3.5,     # 枝上数值大小

  # [配色]
  COLOR_EXP = "#009E73",   # 扩张 (Green)
  COLOR_CON = "#D55E00",   # 收缩 (Orange)
  COLOR_HPD = "#56B4E9"    # 时间置信区间 (Blue)
)

# --- B. 环境加载 ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(aplot)
  library(ape)

  if (!requireNamespace("ggtext", quietly = TRUE)) {
    message(">>> Installing required package: ggtext...")
    install.packages("ggtext", repos = "https://cloud.r-project.org")
  }
  library(ggtext)
})

MAGIC_ROOT <- getwd()
INPUT_PATH <- file.path(MAGIC_ROOT, "input")
OUTPUT_PATH <- file.path(MAGIC_ROOT, "output", CONFIG$OUTPUT_DIR_NAME)

if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH, recursive = TRUE)

FILE_TREE <- file.path(INPUT_PATH, CONFIG$INPUT_TREE_FILE)
FILE_CAFE <- file.path(INPUT_PATH, CONFIG$INPUT_CAFE_FILE)

# --- C. 数据读取与清洗 ---

if (!file.exists(FILE_TREE) || !file.exists(FILE_CAFE)) {
  stop("错误: 输入文件不存在。")
}

message(">>> Reading and patching input tree...")

# UTREE -> TREE
raw_tree_lines <- readLines(FILE_TREE)
if (any(grepl("UTREE", raw_tree_lines))) {
  raw_tree_lines <- gsub("UTREE", "TREE", raw_tree_lines)
}
temp_tree_file <- tempfile(fileext = ".tre")
writeLines(raw_tree_lines, temp_tree_file)

tr <- tryCatch({
  read.beast(temp_tree_file)
}, error = function(e) {
  read.tree(temp_tree_file)
})
unlink(temp_tree_file)

# 0.95HPD -> x0.95HPD（避免以数字开头的列名）
if (inherits(tr, "treedata")) {
  data_cols <- colnames(tr@data)
  hpd_indices <- grep("HPD", data_cols)
  if (length(hpd_indices) > 0) {
    for (i in hpd_indices) {
      if (grepl("^[0-9]", data_cols[i])) {
        colnames(tr@data)[i] <- paste0("x", data_cols[i])
      }
    }
  }
}

message(">>> Reading CAFE results...")
cafe_raw <- read.table(
  FILE_CAFE,
  header       = TRUE,
  comment.char = "",
  check.names  = FALSE,
  stringsAsFactors = FALSE
)
colnames(cafe_raw)[1] <- "Taxon_ID"

clean_and_format_label <- function(raw_names, abbreviate = TRUE) {
  sapply(raw_names, function(x) {
    if (is.na(x) || x == "") return(x)
    clean_x <- str_remove(x, "<\\d+>")
    clean_x <- gsub("'", "", clean_x)
    if (abbreviate) {
      gsub("^([A-Za-z])[^_]+_([a-zA-Z0-9]+).*", "\\1.\\2", clean_x)
    } else {
      gsub("_", " ", clean_x)
    }
  })
}

# --- D. 绘图核心逻辑 (双层隔离) ---

message(">>> Constructing phylogeny plot...")

if (inherits(tr, "treedata")) {
  phylo_obj <- tr@phylo
} else {
  phylo_obj <- tr
}

formatted_tips <- clean_and_format_label(
  phylo_obj$tip.label,
  abbreviate = CONFIG$ABBREVIATE_GENUS
)

if (inherits(tr, "treedata")) {
  tr@phylo$tip.label <- formatted_tips
} else {
  tr$tip.label <- formatted_tips
}
phylo_obj$tip.label <- formatted_tips
tree_tbl_new <- as_tibble(tr)

# 3. 构建映射表：CAFE → 树节点
cafe_clean <- cafe_raw %>%
  mutate(
    Raw_Label = str_remove(Taxon_ID, "<\\d+>"),
    Match_Label = clean_and_format_label(
      Raw_Label,
      abbreviate = CONFIG$ABBREVIATE_GENUS
    ),
    extracted_node_id = as.integer(
      str_extract(Taxon_ID, "(?<=<)\\d+(?=>)")
    ),

    # HTML 彩色标签（无空格版：+541/441）
    Label_HTML = paste0(
      "<span style='color:", CONFIG$COLOR_EXP, "'>+", Increase, "</span>",
      "<span style='color:black'>/</span>",
      "<span style='color:", CONFIG$COLOR_CON, "'>-", Decrease, "</span>"
    )
  )

mapping_data <- cafe_clean %>%
  rowwise() %>%
  mutate(
    # 优先按物种名匹配叶节点，其次按 CAFE 提供的 node id
    node = if (Match_Label %in% formatted_tips) {
      match(Match_Label, formatted_tips)
    } else {
      extracted_node_id
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(node))

# 区分 Tip / Internal
N_tips <- length(phylo_obj$tip.label)

final_mapping <- mapping_data %>%
  mutate(
    Label_Tip = ifelse(node <= N_tips, Label_HTML, NA),
    Label_Int = ifelse(node >  N_tips, Label_HTML, NA)
  ) %>%
  select(node, Label_Tip, Label_Int)

message(paste(
  "    Data split: Tips =", sum(!is.na(final_mapping$Label_Tip)),
  "| Internals =", sum(!is.na(final_mapping$Label_Int))
))

# 4. 计算时间轴（在 Ma 空间生成刻度，再映射回树坐标）
max_age_raw <- max(node.depth.edgelength(phylo_obj))
time_multiplier <- 1
if (max_age_raw < 50) time_multiplier <- 100
age_max <- max_age_raw * time_multiplier   # 根节点年代（Ma）

# 在 “Ma” 空间生成刻度
age_breaks <- pretty(c(0, age_max), n = 8)
age_breaks <- age_breaks[age_breaks >= 0 & age_breaks <= age_max]
age_breaks <- sort(age_breaks, decreasing = TRUE)   # 例如 700,600,...,0

# 将 Ma 刻度映射回树的 x 坐标：root 在 0 Ma，tips 在 0 Ma？（我们让 0 在 tips）
# 这里使用简单线性映射：age = age_max 对应 x=0；age=0 对应 x=max_age_raw
x_breaks <- (age_max - age_breaks) / time_multiplier

# 5. 绘图（不使用 coord_cartesian，避免截断）
p <- ggtree(tr, linewidth = 0.3) +
  theme_tree2() +
  scale_x_continuous(
    breaks = x_breaks,
    labels = age_breaks,
    expand = c(0, 0)
  ) +
  xlab("million years ago") +
  hexpand(0.8) +  # 右侧预留空间（给物种名和 CAFE 数字）
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(
      color  = "black",
      size   = 10,
      family = CONFIG$FONT_FAMILY
    ),
    axis.title.x = element_text(
      color  = "black",
      size   = 11,
      family = CONFIG$FONT_FAMILY,
      hjust  = 1          # 标题靠右
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 60, b = 20, l = 20, unit = "pt")
  )

# 6. HPD 条柱（当前默认关闭）
if (CONFIG$SHOW_HPD_BARS) {
  hpd_col <- grep("HPD", colnames(tree_tbl_new), value = TRUE)[1]
  if (!is.na(hpd_col)) {
    p <- p + geom_range(
      range = hpd_col,
      color = CONFIG$COLOR_HPD,
      alpha = 0.3,
      linewidth = 1
    )
  }
}

# 7. Tip Labels (物种名)
p <- p + geom_tiplab(
  size    = CONFIG$FONT_SIZE_TIP,
  fontface = "italic",
  offset  = 0.01 * max_age_raw,
  family  = CONFIG$FONT_FAMILY,
  align   = TRUE,
  linesize = 0.15,
  linetype = "dotted"
)

# 注入 CAFE 映射数据
p <- p %<+% final_mapping

# 8. 内部节点数字
p <- p + geom_richtext(
  aes(label = Label_Int),
  hjust = 1,
  nudge_x = - (max_age_raw * 0.01),
  vjust = -0.4,
  size   = CONFIG$FONT_SIZE_NUM,
  family = CONFIG$FONT_FAMILY,
  fontface = "plain",   # 不加粗
  fill   = NA,
  label.color   = NA,
  label.padding = unit(0, "lines"),
  na.rm = TRUE
)

# 9. 叶节点数字（物种右侧）
p <- p + geom_richtext(
  aes(label = Label_Tip),
  nudge_x = 0.20 * max_age_raw,
  hjust   = 0,
  size    = CONFIG$FONT_SIZE_NUM,
  family  = CONFIG$FONT_FAMILY,
  fontface = "plain",   # 不加粗
  fill    = NA,
  label.color = NA,
  na.rm   = TRUE
)

# 10. 标题与图例
y_top <- length(phylo_obj$tip.label)

if (CONFIG$SHOW_TITLE) {
  p <- p + ggtitle(CONFIG$TITLE_TEXT) +
    theme(
      plot.title = element_text(
        hjust  = 0.5,
        face   = "bold",
        size   = 14,
        family = CONFIG$FONT_FAMILY
      )
    )
}

p <- p +
  annotate(
    "text",
    x = 0, y = y_top + 1.2,
    label = "Gene Families:",
    hjust = 0,
    fontface = "bold",
    size = 3,
    family = CONFIG$FONT_FAMILY
  ) +
  annotate(
    "text",
    x = 0, y = y_top + 0.6,
    label = "Expansion",
    color = CONFIG$COLOR_EXP,
    hjust = 0,
    size = 3,
    fontface = "bold",
    family = CONFIG$FONT_FAMILY
  ) +
  annotate(
    "text",
    x = 0, y = y_top + 0.0,
    label = "Contraction",
    color = CONFIG$COLOR_CON,
    hjust = 0,
    size = 3,
    fontface = "bold",
    family = CONFIG$FONT_FAMILY
  )

# --- E. 保存 ---

plot_height <- max(6, length(phylo_obj$tip.label) * 0.5)

pdf_file <- file.path(OUTPUT_PATH, paste0(CONFIG$OUTPUT_FILE_NAME, ".pdf"))
ggsave(pdf_file, p, width = 12, height = plot_height, device = cairo_pdf)

png_file <- file.path(OUTPUT_PATH, paste0(CONFIG$OUTPUT_FILE_NAME, ".png"))
ggsave(png_file, p, width = 12, height = plot_height, dpi = 300)

message("========================================================")
message("[SUCCESS] Plots generated successfully!")
message(paste("Output Directory:", OUTPUT_PATH))
message(paste("Files           :", basename(pdf_file), "&", basename(png_file)))
message("========================================================")

