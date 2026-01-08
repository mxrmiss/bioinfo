#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==========================================================================================
# MAGIC PLOTTING KIT - Statistics Plot (v1.3 Legend Fix)
# ==========================================================================================
# [功能] 读取 synteny_summary.tsv，绘制宏观统计图表
# [修复] 1. 修正 Chord 图例位置 (Top -> BottomLeft)
#        2. 修正 Chord 图例物种名格式 (下划线 -> 点)
# ==========================================================================================

# --- 1. 用户配置区 (USER CONFIGURATION) ---

# [输入文件] 默认自动寻找 input/synteny_summary.tsv
INPUT_FILE <- NULL 

# [输出目录设置]
CUSTOM_DIR_NAME <- "synteny_stats"
PREFIX         <- "macro_stats"

# [标签缩写] (应与主图保持一致)
CUSTOM_ABBREV <- NULL 

# [核心控制变量]
SHOW_TITLES <- TRUE
HEATMAP_TITLE <- NULL 
CHORD_TITLE   <- NULL 
SHOW_CHORD_LEGEND <- TRUE

# ==========================================================================================
# ⚠️ 核心逻辑区
# ==========================================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2) 
  library(pheatmap) 
  library(circlize) 
  library(data.table)
  library(stringr)
  library(grid)
})

# --- 1. 环境与数据加载 ---
DIR_ROOT   <- getwd()
DIR_INPUT  <- file.path(DIR_ROOT, "input")

out_dir_name <- if(is.null(CUSTOM_DIR_NAME)) "synteny_stats" else CUSTOM_DIR_NAME
DIR_OUTPUT_SUB <- file.path(DIR_ROOT, "output", out_dir_name)

if (!dir.exists(DIR_OUTPUT_SUB)) {
  dir.create(DIR_OUTPUT_SUB, recursive = TRUE)
  cat(sprintf(">> [Init] 已创建输出目录: output/%s/\n", out_dir_name))
}

cat("======================================================\n")
cat(">> Magic Plotting Kit: Macro Statistics (v1.3)\n")

if (is.null(INPUT_FILE)) {
  f_in <- file.path(DIR_INPUT, "synteny_summary.tsv")
} else {
  f_in <- INPUT_FILE
}

if (!file.exists(f_in)) {
  stop(paste0("[Error] 找不到输入文件: ", f_in, "\n请先运行 Python 脚本生成 summary 数据。"))
}

df <- fread(f_in)
cat(sprintf("   - 读取数据: %s (%d 行)\n", basename(f_in), nrow(df)))

# --- 2. 数据清洗与标签美化 ---
generate_label <- function(species_vec, chr_vec) {
  raw_species <- unique(species_vec)
  sp_map <- list()
  
  for (sp in raw_species) {
    if (!is.null(CUSTOM_ABBREV) && sp %in% names(CUSTOM_ABBREV)) {
      sp_map[[sp]] <- CUSTOM_ABBREV[[sp]]
    } else {
      parts <- unlist(strsplit(sp, "[_\\.]"))
      if (length(parts) >= 2) {
        abbr <- paste0(toupper(substr(parts[1], 1, 1)), tolower(substr(parts[2], 1, 2)))
      } else {
        abbr <- paste0(toupper(substr(sp, 1, 1)), tolower(substr(sp, 2, 3)))
      }
      sp_map[[sp]] <- abbr
    }
  }
  
  labels <- mapply(function(s, c) {
    num <- str_extract(c, "\\d+$")
    if (is.na(num)) num <- c
    paste0(sp_map[[s]], num)
  }, species_vec, chr_vec)
  
  return(labels)
}

df$lbl1 <- generate_label(df$species1, df$chr1)
df$lbl2 <- generate_label(df$species2, df$chr2)

# --- 3. 智能排序 ---
sort_labels <- function(lbls) {
  nums <- as.numeric(str_extract(lbls, "\\d+$"))
  if (all(is.na(nums))) return(sort(unique(lbls)))
  uniq_lbls <- unique(lbls)
  uniq_nums <- as.numeric(str_extract(uniq_lbls, "\\d+$"))
  uniq_lbls[order(uniq_nums)]
}

order1 <- sort_labels(df$lbl1)
order2 <- sort_labels(df$lbl2)


# ==========================================================================================
# 图表 A: 共线性强度热图 (Heatmap)
# ==========================================================================================
cat(">> [Plot 1] 正在绘制热图 (Heatmap)...\n")

h_title <- if (SHOW_TITLES) {
    if (!is.null(HEATMAP_TITLE)) HEATMAP_TITLE else "Synteny Block Counts"
} else {
    NULL 
}

mat_data <- dcast(df, lbl1 ~ lbl2, value.var = "n_blocks", fill = 0)
mat_data <- as.data.frame(mat_data) 
rownames(mat_data) <- mat_data$lbl1
mat_data$lbl1 <- NULL
mat <- as.matrix(mat_data)

valid_r <- order1[order1 %in% rownames(mat)]
valid_c <- order2[order2 %in% colnames(mat)]
mat <- mat[valid_r, valid_c, drop=FALSE]

out_pdf_heat <- file.path(DIR_OUTPUT_SUB, paste0(PREFIX, "_heatmap.pdf"))
out_png_heat <- file.path(DIR_OUTPUT_SUB, paste0(PREFIX, "_heatmap.png"))
heatmap_cols <- colorRampPalette(c("white", "#FFF7BC", "#FEC44F", "#D95F0E", "#993404"))(100)

pheatmap(mat, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         color = heatmap_cols,
         display_numbers = TRUE, number_format = "%.0f", 
         fontsize_number = 8,
         border_color = "grey90",
         main = h_title,
         filename = out_pdf_heat, width = 8, height = 8)

pheatmap(mat, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = heatmap_cols,
         display_numbers = TRUE, number_format = "%.0f",
         fontsize_number = 8,
         border_color = "grey90",
         main = h_title,
         filename = out_png_heat, width = 8, height = 8)

cat(sprintf("   - 生成: %s\n", basename(out_pdf_heat)))

# ==========================================================================================
# 图表 B: 宏观关系弦图 (Chord Diagram)
# ==========================================================================================
cat(">> [Plot 2] 正在绘制弦图 (Chord Diagram)...\n")

out_pdf_chord <- file.path(DIR_OUTPUT_SUB, paste0(PREFIX, "_chord.pdf"))
out_png_chord <- file.path(DIR_OUTPUT_SUB, paste0(PREFIX, "_chord.png"))

chord_df <- df[, c("lbl1", "lbl2", "n_blocks")]

draw_chord <- function() {
  # 1. 标题设置
  c_title <- if (SHOW_TITLES) {
    if (!is.null(CHORD_TITLE)) CHORD_TITLE else "Macro Synteny Relationships"
  } else {
    "" 
  }
  
  # 2. 布局设置
  circos.clear()
  gaps <- c(rep(1, length(order1)-1), 10, rep(1, length(order2)-1), 10)
  
  # 调整边距：增加底部边距，确保 legend 不被裁切
  if (SHOW_CHORD_LEGEND) {
    par(mar = c(4, 1, 3, 1)) # 底部留出更多空间
  } else {
    par(mar = c(1, 1, 1, 1))
  }
  
  circos.par(start.degree = 90, gap.after = gaps)
  
  # 3. 颜色 (物种1=红系，物种2=蓝系)
  sp1_name_raw <- unique(df$species1)
  sp2_name_raw <- unique(df$species2)
  
  # 格式化图例名称 (下划线 -> 点)
  sp1_name_formatted <- gsub("_", ". ", sp1_name_raw)
  sp2_name_formatted <- gsub("_", ". ", sp2_name_raw)
  
  grid_col <- c(
    setNames(colorRampPalette(c("#FEE0D2", "#DE2D26"))(length(order1)), order1),
    setNames(colorRampPalette(c("#DEEBF7", "#3182BD"))(length(order2)), order2)
  )
  
  # 4. 绘图
  chordDiagram(chord_df, 
               order = c(order1, order2), 
               grid.col = grid_col,
               transparency = 0.3,
               annotationTrack = "grid", 
               preAllocateTracks = list(track.height = 0.1))
  
  # 5. 添加标签
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
  }, bg.border = NA)
  
  # 绘制标题
  title(c_title)
  
  # 6. 图例 (移至左下角)
  if (SHOW_CHORD_LEGEND) {
    legend("topright", 
           legend = c(sp1_name_formatted, sp2_name_formatted), # 格式化后的名称
           fill = c(grid_col[order1[1]], grid_col[order2[1]]), # 取每类物种的第一个颜色作为代表
           bty = "n", cex = 0.9,
           text.font = 3, # 斜体
           title = "Species Segments")
  }
}

# 输出 PDF
pdf(out_pdf_chord, width = 8, height = 8)
draw_chord()
dev.off()

# 输出 PNG
png(out_png_chord, width = 2400, height = 2400, res = 300)
draw_chord()
dev.off()

cat(sprintf("   - 生成: %s\n", basename(out_pdf_chord)))
cat(">> [Success] 统计分析完成！请查看 output/synteny_stats/ 目录。\n")
