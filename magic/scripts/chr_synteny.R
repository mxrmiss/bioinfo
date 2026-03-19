#!/usr/bin/env Rscript

# ==========================================================================================
# MAGIC PLOTTING KIT - Synteny Circle Plot (v5.5 Fixed Function)
# ==========================================================================================
# [脚本名称] chr_synteny.R
# [功能] 读取 MCScanX/JCVI 产出的 Block 数据，绘制顶刊级共线性绸带图
# [修复] 修正函数名拼写错误 (circos.genomicRibbon -> circos.genomicLink)
#
# [视觉特性]
# 1. 沉浸式刻度: 色块铺满轨道，刻度尺直接“咬合”在色块内侧边缘，无多余留白
# 2. 自然截断: 刻度只画到整数位 (如 48Mb 画到 40Mb)，余数留白，拒绝臃肿
# 3. 丝滑绸带: 使用 genomicLink 绘制宽带，展示宏观共线性
# 4. 全能配色: 内置 7 种风格，默认 Nature (NPG) 风格
#
# [用法] ./scripts/chr_synteny.R
# ==========================================================================================

# ==========================================================================================
# 1. 用户配置区 (USER CONFIGURATION) - 请按需修改
# ==========================================================================================

# --- [A. 输入/输出设置] ---
# 默认 NULL (自动扫描 input/ 文件夹)
INPUT_GENOME_FILE <- NULL
INPUT_LINKS_FILE  <- NULL
CUSTOM_OUT_DIR    <- NULL   # 自定义输出子文件夹
CUSTOM_FILENAME   <- NULL   # 自定义输出文件名前缀

# --- [B. 显示开关] ---
SHOW_LEGEND <- TRUE   # 是否显示图例
SHOW_TITLE  <- FALSE   # 是否显示标题
SHOW_TICKS  <- TRUE   # 是否显示刻度尺

# --- [C. 标签与命名] ---
# 自动缩写: S_constricta -> Sco
# 如需手动指定: CUSTOM_ABBREV <- c("S_constricta" = "Sco", "S_rivularis" = "Sri")
CUSTOM_ABBREV <- NULL

# --- [D. 配色风格选择] ---
# 1. "nature"      : [默认] NPG 风格。红蓝绿灰主调，稳重优雅，适合大多数期刊。
# 2. "distinct_12" : 12色高反差循环。红蓝绿紫交替，界限分明，防混淆。
# 3. "science"     : Science (AAAS) 风格。深红深蓝高饱和度，对比强烈。
# 4. "viridis"     : Viridis 渐变。黄-绿-紫，色盲友好，现代感强。
# 5. "pastel"      : 糖果色 (Pastel)。色彩柔和，适合做背景。
# 6. "simple"      : 极简双色。物种A全红系，物种B全蓝系。
# 7. "rainbow"     : 经典彩虹色。
COLOR_STYLE     <- "nature"

# --- [E. 绘图细节] ---
LINK_ALPHA      <- 0.5       # 绸带透明度 (0.5 最佳)
GAP_DEGREE      <- 1.5       # [调整] 染色体间隙 (默认 1.5，更紧凑)
BIG_GAP         <- 5         # [调整] 物种间间隙 (默认 5，不用分太开)
LABEL_SIZE      <- 0.8       # 标签字号
AXIS_GAP        <- 10        # 刻度间隔 (Mb)
LABEL_MIN_SIZE_MB <- 10      # 只有长度 ≥ 该阈值 (Mb) 的染色体才绘制标签

# ==========================================================================================
# ⚠️ 核心逻辑区 (CORE LOGIC)
# ==========================================================================================

suppressPackageStartupMessages({
  library(circlize)
  library(RColorBrewer)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(grDevices)
})

# 将标签长度阈值转换为碱基数
LABEL_MIN_SIZE <- LABEL_MIN_SIZE_MB * 1000000

# --- 1. 环境准备 ---
DIR_ROOT   <- getwd()
DIR_INPUT  <- file.path(DIR_ROOT, "input")
DIR_OUTPUT_SUB <- file.path(DIR_ROOT, "output", if(!is.null(CUSTOM_OUT_DIR)) CUSTOM_OUT_DIR else "synteny")
if (!dir.exists(DIR_OUTPUT_SUB)) dir.create(DIR_OUTPUT_SUB, recursive = TRUE)

# --- 2. 智能读取数据 ---
cat("======================================================\n")
cat(">> Magic Plotting Kit: Chromosome Synteny (v5.5 Fixed)\n")

if (is.null(INPUT_GENOME_FILE) || is.null(INPUT_LINKS_FILE)) {
  all_files <- list.files(DIR_INPUT, pattern = "\\.tsv$", full.names = TRUE)
  if (length(all_files) == 0) stop("[错误] input/ 文件夹为空。")

  f_genome <- all_files[grepl("genome_len", basename(all_files))][1]
  f_links  <- all_files[grepl("synteny_blocks", basename(all_files))][1]
  if (is.na(f_links)) f_links <- all_files[grepl("gene_links", basename(all_files))][1]

  if (is.na(f_genome) || is.na(f_links)) stop("[错误] 缺少数据文件。")

  df_chr <- as.data.frame(fread(f_genome))
  df_link <- as.data.frame(fread(f_links))
} else {
  df_chr <- as.data.frame(fread(file.path(DIR_INPUT, INPUT_GENOME_FILE)))
  df_link <- as.data.frame(fread(file.path(DIR_INPUT, INPUT_LINKS_FILE)))
}

req_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
if (all(req_cols %in% colnames(df_link))) {
  df_ribbon <- df_link[, req_cols]
} else {
  df_ribbon <- df_link[, 1:6]
  colnames(df_ribbon) <- req_cols
}

# --- 3. 智能命名处理 ---
raw_species <- unique(df_chr$species)
sp_map <- list()
legend_map <- list()

for (sp in raw_species) {
  if (!is.null(CUSTOM_ABBREV) && sp %in% names(CUSTOM_ABBREV)) {
    abbr <- CUSTOM_ABBREV[[sp]]
  } else {
    parts <- unlist(strsplit(sp, "[_\\.]"))
    if (length(parts) >= 2) {
      abbr <- paste0(toupper(substr(parts[1], 1, 1)), tolower(substr(parts[2], 1, 2)))
    } else {
      abbr <- paste0(toupper(substr(sp, 1, 1)), tolower(substr(sp, 2, 3)))
    }
  }
  sp_map[[sp]] <- abbr
  legend_map[[abbr]] <- gsub("_", ". ", sp)
}

# 按物种 + 染色体名排序，并在每个物种内部重新编号 01,02,... 作为标签编号
df_chr <- df_chr %>%
  group_by(species) %>%
  arrange(as.integer(str_extract(as.character(chr), "\\d+$")), chr, .by_group = TRUE) %>%
  mutate(chr_index = row_number()) %>%
  ungroup()

df_chr$label <- mapply(function(idx, sp, chrname) {
  num <- str_extract(chrname, "\\d+$")
  if (is.na(num) || num == "") {
    num <- sprintf("%02d", idx)
  } else {
    num <- sprintf("%02d", as.integer(num))
  }
  paste0(sp_map[[sp]], num)
}, df_chr$chr_index, df_chr$species, as.character(df_chr$chr))

df_chr$chr <- factor(df_chr$chr, levels = df_chr$chr)

# --- 4. 全能配色引擎 ---
generate_palette <- function(style, n) {
  if (style == "nature") {
    # NPG 风格 (Nature)
    cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148")
    return(colorRampPalette(cols)(n))

  } else if (style == "distinct_12") {
    # 12色高反差
    cols <- c("#E41A1C", "#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00", "#B15928",
              "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99")
    return(rep(cols, length.out = n))

  } else if (style == "science") {
    # AAAS 风格
    cols <- c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021", "#5F559B")
    return(colorRampPalette(cols)(n))

  } else if (style == "viridis") {
    return(hcl.colors(n, palette = "Viridis"))

  } else if (style == "pastel") {
    cols <- brewer.pal(12, "Set3")
    return(rep(cols, length.out = n))

  } else if (style == "simple") {
    sp_counts <- table(df_chr$species)[raw_species]
    col_list <- c()
    main_cols <- c("#EE0000", "#3B4992", "#009900", "#FF9900")
    for (i in 1:length(raw_species)) {
      pal <- colorRampPalette(c("white", main_cols[i]))(sp_counts[i] + 3)[4:(sp_counts[i]+3)]
      col_list <- c(col_list, pal)
    }
    return(col_list)

  } else {
    return(rainbow(n))
  }
}
grid_col <- setNames(generate_palette(COLOR_STYLE, nrow(df_chr)), df_chr$chr)

# 间隙设置
gaps <- rep(GAP_DEGREE, nrow(df_chr))
sp_vec <- as.character(df_chr$species)
change_idx <- which(sp_vec[-length(sp_vec)] != sp_vec[-1])
gaps[change_idx] <- BIG_GAP
gaps[length(gaps)] <- BIG_GAP

# --- 5. 绘图核心函数 ---
draw_core_plot <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.degree = gaps, cell.padding = c(0, 0, 0, 0), track.margin=c(0.01, 0.01))

  circos.initialize(factors = df_chr$chr, xlim = cbind(0, df_chr$size))

  # === 轨道 1: 染色体 + 沉浸式刻度 ===
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr_real_name = CELL_META$sector.index
    chr_label = df_chr$label[df_chr$chr == chr_real_name]
    chr_size  = df_chr$size[df_chr$chr == chr_real_name][1]
    xlim = CELL_META$xlim

    # 1. 染色体色块 (铺满整个轨道，0到1)
    # [集成设计] 不再留白，让色块填满
    circos.rect(xlim[1], 0, xlim[2], 1, col = grid_col[chr_real_name], border = NA)

    # 2. 标签 (弯曲排列) —— 仅对长度 ≥ LABEL_MIN_SIZE 的染色体绘制
    if (!is.na(chr_size) && chr_size >= LABEL_MIN_SIZE) {
      circos.text(mean(xlim), 1.8, chr_label,
                  facing = "bending.outside",
                  niceFacing = TRUE, adj = c(0.5, 0), cex = LABEL_SIZE)
    }

    # 3. 刻度尺 (Integrated / Clean Cut)
    if (SHOW_TICKS) {
      # 刻度画在色块底边 (y=0)
      tick_y_pos <- 0
      tick_step  <- AXIS_GAP * 1000000

      # [自然截断] 计算最大整数刻度
      max_tick_pos <- floor(xlim[2] / tick_step) * tick_step

      if (max_tick_pos > 0) {
        # A. 底线 (Bone): 只画到整数截断点
        # 黑色线条叠加在色块下边缘，形成清晰边界
        circos.lines(c(0, max_tick_pos), c(tick_y_pos, tick_y_pos), col = "black", lwd = 0.5)

        # B. 刻度 (Ticks): 标准整数序列
        major_at = seq(0, max_tick_pos, by = tick_step)

        circos.axis(h = tick_y_pos, major.at = major_at, labels = NULL,
                    major.tick.length = 0.5, lwd = 0.5)
      }
    }
  }, bg.border = NA, track.height = 0.05)

  # === 核心: 绘制绸带 ===
  ribbon_colors <- grid_col[df_ribbon$chr1]

  # [修正点] 使用 genomicLink 绘制宽带 (当输入为区间时，genomicLink = Ribbon)
  circos.genomicLink(
    df_ribbon[, 1:3],
    df_ribbon[, 4:6],
    col = adjustcolor(ribbon_colors, alpha.f = LINK_ALPHA),
    border = NA
  )

  # === 装饰 ===
  if (SHOW_LEGEND) {
    legend_texts <- c()
    for (abbr in names(legend_map)) {
      legend_texts <- c(legend_texts, paste0(abbr, " = ", legend_map[[abbr]]))
    }
    legend("topleft", legend = legend_texts,
           bty = "n", cex = 1.0, text.font = 3,
           title = "Species Legend", title.adj = 0)
  }

  if (SHOW_TITLE) {
    title_text <- paste(unique(legend_map), collapse = " vs ")
    title(title_text, cex.main = 1.2)
  }
}

# --- 6. 输出 ---
prefix <- if(!is.null(CUSTOM_FILENAME)) CUSTOM_FILENAME else paste0("synteny_", COLOR_STYLE)
out_pdf <- file.path(DIR_OUTPUT_SUB, paste0(prefix, ".pdf"))
out_png <- file.path(DIR_OUTPUT_SUB, paste0(prefix, ".png"))

cat(paste0("   [Output] Generating PDF: ", out_pdf, "\n"))
pdf(out_pdf, width = 10, height = 10)
draw_core_plot()
dev.off()

cat(paste0("   [Output] Generating PNG: ", out_png, "\n"))
png(out_png, width = 3000, height = 3000, res = 300)
draw_core_plot()
dev.off()

cat(">> [Success] 绘图完成！\n")

