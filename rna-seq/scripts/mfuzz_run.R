#!/usr/bin/env Rscript

###############################################################################
# mfuzz_run.R
# 功能：
#   1. 读取 gene-level TPM 矩阵和 samples.tsv
#   2. 按 samples.tsv 中 group 列的首次出现顺序，自动构建阶段顺序
#   3. 对同组样本取均值，生成 stage-level TPM 矩阵
#   4. 进行低表达过滤和高变基因筛选
#   5. 执行 Mfuzz 聚类，扫描多个候选 c 值，并输出最终聚类结果
#   6. 导出 cluster membership、各簇基因列表、标准化矩阵、聚类图等
#   7. 导出 mfuzz universe，供后续现有 07/08 富集流程复用
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(Mfuzz)
  library(Biobase)
  library(ggplot2)
  library(grid)
})

###############################################################################
# 一、用户参数区（请按需修改）
###############################################################################

TPM_FILE <- "results/05_matrix/tpms/gene_tpm.tsv"
# 中文说明：输入的 gene-level TPM 矩阵文件
# 要求：第一列为基因 ID，后续各列为样本名

SAMPLES_FILE <- "data/samples.tsv"
# 中文说明：样本信息表
# 要求：至少包含 sample 和 group 两列
# 注意：group 列在文件中的首次出现顺序，将被视为发育时期顺序

OUTDIR <- "results/11_mfuzz"
# 中文说明：Mfuzz 主输出目录

GENEID_COL <- "gene_id"
# 中文说明：表达矩阵中的基因 ID 列名

SAMPLE_COL <- "sample"
# 中文说明：samples.tsv 中的样本名列名

GROUP_COL <- "group"
# 中文说明：samples.tsv 中的分组列名

MIN_TPM <- 1
# 中文说明：阶段平均表达过滤阈值

MIN_STAGE_COUNT <- 2
# 中文说明：至少在多少个阶段达到表达阈值，基因才被保留

VAR_METHOD <- "mad"
# 中文说明：高变基因筛选方法，可选 "mad" 或 "sd"

KEEP_TOP_N <- 4000
# 中文说明：按变化度排序后，保留前多少个高变基因

LOG_PSEUDOCOUNT <- 1
# 中文说明：log2 转换时加的伪计数，即 log2(x + LOG_PSEUDOCOUNT)

CANDIDATE_C <- c(4, 5, 6, 7, 8, 9, 10, 11, 12)
# 中文说明：候选聚类数列表

FINAL_C <- 8
# 中文说明：最终采用的聚类数
# 若手动指定非 NA，则优先使用该值
# 若设为 NA，则脚本会按 Mfuzz 官方路线，
# 结合 Dmin 与 cselection 结果自动推荐 C 值

DMIN_REPEATS <- 5
# 中文说明：Dmin 重复聚类次数
# 建议保持较小整数即可，常用 3

CSELECTION_REPEATS <- 4
# 中文说明：cselection 重复聚类次数
# 这里默认设为 4，便于和候选 C 个数区分，减少结果矩阵方向歧义

CSELECTION_MAX_MEAN_EMPTY <- 0
# 中文说明：自动选 C 时允许的平均空簇数上限
# 官方思路通常优先避免出现空簇，因此默认设为 0

MEMBERSHIP_CORE <- 0.7
# 中文说明：核心基因 membership 阈值

MEMBERSHIP_STRICT <- 0.8
# 中文说明：严格核心基因 membership 阈值

PLOT_WIDTH <- 10
# 中文说明：PDF/PNG 图宽（英寸）

PLOT_HEIGHT <- 7
# 中文说明：PDF/PNG 图高（英寸）

SAVE_PDF <- TRUE
# 中文说明：是否保存 PDF 图

SAVE_PNG <- TRUE
# 中文说明：是否保存 PNG 图

PNG_DPI <- 300
# 中文说明：PNG 输出分辨率

FONT_FAMILY_BASE <- "sans"
# 中文说明：图形基础字体
# 可填常见字体名称，如 Arial、Liberation Sans、DejaVu Sans、sans

FONT_FAMILY_TITLE <- ""
# 中文说明：标题字体；留空则自动继承 FONT_FAMILY_BASE

FONT_FAMILY_AXIS <- ""
# 中文说明：坐标轴字体；留空则自动继承 FONT_FAMILY_BASE

FONT_FAMILY_STRIP <- ""
# 中文说明：分面/簇标题字体；留空则自动继承 FONT_FAMILY_BASE

ENABLE_SHOWTEXT <- FALSE
# 中文说明：是否启用 showtext 接管字体渲染
# 若 Linux + PDF 字体兼容性较差，可尝试改为 TRUE

SHOW_RAW_TRAJECTORIES_IN_PUB <- TRUE
# 中文说明：发表版主图是否叠加少量原始基因轨迹
# 顶刊主图通常建议 FALSE，仅展示均值趋势 + 波动带即可

RAW_TRAJECTORY_MAX_PER_CLUSTER <- 100
# 中文说明：若 SHOW_RAW_TRAJECTORIES_IN_PUB=TRUE，每个 cluster 最多抽样显示多少条原始轨迹

RIBBON_STYLE <- "IQR"
# 中文说明：波动带类型
# 可选：IQR / SD / CI95

RIBBON_ALPHA <- 0.28
# 中文说明：波动带透明度

MEAN_LINE_WIDTH <- 1.0
# 中文说明：主趋势线线宽

RAW_LINE_WIDTH <- 0.25
# 中文说明：原始轨迹线宽（仅在 SHOW_RAW_TRAJECTORIES_IN_PUB=TRUE 时使用）

RAW_LINE_ALPHA <- 0.10
# 中文说明：原始轨迹透明度（仅在 SHOW_RAW_TRAJECTORIES_IN_PUB=TRUE 时使用）

PUB_FILL_COLORS <- c("#F6D6D3", "#D7E7F8", "#F8E5C5", "#D8EEE8", "#EADAF5", "#E2EACC", "#F8D9EA", "#DEE3FB")
# 中文说明：发表版主图各簇波动带填充色，低饱和配色，自动循环使用

PUB_LINE_COLORS <- c("#B55A52", "#4F7FAE", "#AF7A22", "#3E857B", "#8A61B4", "#6D8E4B", "#B16895", "#6776B8")
# 中文说明：发表版主图各簇主趋势线颜色，自动循环使用

RANDOM_SEED <- 12345
# 中文说明：随机种子

VERBOSE <- TRUE
# 中文说明：是否输出详细日志

###############################################################################
# 二、基础函数区
###############################################################################

log_msg <- function(...) {
  if (isTRUE(VERBOSE)) {
    cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
    cat("\n")
  }
}

safe_dir_create <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}


write_tsv <- function(df, file) {
  fwrite(df, file = file, sep = "	", quote = FALSE, na = "NA")
}

resolve_font_family <- function(preferred = "", fallback = c("sans", "Liberation Sans", "DejaVu Sans")) {
  candidates <- unique(c(preferred, fallback))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  if (requireNamespace("systemfonts", quietly = TRUE)) {
    font_tbl <- tryCatch(systemfonts::system_fonts(), error = function(e) NULL)
    if (!is.null(font_tbl) && "family" %in% colnames(font_tbl)) {
      installed <- unique(font_tbl$family)
      for (fam in candidates) {
        if (fam %in% installed) {
          return(fam)
        }
      }
    }
  }

  if ("sans" %in% candidates) {
    return("sans")
  }
  if (length(candidates) > 0) {
    return(candidates[1])
  }
  "sans"
}

open_graphics_device <- function(file, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(file))

  if (ext == "pdf") {
    grDevices::cairo_pdf(
      filename = file,
      width = width,
      height = height,
    )
    return(invisible(NULL))
  }

  if (ext == "png") {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(
        filename = file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        background = "white"
      )
    } else {
      grDevices::png(
        filename = file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        type = "cairo",
        bg = "white"
      )
    }
    return(invisible(NULL))
  }

  stop("暂不支持的输出格式：", ext)
}

init_plot_text_system <- function() {
  if (isTRUE(ENABLE_SHOWTEXT) && requireNamespace("showtext", quietly = TRUE)) {
    showtext::showtext_auto(enable = TRUE)
    showtext::showtext_opts(dpi = PNG_DPI)
  }
}

save_plot_base <- function(plot_fun, out_prefix, width = 10, height = 7, dpi = 300) {
  if (isTRUE(SAVE_PDF)) {
    pdf_file <- paste0(out_prefix, ".pdf")
    open_graphics_device(pdf_file, width = width, height = height, dpi = dpi)
    plot_fun()
    grDevices::dev.off()
  }
  if (isTRUE(SAVE_PNG)) {
    png_file <- paste0(out_prefix, ".png")
    open_graphics_device(png_file, width = width, height = height, dpi = dpi)
    plot_fun()
    grDevices::dev.off()
  }
}


calc_variability <- function(mat, method = "mad") {
  if (method == "mad") {
    return(apply(mat, 1, stats::mad, na.rm = TRUE))
  } else if (method == "sd") {
    return(apply(mat, 1, stats::sd, na.rm = TRUE))
  } else {
    stop("VAR_METHOD 仅支持 'mad' 或 'sd'")
  }
}

theme_pub <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size * 1.15,
        hjust = 0.5,
        colour = "#1F1F1F"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size * 0.92,
        hjust = 0.5,
        colour = "#4D4D4D"
      ),
      axis.title = ggplot2::element_text(
        size = base_size,
        face = "plain",
        colour = "#1F1F1F"
      ),
      axis.text = ggplot2::element_text(
        size = base_size * 0.92,
        colour = "#2F2F2F"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_size,
        colour = "#1F1F1F"
      ),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#333333"),
      legend.position = "none",
      plot.margin = ggplot2::margin(6, 8, 6, 8)
    )
}

get_palette_values <- function(n, fill_values, line_values) {
  list(
    fill = rep(fill_values, length.out = n),
    line = rep(line_values, length.out = n)
  )
}

calc_cluster_band_stats <- function(sub_mat, ribbon_style = "IQR") {
  stage_mean <- colMeans(sub_mat, na.rm = TRUE)

  if (toupper(ribbon_style) == "IQR") {
    lower <- apply(sub_mat, 2, stats::quantile, probs = 0.25, na.rm = TRUE)
    upper <- apply(sub_mat, 2, stats::quantile, probs = 0.75, na.rm = TRUE)
  } else if (toupper(ribbon_style) == "SD") {
    stage_sd <- apply(sub_mat, 2, stats::sd, na.rm = TRUE)
    lower <- stage_mean - stage_sd
    upper <- stage_mean + stage_sd
  } else if (toupper(ribbon_style) == "CI95") {
    stage_sd <- apply(sub_mat, 2, stats::sd, na.rm = TRUE)
    stage_n <- nrow(sub_mat)
    stage_se <- stage_sd / sqrt(max(stage_n, 1))
    lower <- stage_mean - 1.96 * stage_se
    upper <- stage_mean + 1.96 * stage_se
  } else {
    stop("RIBBON_STYLE 仅支持 IQR / SD / CI95")
  }

  list(
    mean = as.numeric(stage_mean),
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}

build_cluster_profile_summary_pub <- function(std_mat, cluster_assign, cluster_summary_df, stage_names) {
  clusters <- sort(unique(cluster_assign))
  out_list <- vector("list", length(clusters))

  for (i in seq_along(clusters)) {
    k <- clusters[i]
    idx <- which(cluster_assign == k)
    sub_mat <- std_mat[idx, , drop = FALSE]
    band_stats <- calc_cluster_band_stats(sub_mat, ribbon_style = RIBBON_STYLE)
    peak_stage <- ""
    if ("stage_peak" %in% colnames(cluster_summary_df)) {
      peak_stage <- cluster_summary_df$stage_peak[match(k, cluster_summary_df$cluster)]
    }

    out_list[[i]] <- data.frame(
      cluster = k,
      stage = stage_names,
      stage_index = seq_along(stage_names),
      mean = band_stats$mean,
      lower = band_stats$lower,
      upper = band_stats$upper,
      n_cluster = nrow(sub_mat),
      peak_stage = peak_stage,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

build_raw_long_df <- function(sub_mat, stage_names, max_n = 80) {
  if (nrow(sub_mat) == 0) {
    return(NULL)
  }

  keep_idx <- seq_len(nrow(sub_mat))
  if (nrow(sub_mat) > max_n) {
    keep_idx <- sort(sample.int(nrow(sub_mat), size = max_n, replace = FALSE))
  }

  sub_use <- sub_mat[keep_idx, , drop = FALSE]
  data.frame(
    gene_id = rep(rownames(sub_use), each = length(stage_names)),
    stage = rep(stage_names, times = nrow(sub_use)),
    stage_index = rep(seq_along(stage_names), times = nrow(sub_use)),
    value = as.vector(t(sub_use)),
    stringsAsFactors = FALSE
  )
}

make_single_cluster_panel <- function(cluster_id,
                                      summary_df,
                                      std_mat,
                                      cluster_assign,
                                      palette,
                                      y_limits) {
  one_df <- summary_df[summary_df$cluster == cluster_id, , drop = FALSE]
  sub_mat <- std_mat[cluster_assign == cluster_id, , drop = FALSE]

  p <- ggplot() +
    geom_ribbon(
      data = one_df,
      aes(x = stage_index, ymin = lower, ymax = upper),
      fill = palette$fill[cluster_id],
      alpha = RIBBON_ALPHA
    )

  if (isTRUE(SHOW_RAW_TRAJECTORIES_IN_PUB)) {
    raw_df <- build_raw_long_df(sub_mat, stage_names = one_df$stage, max_n = RAW_TRAJECTORY_MAX_PER_CLUSTER)
    if (!is.null(raw_df) && nrow(raw_df) > 0) {
      p <- p + geom_line(
        data = raw_df,
        aes(x = stage_index, y = value, group = gene_id),
        colour = palette$line[cluster_id],
        alpha = RAW_LINE_ALPHA,
        linewidth = RAW_LINE_WIDTH,
        lineend = "round"
      )
    }
  }

  p <- p +
    geom_line(
      data = one_df,
      aes(x = stage_index, y = mean, group = 1),
      colour = palette$line[cluster_id],
      linewidth = MEAN_LINE_WIDTH,
      lineend = "round"
    ) +
    geom_point(
      data = one_df,
      aes(x = stage_index, y = mean),
      colour = palette$line[cluster_id],
      fill = "white",
      shape = 21,
      stroke = 0.35,
      size = 1.8
    ) +
    scale_x_continuous(
      breaks = one_df$stage_index,
      labels = one_df$stage,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = sprintf("Cluster %d", cluster_id),
      subtitle = sprintf("n = %d | peak = %s", one_df$n_cluster[1], one_df$peak_stage[1]),
      x = "Time",
      y = "Standardized expression (z-score)"
    ) +
    theme_pub(base_size = 11.5) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, vjust = 0.8),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8))
    )

  p
}

plot_cluster_profiles_pub <- function(std_mat, cluster_assign, cluster_summary_df, outdir, prefix = "cluster") {
  stage_names <- colnames(std_mat)
  summary_df <- build_cluster_profile_summary_pub(
    std_mat = std_mat,
    cluster_assign = cluster_assign,
    cluster_summary_df = cluster_summary_df,
    stage_names = stage_names
  )

  palette <- get_palette_values(
    n = length(unique(cluster_assign)),
    fill_values = PUB_FILL_COLORS,
    line_values = PUB_LINE_COLORS
  )

  y_rng <- range(c(summary_df$lower, summary_df$upper, summary_df$mean), na.rm = TRUE)
  if (isTRUE(SHOW_RAW_TRAJECTORIES_IN_PUB)) {
    y_rng <- range(c(y_rng, std_mat), na.rm = TRUE)
  }
  y_pad <- diff(y_rng) * 0.08
  if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 0.5
  y_limits <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)

  for (k in sort(unique(cluster_assign))) {
    p <- make_single_cluster_panel(
      cluster_id = k,
      summary_df = summary_df,
      std_mat = std_mat,
      cluster_assign = cluster_assign,
      palette = palette,
      y_limits = y_limits
    )
    out_prefix <- file.path(outdir, sprintf("%s%02d_profile", prefix, k))
    save_plot_base(function() print(p), out_prefix, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PNG_DPI)
  }

  write_tsv(summary_df, file.path(outdir, "cluster_profile_summary.tsv"))
}

plot_cluster_overview_pub <- function(std_mat,
                                      cluster_assign,
                                      cluster_summary_df,
                                      out_prefix,
                                      main_title = "Mfuzz final overview") {
  stage_names <- colnames(std_mat)
  cluster_ids <- sort(unique(cluster_assign))
  summary_df <- build_cluster_profile_summary_pub(
    std_mat = std_mat,
    cluster_assign = cluster_assign,
    cluster_summary_df = cluster_summary_df,
    stage_names = stage_names
  )

  palette <- get_palette_values(
    n = length(cluster_ids),
    fill_values = PUB_FILL_COLORS,
    line_values = PUB_LINE_COLORS
  )

  y_rng <- range(c(summary_df$lower, summary_df$upper, summary_df$mean), na.rm = TRUE)
  if (isTRUE(SHOW_RAW_TRAJECTORIES_IN_PUB)) {
    y_rng <- range(c(y_rng, std_mat), na.rm = TRUE)
  }
  y_pad <- diff(y_rng) * 0.08
  if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 0.5
  y_limits <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)

  grob_list <- lapply(cluster_ids, function(k) {
    p <- make_single_cluster_panel(
      cluster_id = k,
      summary_df = summary_df,
      std_mat = std_mat,
      cluster_assign = cluster_assign,
      palette = palette,
      y_limits = y_limits
    )
    ggplotGrob(p)
  })

  ncol_plot <- 2
  nrow_plot <- ceiling(length(grob_list) / ncol_plot)

  plot_fun <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(
      nrow = nrow_plot + 1,
      ncol = ncol_plot,
      heights = unit.c(unit(1.3, "lines"), rep(unit(1, "null"), nrow_plot))
    )
    grid::pushViewport(grid::viewport(layout = lay))

    grid::grid.text(
      label = main_title,
      x = 0.5,
      y = 0.5,
      gp = grid::gpar(
        fontsize = 15,
        fontface = "bold",
        col = "#1F1F1F"
      ),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:ncol_plot)
    )

    for (i in seq_along(grob_list)) {
      rr <- ((i - 1) %/% ncol_plot) + 2
      cc <- ((i - 1) %% ncol_plot) + 1
      grid::pushViewport(grid::viewport(layout.pos.row = rr, layout.pos.col = cc))
      grid::grid.draw(grob_list[[i]])
      grid::upViewport()
    }

    grid::upViewport()
  }

  save_plot_base(
    plot_fun = plot_fun,
    out_prefix = out_prefix,
    width = max(10, PLOT_WIDTH + 1),
    height = max(PLOT_HEIGHT + 2, 3.1 * nrow_plot + 0.8),
    dpi = PNG_DPI
  )

  write_tsv(summary_df, paste0(out_prefix, "_summary.tsv"))
}

plot_heatmap_pub <- function(std_mat, cluster_assign, out_prefix) {
  ord <- order(cluster_assign)
  mat2 <- std_mat[ord, , drop = FALSE]
  cluster_ord <- cluster_assign[ord]

  col_fun <- grDevices::colorRampPalette(c("#3C5488", "#F7F7F7", "#E64B35"))(101)
  zlim <- max(abs(mat2), na.rm = TRUE)
  brks <- seq(-zlim, zlim, length.out = length(col_fun) + 1)

  plot_fun <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    layout(matrix(c(1, 2), nrow = 1), widths = c(6.5, 0.5))
    graphics::par(mar = c(5.2, 5.0, 2.5, 0.8))
    image(
      x = seq_len(ncol(mat2)),
      y = seq_len(nrow(mat2)),
      z = t(mat2[nrow(mat2):1, , drop = FALSE]),
      col = col_fun,
      breaks = brks,
      axes = FALSE,
      xlab = "Stage",
      ylab = "Genes",
      main = "Mfuzz standardized matrix"
    )
    axis(1, at = seq_len(ncol(mat2)), labels = colnames(mat2), las = 1)
    box()
    sep_pos <- which(diff(cluster_ord) != 0)
    if (length(sep_pos) > 0) {
      for (sp in sep_pos) {
        graphics::abline(h = nrow(mat2) - sp + 0.5, col = "#FFFFFF", lwd = 0.8)
      }
    }

    graphics::par(mar = c(5.2, 0.4, 2.5, 3.0))
    image(
      x = 1,
      y = seq(-zlim, zlim, length.out = 101),
      z = matrix(seq(-zlim, zlim, length.out = 101), nrow = 1),
      col = col_fun,
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    axis(4, at = pretty(seq(-zlim, zlim, length.out = 101)))
    mtext("Z-score", side = 4, line = 2.0)
    box()
  }

  save_plot_base(
    plot_fun = plot_fun,
    out_prefix = out_prefix,
    width = max(8, PLOT_WIDTH),
    height = max(10, PLOT_HEIGHT + 2),
    dpi = PNG_DPI
  )
}

get_offdiag_values <- function(mat) {

  if (method == "mad") {
    return(apply(mat, 1, stats::mad, na.rm = TRUE))
  } else if (method == "sd") {
    return(apply(mat, 1, stats::sd, na.rm = TRUE))
  } else {
    stop("VAR_METHOD 仅支持 'mad' 或 'sd'")
  }
}

get_stage_order_from_samples <- function(samples_df, group_col) {
  unique(samples_df[[group_col]])
}

build_stage_mean_matrix <- function(expr_dt, samples_df, geneid_col, sample_col, group_col) {
  sample_to_group <- samples_df[, c(sample_col, group_col), with = FALSE]
  stage_order <- get_stage_order_from_samples(samples_df, group_col)

  expr_mat <- as.data.frame(expr_dt)
  rownames(expr_mat) <- expr_mat[[geneid_col]]
  expr_mat[[geneid_col]] <- NULL

  stage_mean_list <- list()
  for (g in stage_order) {
    ss <- sample_to_group[get(group_col) == g, get(sample_col)]
    ss <- ss[ss %in% colnames(expr_mat)]
    if (length(ss) == 0) {
      stop(sprintf("group '%s' 在表达矩阵中没有对应样本列", g))
    }
    stage_mean_list[[g]] <- rowMeans(expr_mat[, ss, drop = FALSE], na.rm = TRUE)
  }

  stage_mean_mat <- do.call(cbind, stage_mean_list)
  stage_mean_df <- data.frame(
    gene_id = rownames(stage_mean_mat),
    stage_mean_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(stage_mean_df)[1] <- geneid_col
  stage_mean_df
}

make_eset_from_matrix <- function(mat, sample_names) {
  eset <- ExpressionSet(assayData = as.matrix(mat))
  pData(eset) <- data.frame(
    sample = sample_names,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )
  eset
}

extract_cluster_membership <- function(cl_obj, gene_ids) {
  memb_mat <- as.data.frame(cl_obj$membership)
  colnames(memb_mat) <- paste0("cluster_", seq_len(ncol(memb_mat)))

  max_cluster_idx <- apply(memb_mat, 1, which.max)
  max_membership <- apply(memb_mat, 1, max)

  out <- data.frame(
    gene_id = gene_ids,
    cluster = max_cluster_idx,
    membership = max_membership,
    is_core = max_membership >= MEMBERSHIP_CORE,
    is_strict_core = max_membership >= MEMBERSHIP_STRICT,
    stringsAsFactors = FALSE
  )

  list(
    membership_table = out,
    membership_matrix = data.frame(gene_id = gene_ids, memb_mat, check.names = FALSE)
  )
}

get_cluster_profile_summary <- function(std_mat, cluster_assign, stage_names) {
  clusters <- sort(unique(cluster_assign))
  res <- lapply(clusters, function(k) {
    idx <- which(cluster_assign == k)
    sub <- std_mat[idx, , drop = FALSE]
    mean_profile <- colMeans(sub, na.rm = TRUE)
    peak_stage <- stage_names[which.max(mean_profile)]
    trough_stage <- stage_names[which.min(mean_profile)]
    data.frame(
      cluster = k,
      stage_peak = peak_stage,
      stage_trough = trough_stage,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

plot_cluster_profiles_gg <- function(std_mat, cluster_assign, outdir, prefix = "cluster") {
  stage_names <- colnames(std_mat)
  clusters <- sort(unique(cluster_assign))

  for (k in clusters) {
    idx <- which(cluster_assign == k)
    sub <- std_mat[idx, , drop = FALSE]

    long_df <- data.frame(
      gene_id = rep(rownames(sub), times = ncol(sub)),
      stage = rep(stage_names, each = nrow(sub)),
      value = as.vector(t(sub)),
      stringsAsFactors = FALSE
    )

    mean_df <- data.frame(
      stage = stage_names,
      value = colMeans(sub, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    p <- ggplot() +
      geom_line(
        data = long_df,
        aes(x = factor(stage, levels = stage_names), y = value, group = gene_id),
        alpha = 0.15,
        linewidth = 0.3
      ) +
      geom_line(
        data = mean_df,
        aes(x = factor(stage, levels = stage_names), y = value, group = 1),
        linewidth = 1.2
      ) +
      geom_point(
        data = mean_df,
        aes(x = factor(stage, levels = stage_names), y = value),
        size = 2
      ) +
      labs(
        title = sprintf("Cluster %02d profile", k),
        x = "Stage",
        y = "Standardized expression"
      ) +
      theme_classic(base_size = 12)

    out_prefix <- file.path(outdir, sprintf("%s%02d_profile", prefix, k))
    save_plot_base(function() print(p), out_prefix, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PNG_DPI)
  }
}

plot_heatmap_base <- function(std_mat, cluster_assign, out_prefix) {
  ord <- order(cluster_assign)
  mat2 <- std_mat[ord, , drop = FALSE]

  plot_fun <- function() {
    graphics::par(mar = c(6, 6, 3, 2))
    graphics::image(
      x = seq_len(ncol(mat2)),
      y = seq_len(nrow(mat2)),
      z = t(mat2[nrow(mat2):1, , drop = FALSE]),
      axes = FALSE,
      xlab = "Stage",
      ylab = "Gene",
      main = "Mfuzz standardized matrix"
    )
    graphics::axis(1, at = seq_len(ncol(mat2)), labels = colnames(mat2), las = 2)
    graphics::box()
  }

  save_plot_base(plot_fun, out_prefix, width = max(8, PLOT_WIDTH), height = max(10, PLOT_HEIGHT + 2), dpi = PNG_DPI)
}

plot_mfuzz_overview <- function(eset_std, cl_obj, out_prefix, main_title = "Mfuzz clustering overview") {
  ncl <- ncol(cl_obj$membership)
  ncol_plot <- 2
  nrow_plot <- ceiling(ncl / ncol_plot)

  plot_fun <- function() {
    Mfuzz::mfuzz.plot(
      eset_std,
      cl = cl_obj,
      mfrow = c(nrow_plot, ncol_plot),
      time.labels = colnames(exprs(eset_std)),
      new.window = FALSE
    )
    graphics::mtext(main_title, side = 3, line = -2, outer = FALSE, cex = 1.2)
  }

  save_plot_base(
    plot_fun,
    out_prefix,
    width = max(10, PLOT_WIDTH),
    height = max(3 * nrow_plot, PLOT_HEIGHT + 1),
    dpi = PNG_DPI
  )
}

get_offdiag_values <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(numeric(0))
  }
  mat[lower.tri(mat, diag = TRUE)] <- NA_real_
  vals <- as.numeric(mat)
  vals[is.finite(vals)]
}

summarize_overlap <- function(cl_obj) {
  over_mat <- Mfuzz::overlap(cl_obj)
  over_vals <- get_offdiag_values(over_mat)

  data.frame(
    overlap_mean = ifelse(length(over_vals) == 0, NA_real_, mean(over_vals, na.rm = TRUE)),
    overlap_median = ifelse(length(over_vals) == 0, NA_real_, stats::median(over_vals, na.rm = TRUE)),
    overlap_max = ifelse(length(over_vals) == 0, NA_real_, max(over_vals, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

convert_cselection_to_summary <- function(csel_raw, crange) {
  cs <- as.matrix(csel_raw)
  crange_chr <- as.character(crange)

  use_by_col <- FALSE
  use_by_row <- FALSE

  if (!is.null(colnames(cs)) && all(crange_chr %in% colnames(cs))) {
    use_by_col <- TRUE
    cs <- cs[, crange_chr, drop = FALSE]
  } else if (!is.null(rownames(cs)) && all(crange_chr %in% rownames(cs))) {
    use_by_row <- TRUE
    cs <- cs[crange_chr, , drop = FALSE]
  } else if (ncol(cs) == length(crange) && nrow(cs) != length(crange)) {
    use_by_col <- TRUE
  } else if (nrow(cs) == length(crange) && ncol(cs) != length(crange)) {
    use_by_row <- TRUE
  } else if (ncol(cs) == length(crange)) {
    use_by_col <- TRUE
  } else if (nrow(cs) == length(crange)) {
    use_by_row <- TRUE
  } else {
    stop("无法解析 cselection 返回结果的维度，请检查候选 C 与 repeats 设置")
  }

  if (use_by_col) {
    mean_empty <- colMeans(cs, na.rm = TRUE)
    max_empty <- apply(cs, 2, max, na.rm = TRUE)
  } else if (use_by_row) {
    mean_empty <- rowMeans(cs, na.rm = TRUE)
    max_empty <- apply(cs, 1, max, na.rm = TRUE)
  }

  data.frame(
    c_value = crange,
    cselection_mean_empty = as.numeric(mean_empty),
    cselection_max_empty = as.numeric(max_empty),
    stringsAsFactors = FALSE
  )
}

detect_dmin_elbow <- function(crange, dmin_values) {
  crange <- as.numeric(crange)
  dmin_values <- as.numeric(dmin_values)

  ok <- is.finite(crange) & is.finite(dmin_values)
  crange <- crange[ok]
  dmin_values <- dmin_values[ok]

  if (length(crange) == 0) {
    stop("Dmin 结果为空，无法自动选择 C")
  }
  if (length(crange) <= 2) {
    return(crange[1])
  }

  x <- (crange - min(crange)) / (max(crange) - min(crange))
  y <- (dmin_values - min(dmin_values)) / (max(dmin_values) - min(dmin_values))

  if (all(!is.finite(x)) || all(!is.finite(y)) || diff(range(y, na.rm = TRUE)) == 0) {
    return(crange[ceiling(length(crange) / 2)])
  }

  x1 <- x[1]
  y1 <- y[1]
  x2 <- x[length(x)]
  y2 <- y[length(y)]

  numerator <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
  denominator <- sqrt((y2 - y1)^2 + (x2 - x1)^2)

  if (!is.finite(denominator) || denominator == 0) {
    return(crange[ceiling(length(crange) / 2)])
  }

  dist_to_line <- numerator / denominator
  dist_to_line[c(1, length(dist_to_line))] <- -Inf

  crange[which.max(dist_to_line)]
}


plot_dmin_diagnostic <- function(dmin_df, chosen_c, out_prefix) {
  df <- dmin_df[order(dmin_df$c_value), , drop = FALSE]
  p <- ggplot(df, aes(x = c_value, y = dmin_value)) +
    geom_line(linewidth = 0.8, colour = "#4F7FAE", lineend = "round") +
    geom_point(size = 2.2, colour = "#4F7FAE") +
    geom_vline(xintercept = chosen_c, linetype = 2, linewidth = 0.45, colour = "#7A7A7A") +
    geom_point(
      data = df[df$c_value == chosen_c, , drop = FALSE],
      aes(x = c_value, y = dmin_value),
      size = 3.0,
      colour = "#B55A52"
    ) +
    labs(
      title = "Dmin diagnostic plot",
      subtitle = sprintf("Recommended C = %d", chosen_c),
      x = "Candidate C",
      y = "Average minimum centroid distance"
    ) +
    theme_pub(base_size = 12)

  save_plot_base(function() print(p), out_prefix, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PNG_DPI)
}

plot_cselection_diagnostic <- function(csel_df, chosen_c, out_prefix) {
  df <- csel_df[order(csel_df$c_value), , drop = FALSE]
  p <- ggplot(df, aes(x = c_value, y = cselection_mean_empty)) +
    geom_col(width = 0.62, fill = "#D7E7F8", colour = "#4F7FAE", linewidth = 0.35) +
    geom_vline(xintercept = chosen_c, linetype = 2, linewidth = 0.45, colour = "#7A7A7A") +
    geom_point(
      data = df[df$c_value == chosen_c, , drop = FALSE],
      aes(x = c_value, y = cselection_mean_empty),
      size = 2.8,
      colour = "#B55A52"
    ) +
    labs(
      title = "cselection diagnostic plot",
      subtitle = sprintf("Recommended C = %d", chosen_c),
      x = "Candidate C",
      y = "Mean empty clusters"
    ) +
    theme_pub(base_size = 12)

  save_plot_base(function() print(p), out_prefix, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PNG_DPI)
}

plot_cluster_count_barplot <- function(cluster_summary_df, out_prefix) {
  plot_df <- rbind(
    data.frame(
      cluster = factor(cluster_summary_df$cluster, levels = cluster_summary_df$cluster),
      type = "All genes",
      count = cluster_summary_df$n_all,
      stringsAsFactors = FALSE
    ),
    data.frame(
      cluster = factor(cluster_summary_df$cluster, levels = cluster_summary_df$cluster),
      type = "Core genes",
      count = cluster_summary_df$n_core_ge_0.7,
      stringsAsFactors = FALSE
    )
  )

  dodge_pos <- position_dodge(width = 0.72)

  p <- ggplot(plot_df, aes(x = cluster, y = count, fill = type)) +
    geom_col(position = dodge_pos, width = 0.64, colour = "#3A3A3A", linewidth = 0.28) +
    geom_text(
      aes(label = count),
      position = dodge_pos,
      vjust = -0.35,
      size = 3.6,
      colour = "#2B2B2B"
    ) +
    scale_fill_manual(values = c("All genes" = "#95C8F2", "Core genes" = "#F79D93")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      title = "Gene counts per cluster",
      x = "Cluster",
      y = "Number of genes"
    ) +
    theme_pub(base_size = 12)

  save_plot_base(function() print(p), out_prefix, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PNG_DPI)
}


###############################################################################
# 三、目录准备
###############################################################################

set.seed(RANDOM_SEED)

FONT_BASE_RESOLVED <- resolve_font_family(FONT_FAMILY_BASE)
FONT_TITLE_RESOLVED <- resolve_font_family(
  ifelse(nzchar(FONT_FAMILY_TITLE), FONT_FAMILY_TITLE, FONT_BASE_RESOLVED)
)
FONT_AXIS_RESOLVED <- resolve_font_family(
  ifelse(nzchar(FONT_FAMILY_AXIS), FONT_FAMILY_AXIS, FONT_BASE_RESOLVED)
)
FONT_STRIP_RESOLVED <- resolve_font_family(
  ifelse(nzchar(FONT_FAMILY_STRIP), FONT_FAMILY_STRIP, FONT_BASE_RESOLVED)
)

init_plot_text_system()

safe_dir_create(OUTDIR)
DIR_CSCAN <- file.path(OUTDIR, "06_c_scan")
DIR_FINAL <- file.path(OUTDIR, "07_final_clusters")
DIR_LISTS <- file.path(OUTDIR, "08_cluster_gene_lists")
DIR_PLOTS <- file.path(OUTDIR, "09_plots")
DIR_ENRICH <- file.path(OUTDIR, "10_for_enrichment")

safe_dir_create(DIR_CSCAN)
safe_dir_create(DIR_FINAL)
safe_dir_create(DIR_LISTS)
safe_dir_create(DIR_PLOTS)
safe_dir_create(DIR_ENRICH)

###############################################################################
# 四、模块 M1：读取输入并基础对齐
###############################################################################

log_msg("读取 TPM 矩阵：", TPM_FILE)
expr_dt <- fread(TPM_FILE, sep = "\t", header = TRUE, data.table = TRUE)

log_msg("读取样本表：", SAMPLES_FILE)
samples_dt <- fread(SAMPLES_FILE, sep = "\t", header = TRUE, data.table = TRUE)

if (!(GENEID_COL %in% colnames(expr_dt))) {
  stop(sprintf("表达矩阵中未找到基因 ID 列：%s", GENEID_COL))
}
if (!(SAMPLE_COL %in% colnames(samples_dt))) {
  stop(sprintf("samples.tsv 中未找到样本列：%s", SAMPLE_COL))
}
if (!(GROUP_COL %in% colnames(samples_dt))) {
  stop(sprintf("samples.tsv 中未找到分组列：%s", GROUP_COL))
}

samples_dt <- unique(samples_dt[, c(SAMPLE_COL, GROUP_COL), with = FALSE])

expr_sample_cols <- setdiff(colnames(expr_dt), GENEID_COL)
sample_vec <- samples_dt[[SAMPLE_COL]]
group_vec <- samples_dt[[GROUP_COL]]
stage_order <- unique(group_vec)

missing_samples <- setdiff(sample_vec, expr_sample_cols)
if (length(missing_samples) > 0) {
  stop(sprintf(
    "以下 samples.tsv 中的样本，在表达矩阵中找不到：%s",
    paste(missing_samples, collapse = ", ")
  ))
}

expr_dt <- expr_dt[, c(GENEID_COL, sample_vec), with = FALSE]

input_summary <- data.frame(
  item = c("n_genes_raw", "n_samples", "n_groups", "group_order"),
  value = c(
    nrow(expr_dt),
    length(sample_vec),
    length(stage_order),
    paste(stage_order, collapse = " -> ")
  ),
  stringsAsFactors = FALSE
)

group_count_df <- as.data.frame(table(group_vec), stringsAsFactors = FALSE)
colnames(group_count_df) <- c("group", "n_samples")

write_tsv(input_summary, file.path(OUTDIR, "00_input_summary.tsv"))
write_tsv(group_count_df, file.path(OUTDIR, "00_group_sample_counts.tsv"))

log_msg("输入读取完成：", nrow(expr_dt), " 个基因；", length(sample_vec), " 个样本；", length(stage_order), " 个 group")

###############################################################################
# 五、模块 M2：构建 stage-level TPM 矩阵
###############################################################################

log_msg("按 group 聚合同组样本，构建阶段平均 TPM 矩阵")
stage_mean_df <- build_stage_mean_matrix(
  expr_dt = expr_dt,
  samples_df = samples_dt,
  geneid_col = GENEID_COL,
  sample_col = SAMPLE_COL,
  group_col = GROUP_COL
)

write_tsv(stage_mean_df, file.path(OUTDIR, "01_stage_mean_tpm.tsv"))
log_msg("阶段平均 TPM 矩阵已输出：01_stage_mean_tpm.tsv")

###############################################################################
# 六、模块 M3：低表达过滤
###############################################################################

log_msg("进行低表达过滤")
stage_mat <- as.data.frame(stage_mean_df)
rownames(stage_mat) <- stage_mat[[GENEID_COL]]
stage_mat[[GENEID_COL]] <- NULL
stage_mat <- as.matrix(stage_mat)

expr_stage_count <- rowSums(stage_mat >= MIN_TPM, na.rm = TRUE)
keep_expr <- expr_stage_count >= MIN_STAGE_COUNT
stage_expr_filtered <- stage_mat[keep_expr, , drop = FALSE]

expr_stats <- data.frame(
  metric = c("n_input_genes", "n_kept_genes", "n_removed_genes", "MIN_TPM", "MIN_STAGE_COUNT"),
  value = c(
    nrow(stage_mat),
    nrow(stage_expr_filtered),
    nrow(stage_mat) - nrow(stage_expr_filtered),
    MIN_TPM,
    MIN_STAGE_COUNT
  ),
  stringsAsFactors = FALSE
)

stage_expr_filtered_df <- data.frame(
  gene_id = rownames(stage_expr_filtered),
  stage_expr_filtered,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(stage_expr_filtered_df)[1] <- GENEID_COL

write_tsv(stage_expr_filtered_df, file.path(OUTDIR, "02_stage_mean_tpm.filtered_expr.tsv"))
write_tsv(expr_stats, file.path(OUTDIR, "02_filter_expr_stats.tsv"))

log_msg("低表达过滤完成：保留 ", nrow(stage_expr_filtered), " 个基因")

###############################################################################
# 七、模块 M4：高变基因筛选
###############################################################################

log_msg("进行高变基因筛选，方法：", VAR_METHOD)
var_values <- calc_variability(stage_expr_filtered, method = VAR_METHOD)

var_df <- data.frame(
  gene_id = rownames(stage_expr_filtered),
  variability_value = as.numeric(var_values),
  variability_method = VAR_METHOD,
  stringsAsFactors = FALSE
)
colnames(var_df)[1] <- GENEID_COL

var_df <- var_df[order(var_df$variability_value, decreasing = TRUE), , drop = FALSE]
var_df$rank <- seq_len(nrow(var_df))

n_keep_var <- min(KEEP_TOP_N, nrow(var_df))
keep_genes_var <- var_df[[GENEID_COL]][seq_len(n_keep_var)]
stage_var_filtered <- stage_expr_filtered[keep_genes_var, , drop = FALSE]

var_stats <- data.frame(
  metric = c("n_input_genes", "n_kept_genes", "n_removed_genes", "VAR_METHOD", "KEEP_TOP_N"),
  value = c(
    nrow(stage_expr_filtered),
    nrow(stage_var_filtered),
    nrow(stage_expr_filtered) - nrow(stage_var_filtered),
    VAR_METHOD,
    KEEP_TOP_N
  ),
  stringsAsFactors = FALSE
)

stage_var_filtered_df <- data.frame(
  gene_id = rownames(stage_var_filtered),
  stage_var_filtered,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(stage_var_filtered_df)[1] <- GENEID_COL

write_tsv(stage_var_filtered_df, file.path(OUTDIR, "03_stage_mean_tpm.filtered_var.tsv"))
write_tsv(var_df, file.path(OUTDIR, "03_gene_variability.tsv"))
write_tsv(var_stats, file.path(OUTDIR, "03_filter_var_stats.tsv"))

log_msg("高变筛选完成：保留 ", nrow(stage_var_filtered), " 个基因")

###############################################################################
# 八、模块 M5：数据变换与 Mfuzz 输入准备
###############################################################################

log_msg("进行 log2 转换并标准化")
stage_log2 <- log2(stage_var_filtered + LOG_PSEUDOCOUNT)

stage_log2_df <- data.frame(
  gene_id = rownames(stage_log2),
  stage_log2,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(stage_log2_df)[1] <- GENEID_COL
write_tsv(stage_log2_df, file.path(OUTDIR, "04_stage_mean_tpm.log2.tsv"))

eset <- make_eset_from_matrix(stage_log2, sample_names = colnames(stage_log2))
eset_std <- Mfuzz::standardise(eset)
std_mat <- exprs(eset_std)

stage_std_df <- data.frame(
  gene_id = rownames(std_mat),
  std_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(stage_std_df)[1] <- GENEID_COL
write_tsv(stage_std_df, file.path(OUTDIR, "04_stage_mean_tpm.log2.standardised.tsv"))

log_msg("标准化完成")

###############################################################################
# 九、模块 M6：估计 m 值
###############################################################################

log_msg("估计 m 值")
m_est <- Mfuzz::mestimate(eset_std)

m_df <- data.frame(
  metric = c("m_estimate", "n_genes", "n_stages"),
  value = c(m_est, nrow(std_mat), ncol(std_mat)),
  stringsAsFactors = FALSE
)
write_tsv(m_df, file.path(OUTDIR, "05_mestimate.txt"))

log_msg("m 值估计完成：", round(m_est, 4))

###############################################################################
# 十、模块 M7：候选 c 值扫描
###############################################################################

log_msg("开始扫描候选 c 值：", paste(CANDIDATE_C, collapse = ", "))

# 中文说明：
# 这里改为采用 Mfuzz 官方路线：
#   1. Dmin：观察最小质心距离随 c 的变化
#   2. cselection：检测空簇出现情况
#   3. overlap：查看候选聚类下簇间重叠程度（作为官方辅助信息）
# 若 FINAL_C 为 NA，则自动在“无空簇或空簇可接受”的候选范围内，
# 使用 Dmin 曲线拐点来推荐最终 C。

dmin_values <- Mfuzz::Dmin(
  eset_std,
  m = m_est,
  crange = CANDIDATE_C,
  repeats = DMIN_REPEATS,
  visu = FALSE
)

dmin_df <- data.frame(
  c_value = CANDIDATE_C,
  dmin_value = as.numeric(dmin_values),
  stringsAsFactors = FALSE
)
write_tsv(dmin_df, file.path(DIR_CSCAN, "dmin_values.tsv"))

cselection_raw <- Mfuzz::cselection(
  eset_std,
  m = m_est,
  crange = CANDIDATE_C,
  repeats = CSELECTION_REPEATS,
  visu = FALSE
)

cselection_raw_df <- as.data.frame(as.matrix(cselection_raw), stringsAsFactors = FALSE)
write_tsv(cselection_raw_df, file.path(DIR_CSCAN, "cselection_raw.tsv"))

cselection_summary_df <- convert_cselection_to_summary(cselection_raw, CANDIDATE_C)
write_tsv(cselection_summary_df, file.path(DIR_CSCAN, "cselection_summary.tsv"))

cscan_summary_list <- list()
overlap_summary_list <- list()

for (cc in CANDIDATE_C) {
  log_msg("扫描 c = ", cc)
  cl_tmp <- Mfuzz::mfuzz(eset_std, c = cc, m = m_est)

  tmp_membership <- extract_cluster_membership(cl_tmp, gene_ids = rownames(std_mat))
  tmp_tab <- tmp_membership$membership_table

  size_df <- do.call(rbind, lapply(seq_len(cc), function(k) {
    sub <- tmp_tab[tmp_tab$cluster == k, , drop = FALSE]
    data.frame(
      cluster = k,
      n_all = nrow(sub),
      n_core_ge_0.7 = sum(sub$is_core, na.rm = TRUE),
      n_core_ge_0.8 = sum(sub$is_strict_core, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write_tsv(size_df, file.path(DIR_CSCAN, sprintf("cluster_size_c%02d.tsv", cc)))

  tmp_cluster_summary_df <- merge(
    size_df,
    get_cluster_profile_summary(
      std_mat = std_mat,
      cluster_assign = tmp_tab$cluster,
      stage_names = colnames(std_mat)
    ),
    by = "cluster",
    all.x = TRUE
  )

  out_prefix <- file.path(DIR_CSCAN, sprintf("mfuzz_c%02d", cc))
  plot_cluster_overview_pub(
    std_mat = std_mat,
    cluster_assign = tmp_tab$cluster,
    cluster_summary_df = tmp_cluster_summary_df,
    out_prefix = out_prefix,
    main_title = sprintf("Mfuzz overview (c=%d)", cc)
  )

  over_sum <- summarize_overlap(cl_tmp)
  overlap_summary_list[[as.character(cc)]] <- data.frame(
    c_value = cc,
    over_sum,
    stringsAsFactors = FALSE
  )

  cscan_summary_list[[as.character(cc)]] <- data.frame(
    c_value = cc,
    n_clusters = cc,
    min_cluster_size = min(size_df$n_all),
    max_cluster_size = max(size_df$n_all),
    median_cluster_size = stats::median(size_df$n_all),
    stringsAsFactors = FALSE
  )
}

cscan_summary_df <- do.call(rbind, cscan_summary_list)
overlap_summary_df <- do.call(rbind, overlap_summary_list)
write_tsv(overlap_summary_df, file.path(DIR_CSCAN, "overlap_summary.tsv"))

cscan_summary_df <- merge(cscan_summary_df, dmin_df, by = "c_value", all.x = TRUE)
cscan_summary_df <- merge(cscan_summary_df, cselection_summary_df, by = "c_value", all.x = TRUE)
cscan_summary_df <- merge(cscan_summary_df, overlap_summary_df, by = "c_value", all.x = TRUE)
cscan_summary_df <- cscan_summary_df[order(cscan_summary_df$c_value), , drop = FALSE]

write_tsv(cscan_summary_df, file.path(DIR_CSCAN, "c_scan_summary.tsv"))

eligible_c <- cscan_summary_df$c_value[cscan_summary_df$cselection_mean_empty <= CSELECTION_MAX_MEAN_EMPTY]
if (length(eligible_c) == 0) {
  eligible_c <- cscan_summary_df$c_value
  log_msg("所有候选 C 都出现了超过阈值的平均空簇，将退回到全部候选 C 上使用 Dmin 选取")
}

eligible_dmin_df <- dmin_df[dmin_df$c_value %in% eligible_c, , drop = FALSE]
recommended_c <- detect_dmin_elbow(
  crange = eligible_dmin_df$c_value,
  dmin_values = eligible_dmin_df$dmin_value
)

recommended_df <- data.frame(
  item = c(
    "recommended_c",
    "method",
    "DMIN_REPEATS",
    "CSELECTION_REPEATS",
    "CSELECTION_MAX_MEAN_EMPTY",
    "eligible_c_for_dmin"
  ),
  value = c(
    recommended_c,
    "Dmin_elbow_with_cselection_filter",
    DMIN_REPEATS,
    CSELECTION_REPEATS,
    CSELECTION_MAX_MEAN_EMPTY,
    paste(eligible_c, collapse = ",")
  ),
  stringsAsFactors = FALSE
)
write_tsv(recommended_df, file.path(DIR_CSCAN, "c_recommended.tsv"))

plot_dmin_diagnostic(
  dmin_df = dmin_df,
  chosen_c = recommended_c,
  out_prefix = file.path(DIR_CSCAN, "dmin_diagnostic")
)

plot_cselection_diagnostic(
  csel_df = cselection_summary_df,
  chosen_c = recommended_c,
  out_prefix = file.path(DIR_CSCAN, "cselection_diagnostic")
)

log_msg("候选 c 扫描完成；按官方路线自动推荐 C = ", recommended_c)

###############################################################################
# 十一、模块 M8：最终 c 的正式聚类
###############################################################################

if (is.na(FINAL_C)) {
  FINAL_C <- recommended_c
  log_msg("FINAL_C 未指定，按 Dmin + cselection 官方路线自动使用推荐值：", FINAL_C)
} else {
  log_msg("使用用户指定 FINAL_C：", FINAL_C)
}

log_msg("开始正式聚类，FINAL_C = ", FINAL_C)
cl_final <- Mfuzz::mfuzz(eset_std, c = FINAL_C, m = m_est)

final_membership <- extract_cluster_membership(cl_final, gene_ids = rownames(std_mat))
cluster_membership_df <- final_membership$membership_table
cluster_membership_matrix_df <- final_membership$membership_matrix

cluster_profile_summary_df <- get_cluster_profile_summary(
  std_mat = std_mat,
  cluster_assign = cluster_membership_df$cluster,
  stage_names = colnames(std_mat)
)

cluster_count_summary_df <- do.call(rbind, lapply(seq_len(FINAL_C), function(k) {
  sub <- cluster_membership_df[cluster_membership_df$cluster == k, , drop = FALSE]
  data.frame(
    cluster = k,
    n_all = nrow(sub),
    n_core_ge_0.7 = sum(sub$is_core, na.rm = TRUE),
    n_core_ge_0.8 = sum(sub$is_strict_core, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

cluster_summary_df <- merge(cluster_count_summary_df, cluster_profile_summary_df, by = "cluster", all.x = TRUE)

write_tsv(cluster_membership_df, file.path(DIR_FINAL, "cluster_membership.tsv"))
write_tsv(cluster_membership_matrix_df, file.path(DIR_FINAL, "cluster_membership_matrix.tsv"))
write_tsv(cluster_summary_df, file.path(DIR_FINAL, "cluster_summary.tsv"))

final_c_info_df <- data.frame(
  item = c("FINAL_C", "recommended_c_from_scan"),
  value = c(FINAL_C, recommended_c),
  stringsAsFactors = FALSE
)
write_tsv(final_c_info_df, file.path(DIR_FINAL, "final_c_selection_info.tsv"))

log_msg("正式聚类结果已输出")

###############################################################################
# 十二、模块 M9：导出每簇基因列表
###############################################################################

log_msg("导出每簇基因列表")
for (k in seq_len(FINAL_C)) {
  sub <- cluster_membership_df[cluster_membership_df$cluster == k, , drop = FALSE]

  all_df <- data.frame(gene_id = sub$gene_id, stringsAsFactors = FALSE)
  core_df <- data.frame(gene_id = sub$gene_id[sub$is_core], stringsAsFactors = FALSE)
  strict_df <- data.frame(gene_id = sub$gene_id[sub$is_strict_core], stringsAsFactors = FALSE)

  colnames(all_df)[1] <- GENEID_COL
  colnames(core_df)[1] <- GENEID_COL
  colnames(strict_df)[1] <- GENEID_COL

  write_tsv(all_df, file.path(DIR_LISTS, sprintf("cluster%02d_all.tsv", k)))
  write_tsv(core_df, file.path(DIR_LISTS, sprintf("cluster%02d_core.tsv", k)))
  write_tsv(strict_df, file.path(DIR_LISTS, sprintf("cluster%02d_strict.tsv", k)))
}

###############################################################################
# 十三、模块 M10：绘图
###############################################################################

log_msg("输出最终聚类图和每簇趋势图")

plot_cluster_overview_pub(
  std_mat = std_mat,
  cluster_assign = cluster_membership_df$cluster,
  cluster_summary_df = cluster_summary_df,
  out_prefix = file.path(DIR_PLOTS, "mfuzz_final_overview"),
  main_title = sprintf("Mfuzz final overview (c=%d)", FINAL_C)
)

plot_cluster_profiles_pub(
  std_mat = std_mat,
  cluster_assign = cluster_membership_df$cluster,
  cluster_summary_df = cluster_summary_df,
  outdir = DIR_PLOTS,
  prefix = "cluster"
)

plot_heatmap_pub(
  std_mat = std_mat,
  cluster_assign = cluster_membership_df$cluster,
  out_prefix = file.path(DIR_PLOTS, "cluster_heatmap")
)

plot_cluster_count_barplot(
  cluster_summary_df = cluster_summary_df,
  out_prefix = file.path(DIR_PLOTS, "cluster_gene_count_barplot")
)

plot_dmin_diagnostic(
  dmin_df = dmin_df,
  chosen_c = recommended_c,
  out_prefix = file.path(DIR_PLOTS, "dmin_diagnostic")
)

plot_cselection_diagnostic(
  csel_df = cselection_summary_df,
  chosen_c = recommended_c,
  out_prefix = file.path(DIR_PLOTS, "cselection_diagnostic")
)

###############################################################################
# 十四、模块 M11：导出供后续富集使用的 universe
###############################################################################

log_msg("导出 mfuzz universe")
mfuzz_universe_df <- data.frame(
  gene_id = rownames(stage_var_filtered),
  stringsAsFactors = FALSE
)
colnames(mfuzz_universe_df)[1] <- GENEID_COL
write_tsv(mfuzz_universe_df, file.path(DIR_ENRICH, "mfuzz_universe.tsv"))

###############################################################################
# 十五、附加输出：运行参数记录
###############################################################################

param_df <- data.frame(
  parameter = c(
    "TPM_FILE", "SAMPLES_FILE", "OUTDIR",
    "GENEID_COL", "SAMPLE_COL", "GROUP_COL",
    "MIN_TPM", "MIN_STAGE_COUNT",
    "VAR_METHOD", "KEEP_TOP_N",
    "LOG_PSEUDOCOUNT",
    "CANDIDATE_C", "FINAL_C",
    "DMIN_REPEATS", "CSELECTION_REPEATS", "CSELECTION_MAX_MEAN_EMPTY",
    "MEMBERSHIP_CORE", "MEMBERSHIP_STRICT",
    "PLOT_WIDTH", "PLOT_HEIGHT",
    "SAVE_PDF", "SAVE_PNG", "PNG_DPI",
    "FONT_FAMILY_BASE", "FONT_FAMILY_TITLE", "FONT_FAMILY_AXIS", "FONT_FAMILY_STRIP",
    "ENABLE_SHOWTEXT",
    "SHOW_RAW_TRAJECTORIES_IN_PUB", "RAW_TRAJECTORY_MAX_PER_CLUSTER",
    "RIBBON_STYLE", "RIBBON_ALPHA", "MEAN_LINE_WIDTH", "RAW_LINE_WIDTH", "RAW_LINE_ALPHA",
    "RANDOM_SEED"
  ),
  value = c(
    TPM_FILE, SAMPLES_FILE, OUTDIR,
    GENEID_COL, SAMPLE_COL, GROUP_COL,
    as.character(MIN_TPM), as.character(MIN_STAGE_COUNT),
    VAR_METHOD, as.character(KEEP_TOP_N),
    as.character(LOG_PSEUDOCOUNT),
    paste(CANDIDATE_C, collapse = ","),
    as.character(FINAL_C),
    as.character(DMIN_REPEATS), as.character(CSELECTION_REPEATS), as.character(CSELECTION_MAX_MEAN_EMPTY),
    as.character(MEMBERSHIP_CORE), as.character(MEMBERSHIP_STRICT),
    as.character(PLOT_WIDTH), as.character(PLOT_HEIGHT),
    as.character(SAVE_PDF), as.character(SAVE_PNG), as.character(PNG_DPI),
    FONT_FAMILY_BASE, FONT_FAMILY_TITLE, FONT_FAMILY_AXIS, FONT_FAMILY_STRIP,
    as.character(ENABLE_SHOWTEXT),
    as.character(SHOW_RAW_TRAJECTORIES_IN_PUB), as.character(RAW_TRAJECTORY_MAX_PER_CLUSTER),
    RIBBON_STYLE, as.character(RIBBON_ALPHA), as.character(MEAN_LINE_WIDTH), as.character(RAW_LINE_WIDTH), as.character(RAW_LINE_ALPHA),
    as.character(RANDOM_SEED)
  ),
  stringsAsFactors = FALSE
)
write_tsv(param_df, file.path(OUTDIR, "run_params.used.tsv"))

###############################################################################
# 十六、结束
###############################################################################

log_msg("Mfuzz 分析完成")
log_msg("输出目录：", normalizePath(OUTDIR, mustWork = FALSE))
log_msg("推荐查看：06_c_scan/dmin_values.tsv、06_c_scan/cselection_summary.tsv、06_c_scan/c_scan_summary.tsv")
log_msg("建议下一步：使用 08_cluster_gene_lists/ 中的 core 基因列表，结合 10_for_enrichment/mfuzz_universe.tsv，桥接到现有 07/08 富集流程")


