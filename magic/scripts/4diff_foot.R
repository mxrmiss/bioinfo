#!/usr/bin/env Rscript

# ============================================================
# 脚本名称：4diff_foot_shared_axis.R
# 功能：
# 1）读取 GO/KEGG 模块推荐表
# 2）读取四个物种的完整 GO BP / KEGG 富集结果表
# 3）仅保留 plot_recommendation == "Recommended" 的条目
# 4）按 species + term_id/pathway_id 精确匹配完整富集结果
# 5）计算每个物种每个模块的 module score = mean(-log10(p_adjust))
# 6）计算 Others mean = mean(Sr, Tg, Rp)
# 7）绘制 GO 与 KEGG 的 Sc vs Others mean 哑铃图
# 8）导出 PNG / PDF 和多个中间结果表
# ============================================================

# =========================
# 一、顶部参数区
# =========================

# ---- 输入目录与输出目录 ----
input_dir  <- "input"
output_dir <- file.path("output", "4diff_foot")

# ---- 输入文件：模块推荐表 ----
go_recommend_file   <- file.path(input_dir, "GO_BP_4module_plot_recommendation_full_by_species.tsv")
kegg_recommend_file <- file.path(input_dir, "KEGG_4module_plot_recommendation_full_by_species.tsv")

# ---- 输入文件：四个物种 GO 富集结果表 ----
go_result_files <- c(
  Sc = file.path(input_dir, "Sc_GO_BP_by_term_test.tsv"),
  Sr = file.path(input_dir, "Sr_GO_BP_by_term_test.tsv"),
  Tg = file.path(input_dir, "Tg_GO_BP_by_term_test.tsv"),
  Rp = file.path(input_dir, "Rp_GO_BP_by_term_test.tsv")
)

# ---- 输入文件：四个物种 KEGG 富集结果表 ----
kegg_result_files <- c(
  Sc = file.path(input_dir, "Sc_KEGG_by_term_test.tsv"),
  Sr = file.path(input_dir, "Sr_KEGG_by_term_test.tsv"),
  Tg = file.path(input_dir, "Tg_KEGG_by_term_test.tsv"),
  Rp = file.path(input_dir, "Rp_KEGG_by_term_test.tsv")
)

# ---- 允许的物种名称 ----
species_levels <- c("Sc", "Sr", "Tg", "Rp")

# ---- 模块展示顺序 ----
module_order <- c(
  "Neural regulation and excitability",
  "Contractile and cytoskeletal system",
  "Surface-interface remodeling and secretion",
  "Metabolic support"
)

# ---- p_adjust 为 0 时的保护下限 ----
global_min_positive_p <- 1e-300

# ---- 图形输出参数 ----
plot_width  <- 17
plot_height <- 5.2
plot_dpi    <- 600

# ---- 是否让 GO 与 KEGG 共用 x 轴范围 ----
USE_SHARED_X_LIMITS <- TRUE
X_LIMIT_PADDING_RATIO <- 0.08

# ---- x 轴刻度设置：更清爽地只显示 2、4、6 ----
USE_FIXED_X_LIMITS <- TRUE
FIXED_X_LIMITS <- c(1, 6)
X_AXIS_BREAKS <- c(2, 4, 6)

# ---- 配色：与气泡图的御用配色保持一致 ----
color_sc     <- "#F79D93"
color_others <- "#95C8F2"
segment_col  <- "#B8B8B8"

# ---- 哑铃点和连接线设置 ----
point_size <- 4.6
line_width <- 0.75

# ---- 字体大小 ----
base_size        <- 19
plot_title_size  <- 21
axis_text_x_size <- 20
axis_text_y_size <- 22
axis_title_size  <- 20
legend_text_size <- 18

# ---- 面板边框与水平参考线设置 ----
panel_border_width <- 0.35
axis_tick_width    <- 0.35
row_guide_color    <- "#E8E8E8"
row_guide_width    <- 0.28

# =========================
# 二、加载依赖包
# =========================

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "patchwork")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "缺少必要 R 包：", paste(missing_packages, collapse = ", "), "\n",
      "请先安装后再运行，例如：\n",
      "install.packages(c(", paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

# =========================
# 三、基础工具函数
# =========================

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("找不到输入文件：", path), call. = FALSE)
  }
}

check_required_columns <- function(df, required_cols, file_label) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "文件缺少必要列：", file_label, "\n",
        "缺失列为：", paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

normalize_species <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

safe_read_tsv <- function(path, file_label) {
  check_file_exists(path)
  df <- tryCatch(
    {
      readr::read_tsv(
        file = path,
        col_types = cols(.default = col_character()),
        progress = FALSE,
        na = c("", "NA", "NaN", "NULL")
      )
    },
    error = function(e) {
      stop(paste0("读取文件失败：", file_label, "\n原始报错：", e$message), call. = FALSE)
    }
  )
  return(df)
}

sanitize_p_adjust <- function(x, global_floor = 1e-300) {
  x_num <- suppressWarnings(as.numeric(x))
  x_num[!is.finite(x_num)] <- NA_real_

  positive_vals <- x_num[!is.na(x_num) & x_num > 0]
  if (length(positive_vals) > 0) {
    local_floor <- min(positive_vals) / 2
    local_floor <- max(local_floor, global_floor)
  } else {
    local_floor <- global_floor
  }

  x_fixed <- x_num
  x_fixed[!is.na(x_fixed) & x_fixed <= 0] <- local_floor
  return(x_fixed)
}

compute_neglog10_score <- function(p_adjust_vec) {
  p_fixed <- sanitize_p_adjust(p_adjust_vec, global_floor = global_min_positive_p)
  score   <- -log10(p_fixed)
  score[!is.finite(score)] <- NA_real_
  return(score)
}

compute_shared_x_limits <- function(x, padding_ratio = 0.08) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NULL)
  }

  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (x_min == x_max) {
    pad <- max(abs(x_min) * padding_ratio, 0.5)
  } else {
    pad <- (x_max - x_min) * padding_ratio
  }

  lower <- floor((x_min - pad) * 2) / 2
  upper <- ceiling((x_max + pad) * 2) / 2

  if (lower < 0) {
    lower <- 0
  }

  c(lower, upper)
}

msg <- function(...) {
  cat(paste0(..., "\n"))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 四、读取并检查模块推荐表
# =========================

msg("开始读取模块推荐表...")

go_rec <- safe_read_tsv(go_recommend_file, "GO 模块推荐表")
kegg_rec <- safe_read_tsv(kegg_recommend_file, "KEGG 模块推荐表")

check_required_columns(
  go_rec,
  c("species", "term_id", "term_name", "module", "plot_recommendation", "recommendation_reason"),
  "GO 模块推荐表"
)

check_required_columns(
  kegg_rec,
  c("species", "pathway_id", "term_name", "module", "plot_recommendation", "recommendation_reason"),
  "KEGG 模块推荐表"
)

go_rec <- go_rec %>%
  mutate(
    species = normalize_species(species),
    term_id = trimws(term_id),
    term_name = trimws(term_name),
    module = trimws(module),
    plot_recommendation = trimws(plot_recommendation)
  )

kegg_rec <- kegg_rec %>%
  mutate(
    species = normalize_species(species),
    pathway_id = trimws(pathway_id),
    term_name = trimws(term_name),
    module = trimws(module),
    plot_recommendation = trimws(plot_recommendation)
  )

bad_species_go <- setdiff(unique(go_rec$species), species_levels)
bad_species_kegg <- setdiff(unique(kegg_rec$species), species_levels)

if (length(bad_species_go) > 0) {
  stop(paste0("GO 模块推荐表中存在非法 species：", paste(bad_species_go, collapse = ", ")), call. = FALSE)
}
if (length(bad_species_kegg) > 0) {
  stop(paste0("KEGG 模块推荐表中存在非法 species：", paste(bad_species_kegg, collapse = ", ")), call. = FALSE)
}

go_rec_keep <- go_rec %>%
  filter(plot_recommendation == "Recommended") %>%
  filter(module %in% module_order)

kegg_rec_keep <- kegg_rec %>%
  filter(plot_recommendation == "Recommended") %>%
  filter(module %in% module_order)

if (nrow(go_rec_keep) == 0) {
  stop("GO 模块推荐表中没有 plot_recommendation == 'Recommended' 且 module 属于目标四模块的条目。", call. = FALSE)
}
if (nrow(kegg_rec_keep) == 0) {
  stop("KEGG 模块推荐表中没有 plot_recommendation == 'Recommended' 且 module 属于目标四模块的条目。", call. = FALSE)
}

readr::write_tsv(go_rec_keep, file.path(output_dir, "GO_recommended_terms_filtered.tsv"))
readr::write_tsv(kegg_rec_keep, file.path(output_dir, "KEGG_recommended_terms_filtered.tsv"))

# =========================
# 五、读取四个物种的完整富集结果表
# =========================

msg("开始读取四个物种的 GO / KEGG 富集结果表...")

go_required_cols <- c(
  "term_id", "term_name", "test_set", "ontology", "gene_ids", "gene_names",
  "gene_count", "bg_size", "gene_ratio", "bg_ratio", "p_value", "p_adjust",
  "universe_size", "min_gs", "max_gs"
)

go_list <- lapply(names(go_result_files), function(sp) {
  f <- go_result_files[[sp]]
  df <- safe_read_tsv(f, paste0(sp, " GO 富集结果表"))
  check_required_columns(df, go_required_cols, paste0(sp, " GO 富集结果表"))
  df %>%
    mutate(
      species = sp,
      term_id = trimws(term_id),
      term_name = trimws(term_name)
    )
})
names(go_list) <- names(go_result_files)

go_all <- bind_rows(go_list) %>%
  mutate(species = factor(species, levels = species_levels))

kegg_required_cols <- c(
  "pathway_id", "term_name", "test_set", "count_mode", "gene_ids", "gene_names",
  "gene_count", "bg_size", "gene_ratio", "bg_ratio", "p_value", "p_adjust",
  "universe_size", "min_gs", "max_gs", "ontology"
)

kegg_list <- lapply(names(kegg_result_files), function(sp) {
  f <- kegg_result_files[[sp]]
  df <- safe_read_tsv(f, paste0(sp, " KEGG 富集结果表"))
  check_required_columns(df, kegg_required_cols, paste0(sp, " KEGG 富集结果表"))
  df %>%
    mutate(
      species = sp,
      pathway_id = trimws(pathway_id),
      term_name = trimws(term_name)
    )
})
names(kegg_list) <- names(kegg_result_files)

kegg_all <- bind_rows(kegg_list) %>%
  mutate(species = factor(species, levels = species_levels))

# =========================
# 六、按 species + ID 匹配推荐条目与完整富集结果
# =========================

msg("开始匹配推荐条目与完整富集结果...")

go_matched <- go_rec_keep %>%
  left_join(
    go_all,
    by = c("species", "term_id"),
    suffix = c(".rec", ".full")
  ) %>%
  mutate(
    term_name_rec  = term_name.rec,
    term_name_full = term_name.full,
    term_name_match_flag = case_when(
      is.na(term_name_full) ~ "Missing_in_full_result",
      term_name_rec == term_name_full ~ "Exact_match",
      TRUE ~ "ID_matched_but_term_name_diff"
    )
  )

kegg_matched <- kegg_rec_keep %>%
  left_join(
    kegg_all,
    by = c("species", "pathway_id"),
    suffix = c(".rec", ".full")
  ) %>%
  mutate(
    term_name_rec  = term_name.rec,
    term_name_full = term_name.full,
    term_name_match_flag = case_when(
      is.na(term_name_full) ~ "Missing_in_full_result",
      term_name_rec == term_name_full ~ "Exact_match",
      TRUE ~ "ID_matched_but_term_name_diff"
    )
  )

go_match_stat <- go_matched %>%
  count(term_name_match_flag, name = "n")

kegg_match_stat <- kegg_matched %>%
  count(term_name_match_flag, name = "n")

readr::write_tsv(go_match_stat, file.path(output_dir, "GO_match_status_summary.tsv"))
readr::write_tsv(kegg_match_stat, file.path(output_dir, "KEGG_match_status_summary.tsv"))
readr::write_tsv(go_matched, file.path(output_dir, "GO_recommended_terms_matched_full.tsv"))
readr::write_tsv(kegg_matched, file.path(output_dir, "KEGG_recommended_terms_matched_full.tsv"))

go_missing_n <- sum(is.na(go_matched$p_adjust))
kegg_missing_n <- sum(is.na(kegg_matched$p_adjust))

msg("GO 匹配后缺少 p_adjust 的条目数：", go_missing_n)
msg("KEGG 匹配后缺少 p_adjust 的条目数：", kegg_missing_n)

# =========================
# 七、计算 term/pathway 级别分数
# =========================

msg("开始计算 term/pathway 级别的 -log10(p_adjust) ...")

go_term_level <- go_matched %>%
  mutate(
    p_adjust_raw   = p_adjust,
    p_adjust_fixed = sanitize_p_adjust(p_adjust),
    term_score     = compute_neglog10_score(p_adjust)
  ) %>%
  transmute(
    panel = "GO",
    species,
    module,
    term_id,
    term_name_recommend = term_name_rec,
    term_name_full = term_name_full,
    term_name_match_flag,
    plot_recommendation,
    recommendation_reason,
    ontology,
    gene_ids,
    gene_names,
    gene_count,
    bg_size,
    gene_ratio,
    bg_ratio,
    p_value,
    p_adjust_raw,
    p_adjust_fixed,
    term_score
  )

kegg_term_level <- kegg_matched %>%
  mutate(
    p_adjust_raw   = p_adjust,
    p_adjust_fixed = sanitize_p_adjust(p_adjust),
    term_score     = compute_neglog10_score(p_adjust)
  ) %>%
  transmute(
    panel = "KEGG",
    species,
    module,
    pathway_id,
    term_name_recommend = term_name_rec,
    term_name_full = term_name_full,
    term_name_match_flag,
    plot_recommendation,
    recommendation_reason,
    ontology,
    gene_ids,
    gene_names,
    gene_count,
    bg_size,
    gene_ratio,
    bg_ratio,
    p_value,
    p_adjust_raw,
    p_adjust_fixed,
    term_score
  )

readr::write_tsv(go_term_level, file.path(output_dir, "GO_term_level_scores.tsv"))
readr::write_tsv(kegg_term_level, file.path(output_dir, "KEGG_pathway_level_scores.tsv"))

# =========================
# 八、计算每个物种每个模块的 module score
# =========================

msg("开始计算每个物种每个模块的 module score ...")

go_module_species <- go_term_level %>%
  group_by(species, module) %>%
  summarise(
    n_terms_total = n(),
    n_terms_valid = sum(!is.na(term_score)),
    module_score  = ifelse(n_terms_valid > 0, mean(term_score, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    species = species_levels,
    module = module_order,
    fill = list(n_terms_total = 0, n_terms_valid = 0, module_score = NA_real_)
  ) %>%
  mutate(
    species = factor(species, levels = species_levels),
    module = factor(module, levels = module_order)
  ) %>%
  arrange(species, module)

kegg_module_species <- kegg_term_level %>%
  group_by(species, module) %>%
  summarise(
    n_terms_total = n(),
    n_terms_valid = sum(!is.na(term_score)),
    module_score  = ifelse(n_terms_valid > 0, mean(term_score, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    species = species_levels,
    module = module_order,
    fill = list(n_terms_total = 0, n_terms_valid = 0, module_score = NA_real_)
  ) %>%
  mutate(
    species = factor(species, levels = species_levels),
    module = factor(module, levels = module_order)
  ) %>%
  arrange(species, module)

readr::write_tsv(go_module_species, file.path(output_dir, "GO_module_score_by_species.tsv"))
readr::write_tsv(kegg_module_species, file.path(output_dir, "KEGG_module_score_by_species.tsv"))

# =========================
# 九、计算 Others mean
# =========================

msg("开始计算 Others mean ...")

compute_plot_data <- function(module_species_df, panel_name) {
  sc_df <- module_species_df %>%
    filter(species == "Sc") %>%
    transmute(
      panel = panel_name,
      module,
      Sc = module_score,
      Sc_n_terms_total = n_terms_total,
      Sc_n_terms_valid = n_terms_valid
    )

  others_df <- module_species_df %>%
    filter(species %in% c("Sr", "Tg", "Rp")) %>%
    group_by(module) %>%
    summarise(
      Others_mean = ifelse(sum(!is.na(module_score)) > 0, mean(module_score, na.rm = TRUE), NA_real_),
      Others_n_species_with_value = sum(!is.na(module_score)),
      Sr_score = module_score[species == "Sr"][1],
      Tg_score = module_score[species == "Tg"][1],
      Rp_score = module_score[species == "Rp"][1],
      .groups = "drop"
    )

  out <- full_join(sc_df, others_df, by = "module") %>%
    mutate(
      panel = panel_name,
      module = factor(as.character(module), levels = module_order)
    ) %>%
    arrange(module)

  return(out)
}

go_plot_data <- compute_plot_data(go_module_species, "GO")
kegg_plot_data <- compute_plot_data(kegg_module_species, "KEGG")

readr::write_tsv(go_plot_data, file.path(output_dir, "GO_module_summary_for_plot.tsv"))
readr::write_tsv(kegg_plot_data, file.path(output_dir, "KEGG_module_summary_for_plot.tsv"))

final_plot_wide <- bind_rows(go_plot_data, kegg_plot_data) %>%
  arrange(panel, module)

readr::write_tsv(final_plot_wide, file.path(output_dir, "final_plot_data_wide.tsv"))

final_plot_long <- final_plot_wide %>%
  select(panel, module, Sc, Others_mean) %>%
  pivot_longer(
    cols = c("Sc", "Others_mean"),
    names_to = "group",
    values_to = "module_score"
  ) %>%
  mutate(
    group = factor(group, levels = c("Sc", "Others_mean")),
    module = factor(as.character(module), levels = rev(module_order))
  )

readr::write_tsv(final_plot_long, file.path(output_dir, "final_plot_data_long.tsv"))

# =========================
# 十、绘图函数
# =========================

msg("开始绘图...")

shared_x_limits <- NULL
if (USE_SHARED_X_LIMITS) {
  if (USE_FIXED_X_LIMITS) {
    shared_x_limits <- FIXED_X_LIMITS
  } else {
    shared_x_limits <- compute_shared_x_limits(final_plot_long$module_score, X_LIMIT_PADDING_RATIO)
  }
}

make_panel_plot <- function(
  panel_name,
  plot_wide_df,
  plot_long_df,
  x_limits = NULL,
  show_y_axis = TRUE,
  left_margin = 8
) {
  panel_wide <- plot_wide_df %>%
    filter(panel == panel_name) %>%
    mutate(module = factor(as.character(module), levels = rev(module_order)))

  panel_long <- plot_long_df %>%
    filter(panel == panel_name) %>%
    mutate(module = factor(as.character(module), levels = rev(module_order)))

  p <- ggplot() +
    geom_segment(
      data = panel_wide,
      aes(
        x = Others_mean,
        xend = Sc,
        y = module,
        yend = module
      ),
      linewidth = line_width,
      color = segment_col,
      lineend = "round",
      na.rm = TRUE
    ) +
    geom_point(
      data = panel_long %>% filter(group == "Others_mean"),
      aes(x = module_score, y = module, color = group),
      size = point_size,
      na.rm = TRUE
    ) +
    geom_point(
      data = panel_long %>% filter(group == "Sc"),
      aes(x = module_score, y = module, color = group),
      size = point_size,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = c("Sc" = color_sc, "Others_mean" = color_others),
      labels = c("Sc", "Others mean"),
      name = NULL
    ) +
    labs(
      title = ifelse(panel_name == "GO", "GO BP", panel_name),
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = plot_title_size,
        color = "black",
        margin = margin(b = 8)
      ),
      axis.text.y = element_text(
        color = "black",
        size = axis_text_y_size
      ),
      axis.text.x = element_text(
        color = "black",
        size = axis_text_x_size
      ),
      axis.title.x = element_text(
        color = "black",
        size = axis_title_size,
        margin = margin(t = 7)
      ),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(
        linewidth = axis_tick_width,
        color = "black"
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = panel_border_width
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.grid.major.y = element_line(
        color = row_guide_color,
        linewidth = row_guide_width
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(
        color = "black",
        size = legend_text_size
      ),
      legend.key = element_blank(),
      legend.key.width = unit(0.9, "lines"),
      legend.spacing.x = unit(0.4, "lines"),
      plot.margin = margin(t = 6, r = 10, b = 4, l = left_margin)
    )

  if (!show_y_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_line(
          linewidth = axis_tick_width,
          color = "black"
        )
      )
  }

  if (!is.null(x_limits)) {
    p <- p +
      scale_x_continuous(
        limits = x_limits,
        breaks = X_AXIS_BREAKS
      )
  }

  return(p)
}

p_go <- make_panel_plot(
  panel_name = "GO",
  plot_wide_df = final_plot_wide,
  plot_long_df = final_plot_long,
  x_limits = shared_x_limits,
  show_y_axis = TRUE,
  left_margin = 8
)

p_kegg <- make_panel_plot(
  panel_name = "KEGG",
  plot_wide_df = final_plot_wide,
  plot_long_df = final_plot_long,
  x_limits = shared_x_limits,
  show_y_axis = FALSE,
  left_margin = 3
)

x_axis_label_plot <- ggplot() +
  annotate(
    "text",
    x = 0.5,
    y = 0.5,
    label = "Module score",
    size = axis_title_size / ggplot2::.pt,
    color = "black"
  ) +
  theme_void()

panel_plot <- p_go + p_kegg +
  plot_layout(ncol = 2, guides = "collect", widths = c(1, 1)) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = legend_text_size, color = "black")
  )

combined_plot <- panel_plot / x_axis_label_plot +
  plot_layout(heights = c(1, 0.09))

# =========================
# 十一、保存图形
# =========================

msg("开始导出图形文件...")

png_file <- file.path(output_dir, "Sc_vs_others_mean_4module_GO_KEGG.png")
pdf_file <- file.path(output_dir, "Sc_vs_others_mean_4module_GO_KEGG.pdf")

ggsave(
  filename = png_file,
  plot = combined_plot,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  units = "in",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = combined_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# =========================
# 十二、导出额外检查表
# =========================

msg("导出额外检查表...")

go_module_term_map <- go_term_level %>%
  arrange(species, module, desc(term_score), term_id)

kegg_module_term_map <- kegg_term_level %>%
  arrange(species, module, desc(term_score), pathway_id)

readr::write_tsv(go_module_term_map, file.path(output_dir, "GO_module_term_membership.tsv"))
readr::write_tsv(kegg_module_term_map, file.path(output_dir, "KEGG_module_pathway_membership.tsv"))

run_summary <- tibble::tibble(
  item = c(
    "GO 推荐条目数",
    "KEGG 推荐条目数",
    "GO 匹配后条目数",
    "KEGG 匹配后条目数",
    "GO 缺少 p_adjust 条目数",
    "KEGG 缺少 p_adjust 条目数",
    "输出目录"
  ),
  value = c(
    nrow(go_rec_keep),
    nrow(kegg_rec_keep),
    nrow(go_matched),
    nrow(kegg_matched),
    go_missing_n,
    kegg_missing_n,
    output_dir
  )
)

readr::write_tsv(run_summary, file.path(output_dir, "run_summary.tsv"))

# =========================
# 十三、运行结束提示
# =========================

msg("全部完成。")
msg("输出目录：", output_dir)
msg("主要图文件：")
msg("  - ", png_file)
msg("  - ", pdf_file)
msg("主要数据文件：")
msg("  - GO_module_score_by_species.tsv")
msg("  - KEGG_module_score_by_species.tsv")
msg("  - GO_module_summary_for_plot.tsv")
msg("  - KEGG_module_summary_for_plot.tsv")
msg("  - final_plot_data_wide.tsv")
msg("  - final_plot_data_long.tsv")
msg("  - GO_module_term_membership.tsv")
msg("  - KEGG_module_pathway_membership.tsv")

