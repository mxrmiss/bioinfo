#!/usr/bin/env Rscript

# ============================================================
# 脚本名称：plot_4diff_foot.R
# 功能：
# 1）读取 GO/KEGG 模块推荐表
# 2）读取四个物种的完整 GO BP / KEGG 富集结果表
# 3）仅保留 plot_recommendation == "Recommended" 的条目
# 4）按 species + term_id/pathway_id 精确匹配完整富集结果
# 5）计算每个物种每个模块的 module score = mean(-log10(p_adjust))
# 6）计算 Others mean = mean(Sr, Tg, Rp)
# 7）绘制左右并排双 panel 图（左 GO，右 KEGG）
# 8）导出 PNG / PDF 和多个中间结果表
#
# 运行方式：
# 在项目根目录 ~/project/magic 下运行
# Rscript scripts/plot_4diff_foot.R
#
# 注意：
# - 所有路径均使用相对路径
# - 所有参数集中在脚本顶部
# - 输出目录：output/4diff_foot/
# ============================================================


# =========================
# 一、顶部参数区
# =========================

# ---- 输入目录与输出目录（相对于项目根目录） ----
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

# ---- 必须使用的 4 个模块名称（严格按要求） ----
module_order <- c(
  "Neural regulation and excitability",
  "Contractile and cytoskeletal system",
  "Surface-interface remodeling and secretion",
  "Metabolic support"
)

# ---- 模块得分中 p_adjust 为 0 时的保护下限 ----
# 说明：
# 当 p_adjust == 0 时，-log10(0) 会变成 Inf，因此需要替换成一个极小正数
# 这里设置一个理论最小值，同时在实际计算时也会优先使用“数据中最小正值的一半”
global_min_positive_p <- 1e-300

# ---- 图形输出参数 ----
plot_width  <- 17
plot_height <- 4.7
plot_dpi    <- 600

# ---- 颜色设置（简洁、适合论文）----
color_sc     <- "#C44E52"   # Sc
color_others <- "#4C72B0"   # Others mean
segment_col  <- "#7A7A7A"

# ---- 点大小与线宽 ----
point_size   <- 3.2
line_width   <- 0.65

# ---- 字体大小 ----
base_size <- 21


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

# ---- 路径存在性检查 ----
check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("找不到输入文件：", path), call. = FALSE)
  }
}

# ---- 检查必要列是否存在 ----
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

# ---- 标准化物种名 ----
normalize_species <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

# ---- 安全读取 TSV ----
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

# ---- 将 p_adjust 转成安全的可计算数值 ----
# 逻辑：
# 1）先转成数值
# 2）保留有限值
# 3）如果 p_adjust == 0，则替换为“当前向量中最小正数的一半”
# 4）如果没有任何正数，则使用 global_min_positive_p
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

# ---- 根据 p_adjust 计算 term/pathway 级别得分 ----
# score_term = -log10(p_adjust_fixed)
compute_neglog10_score <- function(p_adjust_vec) {
  p_fixed <- sanitize_p_adjust(p_adjust_vec, global_floor = global_min_positive_p)
  score   <- -log10(p_fixed)
  score[!is.finite(score)] <- NA_real_
  return(score)
}

# ---- 创建目录（若不存在则自动创建） ----
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 输出消息函数 ----
msg <- function(...) {
  cat(paste0(..., "\n"))
}


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

# ---- 检查物种名 ----
bad_species_go <- setdiff(unique(go_rec$species), species_levels)
bad_species_kegg <- setdiff(unique(kegg_rec$species), species_levels)

if (length(bad_species_go) > 0) {
  stop(paste0("GO 模块推荐表中存在非法 species：", paste(bad_species_go, collapse = ", ")), call. = FALSE)
}
if (length(bad_species_kegg) > 0) {
  stop(paste0("KEGG 模块推荐表中存在非法 species：", paste(bad_species_kegg, collapse = ", ")), call. = FALSE)
}

# ---- 仅保留推荐条目 ----
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

# 导出保留后的推荐表，方便检查
readr::write_tsv(go_rec_keep, file.path(output_dir, "GO_recommended_terms_filtered.tsv"))
readr::write_tsv(kegg_rec_keep, file.path(output_dir, "KEGG_recommended_terms_filtered.tsv"))


# =========================
# 五、读取四个物种的完整富集结果表
# =========================

msg("开始读取四个物种的 GO / KEGG 富集结果表...")

# ---- 读取 GO ----
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

# ---- 读取 KEGG ----
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

# ---- GO：用 species + term_id 优先精确匹配 ----
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

# ---- KEGG：用 species + pathway_id 优先精确匹配 ----
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

# ---- 输出匹配状态汇总，方便检查 ----
go_match_stat <- go_matched %>%
  count(term_name_match_flag, name = "n")

kegg_match_stat <- kegg_matched %>%
  count(term_name_match_flag, name = "n")

readr::write_tsv(go_match_stat, file.path(output_dir, "GO_match_status_summary.tsv"))
readr::write_tsv(kegg_match_stat, file.path(output_dir, "KEGG_match_status_summary.tsv"))

# ---- 输出完整匹配明细 ----
readr::write_tsv(go_matched, file.path(output_dir, "GO_recommended_terms_matched_full.tsv"))
readr::write_tsv(kegg_matched, file.path(output_dir, "KEGG_recommended_terms_matched_full.tsv"))

# ---- 如果某些条目在完整富集结果表中找不到，给出提示，但不直接崩溃 ----
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

# 导出 term/pathway 级别分数
readr::write_tsv(go_term_level, file.path(output_dir, "GO_term_level_scores.tsv"))
readr::write_tsv(kegg_term_level, file.path(output_dir, "KEGG_pathway_level_scores.tsv"))


# =========================
# 八、计算每个物种每个模块的 module score
# =========================

msg("开始计算每个物种每个模块的 module score ...")

# 说明：
# module score = mean(-log10(p_adjust_fixed))
# 若某模块下没有任何可用 term_score，则记为 NA，不报错

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

# 导出每个物种每个模块的 score 汇总表
readr::write_tsv(go_module_species, file.path(output_dir, "GO_module_score_by_species.tsv"))
readr::write_tsv(kegg_module_species, file.path(output_dir, "KEGG_module_score_by_species.tsv"))


# =========================
# 九、计算 Others mean（Sr + Tg + Rp）
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

# 导出 GO/KEGG 模块汇总表
readr::write_tsv(go_plot_data, file.path(output_dir, "GO_module_summary_for_plot.tsv"))
readr::write_tsv(kegg_plot_data, file.path(output_dir, "KEGG_module_summary_for_plot.tsv"))

# 合并成最终作图表
final_plot_wide <- bind_rows(go_plot_data, kegg_plot_data) %>%
  arrange(panel, module)

readr::write_tsv(final_plot_wide, file.path(output_dir, "final_plot_data_wide.tsv"))

# 转成长表，方便画点
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

make_panel_plot <- function(panel_name, plot_wide_df, plot_long_df) {

  # 当前 panel 的宽表数据，用于画连接线
  panel_wide <- plot_wide_df %>%
    filter(panel == panel_name) %>%
    mutate(module = factor(as.character(module), levels = rev(module_order)))

  # 当前 panel 的长表数据，用于画点
  panel_long <- plot_long_df %>%
    filter(panel == panel_name) %>%
    mutate(module = factor(as.character(module), levels = rev(module_order)))

  p <- ggplot() +
    # 先画连接线
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
      na.rm = TRUE
    ) +
    # 再画点
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
      title = panel_name,
      x = "Module score",
      y = NULL
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(color = "black"),
	  axis.line.x = element_line(linewidth = 0.5, color = "black"),
	  axis.line.y = element_line(linewidth = 0.5, color = "black"),
	  axis.ticks  = element_line(linewidth = 0.4, color = "black")
    )

  return(p)
}

p_go   <- make_panel_plot("GO", final_plot_wide, final_plot_long)
p_kegg <- make_panel_plot("KEGG", final_plot_wide, final_plot_long)

combined_plot <- p_go + p_kegg + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")


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

# GO：每个模块由哪些 term 构成
go_module_term_map <- go_term_level %>%
  arrange(species, module, desc(term_score), term_id)

# KEGG：每个模块由哪些 pathway 构成
kegg_module_term_map <- kegg_term_level %>%
  arrange(species, module, desc(term_score), pathway_id)

readr::write_tsv(go_module_term_map, file.path(output_dir, "GO_module_term_membership.tsv"))
readr::write_tsv(kegg_module_term_map, file.path(output_dir, "KEGG_module_pathway_membership.tsv"))

# 最终汇总说明表
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
