################################################################################
#
#           通用富集分析脚本 (最终版 - 智能检查)
#
# 功能:
#   - 智能检查包依赖，只在需要时安装，避免重复警告。
#   - 兼容模式与非模式生物 (通过 ANALYSIS_MODE 切换)。
#   - KEGG分析采用在线查询，兼容性最强。
#   - 输出标准化的结果表格和多种可视化图表。
#
# 作者: [你的名字] (由Gemini最终优化)
# 日期: 2025-08-22
#
################################################################################

#==============================================================================
# 步骤 0: 参数配置与环境准备
#==============================================================================

# --- 1. [核心开关] 选择您的分析模式 ---
# "model"    : 适用于有Bioconductor OrgDb包的物种 (如 Human, Mouse, Fly等)。
# "non_model": 适用于任何其他物种，需要您提供注释文件。
ANALYSIS_MODE <- "non_model"  # <-- 在这里切换 "model" 或 "non_model"


# --- 2. [模式生物] 配置区 (仅在 ANALYSIS_MODE = "model" 时填写) ---
ORGANISM_SCIDB  <- "hsa"
ORGANISM_DB     <- "org.Hs.eg.db"
KEY_TYPE        <- "ENSEMBL"
DEG_GENE_FILE   <- "DEGs_gene_ids.txt"
ALL_GENE_FILE   <- "all_genes.txt"


# --- 3. [非模式生物] 配置区 (仅在 ANALYSIS_MODE = "non_model" 时填写) ---
INPUT_DIR      <- "inmaterial"
GO_ANNO_FILE   <- "parsed_gene_go_annotation_updated.tsv"
KEGG_ANNO_FILE <- "parsed_gene_kegg_annotation_updated.tsv"


# --- 4. [通用] 配置区 ---
OUTPUT_DIR <- "oumaterial"
P_ADJUST_METHOD <- "BH"
P_VALUE_CUTOFF  <- 0.05
Q_VALUE_CUTOFF  <- 0.2
SHOW_CATEGORY_NUM <- 20


# --- 5. [全新] 智能化的包管理 ---
# 确保BiocManager已安装
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 定义所有需要的包
base_packages <- c("clusterProfiler", "enrichplot", "dplyr", "ggplot2", "data.table", "pacman")
if (ANALYSIS_MODE == "model") {
  required_packages <- c(base_packages, ORGANISM_DB)
} else {
  required_packages <- c(base_packages, "KEGGREST")
}

# 找出尚未安装的包
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

# 如果有缺失的包，则只安装这些缺失的包
if (length(missing_packages) > 0) {
  cat("发现以下缺失的包:", paste(missing_packages, collapse = ", "), "\n")
  cat("正在进行安装...\n")
  BiocManager::install(missing_packages, ask = FALSE)
} else {
  cat("所有必需的R包均已安装。\n")
}

# 使用pacman加载所有包 (现在会非常安静)
pacman::p_load(char = required_packages)

options(timeout = 600)
cat("所有必要的R包已成功加载。\n")


#==============================================================================
# 步骤 1: 定义核心函数 (无修改)
#==============================================================================

load_gene_list <- function(file_path) {
  if (!file.exists(file_path)) {stop(paste("错误: 基因列表文件不存在:", file_path))}
  cat(paste0("  - 正在加载基因列表: ", basename(file_path), "\n"))
  return(fread(file_path, header = FALSE)$V1)
}

load_annotation_file <- function(file_path) {
  if (!file.exists(file_path)) {stop(paste("错误: 注释文件不存在:", file_path))}
  cat(paste0("  - 正在加载注释文件: ", basename(file_path), "\n"))
  return(as.data.frame(fread(file_path, header = TRUE)))
}

save_enrichment_results <- function(enrichment_type, enrich_result, output_dir) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    cat(paste0("警告: ", enrichment_type, " 分析无显著结果，跳过保存。\n")); return()
  }
  cat(paste0("  - 正在保存 ", enrichment_type, " 结果...\n"))
  write.table(as.data.frame(enrich_result), file.path(output_dir, paste0(enrichment_type, "_results.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  p_bar <- barplot(enrich_result, showCategory = SHOW_CATEGORY_NUM)
  ggsave(file.path(output_dir, paste0(enrichment_type, "_barplot.png")), p_bar, width = 12, height = 8, dpi = 300)
  p_dot <- dotplot(enrich_result, showCategory = SHOW_CATEGORY_NUM)
  ggsave(file.path(output_dir, paste0(enrichment_type, "_dotplot.png")), p_dot, width = 12, height = 8, dpi = 300)
  show_net_num <- min(nrow(enrich_result), 5)
  if (show_net_num > 0) {
    p_cnet <- cnetplot(enrich_result, showCategory = show_net_num, layout = "fr")
    ggsave(file.path(output_dir, paste0(enrichment_type, "_cnetplot.png")), p_cnet, width = 10, height = 10, dpi = 300)
  }
  cat(paste0("    - 结果已保存至 '", output_dir, "' 文件夹。\n"))
}

convert_gene_ids <- function(gene_list, target_id, org_db, key_type) {
  cat(paste0("  - 正在将 ", key_type, " 转换为 ", target_id, "...\n"))
  conversion_result <- bitr(gene_list, fromType = key_type, toType = target_id, OrgDb = org_db)
  original_count <- length(unique(gene_list))
  converted_count <- length(unique(conversion_result[[target_id]]))
  cat(sprintf("    - %d 个原始ID成功转换为 %d 个目标ID (转换率: %.1f%%)\n", 
              original_count, converted_count, 100 * converted_count / original_count))
  return(conversion_result)
}


#==============================================================================
# 主程序: 根据ANALYSIS_MODE执行不同工作流 (无修改)
#==============================================================================

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

if (ANALYSIS_MODE == "model") {
  
  cat("\n--- [模式生物模式] 开始分析 ---\n")
  deg_genes_original <- load_gene_list(DEG_GENE_FILE)
  all_genes_original <- load_gene_list(ALL_GENE_FILE)
  
  deg_genes_entrez <- convert_gene_ids(deg_genes_original, "ENTREZID", ORGANISM_DB, KEY_TYPE)$ENTREZID
  all_genes_entrez <- convert_gene_ids(all_genes_original, "ENTREZID", ORGANISM_DB, KEY_TYPE)$ENTREZID
  
  if(length(deg_genes_entrez) == 0) {stop("错误: 没有任何差异基因成功转换ID，无法继续分析。")}
  
  cat("\n--- 正在执行 GO 富集分析 ---\n")
  ego_result <- enrichGO(gene = deg_genes_entrez, universe = all_genes_entrez, OrgDb = get(ORGANISM_DB), keyType = "ENTREZID", ont = "ALL", pAdjustMethod = P_ADJUST_METHOD, pvalueCutoff = P_VALUE_CUTOFF, qvalueCutoff = Q_VALUE_CUTOFF)
  save_enrichment_results("GO", ego_result, OUTPUT_DIR)
  
  cat("\n--- 正在执行 KEGG 富集分析 ---\n")
  ekk_result <- enrichKEGG(gene = deg_genes_entrez, universe = all_genes_entrez, organism = ORGANISM_SCIDB, pAdjustMethod = P_ADJUST_METHOD, pvalueCutoff = P_VALUE_CUTOFF, qvalueCutoff = Q_VALUE_CUTOFF)
  save_enrichment_results("KEGG", ekk_result, OUTPUT_DIR)
  
} else if (ANALYSIS_MODE == "non_model") {
  
  cat("\n--- [非模式生物模式] 开始分析 ---\n")
  gene_to_go   <- load_annotation_file(file.path(INPUT_DIR, GO_ANNO_FILE))
  gene_to_kegg <- load_annotation_file(file.path(INPUT_DIR, KEGG_ANNO_FILE))
  deg_genes    <- load_gene_list(file.path(INPUT_DIR, DEG_GEN_FILE))
  all_genes    <- load_gene_list(file.path(INPUT_DIR, ALL_GEN_FILE))
  
  cat("\n--- 正在执行 GO 富集分析 ---\n")
  go_term2gene <- gene_to_go %>% dplyr::select(GO_ID, GeneID) %>% filter(GO_ID != "" & !is.na(GO_ID))
  ego_result <- enricher(gene = deg_genes, TERM2GENE = go_term2gene, universe = all_genes, pAdjustMethod = P_ADJUST_METHOD, pvalueCutoff = P_VALUE_CUTOFF, qvalueCutoff = Q_VALUE_CUTOFF)
  save_enrichment_results("GO", ego_result, OUTPUT_DIR)
  
  cat("\n--- 正在执行 KEGG 富集分析 (在线模式) ---\n")
  ekk_result <- NULL
  tryCatch({
    cat("  - 正在从KEGG官网实时下载通路注释...\n")
    kegg_term2gene <- keggLink("ko", "pathway") %>% data.frame(pathway = gsub("path:", "", names(.)), ko = gsub("ko:", "", .)) %>% dplyr::select(pathway, ko)
    kegg_term2name <- keggList("pathway") %>% data.frame(pathway = gsub("path:ko", "ko", names(.)), description = .)
    de_ko <- gene_to_kegg %>% filter(GeneID %in% deg_genes, !is.na(KO_ID), KO_ID != "") %>% pull(KO_ID) %>% unique()
    universe_ko <- gene_to_kegg %>% filter(GeneID %in% all_genes, !is.na(KO_ID), KO_ID != "") %>% pull(KO_ID) %>% unique()
    if (length(de_ko) > 0) {
      cat("  - 开始在线富集分析...\n")
      ekk_result <- enricher(gene = de_ko, universe = universe_ko, TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name, pAdjustMethod = P_ADJUST_METHOD, pvalueCutoff = P_VALUE_CUTOFF, qvalueCutoff = Q_VALUE_CUTOFF)
    } else {cat("警告: 差异基因未能映射到任何KO条目，跳过KEGG分析。\n")}
  }, error = function(e) {cat("错误: KEGG在线分析失败！信息:", e$message, "\n")})
  
  save_enrichment_results("KEGG", ekk_result, OUTPUT_DIR)
  
} else {
  stop("错误: 无效的 ANALYSIS_MODE。请选择 'model' 或 'non_model'。")
}

#==============================================================================
# 脚本结束
#==============================================================================
cat("\n分析全部完成。请检查 '", OUTPUT_DIR, "' 文件夹获取结果。\n")
################################################################################