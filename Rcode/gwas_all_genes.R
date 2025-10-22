################################################################################
#
#       脚本: 从GFF3文件生成背景基因列表 (all_genes.txt) - 灵活配置版
#
# 目的:
#   - 允许用户通过修改配置变量，灵活地从GFF3文件中提取基因ID。
#   - 从指定的输入文件夹读取文件，并将结果保存到指定的输出文件夹。
#
# 使用前须知:
#   - !!! 最重要 !!! 请根据您的GFF3文件格式，正确修改下方的【用户配置区】。
#   - GFF3文件应放在 'inmaterial' 文件夹内。
#
################################################################################

#==============================================================================
# 步骤 0: 安装和加载必要的R包
#==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(rtracklayer)
library(dplyr)

cat("所有必要的R包已成功加载！\n")


#==============================================================================
#
#                      ----------- 用户配置区 -----------
#
# !!! 请根据您的GFF3文件，仔细检查并修改以下两个变量 !!!
#==============================================================================

# 1. 您想从哪种特征类型中提取基因ID？
#    打开GFF3文件，查看第3列。常见选项有 "gene", "mRNA", "CDS" 等。
FEATURE_TYPE_TO_EXTRACT <- "gene"

# 2. 包含您想要的基因ID的属性名叫什么？
#    打开GFF3文件，查看第9列。查看您的基因ID是跟在哪个标签后面的。
#    常见选项有 "ID", "Name", "Parent", "gene_id" 等。
GENE_ID_ATTRIBUTE <- "ID"

#==============================================================================


#==============================================================================
# 步骤 1: 定义文件路径并读取GFF3文件
#==============================================================================
input_file_path <- "inmaterial/Sco_annotation.gff3"
output_dir <- "oumaterial"
output_file_path <- file.path(output_dir, "all_genes.txt")

cat(paste("准备读取输入文件:", input_file_path, "\n"))

if (!file.exists(input_file_path)) {
  stop(paste("错误：找不到输入文件 '", input_file_path, "'。"))
}

gff_data <- rtracklayer::import(input_file_path)
gff_df <- as.data.frame(gff_data)

cat("GFF3文件读取成功！\n")


#==============================================================================
# 步骤 2: 根据您的配置，灵活地筛选和提取基因ID
#==============================================================================
cat(paste("正在从特征 '", FEATURE_TYPE_TO_EXTRACT, "' 的属性 '", GENE_ID_ATTRIBUTE, "' 中提取ID...\n", sep=""))

# 检查配置的属性列是否存在于数据中
if (!(GENE_ID_ATTRIBUTE %in% colnames(gff_df))) {
  stop(paste("错误：在GFF文件中找不到名为 '", GENE_ID_ATTRIBUTE, "' 的属性列。\n请检查您的配置，并确保它与GFF文件第9列中的标签完全匹配。"))
}

all_genes_list <- gff_df %>%
  # 第一步：根据用户配置的特征类型进行筛选
  filter(type == FEATURE_TYPE_TO_EXTRACT) %>%
  
  # 第二步：根据用户配置的属性名提取ID
  # `!!sym()` 是一个特殊的技巧，它允许我们在dplyr函数中使用字符串变量来指定列名
  pull(!!sym(GENE_ID_ATTRIBUTE)) %>%
  
  # 第三步：确保所有ID都是唯一的
  # GFF的属性列有时是列表，先转换为字符，再取唯一值
  as.character() %>%
  unique()

# 检查是否成功提取到ID
if (length(all_genes_list) == 0) {
  warning("警告：根据您的配置，没有提取到任何基因ID。请仔细检查您的GFF文件和【用户配置区】的设置是否正确。")
} else {
  cat(paste("成功提取到", length(all_genes_list), "个唯一的基因ID。\n"))
  cat("基因ID预览：\n")
  print(head(all_genes_list))
}


#==============================================================================
# 步骤 3: 创建输出目录并保存文件
#==============================================================================
cat(paste("\n准备将基因ID列表写入文件:", output_file_path, "...\n"))

if (!dir.exists(output_dir)) {
  cat(paste("输出文件夹 '", output_dir, "' 不存在，正在创建...\n", sep=""))
  dir.create(output_dir)
}

write.table(all_genes_list,
            file = output_file_path,
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

cat("----------------------------------------\n")
cat(paste("脚本运行成功！背景基因文件已在 '", output_dir, "' 文件夹中生成。\n", sep=""))