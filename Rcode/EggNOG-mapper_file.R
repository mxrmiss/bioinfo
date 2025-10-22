################################################################################
#
#           脚本 1: 从EggNOG-mapper结果中生成注释文件 (修改版)
#
# 目的:
#   - 读取 EggNOG-mapper 的 '.emapper.annotations' 输出文件。
#   - 提取基因ID与GO term的对应关系。
#   - 提取基因ID与KEGG KO的对应关系。
#   - 将提取出的干净数据保存为两个TSV文件，供后续分析使用。
#
# 使用方法:
#   1. 在此脚本的同级目录下创建名为 'inmaterial' 的文件夹，
#      并将您的 '.emapper.annotations' 文件放入其中。
#   2. 在下方【用户配置区】修改输入文件名。
#   3. 运行此脚本。成功后，会在脚本同级目录下自动创建 'oumaterial' 文件夹，
#      并生成 'parsed_gene_go_annotation.tsv' 和
#      'parsed_gene_kegg_annotation.tsv' 两个文件。
#
################################################################################

#==============================================================================
# 步骤 0: 安装和加载必要的R包
#==============================================================================
packages_to_install <- c("dplyr", "tidyr", "stringr")
for(pkg in packages_to_install){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(dplyr)
library(tidyr)
library(stringr)

cat("所有必要的R包已成功加载！\n")

#==============================================================================
#                         ----------- 用户配置区 -----------
#==============================================================================
# 请将这里的文件名修改为您从EggNOG-mapper下载的annotations文件名
EGGNOG_ANNOTATION_FILE_BASENAME <- "MM_2z45cx4s.emapper.annotations.tsv"

# --- 目录设置 ---
# 输入文件将在与此脚本同级的 'inmaterial' 文件夹中寻找
INPUT_DIR <- "inmaterial"
# 输出文件将保存在与此脚本同级的 'oumaterial' 文件夹中
OUTPUT_DIR <- "oumaterial"

# 构建完整输入文件路径
EGGNOG_ANNOTATION_FILE <- file.path(INPUT_DIR, EGGNOG_ANNOTATION_FILE_BASENAME)

#==============================================================================
# 步骤 1: 读取并解析EggNOG-mapper的输出文件
#==============================================================================
cat("正在读取并解析EggNOG-mapper文件...\n")

# 检查输入文件是否存在
if (!file.exists(EGGNOG_ANNOTATION_FILE)) {
  stop(paste("错误：输入文件不存在，请检查路径:", EGGNOG_ANNOTATION_FILE))
}

# 读取文件，EggNOG的输出文件可能会有很多'#'开头的注释行，我们跳过它们
eggnog_data <- read.delim(EGGNOG_ANNOTATION_FILE,
                          header = FALSE,          # 文件没有标准的表头
                          comment.char = "#",      # 跳过以'#'开头的行
                          stringsAsFactors = FALSE,
                          sep = "\t")              # 明确指定制表符为分隔符

# 为数据框添加列名，根据EggNOG v2的输出格式
# 第1列是基因ID, 第6列是GOs, 第8列是KEGG_ko (请根据您的文件版本核对)
# 为了代码的稳健性，我们只选择需要的列进行命名
colnames(eggnog_data)[c(1, 6, 8)] <- c("GeneID", "GOs", "KEGG_ko")


# 1.1 提取并处理GO注释
cat("正在提取GO注释...\n")
gene_to_go <- eggnog_data %>%
  dplyr::select(GeneID, GOs) %>%
  filter(GOs != "" & !is.na(GOs)) %>% # 过滤掉没有GO注释的行
  separate_rows(GOs, sep = ",") %>%   # 按逗号拆分GOs
  rename(GO_ID = GOs) %>%
  distinct()

# 1.2 提取并处理KEGG注释
cat("正在提取KEGG注释...\n")
gene_to_kegg <- eggnog_data %>%
  dplyr::select(GeneID, KEGG_ko) %>%
  filter(KEGG_ko != "" & !is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  # EggNOG的输出是 "ko:K12345", 我们需要去掉 "ko:"
  mutate(KO_ID = str_remove(KEGG_ko, "ko:")) %>%
  dplyr::select(GeneID, KO_ID) %>%
  distinct()

#==============================================================================
# 步骤 2: 保存处理好的注释文件
#==============================================================================
# 检查输出目录是否存在，如果不存在则创建
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(paste("输出目录 '", OUTPUT_DIR, "' 已创建。\n", sep=""))
}

# 构建完整的输出文件路径
output_go_file <- file.path(OUTPUT_DIR, "parsed_gene_go_annotation.tsv")
write.table(gene_to_go, output_go_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("GO注释文件已成功保存为:", output_go_file, "\n"))

output_kegg_file <- file.path(OUTPUT_DIR, "parsed_gene_kegg_annotation.tsv")
write.table(gene_to_kegg, output_kegg_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("KEGG注释文件已成功保存为:", output_kegg_file, "\n"))

cat("脚本 '01_Parse_EggNOG_to_Annotation.R' 运行完毕！\n")