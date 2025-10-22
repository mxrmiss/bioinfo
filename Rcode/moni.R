################################################################################
#
#       脚本: 从GEO下载的数据生成 all_genes.txt 和 DEG_genes.txt
#
# 目的:
#   1. 读取从GEO下载的基因表达矩阵文件。
#   2. 使用DESeq2包进行标准的差异表达分析。
#   3. 从分析结果中，提取所有参与分析的基因作为背景基因集 (all_genes.txt)。
#   4. 根据p值和表达倍数变化的阈值，筛选出差异表达基因 (DEG_genes.txt)。
#
################################################################################

#==============================================================================
# 步骤 0: 安装和加载必要的R包
#================================################==============================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(DESeq2)
library(dplyr)

cat("所有必要的R包已成功加载！\n")


#==============================================================================
# 步骤 1: 读取并准备数据
#==============================================================================
cat("--- 正在读取从GEO下载的表达矩阵 ---\n")

# 输入文件名
count_file <- "GSE125511_processed_data.txt.gz"

# 读取数据，第一列是基因名，所以我们将其设为行名
count_data <- read.delim(count_file, header = TRUE, row.names = 1, check.names = FALSE)

# 查看数据前几行
print(head(count_data))

# 根据GEO页面的样本信息，创建样本信息表 (colData)
# 前3个样本是对照组(CK)，后3个样本是酸化处理组(T)
sample_info <- data.frame(
  condition = factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")),
  row.names = colnames(count_data)
)

cat("\n样本信息表:\n")
print(sample_info)


#==============================================================================
# 步骤 2: 运行DESeq2差异表达分析
#==============================================================================
cat("\n--- 正在运行DESeq2进行差异表达分析... ---\n")

# 构建DESeq2对象
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

# 过滤掉在所有样本中表达量都极低的基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 运行差异分析
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "Treatment", "Control"))

# 将结果转换为数据框并查看
res_df <- as.data.frame(res)
cat("\n差异表达分析结果预览:\n")
print(head(res_df))


#==============================================================================
# 步骤 3: 生成 all_genes.txt 和 DEG_genes.txt
#==============================================================================
cat("\n--- 正在生成最终的基因列表文件 ---\n")

# 3.1 生成 all_genes.txt
# 背景基因集应该是所有经过低表达过滤后、实际参与了差异分析的基因
# 我们需要去除那些因为表达量太低或包含异常值而没有p值的基因
all_genes <- res_df %>%
  na.omit() %>%       # 去除含有NA值的行 (通常是p值为NA的)
  rownames()          # 提取行名，即基因ID

cat(paste("共计", length(all_genes), "个基因被用作背景基因集 (all_genes.txt)\n"))


# 3.2 生成 DEG_genes.txt
# 设定筛选阈值，例如：调整后的p值 < 0.05 并且 |log2FC| > 1
padj_threshold <- 0.05
logfc_threshold <- 1

DEG_genes <- res_df %>%
  na.omit() %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > logfc_threshold) %>%
  rownames()

cat(paste("根据阈值 (padj <", padj_threshold, ", |log2FC| >", logfc_threshold, ")，筛选出", length(DEG_genes), "个差异表达基因 (DEG_genes.txt)\n"))


# 3.3 保存文件
write.table(all_genes, "all_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(DEG_genes, "DEG_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\n成功生成 'all_genes.txt' 和 'DEG_genes.txt' 文件！\n")
cat("您现在拥有了进行富集分析所需的全部四个文件。\n")