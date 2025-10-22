#!/usr/bin/env Rscript
################################################################################
#  脚本 1: 从 EggNOG-mapper 结果中生成注释文件（同级目录版本）
#
# 变更说明（相对旧版）：
#   - 取消 inmaterial/oumaterial 子文件夹约定；
#   - 自动定位“脚本自身所在目录”为输入/输出目录（同级目录）；
#   - 无论在 Rscript / RStudio / source() 下运行，均尽量稳健获取脚本目录。
################################################################################

# ==============================================================================
# 步骤 0: 安装和加载必要的 R 包
# ==============================================================================
packages_to_install <- c("dplyr", "tidyr", "stringr")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(dplyr)
library(tidyr)
library(stringr)

cat("所有必要的 R 包已成功加载！\n")

# ==============================================================================
# 工具函数：稳健地获取脚本所在目录（支持 Rscript / RStudio / source）
# ==============================================================================
get_script_dir <- function() {
  # 1) Rscript 方式：从命令行参数 --file=... 抓路径
  args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("^--file=", args)
  if (length(file_idx) > 0) {
    f <- sub("^--file=", "", args[file_idx[1]])
    return(dirname(normalizePath(f)))
  }
  # 2) 被 source() 调用：从调用栈获取
  if (!is.null(sys.frames()) && length(sys.frames()) > 0) {
    frame_files <- lapply(sys.frames(), function(x) x$ofile)
    frame_files <- Filter(Negate(is.null), frame_files)
    if (length(frame_files) > 0) {
      return(dirname(normalizePath(frame_files[[length(frame_files)]])))
    }
  }
  # 3) RStudio 运行：尝试使用 rstudioapi
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  # 4) 兜底：当前工作目录（若无法探测脚本路径）
  getwd()
}

# ==============================================================================
#                         ----------- 用户配置区 -----------
# ==============================================================================
# 请将这里的文件名改为您的 EggNOG-mapper 下载文件名（与本脚本同一目录）
EGGNOG_ANNOTATION_FILE_BASENAME <- "MM_2z45cx4s.emapper.annotations.tsv"

# ==============================================================================
# 步骤 1: 组装“同级目录”路径并读取文件
# ==============================================================================
SCRIPT_DIR <- get_script_dir()
INPUT_FILE <- file.path(SCRIPT_DIR, EGGNOG_ANNOTATION_FILE_BASENAME)

cat("脚本所在目录：", SCRIPT_DIR, "\n", sep = "")
cat("将从同级目录读取：", INPUT_FILE, "\n", sep = "")

if (!file.exists(INPUT_FILE)) {
  stop(paste0("错误：未在脚本同级目录找到输入文件：\n  ", INPUT_FILE,
              "\n请确认文件名与位置正确。"))
}

cat("正在读取并解析 EggNOG-mapper 文件...\n")

# 注意：EggNOG 输出常含以 "#" 开头的注释行，这里通过 comment.char 跳过
eggnog_data <- read.delim(
  INPUT_FILE,
  header = FALSE,          # 通常无标准表头
  comment.char = "#",      # 跳过注释
  stringsAsFactors = FALSE,
  sep = "\t"
)

# 为数据框添加必要列名（按 EggNOG v2 常见列位：1=GeneID, 6=GOs, 8=KEGG_ko）
# 若列不足会抛错，避免误解析
if (ncol(eggnog_data) < 8) {
  stop(paste0("错误：检测到输入文件列数不足（ncol=", ncol(eggnog_data),
              "），请核对是否为完整的 .emapper.annotations 文件。"))
}
colnames(eggnog_data)[c(1, 6, 8)] <- c("GeneID", "GOs", "KEGG_ko")

# ==============================================================================
# 步骤 2: 提取 GO 与 KEGG KO 对应关系
# ==============================================================================
cat("正在提取 GO 注释...\n")
gene_to_go <- eggnog_data %>%
  dplyr::select(GeneID, GOs) %>%
  filter(!is.na(GOs) & GOs != "") %>%
  separate_rows(GOs, sep = ",") %>%
  rename(GO_ID = GOs) %>%
  distinct()

cat("正在提取 KEGG 注释...\n")
gene_to_kegg <- eggnog_data %>%
  dplyr::select(GeneID, KEGG_ko) %>%
  filter(!is.na(KEGG_ko) & KEGG_ko != "") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KO_ID = str_remove(KEGG_ko, "^ko:")) %>%  # 去掉前缀 "ko:"
  dplyr::select(GeneID, KO_ID) %>%
  distinct()

cat(sprintf("摘要：GO 关联 %d 条，KEGG 关联 %d 条。\n",
            nrow(gene_to_go), nrow(gene_to_kegg)))

# ==============================================================================
# 步骤 3: 结果写回“脚本同级目录”
# ==============================================================================
out_go   <- file.path(SCRIPT_DIR, "parsed_gene_go_annotation.tsv")
out_kegg <- file.path(SCRIPT_DIR, "parsed_gene_kegg_annotation.tsv")

write.table(gene_to_go, out_go, sep = "\t", row.names = FALSE, quote = FALSE)
cat("GO 注释文件已保存：", out_go, "\n", sep = "")

write.table(gene_to_kegg, out_kegg, sep = "\t", row.names = FALSE, quote = FALSE)
cat("KEGG 注释文件已保存：", out_kegg, "\n", sep = "")

cat("脚本 '01_Parse_EggNOG_to_Annotation_same_dir.R' 运行完毕！\n")
