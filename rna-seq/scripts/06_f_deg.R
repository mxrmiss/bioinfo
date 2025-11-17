#!/usr/bin/env Rscript
# 06_f_deg.R — DESeq2 差异表达（严格列名版：samples=sample,group；contrasts=contrast,case,control）

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(yaml)
  library(stringr)
})

# 小工具
trim <- function(x) { x <- as.character(x); trimws(x) }
getp <- function(lst, ..., default=NULL) {
  ref <- lst; keys <- list(...)
  for (k in keys) {
    if (is.null(ref)) return(default)
    if (!(is.list(ref) || is.environment(ref))) return(default)
    v <- ref[[as.character(k)]]
    if (is.null(v)) return(default)
    ref <- v
  }
  ref
}
norm_path <- function(x, fallback) {
  if (is.null(x) || length(x)==0) return(as.character(fallback))
  x <- as.character(x)[1]; if (!nzchar(x)) return(as.character(fallback)); x
}
norm_num <- function(x, fallback) {
  if (is.null(x) || length(x)==0) return(as.numeric(fallback))
  y <- suppressWarnings(as.numeric(x[1])); if (is.na(y)) return(as.numeric(fallback)); y
}

# 读配置（严格不推断）
cfg_path <- "config.yaml"
if (!file.exists(cfg_path)) stop("缺少 config.yaml（请在项目根运行）")
cfg <- tryCatch(yaml::read_yaml(cfg_path), error=function(e) NULL)
if (is.null(cfg) || !(is.list(cfg) || is.environment(cfg))) stop("config.yaml 解析失败")

matrix_dir  <- norm_path(getp(cfg,"paths","matrix_dir"),  "results/matrix")
deg_dir     <- norm_path(getp(cfg,"paths","deg_dir"),     "results/deg")
samples_tsv <- norm_path(getp(cfg,"tables","samples"),    "data/samples.tsv")
contr_tsv   <- norm_path(getp(cfg,"tables","contrasts"),  "data/contrasts.tsv")
alpha       <- norm_num(getp(cfg,"deg","alpha", 0.05), 0.05)
lfc_thresh  <- norm_num(getp(cfg,"deg","lfc",   1.0),  1.0)

message("[INFO] Resolved paths/params:")
message("  - matrix_dir = ", matrix_dir)
message("  - deg_dir    = ", deg_dir)
message("  - samples    = ", samples_tsv)
message("  - contrasts  = ", contr_tsv)
message("  - alpha      = ", alpha)
message("  - lfc        = ", lfc_thresh)

ok1 <- try(dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
if (inherits(ok1,"try-error")) stop("[ERR] 无法创建 deg_dir：", deg_dir)
ok2 <- try(dir.create(file.path(deg_dir,"tables"), recursive = TRUE, showWarnings = FALSE), silent = TRUE)
if (inherits(ok2,"try-error")) stop("[ERR] 无法创建结果子目录：", file.path(deg_dir,"tables"))

# 读计数并 round
counts_fp <- file.path(matrix_dir, "gene_counts.tsv")
if (!file.exists(counts_fp)) stop("[ERR] 缺少计数矩阵：", counts_fp)
counts <- readr::read_tsv(counts_fp, show_col_types = FALSE)
if (!"gene_id" %in% names(counts)) stop("[ERR] gene_counts.tsv 需包含列：gene_id")

gene_ids <- counts$gene_id
mat <- counts %>% select(-gene_id) %>% as.matrix()
mode(mat) <- "numeric"; mat <- round(mat, 0)
rownames(mat) <- gene_ids

# 读 samples.tsv（严格列名）
if (!file.exists(samples_tsv)) stop("[ERR] 缺少样本表：", samples_tsv)
samples <- readr::read_tsv(samples_tsv, show_col_types = FALSE)
need_s <- c("sample","group")
miss_s <- setdiff(need_s, names(samples))
if (length(miss_s)) stop("[ERR] samples.tsv 缺少必需列：", paste(miss_s, collapse=", "),
                         "（必须严格为列名：sample, group）")
samples$sample <- trim(samples$sample)
samples$group  <- trim(samples$group)

# 对齐顺序
if (!all(colnames(mat) %in% samples$sample)) {
  miss <- setdiff(colnames(mat), samples$sample)
  stop("[ERR] 样本表缺少以下样本：", paste(miss, collapse=", "))
}
samples <- samples[match(colnames(mat), samples$sample), , drop=FALSE]
grp <- factor(samples$group)
coldata <- data.frame(row.names = colnames(mat), group = grp, check.names = FALSE)

# 构建与拟合
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ group)
dds <- DESeq(dds)

# 读 contrasts.tsv（严格列名）
if (!file.exists(contr_tsv)) stop("[ERR] 缺少对比表：", contr_tsv)
contr <- readr::read_tsv(contr_tsv, show_col_types = FALSE)
need_c <- c("contrast","case","control")
miss_c <- setdiff(need_c, names(contr))
if (length(miss_c)) stop("[ERR] contrasts.tsv 缺少必需列：", paste(miss_c, collapse=", "),
                          "（必须严格为列名：contrast, case, control）")
contr$contrast <- trim(contr$contrast)
contr$case     <- trim(contr$case)
contr$control  <- trim(contr$control)

# 校验 group 水平
levs <- levels(grp)
bad <- !(contr$case %in% levs & contr$control %in% levs)
if (any(bad)) {
  bad_rows <- paste0(contr$contrast[bad], " (", contr$case[bad], " vs ", contr$control[bad], ")")
  stop("[ERR] contrasts.tsv 中存在 group 不在样本分组水平内的对比：\n - ",
       paste(bad_rows, collapse="\n - "),
       "\n样本分组水平：", paste(levs, collapse=", "))
}

# 逐对比输出
for (i in seq_len(nrow(contr))) {
  cname <- contr$contrast[i]; g1 <- contr$case[i]; g2 <- contr$control[i]
  message("[INFO] 运行对比：", cname, " (", g1, " vs ", g2, ")")
  res <- results(dds, contrast = c("group", g1, g2), alpha = alpha)
  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      padj = ifelse(is.na(padj), 1, padj),
      sig  = ifelse(padj <= alpha & abs(log2FoldChange) >= lfc_thresh,
                    ifelse(log2FoldChange > 0, "Up", "Down"),
                    "NS")
    )
  out_tsv <- file.path(deg_dir, "tables", paste0(cname, ".deseq2.tsv"))
  readr::write_tsv(res_tbl, out_tsv)
  message("[OK] 写出：", out_tsv)
}

# 汇总
summary_fp <- file.path(deg_dir, "DE_summary.tsv")
summ <- list.files(file.path(deg_dir, "tables"), pattern="\\.deseq2\\.tsv$", full.names=TRUE) %>%
  lapply(function(f){
    tb <- readr::read_tsv(f, show_col_types = FALSE)
    tibble(
      contrast = sub("\\.deseq2\\.tsv$", "", basename(f)),
      n_up   = sum(tb$sig=="Up",   na.rm=TRUE),
      n_down = sum(tb$sig=="Down", na.rm=TRUE),
      n_sig  = sum(tb$sig %in% c("Up","Down"), na.rm=TRUE)
    )
  }) %>% bind_rows()
readr::write_tsv(summ, summary_fp)
message("DESeq2 完成 ✅ 结果目录：", deg_dir)

