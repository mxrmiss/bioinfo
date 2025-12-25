# -*- coding: utf-8 -*-
# ess_report.R —— MCMCTree MCMC 轨迹 ESS 报告（Tracer风 IPS + coda）
# 硬编码路径 + 去除所有转义引号；不读配置、不读环境变量。

## ========= 硬编码路径（需要改目录时，只改这里）=========
WORK_DIR  <- "results/06_cafe/mcmctree"
QC_SUBDIR <- "qc_report"
## ================================================

MCMC_PATH <- file.path(WORK_DIR, "mcmc.txt")
OUT_DIR   <- file.path(WORK_DIR, QC_SUBDIR)
suppressWarnings(dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE))

tsv_path <- file.path(OUT_DIR, "ess.tsv")
md_path  <- file.path(OUT_DIR, "ess_summary.md")
rec_path <- file.path(OUT_DIR, "ess_recommendation.txt")

TARGET_ESS <- 250
TOP_K      <- 15
KEEP_REGEX <- "^(t_|sigma2|rgene|rates|kappa|alpha|lambda|mu)"

if (!file.exists(MCMC_PATH)) {
  stop(paste("找不到 mcmc.txt：", MCMC_PATH))
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("需要 R 包：data.table（install.packages('data.table')）")
}

# 读取与列筛选
x <- data.table::fread(MCMC_PATH, data.table = FALSE)
num_cols <- sapply(x, is.numeric)
x <- x[, num_cols, drop = FALSE]
keep <- grep(KEEP_REGEX, colnames(x), value = TRUE)
if (length(keep) == 0) {
  stop("没有匹配到可评估的参数列；请检查 KEEP_REGEX 与 mcmc.txt 列名")
}
M <- as.matrix(x[, keep, drop = FALSE])

# ESS 计算
ess_ips_one <- function(v){
  v <- as.numeric(v); v <- v[is.finite(v)]
  n <- length(v); if (n < 50) return(NA_real_)
  ac <- stats::acf(v, plot = FALSE, lag.max = min(1000, n - 1), demean = TRUE)$acf[-1]
  s <- 0
  for (k in seq(1, length(ac), by = 2)) {
    a1 <- ac[k]; a2 <- if (k + 1 <= length(ac)) ac[k + 1] else 0
    pair <- a1 + a2
    if (!is.finite(pair) || pair <= 0) break
    s <- s + pair
  }
  n / (1 + 2 * s)
}
calc_ess_ips  <- function(M) apply(M, 2, ess_ips_one)
calc_ess_coda <- function(M){
  if (!requireNamespace("coda", quietly = TRUE)) return(rep(NA_real_, ncol(M)))
  as.numeric(coda::effectiveSize(coda::as.mcmc(M)))
}

ESS_ips  <- suppressWarnings(calc_ess_ips(M))
ESS_coda <- suppressWarnings(calc_ess_coda(M))

# 结果表（列名修复）
res <- data.frame(
  param     = colnames(M),
  ESS       = as.numeric(ESS_ips),
  ESS_ips   = as.numeric(ESS_ips),
  ESS_coda  = as.numeric(ESS_coda),
  delta_pct = as.numeric(100 * (ESS_coda - ESS_ips) / ESS_ips),
  stringsAsFactors = FALSE
)
res <- res[order(res$ESS), ]
utils::write.table(res, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

# 概要统计与推荐
minESS <- suppressWarnings(min(res$ESS, na.rm = TRUE))
medESS <- suppressWarnings(median(res$ESS, na.rm = TRUE))
p25ESS <- suppressWarnings(as.numeric(stats::quantile(res$ESS, 0.25, na.rm = TRUE)))
n_bad  <- sum(is.finite(res$ESS) & res$ESS < TARGET_ESS)
mult   <- if (is.finite(minESS) && minESS > 0) ceiling(TARGET_ESS / minESS) else 1L
writeLines(as.character(mult), rec_path)

# 表格工具（不再使用任何转义引号）
mk_table <- function(df){
  hdr <- "| rank | param | ESS |\n|---:|:------|---:|\n"
  if (nrow(df) == 0) return(paste0(hdr, "| 1 | (none) | NA |\n"))
  rows <- paste0(
    "| ", seq_len(nrow(df)),
    " | ", df$param,
    " | ", sprintf("%.2f", df$ESS),
    " |",
    collapse = "\n"
  )
  paste0(hdr, rows, "\n")
}
mk_table_delta <- function(df){
  hdr <- "| rank | param | ESS_ips | ESS_coda | Δ% |\n|---:|:------|-------:|--------:|---:|\n"
  if (nrow(df) == 0) return(paste0(hdr, "| 1 | (none) | NA | NA | NA |\n"))
  rows <- paste0(
    "| ", seq_len(nrow(df)),
    " | ", df$param,
    " | ", sprintf("%.2f", df$ESS_ips),
    " | ", sprintf("%.2f", df$ESS_coda),
    " | ", sprintf("%+.1f", df$delta_pct),
    " |",
    collapse = "\n"
  )
  paste0(hdr, rows, "\n")
}

top_bad <- head(res[, c("param","ESS")], min(nrow(res), TOP_K))
delta_score <- ifelse(is.finite(res$delta_pct), abs(res$delta_pct), -Inf)
delta_ord   <- order(delta_score, decreasing = TRUE)
top_delta   <- head(res[delta_ord, c("param","ESS_ips","ESS_coda","delta_pct")],
                    min(nrow(res), TOP_K))

# 写 MD
lines <- c(
  "# ESS Report (IPS + coda)",
  "",
  paste0("- Hard-wired work_dir: `", WORK_DIR, "`"),
  paste0("- File: `", MCMC_PATH, "`"),
  paste0("- Parameters evaluated: **", nrow(res), "**"),
  paste0("- Target ESS: **", TARGET_ESS, "**"),
  "",
  "## Overview",
  paste0("- Min ESS (IPS): **", sprintf("%.2f", minESS), "**"),
  paste0("- Median ESS (IPS): **", sprintf("%.2f", medESS), "**"),
  paste0("- 25% quantile (IPS): **", sprintf("%.2f", p25ESS), "**"),
  paste0("- Count < Target (IPS): **", n_bad, "**"),
  "",
  "## Worst parameters (Top-K by ESS, IPS)",
  mk_table(top_bad),
  "## Method discrepancy (Top-K by |Δ%|)",
  mk_table_delta(top_delta),
  "## Recommendation",
  paste0("- Suggested `nsample` multiplier: **×", mult, "**  (using IPS min ESS)"),
  "- Prefer two independent chains; IPS (Tracer风) 达标 + trace 平稳更审稿友好。",
  "",
  "_Artifacts_: `ess.tsv`（含 ESS 明细）, `ess_recommendation.txt`。"
)
writeLines(paste(lines, collapse = "\n"), md_path)

cat("[OK] Hard-wired work_dir: ", WORK_DIR, "\n", sep = "")
cat("[OK] ESS report written:\n", md_path, "\n", tsv_path, "\n", rec_path, "\n", sep = "")

