# -*- coding: utf-8 -*-
MCMC_PATH   <- "mcmctree/mcmc.txt"
OUT_DIR     <- "mcmctree/qc_report"
TARGET_ESS  <- 250
TOP_K       <- 15
KEEP_REGEX  <- "^(t_|sigma2|rgene|rates|kappa|alpha|lambda|mu)"

suppressWarnings(dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE))
tsv_path  <- file.path(OUT_DIR, "ess.tsv")
md_path   <- file.path(OUT_DIR, "ess_summary.md")
rec_path  <- file.path(OUT_DIR, "ess_recommendation.txt")

read_dt <- function(p){
  if (!requireNamespace("data.table", quietly = TRUE)) stop("需要 data.table")
  data.table::fread(p, data.table = FALSE)
}

ess_ips <- function(v){
  v <- as.numeric(v); v <- v[is.finite(v)]
  n <- length(v); if (n < 10) return(NA_real_)
  ac <- acf(v, plot = FALSE, lag.max = min(1000, n - 1), demean = TRUE)$acf[-1]
  kmax <- length(ac) %/% 2; s <- 0
  for (k in seq_len(kmax)) {
    pair <- ac[2*k - 1] + ac[2*k]
    if (!is.finite(pair) || pair <= 0) break
    s <- s + pair
  }
  ess <- n / (1 + 2*s)
  max(1, ess)
}

calc_ess <- function(M){
  if (requireNamespace("coda", quietly = TRUE)) {
    as.numeric(coda::effectiveSize(coda::as.mcmc(M)))
  } else {
    apply(M, 2, ess_ips)
  }
}

x <- read_dt(MCMC_PATH)
num_cols <- sapply(x, is.numeric)
x <- x[, num_cols, drop = FALSE]
keep <- grep(KEEP_REGEX, colnames(x), value = TRUE)
if (length(keep) == 0) stop("没有匹配到可评估的参数列；请检查 KEEP_REGEX")
M <- as.matrix(x[, keep, drop = FALSE])

ess <- calc_ess(M)
res <- data.frame(param = colnames(M), ESS = as.numeric(ess), row.names = NULL)
res <- res[order(res$ESS), ]
utils::write.table(res, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

minESS <- suppressWarnings(min(res$ESS, na.rm = TRUE))
medESS <- suppressWarnings(median(res$ESS, na.rm = TRUE))
p25ESS <- suppressWarnings(as.numeric(stats::quantile(res$ESS, 0.25, na.rm = TRUE)))
n_bad  <- sum(is.finite(res$ESS) & res$ESS < TARGET_ESS)
mult   <- if (is.finite(minESS) && minESS > 0) ceiling(TARGET_ESS / minESS) else 1L
writeLines(as.character(mult), rec_path)

mk_table <- function(df){
  hdr <- "| rank | param | ESS |\n|---:|:------|---:|\n"
  if (nrow(df) == 0) return(paste0(hdr, "| 1 | (none) | NA |\n"))
  rows <- paste0("| ", seq_len(nrow(df)), " | ", df$param, " | ", sprintf("%.2f", df$ESS), " |", collapse = "\n")
  paste0(hdr, rows, "\n")
}

top_bad <- head(res, min(nrow(res), TOP_K))

lines <- c(
  "# ESS Report",
  "",
  paste0("- File: `", MCMC_PATH, "`"),
  paste0("- Parameters evaluated: **", nrow(res), "**"),
  paste0("- Target ESS: **", TARGET_ESS, "**"),
  "",
  "## Overview",
  paste0("- Min ESS: **", sprintf("%.2f", minESS), "**"),
  paste0("- Median ESS: **", sprintf("%.2f", medESS), "**"),
  paste0("- 25% quantile: **", sprintf("%.2f", p25ESS), "**"),
  paste0("- Count < Target: **", n_bad, "**"),
  "",
  "## Worst parameters (Top-K by ESS)",
  mk_table(top_bad),
  "## Recommendation",
  paste0("- Suggested `nsample` multiplier: **×", mult, "**"),
  "- If many parameters < Target, consider running a second chain (different seed) and tuning `finetune` to keep acceptance in ~0.2–0.4.",
  "- Do not rely on larger `sampfreq` (thinning) to increase ESS; it only sparsifies logs.",
  "",
  "_Artifacts_: `ess.tsv`, `ess_recommendation.txt`."
)

writeLines(paste(lines, collapse = "\n"), md_path)
cat("[OK] ESS report written:\n", md_path, "\n", tsv_path, "\n", rec_path, "\n", sep = "")
