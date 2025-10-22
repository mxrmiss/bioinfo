#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(readr) })

outdir <- "results/matrix"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

files <- Sys.glob("results/quant/*/quant.genes.sf")
stopifnot(length(files) > 0)
samples <- basename(dirname(files))
names(files) <- samples

gene_union <- NULL
counts_list <- list()
tpm_list <- list()

for (s in samples) {
  df <- readr::read_tsv(files[s], show_col_types=FALSE)
  # Name Length EffectiveLength TPM NumReads
  counts_list[[s]] <- setNames(df$NumReads, df$Name)
  tpm_list[[s]]    <- setNames(df$TPM,      df$Name)
  gene_union <- union(gene_union, df$Name)
}

mat_counts <- do.call(cbind, lapply(samples, function(s) counts_list[[s]][gene_union]))
mat_tpm    <- do.call(cbind, lapply(samples, function(s) tpm_list[[s]][gene_union]))
colnames(mat_counts) <- samples; rownames(mat_counts) <- gene_union
colnames(mat_tpm)    <- samples; rownames(mat_tpm)    <- gene_union

mat_counts[is.na(mat_counts)] <- 0
mat_tpm[is.na(mat_tpm)] <- 0

readr::write_tsv(data.frame(gene_id=rownames(mat_counts), mat_counts, check.names=FALSE),
                 file.path(outdir,"gene_counts.tsv"))
readr::write_tsv(data.frame(gene_id=rownames(mat_tpm),    mat_tpm,    check.names=FALSE),
                 file.path(outdir,"gene_tpm.tsv"))

cat("✅ 合并完成：", file.path(outdir,"gene_counts.tsv"), "与", file.path(outdir,"gene_tpm.tsv"), "\n")
