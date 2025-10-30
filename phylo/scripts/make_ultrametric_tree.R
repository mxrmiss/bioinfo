#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 等时树生成：仅依赖 ape；内置简易 midpoint rooting；强制分支长度≥eps

suppressPackageStartupMessages(library(ape))

# ---- 轻量参数解析（不用 optparse） ----
args <- commandArgs(trailingOnly = TRUE)
getArg <- function(key, default=NULL){
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i+1] else default
}
fin  <- getArg("--in")
fout <- getArg("--out")
lam  <- as.numeric(getArg("--lambda", "1.0"))

if (is.null(fin) || is.null(fout)) {
  stop("用法: Rscript make_ultrametric_tree.R --in <in.nwk> --out <out.nwk> [--lambda 1.0]")
}

# ---- 简易 midpoint rooting（不依赖 phytools）----
midpoint_root_simple <- function(tr){
  D  <- cophenetic(tr)
  ij <- which(D == max(D), arr.ind = TRUE)[1, ]
  tip1 <- rownames(D)[ij[1]]
  tip2 <- rownames(D)[ij[2]]
  root(tr, outgroup = c(tip1, tip2), resolve.root = TRUE)
}

# ---- 主流程 ----
tr <- read.tree(fin)
if (!is.rooted(tr)) {
  ok <- TRUE
  tryCatch({ tr <- midpoint_root_simple(tr) }, error=function(e) ok<<-FALSE)
  if (!ok || !is.rooted(tr)) {
    tr <- root(tr, outgroup = tr$tip.label[1], resolve.root = TRUE)
  }
}

# 等时化（惩罚似然；lambda 可从命令行传入）
utr <- chronos(tr, lambda = lam, quiet = TRUE)

# ---- 关键修复：把过小/非正的边长抬到极小正数 eps ----
eps <- 1e-5  # 若 CAFÉ 仍提示无效分支长度，可把 eps 提到 1e-4
bl  <- utr$edge.length
bl[is.na(bl) | bl < eps] <- eps
utr$edge.length <- bl

# 写出，保留足够小数位
write.tree(utr, file = fout, digits = 10)

