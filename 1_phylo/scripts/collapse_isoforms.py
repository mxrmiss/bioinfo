#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collapse_isoforms.py
=====================
作用：
  - 从 trim_norm/*.faa 文件中，按物种 ID 折叠 isoforms
  - 每个物种在每个 OG (locus) 里只保留一条最长序列
  - 输出到指定目录，供拼接 supermatrix 使用

用法：
  python collapse_isoforms.py --indir results/trim_norm --outdir results/trim_norm_collapse
"""

import os
import argparse
from pathlib import Path
from Bio import SeqIO

def collapse_faa(infile, outfile):
    records = list(SeqIO.parse(infile, "fasta"))
    best = {}  # {speciesID: record}

    for rec in records:
        # header 格式是 SpeciesID|SeqID
        species = rec.id.split("|")[0]
        if species not in best or len(rec.seq) > len(best[species].seq):
            best[species] = rec

    SeqIO.write(best.values(), outfile, "fasta")

def main():
    ap = argparse.ArgumentParser(description="Collapse isoforms per species (keep longest).")
    ap.add_argument("--indir", required=True, help="输入目录 (trim_norm/*.faa)")
    ap.add_argument("--outdir", required=True, help="输出目录")
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n = 0
    for f in indir.glob("*.faa"):
        outfile = outdir / f.name
        collapse_faa(f, outfile)
        n += 1
    print(f"[OK] 处理完成：共 {n} 个文件写入 {outdir}")

if __name__ == "__main__":
    main()

