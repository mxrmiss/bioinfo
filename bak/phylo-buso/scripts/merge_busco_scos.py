#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_busco_scos.py（精简稳妥版）
- 从多个 BUSCO 运行目录中收集“单拷贝”蛋白序列，合并为“每个 BUSCO 一份多物种 FASTA”。
- 兼容 v5/v6 常见布局与 .gz。
- 仅需三个参数：--busco_runs --outdir --mincov
- 额外输出 _coverage.tsv 便于诊断（可选但很轻量）。

用法：
  python scripts/merge_busco_scos.py \
      --busco_runs results/busco \
      --outdir     results/orthogroups \
      --mincov     6
"""

import os
import re
import sys
import glob
import gzip
import argparse
from collections import defaultdict, OrderedDict

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode, encoding=None if "b" in mode else "utf-8", errors="ignore")
    return open(path, mode=mode, encoding=None if "b" in mode else "utf-8", errors="ignore")

def read_fasta(path):
    recs, hdr, seq = [], None, []
    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if hdr is not None:
            recs.append((hdr, "".join(seq)))
    return recs

def looks_like_protein(seq):
    seq = seq.upper()
    aa_chars = set("ACDEFGHIKLMNPQRSTVWYBXZ*")
    dna_chars = set("ACGTN")
    non_dna = sum(1 for c in seq if c not in dna_chars and c.isalpha())
    return non_dna >= max(1, int(0.05 * len(seq)))

def find_single_copy_dirs(run_root):
    # 常见结构：busco_sequences/single_copy_busco_sequences/ 或 translated_proteins/single_copy_busco_sequences/
    candidates = []
    candidates += glob.glob(os.path.join(run_root, "busco_sequences", "single_copy_busco_sequences"))
    candidates += glob.glob(os.path.join(run_root, "translated_proteins", "single_copy_busco_sequences"))
    candidates += glob.glob(os.path.join(run_root, "**", "single_copy_busco_sequences"), recursive=True)
    seen, uniq = set(), []
    for d in candidates:
        if os.path.isdir(d) and d not in seen:
            uniq.append(d)
            seen.add(d)
    return uniq

def guess_species_id_from_dir(path):
    m = re.search(r"/([^/]+)_busco(?:/|$)", path.replace("\\", "/"))
    if m:
        return m.group(1)
    return os.path.basename(os.path.dirname(path))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--busco_runs", required=True, help="BUSCO 结果父目录（包含 <sid>_busco 子目录）")
    ap.add_argument("--outdir",     required=True, help="输出 orthogroups 目录")
    ap.add_argument("--mincov",     type=int, required=True, help="每 BUSCO 至少覆盖的物种数")
    args = ap.parse_args()

    busco_root = os.path.abspath(args.busco_runs)
    outdir     = os.path.abspath(args.outdir)
    mincov     = args.mincov
    os.makedirs(outdir, exist_ok=True)

    # 1) 列出各物种最近一次 run_* 目录
    run_dirs = []
    for sid_dir in glob.glob(os.path.join(busco_root, "*_busco")):
        runs = sorted(glob.glob(os.path.join(sid_dir, "run_*")), key=os.path.getmtime)
        if runs:
            run_dirs.append(runs[-1])
    if not run_dirs:
        print(f"[ERROR] 未在 {busco_root} 下找到任何 run_* 目录。", file=sys.stderr)
        sys.exit(2)

    # 2) 收集 BUSCO 单拷贝蛋白
    combined = defaultdict(OrderedDict)  # BUSCO_ID -> OrderedDict(species_id -> (hdr, seq))
    for run in run_dirs:
        species_id = guess_species_id_from_dir(run)
        picked = None
        for d in find_single_copy_dirs(run):
            files = sum([glob.glob(os.path.join(d, f"*{ext}"))
                         for ext in (".faa", ".fa", ".fasta", ".faa.gz", ".fa.gz", ".fasta.gz")], [])
            if files:
                picked = d
                break
        if not picked:
            print(f"[WARN] {species_id}: 未找到单拷贝蛋白目录，跳过", file=sys.stderr)
            continue

        for fpath in sum([glob.glob(os.path.join(picked, f"*{ext}"))
                          for ext in (".faa", ".fa", ".fasta", ".faa.gz", ".fa.gz", ".fasta.gz")], []):
            recs = read_fasta(fpath)
            if not recs:
                continue
            busco_id = os.path.splitext(os.path.basename(fpath))[0]
            busco_id = re.sub(r"\.(fa|faa|fasta)$", "", busco_id, flags=re.IGNORECASE)
            busco_id = re.sub(r"\.(gz)$", "", busco_id, flags=re.IGNORECASE)

            hdr, seq = recs[0]
            if not seq or not looks_like_protein(seq):
                continue
            gid = hdr.split()[0]
            new_hdr = f"{species_id}|{gid}"
            if species_id not in combined[busco_id]:
                combined[busco_id][species_id] = (new_hdr, seq)

    # 3) 覆盖度过滤与输出
    kept, dropped = 0, 0
    cov_hist = defaultdict(int)
    for busco_id, spmap in combined.items():
        cov = len(spmap)
        cov_hist[cov] += 1
        if cov >= mincov:
            out_fa = os.path.join(outdir, f"{busco_id}.faa")
            with open(out_fa, "w", encoding="utf-8") as out:
                for _, (hdr, seq) in spmap.items():
                    out.write(f">{hdr}\n")
                    for i in range(0, len(seq), 60):
                        out.write(seq[i:i+60] + "\n")
            kept += 1
        else:
            dropped += 1

    # 4) 覆盖度统计（轻量输出）
    cov_path = os.path.join(outdir, "_coverage.tsv")
    with open(cov_path, "w", encoding="utf-8") as fo:
        fo.write("coverage\tcount\n")
        for k in sorted(cov_hist):
            fo.write(f"{k}\t{cov_hist[k]}\n")

    print(f"[汇总] 通过覆盖度(≥{mincov})的 BUSCO 位点数: {kept}", file=sys.stderr)
    print(f"[提示] 覆盖度直方图见: {cov_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
