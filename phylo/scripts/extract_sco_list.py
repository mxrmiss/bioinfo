#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 Orthogroups.tsv 提取第一列（去表头）为 OG 列表，严格 TSV 解析。
"""
import argparse, sys, csv
from pathlib import Path

def log(m): sys.stderr.write(m.rstrip()+"\n")

def read_results_path(pf: Path) -> Path:
    return Path(pf.read_text(encoding="utf-8").strip())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--results-file", required=True)
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    resdir = read_results_path(Path(args.results_file))
    tsv = resdir/"Orthogroups/Orthogroups.tsv"
    if not tsv.exists(): log(f"[ERR] 缺少 Orthogroups.tsv: {tsv}"); return 3
    ogs=[]
    with tsv.open("r", encoding="utf-8", errors="replace", newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        header=True
        for row in rdr:
            if header: header=False; continue
            if not row: continue
            og=row[0].strip()
            if og: ogs.append(og)
    outp=Path(args.out); outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text("\n".join(ogs)+"\n", encoding="utf-8")
    log(f"[OK] OG count: {len(ogs)} → {outp}")
    return 0

if __name__=="__main__":
    sys.exit(main())

