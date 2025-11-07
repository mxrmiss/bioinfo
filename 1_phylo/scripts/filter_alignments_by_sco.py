#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_alignments_by_sco.py
功能：根据严格SCO白名单，从 collapse_to_species 目录筛出存在的对齐文件，
在一个新目录中创建指向它们的软链接，供后续拼接仅扫描 SCO 位点。
"""

import argparse
import sys
from pathlib import Path
import os

def parse_args():
    ap = argparse.ArgumentParser(description="Link alignments by SCO keep list")
    ap.add_argument("--indir", required=True,
                    help="collapse_to_species 目录，如 results/trim_norm_collapse_species")
    ap.add_argument("--keep-list", required=True,
                    help="严格SCO白名单文件，每行一个OG，如 results/orthogroups/sco.keep.list")
    ap.add_argument("--outdir", required=True,
                    help="输出目录，如 results/trim_norm_sco_species")
    return ap.parse_args()

def main():
    args = parse_args()
    indir = Path(args.indir)
    keep = Path(args.keep_list)
    outdir = Path(args.outdir)

    if not indir.is_dir():
        print(f"[ERR] 输入目录不存在: {indir}", file=sys.stderr)
        sys.exit(2)
    if not (keep.is_file() and keep.stat().st_size > 0):
        print(f"[ERR] 白名单文件不存在或为空: {keep}", file=sys.stderr)
        sys.exit(3)

    outdir.mkdir(parents=True, exist_ok=True)

    with keep.open() as fh:
        ogs = [ln.strip() for ln in fh if ln.strip()]

    kept = 0
    missed = 0
    total = len(ogs)

    for og in ogs:
        src = indir / f"{og}.trim.faa"
        dst = outdir / f"{og}.trim.faa"
        if src.is_file() and src.stat().st_size > 0:
            # 绝对路径软链接，更稳妥
            if dst.exists() or dst.is_symlink():
                try:
                    dst.unlink()
                except Exception:
                    pass
            os.symlink(src.resolve(), dst)
            kept += 1
        else:
            # 源缺失/为空，跳过并统计
            if dst.is_symlink():
                try:
                    dst.unlink()
                except Exception:
                    pass
            missed += 1

    print(f"[INFO] 严格SCO白名单总计: {total} | 成功链接: {kept} | 缺失/空: {missed}", file=sys.stderr)

if __name__ == "__main__":
    main()

