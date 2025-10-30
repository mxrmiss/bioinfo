#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_sco_keep_list.py
功能：生成“严格单拷贝（SCO）”白名单（每行一个 OG ID）
优先使用 OrthoFinder 自带的 Orthogroups_SingleCopyOrthologues.txt；
若缺失则从 Orthogroups.tsv 计算（每物种恰好1条且非空）。
"""

import argparse
import sys
from pathlib import Path

def parse_args():
    ap = argparse.ArgumentParser(description="Build strict SCO keep list")
    ap.add_argument("--results-file", required=True,
                    help="results/orthogroups/RESULTS_PATH.txt，内容是一行 OrthoFinder Results_* 目录的绝对路径")
    ap.add_argument("--out", required=True,
                    help="输出文件路径，如 results/orthogroups/sco.keep.list")
    return ap.parse_args()

def main():
    args = parse_args()
    resdir_file = Path(args.results_file)
    if not resdir_file.is_file():
        print(f"[ERR] RESULTS_PATH.txt 不存在: {resdir_file}", file=sys.stderr)
        sys.exit(2)

    resdir = Path(resdir_file.read_text().strip())
    if not resdir.exists():
        print(f"[ERR] OrthoFinder 结果目录不存在: {resdir}", file=sys.stderr)
        sys.exit(3)

    sco_txt = resdir / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt"
    tsv = resdir / "Orthogroups" / "Orthogroups.tsv"

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    ogs = []

    if sco_txt.is_file() and sco_txt.stat().st_size > 0:
        # 直接读取官方 SCO 清单（形如：OG0000123: geneA ...）
        with sco_txt.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                og = line.split(":", 1)[0].strip()
                if og:
                    ogs.append(og)
        print(f"[INFO] 使用 OrthoFinder 严格SCO清单: {len(ogs)}", file=sys.stderr)
    else:
        # 兜底：从 Orthogroups.tsv 计算严格 SCO（每物种恰好1条且非空）
        if not (tsv.is_file() and tsv.stat().st_size > 0):
            print(f"[ERR] 缺少 Orthogroups.tsv: {tsv}", file=sys.stderr)
            sys.exit(4)
        with tsv.open() as fh:
            header = fh.readline()
            if not header:
                print("[ERR] Orthogroups.tsv 为空", file=sys.stderr)
                sys.exit(5)
            ncol = len(header.rstrip("\n").split("\t"))
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                # 列1是 OG，2..n 是各物种的基因列表（逗号分隔）
                if len(parts) < ncol:
                    # 宽松处理：列数不足直接跳过
                    continue
                ok = True
                for cell in parts[1:]:
                    cell = cell.strip()
                    if cell == "" or "," in cell:
                        ok = False
                        break
                if ok:
                    ogs.append(parts[0].strip())
        print(f"[INFO] 自 Orthogroups.tsv 计算严格SCO清单: {len(ogs)}", file=sys.stderr)

    # 去重并写出
    seen = set()
    with out.open("w") as fo:
        for og in ogs:
            if og and og not in seen:
                fo.write(f"{og}\n")
                seen.add(og)

    if out.stat().st_size == 0:
        print("[WARN] 生成的严格SCO清单为空，请检查数据质量/物种列是否缺失", file=sys.stderr)

if __name__ == "__main__":
    main()

