#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import csv
import shutil
import sys
from pathlib import Path
from typing import Dict, Any, List

# 让“项目根目录”进入 sys.path，保证可以 import scripts._common
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts._common import (
    find_project_root,
    read_yaml_config,
    parse_rel_from_01_db,
    get_species_list,
    ensure_dir,
)

def read_kv_tsv(p: Path) -> Dict[str, str]:
    """
    读取形如：
      KEY\tVALUE
    的 TSV（无表头），返回 dict。文件不存在则返回空 dict。
    """
    d: Dict[str, str] = {}
    if not p.is_file():
        return d
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t", 1)
        if len(parts) == 2:
            d[parts[0]] = parts[1]
    return d

def copy_species_annotation(annot_src: Path, outdir: Path, species_id: str) -> Path:
    """
    把每个物种最终注释表（04_annot产物）原样复制到 99_reports 目录。
    输出命名：
      {outdir}/{species_id}.swissprot_annotation.tsv
    """
    dst = outdir / f"{species_id}.swissprot_annotation.tsv"
    ensure_dir(outdir)

    if not annot_src.is_file() or annot_src.stat().st_size == 0:
        # 不强制报错，避免因为某个物种没产物导致 99 步骤整体失败
        # 但会在屏幕提示，方便你排查上游是否跑完
        print(f"[99_reports] WARN missing/empty annot: {annot_src}", flush=True)
        return dst

    # 原样复制（保留时间戳/权限等元信息）
    shutil.copy2(annot_src.as_posix(), dst.as_posix())
    return dst

def main() -> None:
    root = find_project_root()
    os.chdir(root)

    cfg = read_yaml_config(root)
    rel = parse_rel_from_01_db(root)
    species = get_species_list(cfg)

    outdir = root / "results" / "99_reports" / rel
    ensure_dir(outdir)

    out_species = outdir / "species_summary.tsv"
    out_long = outdir / "all_annotations.long.tsv"

    print(f"[99_reports] REL={rel}", flush=True)
    print(f"[99_reports] N_species={len(species)}", flush=True)
    print(f"[99_reports] OUT={out_species}", flush=True)
    print(f"[99_reports] OUT_LONG={out_long}", flush=True)
    print(f"[99_reports] COPY per-species annot -> {outdir}", flush=True)

    # -----------------------------
    # 1) 生成 species_summary.tsv
    # -----------------------------
    header = [
        "species_id", "REL",
        "N_QUERY", "N_FILTERED_HITS", "N_BESTHIT", "HIT_RATE",
        "N_annot", "NO_OS", "NO_GN", "NO_OS_rate", "NO_GN_rate",
    ]
    rows: List[List[str]] = []

    # -----------------------------
    # 2) 生成 all_annotations.long.tsv
    # -----------------------------
    long_header: List[str] = []
    long_rows: List[List[str]] = []

    for sid in species:
        print(f"[99_reports] collecting {sid}", flush=True)

        summ = read_kv_tsv(root / "results" / "03_filter" / rel / sid / "summary.tsv")
        qc = read_kv_tsv(root / "results" / "04_annot" / rel / sid / "qc_missing_fields.tsv")

        n_query = summ.get("N_QUERY", "0")
        n_filt = summ.get("N_FILTERED_HITS", "0")
        n_best = summ.get("N_BESTHIT", "0")
        hit_rate = summ.get("HIT_RATE", "0")

        n_annot = qc.get("N_annot", "0")
        no_os = qc.get("NO_OS", "0")
        no_gn = qc.get("NO_GN", "0")
        no_os_rate = qc.get("NO_OS_rate", "0")
        no_gn_rate = qc.get("NO_GN_rate", "0")

        rows.append([
            sid, rel,
            n_query, n_filt, n_best, hit_rate,
            n_annot, no_os, no_gn, no_os_rate, no_gn_rate,
        ])

        # 物种注释表路径（04_annot产物）
        annot = root / "results" / "04_annot" / rel / sid / "swissprot_annotation.tsv"

        # 新增：复制到 99_reports/{REL}/ 下
        copied = copy_species_annotation(annot, outdir, sid)
        if copied.is_file() and copied.stat().st_size > 0:
            print(f"[99_reports] copied -> {copied}", flush=True)

        # 原有逻辑：用于拼 all_annotations.long.tsv
        if annot.is_file() and annot.stat().st_size > 0:
            with annot.open("r", encoding="utf-8", errors="replace") as fin:
                reader = csv.reader(fin, delimiter="\t")
                hdr = next(reader, None)
                if hdr is None:
                    continue
                if not long_header:
                    long_header = ["species_id"] + hdr
                for r in reader:
                    if not r:
                        continue
                    long_rows.append([sid] + r)

    # 写 species_summary.tsv
    with out_species.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(r)

    # 写 all_annotations.long.tsv
    with out_long.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        if not long_header:
            long_header = ["species_id"] + [
                "protein_id", "swissprot_accession", "entry_name", "protein_name", "organism", "gene_name",
                "evalue", "bitscore", "pident", "qcov", "scov",
            ]
        w.writerow(long_header)
        for r in long_rows:
            w.writerow(r)

    print(f"[99_reports] DONE species_summary={out_species} all_annotations={out_long}", flush=True)

if __name__ == "__main__":
    main()

