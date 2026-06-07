#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
99_reports.py
================================================================================
[定位（皇上已确认）]
- 05_go.py 负责生成“最终注释大表”：results/05_go/{REL}/{species}/annot.tsv
- 99_reports.py 只做复制整理 + 汇总（不生成、不加工 annot.tsv 内容）

[功能]
1) 复制每个物种的最终注释表（05 产物）到 results/99_reports/{REL}/
   - {species}.annot.tsv
2) 生成 QC 汇总：species_summary.tsv
   - 汇总 03_filter / 04_annot / 05_go 的关键统计
3) 可选：拼接所有物种注释表为长表 all_annotations.long.tsv
   - 由 config.yaml: reports.make_long_annotation 控制
4) 可选：把 99 产出的“注释 TSV”（每物种 + 可选 long）额外导出为 XLSX
   - 由 config.yaml: reports.export_annotation_xlsx 控制
   - 注意 Excel 单 sheet 行数上限 ~1,048,576，超限会跳过并 WARN

[配置开关（config.yaml）]
reports:
  enabled: true
  make_long_annotation: true
  export_annotation_xlsx: false

[输出目录]
results/99_reports/{REL}/
  - species_summary.tsv
  - {species}.annot.tsv
  - (可选) all_annotations.long.tsv
  - (可选) {species}.annot.xlsx
  - (可选) all_annotations.long.xlsx
================================================================================
"""

from __future__ import annotations

import os
import csv
import shutil
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple


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
    eprint,
)

EXCEL_MAX_ROWS = 1048576


def _cfg_get_reports(cfg: Dict[str, Any]) -> Dict[str, Any]:
    d = cfg.get("reports", {})
    return d if isinstance(d, dict) else {}


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


def read_go_qc_tsv(p: Path) -> Dict[str, str]:
    """
    读取 05_go 的 qc.tsv（有表头，单行数据）：
    species_id  REL  N_total  N_has_go  N_no_go  HAS_GO_rate
    """
    d: Dict[str, str] = {}
    if not p.is_file() or p.stat().st_size == 0:
        return d
    with p.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        line = f.readline()
        if not header or not line:
            return d
        h = header.rstrip("\n").split("\t")
        v = line.rstrip("\n").split("\t")
        for i, k in enumerate(h):
            if i < len(v):
                d[k] = v[i]
    return d


def is_nonempty_file(p: Path) -> bool:
    return p.is_file() and p.stat().st_size > 0


def copy_file_if_ok(src: Path, dst: Path) -> bool:
    """
    原样复制，成功返回 True；若 src 缺失/为空，WARN 并返回 False。
    """
    ensure_dir(dst.parent)
    if not is_nonempty_file(src):
        print(f"[99_reports] WARN missing/empty file (skip copy): {src}", flush=True)
        return False
    shutil.copy2(src.as_posix(), dst.as_posix())
    return True


def _tsv_to_xlsx_stream(tsv_path: Path, xlsx_path: Path, sheet_name: str = "Sheet1") -> Tuple[bool, int]:
    """
    把 TSV 流式写成 XLSX（单 sheet）。
    返回 (ok, n_rows_written)。
    若 TSV 行数 > Excel 上限，返回 (False, n_rows) 并不写 xlsx。
    """
    try:
        from openpyxl import Workbook
    except Exception as ex:
        raise RuntimeError("export_annotation_xlsx=true but openpyxl is not installed. Please install openpyxl.") from ex

    if not is_nonempty_file(tsv_path):
        print(f"[99_reports] WARN missing/empty TSV (skip xlsx): {tsv_path}", flush=True)
        return (False, 0)

    # 先快速估算行数：逐行计数（不加载内存）
    n_lines = 0
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        for _ in f:
            n_lines += 1
            if n_lines > EXCEL_MAX_ROWS:
                print(
                    f"[99_reports] WARN TSV rows exceed Excel limit ({EXCEL_MAX_ROWS}): {tsv_path} "
                    f"(rows>={n_lines}), skip xlsx export",
                    flush=True
                )
                return (False, n_lines)

    ensure_dir(xlsx_path.parent)
    wb = Workbook(write_only=True)
    ws = wb.create_sheet(title=sheet_name)

    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            row = line.rstrip("\n").split("\t")
            ws.append(row)

    # openpyxl 的 write_only 模式会自动处理写出
    wb.save(xlsx_path.as_posix())
    return (True, n_lines)


def make_long_annotation(species: List[str], per_species_tsv: Dict[str, Path], out_long: Path) -> None:
    """
    拼接所有物种的 {sid}.annot.tsv 成长表：
      species_id + 原表所有列
    不做任何加工，不补列；要求所有参与拼接的表头一致，否则跳过不一致的物种。
    """
    ensure_dir(out_long.parent)

    std_header: Optional[List[str]] = None
    n_written = 0
    n_species_used = 0

    with out_long.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout, delimiter="\t", lineterminator="\n")

        for sid in species:
            tsv = per_species_tsv.get(sid)
            if tsv is None or (not is_nonempty_file(tsv)):
                print(f"[99_reports] WARN no per-species annot to join (skip): {sid}", flush=True)
                continue

            with tsv.open("r", encoding="utf-8", errors="replace") as fin:
                hdr_line = fin.readline()
                if not hdr_line:
                    print(f"[99_reports] WARN empty header (skip): {tsv}", flush=True)
                    continue
                hdr = hdr_line.rstrip("\n").split("\t")

                if std_header is None:
                    std_header = hdr
                    w.writerow(["species_id"] + std_header)
                else:
                    if hdr != std_header:
                        print(f"[99_reports] WARN header mismatch (skip species): {sid}", flush=True)
                        print(f"[99_reports] WARN expected header from first species, got different header in: {tsv}", flush=True)
                        continue

                # 写数据行
                for line in fin:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    w.writerow([sid] + parts)
                    n_written += 1

            n_species_used += 1
            print(f"[99_reports] join progress: used_species={n_species_used} rows={n_written:,}", flush=True)

        # 若一个物种都没拼进去，也要写一个最小表头，避免下游找不到文件
        if std_header is None:
            w.writerow(["species_id"])
            print(f"[99_reports] WARN no species joined; wrote minimal header only: {out_long}", flush=True)

    print(f"[99_reports] DONE join: out={out_long} used_species={n_species_used} rows={n_written:,}", flush=True)


def main() -> None:
    root = find_project_root()
    os.chdir(root)

    cfg = read_yaml_config(root)
    rel = parse_rel_from_01_db(root)
    species = get_species_list(cfg)

    rcfg = _cfg_get_reports(cfg)
    enabled = bool(rcfg.get("enabled", True))
    make_long = bool(rcfg.get("make_long_annotation", True))
    export_xlsx = bool(rcfg.get("export_annotation_xlsx", False))

    if not enabled:
        print("[99_reports] reports.enabled=false -> skip", flush=True)
        return

    outdir = root / "results" / "99_reports" / rel
    ensure_dir(outdir)

    out_species = outdir / "species_summary.tsv"
    out_long = outdir / "all_annotations.long.tsv"

    print(f"[99_reports] REL={rel}", flush=True)
    print(f"[99_reports] N_species={len(species)}", flush=True)
    print(f"[99_reports] OUTDIR={outdir}", flush=True)
    print(f"[99_reports] make_long_annotation={make_long}", flush=True)
    print(f"[99_reports] export_annotation_xlsx={export_xlsx}", flush=True)

    # -----------------------------
    # 1) 复制每物种最终注释表（05 产物）
    # -----------------------------
    per_species_annot: Dict[str, Path] = {}
    for sid in species:
        src = root / "results" / "05_go" / rel / sid / "annot.tsv"
        dst = outdir / f"{sid}.annot.tsv"
        ok = copy_file_if_ok(src, dst)
        if ok:
            per_species_annot[sid] = dst
            print(f"[99_reports] copied annot -> {dst}", flush=True)

    # -----------------------------
    # 2) 生成 species_summary.tsv（QC 汇总）
    # -----------------------------
    header = [
        "species_id", "REL",
        "N_QUERY", "N_FILTERED_HITS", "N_BESTHIT", "HIT_RATE",
        "N_annot", "NO_OS", "NO_GN", "NO_OS_rate", "NO_GN_rate",
        "N_HAS_GO", "N_NO_GO", "HAS_GO_rate",
    ]

    rows: List[List[str]] = []

    for sid in species:
        # 03_filter summary
        summ = read_kv_tsv(root / "results" / "03_filter" / rel / sid / "summary.tsv")
        # 04_annot qc
        qc4 = read_kv_tsv(root / "results" / "04_annot" / rel / sid / "qc_missing_fields.tsv")
        # 05_go qc
        qc5 = read_go_qc_tsv(root / "results" / "05_go" / rel / sid / "qc.tsv")

        n_query = summ.get("N_QUERY", "0")
        n_filt = summ.get("N_FILTERED_HITS", "0")
        n_best = summ.get("N_BESTHIT", "0")
        hit_rate = summ.get("HIT_RATE", "0")

        n_annot = qc4.get("N_annot", "0")
        no_os = qc4.get("NO_OS", "0")
        no_gn = qc4.get("NO_GN", "0")
        no_os_rate = qc4.get("NO_OS_rate", "0")
        no_gn_rate = qc4.get("NO_GN_rate", "0")

        # 05_go qc.tsv
        # 若缺失，默认 NA 更诚实（区分“确实为0”与“没跑到05”）
        n_has_go = qc5.get("N_has_go", "NA")
        n_no_go = qc5.get("N_no_go", "NA")
        has_go_rate = qc5.get("HAS_GO_rate", "NA")

        rows.append([
            sid, rel,
            n_query, n_filt, n_best, hit_rate,
            n_annot, no_os, no_gn, no_os_rate, no_gn_rate,
            n_has_go, n_no_go, has_go_rate,
        ])

    with out_species.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(r)

    print(f"[99_reports] wrote species_summary -> {out_species}", flush=True)

    # -----------------------------
    # 3) 可选：拼接 long 注释表
    # -----------------------------
    if make_long:
        print(f"[99_reports] step: make long annotation -> {out_long}", flush=True)
        make_long_annotation(species, per_species_annot, out_long)
    else:
        print("[99_reports] make_long_annotation=false -> skip long annotation", flush=True)

    # -----------------------------
    # 4) 可选：导出注释 TSV 为 XLSX
    # -----------------------------
    if export_xlsx:
        print("[99_reports] step: export annotation TSV -> XLSX", flush=True)

        # per-species
        for sid, tsv in per_species_annot.items():
            xlsx = tsv.with_suffix(".xlsx")
            ok, n_rows = _tsv_to_xlsx_stream(tsv, xlsx, sheet_name="annot")
            if ok:
                print(f"[99_reports] xlsx OK: {xlsx} (rows={n_rows:,})", flush=True)

        # long (only if exists and non-empty and make_long enabled)
        if make_long and is_nonempty_file(out_long):
            long_xlsx = out_long.with_suffix(".xlsx")
            ok, n_rows = _tsv_to_xlsx_stream(out_long, long_xlsx, sheet_name="long")
            if ok:
                print(f"[99_reports] xlsx OK: {long_xlsx} (rows={n_rows:,})", flush=True)

    else:
        print("[99_reports] export_annotation_xlsx=false -> skip xlsx export", flush=True)

    print("[99_reports] ALL DONE", flush=True)


if __name__ == "__main__":
    main()

