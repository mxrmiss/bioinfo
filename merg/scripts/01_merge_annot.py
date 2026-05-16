#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_merge_annot.py
功能：合并 InterProScan wide 表 与 Swiss-Prot best-hit 注释表，生成单物种“总注释表”

皇上指定规则：
  - protein_id 全量保留：两表没匹配到的行也要进入总表（FULL OUTER JOIN）
  - 合并后不包含 GO 列（InterPro 的 GO；Swiss-Prot 的 go_ids / n_go）
  - 保留 gene_name；不保留 organism 和 entry_name
  - 同时输出 TSV + XLSX（不加开关）
  - 流式屏幕输出 + 日志文件
  - 输出命名：从 Swiss-Prot 输入文件名提取 Genus_species，输出为 Genus_species_merged.*
  - 路径全部使用相对路径（相对于 anno/merg 运行目录）
"""

from __future__ import annotations

import os
import sys
import time
import datetime
from typing import List

import pandas as pd


# =============================================================================
# 【皇上手填参数区】不支持命令行参数，全部在这里改（相对 anno/merg）
# =============================================================================

# 1) InterPro wide 表（脚本会自动丢弃 GO 列）
INTERPRO_WIDE_TSV = "../interpro/results/02_wide/5.76-107.0.core/Sinonovacula_rivularis.faa.wide.tsv"

# 2) Swiss-Prot 单物种注释表（脚本会自动丢弃 organism/entry_name/go_ids/n_go）
SWISSPROT_ANNOT_TSV = "../swissprot/results/99_reports/2025_04/Sinonovacula_rivularis.annot.tsv"

# 3) 产物根目录（相对 anno/merg）
#    会自动创建一个 01_ 前缀目录：results/01_merge_annot/
RESULTS_ROOT = "results"
RESULTS_SUBDIR_01 = "01_merge_annot"

# 4) 读取参数（一般不用改）
TSV_SEP = "\t"
ENCODING = "utf-8"
LOW_MEMORY = False

# 5) 合并后列顺序（严格按皇上保留/删除规则）
FINAL_COL_ORDER = [
    "protein_id",
    "swissprot_accession",
    "protein_name",
    "gene_name",
    "evalue",
    "bitscore",
    "pident",
    "qcov",
    "scov",
    "Pfam",
    "SMART",
    "ProSitePatterns",
    "ProSiteProfiles",
    "SUPERFAMILY",
    "PANTHER",
    "Gene3D",
]


# =============================================================================
# 流式输出 + 日志
# =============================================================================

def now_ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class Logger:
    def __init__(self, log_path: str):
        self.log_path = log_path
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        self.fh = open(log_path, "w", encoding="utf-8")

    def close(self):
        try:
            self.fh.close()
        except Exception:
            pass

    def log(self, msg: str):
        line = f"[{now_ts()}] {msg}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()


def die(logger: Logger, msg: str, code: int = 1):
    logger.log("ERROR: " + msg)
    logger.close()
    sys.exit(code)


def infer_prefix_from_swissprot_path(path: str) -> str:
    """
    从 Swiss-Prot 输入文件名推断输出前缀：
      - basename 如：Sinonovacula_constricta.annot.tsv
      - 取 '_' 分割的前两个 token：Sinonovacula + constricta -> Sinonovacula_constricta
    兜底：
      - 如果没有 '_' 或 token 不足 2 个，则取第一个 '.' 之前的部分
    """
    base = os.path.basename(path).strip()
    parts = [p for p in base.split("_") if p]
    if len(parts) >= 2:
        genus = parts[0].split(".", 1)[0]
        species = parts[1].split(".", 1)[0]
        if genus and species:
            return f"{genus}_{species}"
    stem = base.split(".", 1)[0]
    return stem if stem else "output"


# =============================================================================
# 读表 + 清洗
# =============================================================================

def read_tsv_required(logger: Logger, path: str, name: str) -> pd.DataFrame:
    if not os.path.exists(path):
        die(logger, f"{name} file not found: {path}")
    if os.path.getsize(path) == 0:
        die(logger, f"{name} file is empty: {path}")

    logger.log(f"Read {name}: {path}")
    t0 = time.time()
    df = pd.read_csv(
        path,
        sep=TSV_SEP,
        dtype=str,
        encoding=ENCODING,
        low_memory=LOW_MEMORY,
        keep_default_na=False,
    )
    logger.log(f"{name} rows={len(df):,}, cols={len(df.columns)} ({time.time() - t0:.2f}s)")
    return df


def normalize_protein_id(df: pd.DataFrame, col: str = "protein_id") -> pd.DataFrame:
    if col in df.columns:
        df[col] = df[col].astype(str).str.strip()
    return df


def drop_cols_if_exist(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    drop = [c for c in cols if c in df.columns]
    if drop:
        df = df.drop(columns=drop)
    return df


def ensure_has_column(df: pd.DataFrame, col: str, default: str = "") -> pd.DataFrame:
    if col not in df.columns:
        df[col] = default
    return df


# =============================================================================
# 去重策略
# =============================================================================

def dedup_swissprot(logger: Logger, df: pd.DataFrame) -> pd.DataFrame:
    """
    Swiss-Prot 表理论上一蛋白一行；若出现重复：
      - 按 bitscore 取最大的一行
    """
    if "protein_id" not in df.columns:
        die(logger, "Swiss-Prot table missing required column: protein_id")

    if df["protein_id"].duplicated().any():
        logger.log("Swiss-Prot duplicated protein_id detected -> deduplicate by max(bitscore)")
        bs = pd.to_numeric(df.get("bitscore", ""), errors="coerce")
        df = df.assign(_bitscore_num=bs)
        df = df.sort_values(by=["protein_id", "_bitscore_num"], ascending=[True, False])
        df = df.drop_duplicates(subset=["protein_id"], keep="first")
        df = df.drop(columns=["_bitscore_num"])

    return df


def merge_interpro_duplicates(logger: Logger, df: pd.DataFrame) -> pd.DataFrame:
    """
    InterPro wide 表如果出现重复 protein_id：
      - 对每个库列，把非空值按 ';' 拆条目去重后再拼回
    """
    if "protein_id" not in df.columns:
        die(logger, "InterPro wide table missing required column: protein_id")

    if not df["protein_id"].duplicated().any():
        return df

    logger.log("InterPro wide duplicated protein_id detected -> merge rows by unique entries per column")

    cols = [c for c in df.columns if c != "protein_id"]

    def merge_cell(values: List[str]) -> str:
        items: List[str] = []
        for v in values:
            s = str(v).strip() if v is not None else ""
            if not s:
                continue
            parts = [p.strip() for p in s.split(";") if p.strip()]
            items.extend(parts)

        seen = set()
        uniq: List[str] = []
        for it in items:
            if it not in seen:
                seen.add(it)
                uniq.append(it)
        return "; ".join(uniq)

    grouped = []
    for pid, sub in df.groupby("protein_id", sort=False):
        row = {"protein_id": pid}
        for c in cols:
            row[c] = merge_cell(sub[c].tolist())
        grouped.append(row)

    return pd.DataFrame(grouped)


# =============================================================================
# 主流程
# =============================================================================

def main():
    out_prefix = infer_prefix_from_swissprot_path(SWISSPROT_ANNOT_TSV)

    outdir = os.path.join(RESULTS_ROOT, RESULTS_SUBDIR_01)
    os.makedirs(outdir, exist_ok=True)

    out_tsv = os.path.join(outdir, f"{out_prefix}_merged.tsv")
    out_xlsx = os.path.join(outdir, f"{out_prefix}_merged.xlsx")
    log_path = os.path.join(outdir, f"{out_prefix}_merged.log")

    logger = Logger(log_path)

    logger.log("Start merge annotation")
    logger.log(f"RunPWD={os.getcwd()}")
    logger.log(f"InterPro={INTERPRO_WIDE_TSV}")
    logger.log(f"SwissProt={SWISSPROT_ANNOT_TSV}")
    logger.log(f"OUTDIR={outdir}")
    logger.log(f"OUT_TSV={out_tsv}")
    logger.log(f"OUT_XLSX={out_xlsx}")
    logger.log(f"LOG={log_path}")

    df_ipr = read_tsv_required(logger, INTERPRO_WIDE_TSV, "InterPro wide")
    df_sp = read_tsv_required(logger, SWISSPROT_ANNOT_TSV, "Swiss-Prot annot")

    df_ipr = normalize_protein_id(df_ipr, "protein_id")
    df_sp = normalize_protein_id(df_sp, "protein_id")

    df_ipr = drop_cols_if_exist(df_ipr, ["GO"])
    df_sp = drop_cols_if_exist(df_sp, ["go_ids", "n_go", "organism", "entry_name"])

    df_ipr = ensure_has_column(df_ipr, "protein_id")
    df_sp = ensure_has_column(df_sp, "protein_id")

    df_sp = dedup_swissprot(logger, df_sp)
    df_ipr = merge_interpro_duplicates(logger, df_ipr)

    set_ipr = set(df_ipr["protein_id"].tolist())
    set_sp = set(df_sp["protein_id"].tolist())
    n_ipr = len(set_ipr)
    n_sp = len(set_sp)
    n_overlap = len(set_ipr & set_sp)
    n_total_expected = n_ipr + n_sp - n_overlap

    logger.log(f"N_interpro_unique={n_ipr:,}")
    logger.log(f"N_swissprot_unique={n_sp:,}")
    logger.log(f"N_overlap={n_overlap:,}")
    logger.log(f"N_total_expected(full_outer)={n_total_expected:,}")

    logger.log("Do FULL OUTER JOIN on protein_id")
    t0 = time.time()
    merged = pd.merge(df_sp, df_ipr, on="protein_id", how="outer")
    logger.log(f"Merged rows={len(merged):,} ({time.time() - t0:.2f}s)")

    if len(merged) != n_total_expected:
        logger.log("WARNING: merged row count != expected (possible duplicate ids or hidden whitespace)")
        logger.log(f"merged={len(merged):,}, expected={n_total_expected:,}")

    for c in FINAL_COL_ORDER:
        if c not in merged.columns:
            merged[c] = ""

    merged = merged[FINAL_COL_ORDER].copy()
    merged = merged.sort_values(by=["protein_id"], ascending=True)

    logger.log(f"Write TSV: {out_tsv}")
    merged.to_csv(out_tsv, sep="\t", index=False, encoding="utf-8")

    logger.log(f"Write XLSX: {out_xlsx}")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        merged.to_excel(writer, sheet_name="merged", index=False)

    logger.log("Done")
    logger.close()


if __name__ == "__main__":
    main()

