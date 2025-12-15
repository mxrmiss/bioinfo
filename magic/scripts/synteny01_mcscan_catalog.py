#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny01_mcscan_catalog.py
生成 JCVI karyotype 所需的 per-species BED（geneorder bed）与 gene 索引

硬合同：
- 读取 raw_data/synteny_species_meta.tsv（仅两列：species_id, group；顺序即全流程顺序）
- 读取 Step00 的 chr_rename_<species>.tsv，把 seqid_raw 映射为 Chr01..
- 从 raw_data/gff/<species>.* 提取基因坐标，输出：
    output/synteny_01_mcscan_catalog/<species>.geneorder.bed     (BED4，无表头)
    output/synteny_01_mcscan_catalog/geneorder_index_<species>.tsv
- gene_id 抽取/标准化规则与全流程一致（见 common utils）
- 不生成任何哨兵文件
- 脚本不接受命令行参数
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional


# -*- coding: utf-8 -*-
"""
公共工具函数（仅供本脚本使用；不依赖外部第三方库）。
"""

import os
import sys
import re
import csv
import gzip
import time
import shutil
import traceback
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterable


def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


class Logger:
    def __init__(self, step_name: str, log_file: Path):
        self.step_name = step_name
        self.log_file = log_file
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

    def _write(self, level: str, msg: str) -> None:
        line = f"{now_ts()} [{level}] {msg}"
        print(line, flush=True)
        with self.log_file.open("a", encoding="utf-8") as fw:
            fw.write(line + "\n")

    def info(self, msg: str) -> None:
        self._write("INFO", msg)

    def warn(self, msg: str) -> None:
        self._write("WARN", msg)

    def error(self, msg: str) -> None:
        self._write("ERROR", msg)


def abort(logger: Logger, msg: str, code: int = 1) -> None:
    logger.error(msg)
    raise SystemExit(code)


def clean_dir(path: Path, clean: bool, logger: Logger) -> None:
    if clean and path.exists():
        logger.info(f"Clean output dir: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def open_text_auto(path: Path):
    """
    支持普通文本与 .gz。
    返回一个 text mode file handle（utf-8，遇到奇怪字符直接替换）。
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def detect_first_existing(base_dir: Path, candidates: List[str]) -> Optional[Path]:
    for name in candidates:
        p = base_dir / name
        if p.exists():
            return p
    return None


def read_meta(meta_tsv: Path) -> Tuple[List[str], Dict[str, str]]:
    """
    meta 表硬合同：两列（tab 分隔）
      species_id, group
    """
    species: List[str] = []
    group: Dict[str, str] = {}
    with open_text_auto(meta_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != ["species_id", "group"]:
            raise ValueError(f"meta header must be: species_id<TAB>group, got: {header}")
        for ln, line in enumerate(fr, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"meta line {ln} must have 2 columns, got {len(parts)}: {line}")
            sid, grp = parts
            sid = sid.strip()
            grp = grp.strip()
            if not sid:
                raise ValueError(f"meta line {ln} species_id empty")
            species.append(sid)
            group[sid] = grp
    if len(species) < 2:
        raise ValueError("meta must contain at least 2 species")
    return species, group


def parse_species_display_label(species_id: str) -> str:
    """
    由 Genus_species_xxx 生成显示名：G. species
    """
    parts = species_id.split("_")
    if len(parts) < 2:
        return species_id
    genus = parts[0]
    sp = parts[1]
    if not genus:
        return species_id
    return f"{genus[0].upper()}. {sp.lower()}"


_HEX_RE = re.compile(r"^#[0-9a-fA-F]{6}\*")


def strip_color_prefix(gene_id: str) -> str:
    """
    如果 gene_id 是 #RRGGBB*xxx，返回 xxx，否则原样返回。
    """
    if _HEX_RE.match(gene_id):
        return gene_id.split("*", 1)[1]
    return gene_id


def normalize_gene_id(raw: str) -> str:
    """
    ID 标准化硬合同：
      - 去引号/空白
      - 若含 |，取最后一个 | 后
      - 去 rna- 前缀
    """
    x = raw.strip().strip('"').strip("'").strip()
    if "|" in x:
        x = x.split("|")[-1]
    if x.startswith("rna-"):
        x = x[4:]
    return x


def parse_gff_attributes(attr: str) -> Dict[str, str]:
    """
    GFF3 attributes 解析为 dict；允许出现无 '=' 的碎片但会忽略。
    """
    d: Dict[str, str] = {}
    for item in attr.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        k, v = item.split("=", 1)
        d[k.strip()] = v.strip()
    return d


def extract_gene_id_from_gff_attr(attr: str) -> str:
    """
    gene_id 抽取优先级：
      transcript_id > orig_transcript_id > ID > Parent
    """
    d = parse_gff_attributes(attr)
    for k in ("transcript_id", "orig_transcript_id", "ID", "Parent"):
        if k in d and d[k]:
            return normalize_gene_id(d[k])
    # 没有任何键，直接返回整段（仍做 normalize）
    return normalize_gene_id(attr)


def safe_int(x: str, default: int = 0) -> int:
    try:
        return int(x)
    except Exception:
        return default


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fw:
        w = csv.writer(fw, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


# ============================
# 用户参数区（手填；不走命令行）# ============================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
META_TSV = RAW_DATA_DIR / "synteny_species_meta.tsv"
GFF_DIR = RAW_DATA_DIR / "gff"

STEP00_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"

OUTPUT_DIR = PROJECT_ROOT / "output" / "synteny_01_mcscan_catalog"
LOG_DIR = PROJECT_ROOT / "logs"
CLEAN_OUTPUT = True

# 从 GFF 中提取哪些 feature 作为“基因坐标”
# 规则：优先使用 mRNA/transcript；若文件中没有 mRNA/transcript，则使用 gene
PREFERRED_FEATURES = ("mRNA", "transcript")


def load_chr_rename(chr_rename_tsv: Path) -> Dict[str, str]:
    """
    seqid_raw -> seqid_renamed（仅保留 is_chromosome=yes）
    """
    mp: Dict[str, str] = {}
    with open_text_auto(chr_rename_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        need = ["species_id", "seqid_raw", "seqid_renamed", "rank", "length_bp", "is_chromosome"]
        if header != need:
            raise ValueError(f"Bad header in {chr_rename_tsv.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                continue
            seqid_raw = parts[1]
            seqid_renamed = parts[2]
            is_chr = parts[5]
            if is_chr == "yes" and seqid_renamed != "NA":
                mp[seqid_raw] = seqid_renamed
    if not mp:
        raise ValueError(f"No chromosomes kept in {chr_rename_tsv.name}")
    return mp


def detect_gff(sid: str) -> Path:
    p = detect_first_existing(GFF_DIR, [
        f"{sid}.gff3", f"{sid}.gff",
        f"{sid}.gff3.gz", f"{sid}.gff.gz",
    ])
    if p is None:
        raise FileNotFoundError(f"Missing GFF for {sid} under raw_data/gff/")
    return p


def scan_gff_for_feature_types(gff_path: Path) -> Tuple[bool, bool]:
    """
    返回：(has_mrna_or_transcript, has_gene)
    """
    has_m = False
    has_g = False
    with open_text_auto(gff_path) as fr:
        for line in fr:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            ftype = parts[2].lower()
            if ftype in ("mrna", "transcript"):
                has_m = True
            elif ftype == "gene":
                has_g = True
            if has_m and has_g:
                break
    return has_m, has_g


def main() -> None:
    logger = Logger("synteny01", LOG_DIR / "synteny01_mcscan_catalog.log")
    t0 = time.time()

    logger.info("========== synteny01 — geneorder BED ==========")
    logger.info(f"PROJECT_ROOT={PROJECT_ROOT}")
    logger.info(f"OUTPUT_DIR={OUTPUT_DIR}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"PREFERRED_FEATURES={','.join(PREFERRED_FEATURES)}")

    if not META_TSV.exists():
        abort(logger, f"Missing meta: {META_TSV}")
    if not STEP00_DIR.exists():
        abort(logger, f"Missing step00 output dir: {STEP00_DIR}")

    species, _group = read_meta(META_TSV)
    logger.info(f"META species count={len(species)}")

    clean_dir(OUTPUT_DIR, CLEAN_OUTPUT, logger)

    for sid in species:
        chr_rename = STEP00_DIR / f"chr_rename_{sid}.tsv"
        if not chr_rename.exists():
            abort(logger, f"[{sid}] Missing step00 file: {chr_rename}")
        mp = load_chr_rename(chr_rename)

        gff = detect_gff(sid)
        has_m, has_g = scan_gff_for_feature_types(gff)
        use_mrna = has_m
        fset = set(x.lower() for x in (PREFERRED_FEATURES if use_mrna else ("gene",)))

        out_bed = OUTPUT_DIR / f"{sid}.geneorder.bed"
        out_idx = OUTPUT_DIR / f"geneorder_index_{sid}.tsv"

        bed_lines: List[List[str]] = []
        idx_rows: List[List[str]] = []

        kept_gene = 0
        dropped_non_chr = 0
        dropped_no_id = 0

        with open_text_auto(gff) as fr:
            for line in fr:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seqid_raw = parts[0]
                ftype = parts[2].lower()
                if ftype not in fset:
                    continue
                if seqid_raw not in mp:
                    dropped_non_chr += 1
                    continue

                start = safe_int(parts[3], 0)
                end = safe_int(parts[4], 0)
                strand = parts[6] if len(parts) > 6 else "."
                attr = parts[8]

                gene_id = extract_gene_id_from_gff_attr(attr)
                gene_id = normalize_gene_id(gene_id)
                if not gene_id:
                    dropped_no_id += 1
                    continue

                chr_name = mp[seqid_raw]
                # BED 坐标：0-based start，1-based end（与 jcvi.formats.gff bed 的习惯一致）
                bed_start = max(0, start - 1)
                bed_end = max(bed_start + 1, end)

                bed_lines.append([chr_name, str(bed_start), str(bed_end), gene_id])
                idx_rows.append([gene_id, chr_name, str(bed_start), str(bed_end), strand])
                kept_gene += 1

        if kept_gene == 0:
            abort(logger, f"[{sid}] geneorder is empty (no features extracted). Check GFF and chr_rename.")

        # 排序：Chr01.. + start
        def chr_key(c: str) -> Tuple[int, str]:
            m = re.match(r"^Chr(\d+)$", c)
            if m:
                return (int(m.group(1)), c)
            return (10**9, c)

        bed_lines.sort(key=lambda r: (chr_key(r[0]), int(r[1]), r[3]))
        idx_rows.sort(key=lambda r: (chr_key(r[1]), int(r[2]), r[0]))

        # 写文件
        out_bed.parent.mkdir(parents=True, exist_ok=True)
        with out_bed.open("w", encoding="utf-8", newline="") as fw:
            for r in bed_lines:
                fw.write("\t".join(r) + "\n")

        write_tsv(out_idx, ["gene_id", "chr", "start", "end", "strand"], idx_rows)

        logger.info(
            f"[{sid}] gff={gff.name} use={'mRNA/transcript' if use_mrna else 'gene'} "
            f"kept_gene={kept_gene} dropped_non_chr={dropped_non_chr} dropped_no_id={dropped_no_id}"
        )

    logger.info(f"Done. runtime_sec={int(time.time()-t0)}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        log_path = Path(__file__).resolve().parents[1] / "logs" / "synteny01_mcscan_catalog.log"
        lg = Logger("synteny01", log_path)
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise
