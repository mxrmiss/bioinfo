#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_genespace_prepare_inputs.py
—— 为 GENESPACE 量身定制：从 phylo/data/annotation + phylo/data/proteomes 生成可直接使用的 bed + peptide

核心目标（GENESPACE 输入契约）：
  results/genespace/peptide/<GenomeID>.fa
  results/genespace/bed/<GenomeID>.bed    (chr, start0, end, name)

关键规则：
  1) 物种集合以 data/proteomes 中的蛋白文件为准（从文件名提取 GenomeID）；
     - 蛋白多于注释（某蛋白物种找不到对应 GFF）=> 报错终止
     - 注释多于蛋白（annotation 目录里有多余物种）=> 忽略，不报错
  2) BED 坐标来源：沿用原 synteny01 逻辑，只使用 mRNA/transcript 层，写 BED 0-based start (start-1)
  3) ID 匹配：
     - 从 GFF 第 9 列提取“core id”（按你原脚本的 TRANSCRIPT_ATTR_KEYS 优先级；否则退回 ID/Parent）
     - 从蛋白 FASTA header 取完整主键 full_id（> 后第一个 token，通常无空格），其 core 为最后一个 '|' 后的片段
     - 用 core 对齐：GFF_core -> full_protein_id
     - 写 BED 第 4 列时写 full_protein_id（例如 Species|XM_...），保证 GENESPACE 可读
  4) “匹配到但不唯一”：同一 core 在蛋白中对应多个不同 full_id => 报错终止并列出冲突

运行方式：
  在 phylo/ 目录下运行：
    python scripts/12_genespace_prepare_inputs.py

不接受命令行参数：所有参数集中在脚本顶部“参数区”。
"""

from __future__ import annotations

import sys
import gzip
import shutil
import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional


# =========================
# 参数区（皇上可在此修改）
# =========================

# 项目根目录（默认：脚本所在目录的上一级 phylo/）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 输入目录（按皇上已确定的结构）
ANNOT_DIR = PROJECT_ROOT / "data" / "annotation"
PROTEOME_DIR = PROJECT_ROOT / "data" / "proteomes"

# 输出目录（GENESPACE working directory 的核心子目录）
OUT_ROOT = PROJECT_ROOT / "results" / "genespace"
OUT_BED_DIR = OUT_ROOT / "bed"
OUT_PEPTIDE_DIR = OUT_ROOT / "peptide"
OUT_LOG_DIR = OUT_ROOT / "logs"

# 蛋白文件后缀兼容（以此目录中出现的文件作为物种集合）
PROTEOME_SUFFIXES = [
    ".faa",
    ".fa",
    ".fasta",
    ".fas",
    ".faa.gz",
    ".fa.gz",
    ".fasta.gz",
    ".fas.gz",
]

# 注释文件后缀兼容（每物种仅一个文件，后缀可不同）
GFF_SUFFIXES = [
    ".gff3",
    ".gff",
    ".gff3.gz",
    ".gff.gz",
]

# 仅解析这些 feature 类型（沿用原 synteny01 的“只取 mRNA/transcript”策略）
GFF_FEATURE_TYPES = {"mRNA", "transcript"}

# 可选：仅保留满足正则的 seqid（不想过滤就留空字符串）
# 例：只保留 Chr01..Chr99 => r"^Chr\d+$"
KEEP_SEQID_REGEX = ""  # 默认不过滤

# GFF 第 9 列中用于提取“转录本 ID”的字段优先级（沿用原脚本）
TRANSCRIPT_ATTR_KEYS: List[str] = [
    "transcript_id",
    "orig_transcript_id",
]

# 日志等级
LOG_LEVEL = "INFO"


# =========================
# 日志与目录
# =========================

def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "12_genespace_prepare_inputs.log"

    logger = logging.getLogger("genespace12")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清空旧 handler，避免重复打印
    for h in list(logger.handlers):
        logger.removeHandler(h)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    ch = logging.StreamHandler(sys.stdout)

    fh.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))
    ch.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)

    logger.info("========== 12_genespace_prepare_inputs — 生成 GENESPACE bed + peptide ==========")
    logger.info("PROJECT_ROOT   = %s", PROJECT_ROOT)
    logger.info("ANNOT_DIR      = %s", ANNOT_DIR)
    logger.info("PROTEOME_DIR   = %s", PROTEOME_DIR)
    logger.info("OUT_ROOT       = %s", OUT_ROOT)
    logger.info("KEEP_SEQID_REGEX = %s", KEEP_SEQID_REGEX if KEEP_SEQID_REGEX else "(none)")
    return logger


def reset_output_dirs() -> None:
    """每次运行都删除旧 OUT_ROOT 并重建子目录。"""
    if OUT_ROOT.exists():
        shutil.rmtree(OUT_ROOT)
    OUT_BED_DIR.mkdir(parents=True, exist_ok=True)
    OUT_PEPTIDE_DIR.mkdir(parents=True, exist_ok=True)
    OUT_LOG_DIR.mkdir(parents=True, exist_ok=True)


# =========================
# 通用工具：文件名/后缀/打开
# =========================

def is_gz(path: Path) -> bool:
    return path.name.endswith(".gz")


def open_text(path: Path):
    """以文本方式打开（自动兼容 gz / 非 gz）。"""
    if is_gz(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def copy_to_fa(src: Path, dst: Path, logger: logging.Logger) -> None:
    """复制蛋白文件到 dst（.fa），支持 src 为 gz 或非 gz。"""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if is_gz(src):
        # 解压后写入
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(src, dst)
    logger.info("复制 peptide：%s -> %s", src.name, dst)


def strip_known_suffix(filename: str, suffixes: List[str]) -> Optional[str]:
    """
    从文件名中剥离已知后缀，返回 basename（物种名）。
    若不匹配任何后缀则返回 None。
    """
    for suf in sorted(suffixes, key=len, reverse=True):
        if filename.endswith(suf):
            return filename[: -len(suf)]
    return None


# =========================
# ID 规则（沿用原 synteny01 思路）
# =========================

def normalize_id_core(raw_id: str) -> str:
    """
    标准化 core：
      1) 去空白/引号
      2) 若包含 '|'，取最后一个 '|' 后的片段
    """
    val = raw_id.strip().strip('"').strip()
    if "|" in val:
        val = val.rsplit("|", 1)[1]
    return val


def extract_id_from_attr(attr: str) -> Optional[str]:
    """
    从 GFF 第 9 列提取 core id（沿用原脚本逻辑）：
      1) 优先 TRANSCRIPT_ATTR_KEYS（transcript_id / orig_transcript_id）
      2) 否则退回到 ID= 或 Parent=；若以 rna- 开头则去掉
      3) 统一做 normalize_id_core（取最后一个 | 后片段）
    """
    kv_list: List[Tuple[str, str]] = []
    for field in attr.split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        else:
            parts = field.split()
            if len(parts) >= 2:
                k, v = parts[0], parts[-1]
            elif len(parts) == 1:
                k, v = parts[0], ""
            else:
                continue
        k = k.strip()
        v = v.strip().strip('"').strip()
        if k:
            kv_list.append((k, v))

    # 1) 优先 keys
    for key in TRANSCRIPT_ATTR_KEYS:
        for k, v in kv_list:
            if k == key and v:
                return normalize_id_core(v)

    # 2) fallback ID/Parent
    id_val: Optional[str] = None
    for k, v in kv_list:
        if k == "ID" and v:
            id_val = v
            break
    if id_val is None:
        for k, v in kv_list:
            if k == "Parent" and v:
                id_val = v
                break
    if not id_val:
        return None

    if id_val.startswith("rna-"):
        id_val = id_val[len("rna-") :]

    return normalize_id_core(id_val)


def parse_fasta_core_map(fa_path: Path, logger: logging.Logger) -> Dict[str, str]:
    """
    读取蛋白 FASTA，建立 core -> full_id 的映射。
    full_id 取 '>' 后第一个 token（遇到空格即截断）。
    core 为 full_id 最后一个 '|' 后片段。
    若同一 core 映射到多个不同 full_id，直接报错终止（按皇上要求）。
    """
    core2full: Dict[str, str] = {}
    conflicts: Dict[str, List[str]] = {}

    n_headers = 0
    with open_text(fa_path) as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                n_headers += 1
                full = line[1:].strip().split()[0]  # 主键：第一个 token
                core = normalize_id_core(full)
                if core in core2full and core2full[core] != full:
                    conflicts.setdefault(core, sorted({core2full[core], full}))
                else:
                    core2full[core] = full

    if conflicts:
        logger.error("检测到同一 core 对应多个 peptide ID（不唯一匹配），将终止。冲突如下：")
        for core, ids in conflicts.items():
            logger.error("  core=%s  peptide_ids=%s", core, ", ".join(ids))
        sys.exit(20)

    logger.info("解析 peptide：%s  headers=%d  core2full=%d",
                fa_path.name, n_headers, len(core2full))
    return core2full


# =========================
# 查找输入文件
# =========================

def list_proteomes(logger: logging.Logger) -> Dict[str, Path]:
    """
    扫描 PROTEOME_DIR，返回 genomeID -> proteome_path。
    genomeID 从文件名剥离 PROTEOME_SUFFIXES 得到。
    """
    if not PROTEOME_DIR.exists():
        logger.error("蛋白目录不存在：%s", PROTEOME_DIR)
        sys.exit(1)

    mapping: Dict[str, Path] = {}
    for p in sorted(PROTEOME_DIR.iterdir()):
        if not p.is_file():
            continue
        genome = strip_known_suffix(p.name, PROTEOME_SUFFIXES)
        if genome is None:
            continue
        # 若同一 genome 出现多个后缀文件，优先选“更长后缀”（例如 .faa.gz 比 .faa）
        if genome in mapping:
            old = mapping[genome]
            if len(p.name) > len(old.name):
                mapping[genome] = p
        else:
            mapping[genome] = p

    if not mapping:
        logger.error("未在 %s 中发现任何蛋白文件（支持后缀：%s）",
                     PROTEOME_DIR, ", ".join(PROTEOME_SUFFIXES))
        sys.exit(2)

    logger.info("从 proteomes 目录识别到 %d 个物种：%s",
                len(mapping), ", ".join(list(mapping.keys())[:10]) + (" ..." if len(mapping) > 10 else ""))
    return mapping


def find_gff_for_species(genome_id: str, logger: logging.Logger) -> Path:
    """
    在 ANNOT_DIR 中查找指定物种的注释文件（仅一个，但后缀可能不同）。
    """
    if not ANNOT_DIR.exists():
        logger.error("注释目录不存在：%s", ANNOT_DIR)
        sys.exit(1)

    for suf in GFF_SUFFIXES:
        p = ANNOT_DIR / f"{genome_id}{suf}"
        if p.exists():
            logger.info("物种 %s 使用注释文件：%s", genome_id, p.name)
            return p

    logger.error("物种 %s 在 %s 中未找到注释文件（支持后缀：%s）",
                 genome_id, ANNOT_DIR, ", ".join(GFF_SUFFIXES))
    raise FileNotFoundError(genome_id)


# =========================
# GFF -> BED（写入 full peptide ID）
# =========================

def parse_gff_to_bed_intervals(
    gff_path: Path,
    core2full: Dict[str, str],
    logger: logging.Logger,
) -> Tuple[List[Tuple[str, int, int, str]], Dict[str, int]]:
    """
    解析 GFF，输出满足：
      - feature_type in {"mRNA","transcript"}
      - seqid 可选正则过滤
      - GFF_core 能在 core2full 中找到对应 full peptide id
    返回：
      intervals: [(chr, start0, end, full_id), ...]
      stats: 统计信息（便于日志输出）
    """
    intervals: List[Tuple[str, int, int, str]] = []
    stats = {
        "seen_feature": 0,
        "kept": 0,
        "skip_no_core": 0,
        "skip_no_peptide": 0,
        "skip_bad_coord": 0,
        "skip_seqid_filter": 0,
        "dup_full_id": 0,
    }

    seqid_re = re.compile(KEEP_SEQID_REGEX) if KEEP_SEQID_REGEX else None
    seen_full: set[str] = set()

    with open_text(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqid = parts[0]
            feature_type = parts[2]
            start_str = parts[3]
            end_str = parts[4]
            attr = parts[8]

            if feature_type not in GFF_FEATURE_TYPES:
                continue

            stats["seen_feature"] += 1

            if seqid_re and (seqid_re.search(seqid) is None):
                stats["skip_seqid_filter"] += 1
                continue

            core = extract_id_from_attr(attr)
            if not core:
                stats["skip_no_core"] += 1
                continue

            full_id = core2full.get(core)
            if not full_id:
                # 注释多于蛋白很正常：跳过即可
                stats["skip_no_peptide"] += 1
                continue

            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                stats["skip_bad_coord"] += 1
                continue

            # 防御性：若 start/end 颠倒则纠正
            if start <= 0 or end <= 0:
                stats["skip_bad_coord"] += 1
                continue
            if start > end:
                start, end = end, start

            start0 = start - 1  # BED 0-based
            if start0 < 0:
                start0 = 0

            # 防止重复 full_id（少数 GFF 会重复写 mRNA 行）
            if full_id in seen_full:
                stats["dup_full_id"] += 1
                continue
            seen_full.add(full_id)

            intervals.append((seqid, start0, end, full_id))
            stats["kept"] += 1

    # 可选排序：按 chr、start0 排序，便于阅读/调试（GENESPACE 不强制）
    intervals.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    return intervals, stats


def write_bed(genome_id: str, intervals: List[Tuple[str, int, int, str]], logger: logging.Logger) -> Path:
    out_path = OUT_BED_DIR / f"{genome_id}.bed"
    with out_path.open("w", encoding="utf-8") as f:
        for chr_name, start0, end, full_id in intervals:
            f.write(f"{chr_name}\t{start0}\t{end}\t{full_id}\n")
    logger.info("写出 BED：%s（%d lines）", out_path, len(intervals))
    return out_path


# =========================
# 主流程
# =========================

def main() -> None:
    reset_output_dirs()
    logger = setup_logging(OUT_LOG_DIR)

    # 1) 以 proteomes 目录定义物种集合
    proteomes = list_proteomes(logger)
    genome_ids = sorted(proteomes.keys())
    logger.info("将处理 %d 个物种（以 proteomes 为准）。", len(genome_ids))

    # 2) 检查每个物种必须存在注释文件（蛋白多于注释 => 报错终止）
    missing_annot: List[str] = []
    gff_paths: Dict[str, Path] = {}
    for gid in genome_ids:
        try:
            gff_paths[gid] = find_gff_for_species(gid, logger)
        except FileNotFoundError:
            missing_annot.append(gid)

    if missing_annot:
        logger.error("检测到蛋白目录中的以下物种找不到对应注释文件，将终止：")
        for gid in missing_annot:
            logger.error("  - %s", gid)
        logger.error("请在 %s 中为上述物种提供 {gid}.gff/.gff3(.gz) 之一。", ANNOT_DIR)
        sys.exit(10)

    # 3) 复制 peptide 到 results/genespace/peptide，并为每个物种建立 core->full 映射
    core_maps: Dict[str, Dict[str, str]] = {}
    for gid in genome_ids:
        src = proteomes[gid]
        dst = OUT_PEPTIDE_DIR / f"{gid}.fa"
        copy_to_fa(src, dst, logger)
        core_maps[gid] = parse_fasta_core_map(dst, logger)

    # 4) 逐物种解析 GFF -> BED（第 4 列写 full peptide id）
    total_bed_lines = 0
    for gid in genome_ids:
        logger.info("====== 处理物种：%s ======", gid)
        gff = gff_paths[gid]
        core2full = core_maps[gid]

        intervals, stats = parse_gff_to_bed_intervals(gff, core2full, logger)
        write_bed(gid, intervals, logger)
        total_bed_lines += len(intervals)

        logger.info(
            "统计[%s] seen_feature=%d kept=%d skip_no_core=%d skip_no_peptide=%d skip_bad_coord=%d skip_seqid_filter=%d dup_full_id=%d",
            gid,
            stats["seen_feature"],
            stats["kept"],
            stats["skip_no_core"],
            stats["skip_no_peptide"],
            stats["skip_bad_coord"],
            stats["skip_seqid_filter"],
            stats["dup_full_id"],
        )

        if len(intervals) == 0:
            logger.warning("物种 %s 产出 BED 为 0 行：通常意味着 GFF 的 core 无法在 peptide 中命中（请检查 header/core 规则）。", gid)

    logger.info("========== 完成 ==========")
    logger.info("GENESPACE 输入已生成：")
    logger.info("  peptide: %s（%d files）", OUT_PEPTIDE_DIR, len(genome_ids))
    logger.info("  bed    : %s（%d files, total lines=%d）", OUT_BED_DIR, len(genome_ids), total_bed_lines)
    logger.info("下一步可在 R 中对 wd=%s 运行 init_genespace + run_genespace。", OUT_ROOT)


if __name__ == "__main__":
    main()

