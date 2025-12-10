#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_00_build_chr_rename_table.py
—— 染色体长度统计 + 自动筛选大染色体并重命名为 Chr1..N（v2 · 16 物种版）

职责概述：
  1) 从 raw_data/synteny_species_meta.tsv 读取物种列表（只要有 species_id 就行）；
  2) 对每个物种：
       - 自动在 raw_data/gff/ 下寻找对应的 GFF 文件
         （支持 .gff / .gff3 / .gff.gz / .gff3.gz 后缀）；
       - 基于 GFF 中 feature 的坐标估计每个 seqid 的长度；
       - 按长度筛选出“主染色体”（长度 ≥ 阈值），其余视为小片段；
       - 对主染色体按长度从大到小排序，重命名为 Chr1, Chr2, ...；
       - 输出每物种的 chr 长度表 + 重命名表；
  3) 汇总所有物种到一个 summary 表；
  4) 全程输出详细日志到屏幕与日志文件。

重要说明：
  - 本脚本完全不涉及进化树信息；
  - synteny_species_meta.tsv 只需要最简结构：
        species_id  group  is_reference
    本脚本当前只用到 species_id 列，其余列后续脚本使用。

使用前准备：
  - 保证本脚本位于 magic/scripts/ 目录下；
  - 物种 meta 表位于 magic/raw_data/synteny_species_meta.tsv；
  - GFF 文件位于 magic/raw_data/gff/，文件名以 species_id 为前缀；
  - 若目录结构不同，可在脚本顶部参数区自行修改。
"""

from __future__ import annotations

import os
import sys
import csv
import gzip
import shutil
import logging
from pathlib import Path
from typing import Dict, Tuple, List


# =========================
# 参数区（皇上可在此修改）
# =========================

# 项目根目录（默认自动推断为脚本所在目录的上一级：magic/）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 原始数据与输出目录
RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
GFF_DIR = RAW_DATA_DIR / "gff"
SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_00_chr_rename"

# 染色体长度阈值（bp），默认 10 Mb，低于此长度的序列不会被视为“主染色体”
CHR_LENGTH_THRESHOLD_BP = 10_000_000

# 是否基于 GFF 估计长度（推荐 True）
USE_GFF_FOR_LENGTH = True

# 每个物种最多保留的主染色体数量（None 表示不限制）
MAX_CHR_TO_KEEP: int | None = None

# 日志等级：DEBUG / INFO / WARNING / ERROR
LOG_LEVEL = "INFO"


# =========================
# 工具函数
# =========================

def setup_logging(log_dir: Path) -> logging.Logger:
    """
    初始化日志系统：同时输出到屏幕与日志文件。
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny_00_build_chr_rename_table.log"

    logger = logging.getLogger("synteny_chr_rename")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清空旧 handler，避免重复输出
    for handler in list(logger.handlers):
        logger.removeHandler(handler)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    logger.info("========== synteny_00 — 染色体长度统计 + 重命名 ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("RAW_DATA_DIR = %s", RAW_DATA_DIR)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info("CHR_LENGTH_THRESHOLD_BP = %d", CHR_LENGTH_THRESHOLD_BP)
    logger.info("USE_GFF_FOR_LENGTH      = %s", USE_GFF_FOR_LENGTH)
    logger.info("MAX_CHR_TO_KEEP         = %s", str(MAX_CHR_TO_KEEP))

    return logger


def clean_output_root(output_root: Path, logger: logging.Logger) -> None:
    """
    删除旧的输出目录并重建。
    """
    if output_root.exists():
        logger.info("删除旧的输出目录：%s", output_root)
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)


def load_species_meta(meta_file: Path, logger: logging.Logger) -> List[str]:
    """
    读取 synteny_species_meta.tsv，返回需要处理的 species_id 列表。

    约定：
      - 至少包含一列 species_id；
      - group / is_reference 等列若不存在，本脚本也不会报错。
    """
    if not meta_file.exists():
        logger.error("物种 meta 文件不存在：%s", meta_file)
        sys.exit(1)

    logger.info("读取物种 meta 表：%s", meta_file)
    species_ids: List[str] = []
    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("synteny_species_meta.tsv 必须包含列：species_id")
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                # 空行直接跳过
                continue
            species_ids.append(sid)

    logger.info("总计 %d 个物种将用于 synteny_00（按表格行顺序处理）。", len(species_ids))
    return species_ids


def find_gff_for_species(species_id: str, logger: logging.Logger) -> Path | None:
    """
    根据 species_id 在 GFF_DIR 中尝试匹配 GFF 文件，支持多种后缀：
      - .gff3
      - .gff
      - .gff3.gz
      - .gff.gz

    优先级顺序：
      .gff3 > .gff > .gff3.gz > .gff.gz
    """
    candidates = [
        GFF_DIR / f"{species_id}.gff3",
        GFF_DIR / f"{species_id}.gff",
        GFF_DIR / f"{species_id}.gff3.gz",
        GFF_DIR / f"{species_id}.gff.gz",
    ]
    for path in candidates:
        if path.exists():
            logger.info("物种 %s 使用 GFF 文件：%s", species_id, path)
            return path

    # 宽松匹配：允许目录里已经存在完全相同前缀+后缀组合
    pattern_list = [
        f"{species_id}.gff3",
        f"{species_id}.gff",
        f"{species_id}.gff3.gz",
        f"{species_id}.gff.gz",
    ]
    try:
        for fname in os.listdir(GFF_DIR):
            for pat in pattern_list:
                if fname == pat:
                    path = GFF_DIR / fname
                    logger.info("物种 %s 使用 GFF 文件（宽松匹配）：%s", species_id, path)
                    return path
    except FileNotFoundError:
        logger.error("GFF 目录不存在：%s", GFF_DIR)
        return None

    logger.error(
        "未能在 %s 下找到物种 %s 对应的 GFF 文件（尝试后缀：.gff/.gff3/.gff.gz/.gff3.gz）",
        GFF_DIR,
        species_id,
    )
    return None


def estimate_lengths_from_gff(gff_path: Path, logger: logging.Logger) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    基于 GFF 估计每个 seqid 的长度（取该 seqid 上所有 feature 的最大 end 值），
    同时统计每个 seqid 上 feature 的计数（近似 n_genes）。

    支持普通文本 GFF 和 .gz 压缩 GFF。
    """
    lengths: Dict[str, int] = {}
    counts: Dict[str, int] = {}

    if not gff_path.exists():
        logger.error("GFF 文件不存在：%s", gff_path)
        raise FileNotFoundError(str(gff_path))

    is_gz = gff_path.suffix == ".gz"

    open_func = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    with open_func(gff_path, mode, encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            seqid = parts[0]
            try:
                end = int(parts[4])
            except ValueError:
                continue

            prev_len = lengths.get(seqid, 0)
            if end > prev_len:
                lengths[seqid] = end
            counts[seqid] = counts.get(seqid, 0) + 1

    return lengths, counts


def select_and_rename_chromosomes(
    seq_lengths: Dict[str, int],
    seq_counts: Dict[str, int],
    logger: logging.Logger,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """
    根据长度阈值筛选“主染色体”并按长度排序重命名为 Chr1..N。

    返回两个列表：
      1) chr_lengths_records：每个 seqid 的长度和 n_genes 信息；
      2) chr_rename_records：每个 seqid 对应的 is_chromosome 和 new_chr_name 等信息。
    """
    chr_lengths_records: List[Dict[str, object]] = []
    for seqid, length in seq_lengths.items():
        chr_lengths_records.append(
            {
                "seqid_raw": seqid,
                "length_bp": int(length),
                "n_genes": int(seq_counts.get(seqid, 0)),
            }
        )

    # 按长度由大到小排序
    chr_lengths_records.sort(key=lambda r: r["length_bp"], reverse=True)

    # 选出满足阈值的候选“主染色体”
    candidates = [r for r in chr_lengths_records if r["length_bp"] >= CHR_LENGTH_THRESHOLD_BP]

    if MAX_CHR_TO_KEEP is not None and MAX_CHR_TO_KEEP > 0:
        candidates = candidates[:MAX_CHR_TO_KEEP]

    chr_rename_records: List[Dict[str, object]] = []
    seqid_to_new: Dict[str, str] = {}
    rank = 0
    for rec in candidates:
        rank += 1
        new_name = f"Chr{rank}"
        seqid_raw = rec["seqid_raw"]
        seqid_to_new[seqid_raw] = new_name

    # 补充所有 seqid 的信息
    for rec in chr_lengths_records:
        seqid_raw = rec["seqid_raw"]
        length_bp = rec["length_bp"]
        n_genes = rec["n_genes"]
        if seqid_raw in seqid_to_new:
            is_chr = "yes"
            new_name = seqid_to_new[seqid_raw]
            rank_val = int(new_name.replace("Chr", "")) if new_name.startswith("Chr") else None
        else:
            is_chr = "no"
            new_name = ""
            rank_val = None

        chr_rename_records.append(
            {
                "seqid_raw": seqid_raw,
                "is_chromosome": is_chr,
                "length_bp": int(length_bp),
                "new_chr_name": new_name,
                "rank": rank_val,
                "n_genes": int(n_genes),
            }
        )

    return chr_lengths_records, chr_rename_records


def write_tsv(path: Path, fieldnames: List[str], records: List[Dict[str, object]]) -> None:
    """
    将记录列表写出为 TSV 文件。
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


# =========================
# 主流程
# =========================

def main() -> None:
    logger = setup_logging(OUTPUT_ROOT / "logs")
    clean_output_root(OUTPUT_ROOT, logger)

    species_ids = load_species_meta(SPECIES_META_FILE, logger)
    if not species_ids:
        logger.error("synteny_species_meta.tsv 中没有任何有效的 species_id，无法继续。")
        sys.exit(1)

    summary_records: List[Dict[str, object]] = []

    for species_id in species_ids:
        logger.info("------ 处理物种：%s ------", species_id)
        gff_path = find_gff_for_species(species_id, logger)
        if gff_path is None:
            logger.error("跳过物种 %s，因为未找到对应 GFF 文件。", species_id)
            continue

        try:
            if USE_GFF_FOR_LENGTH:
                seq_lengths, seq_counts = estimate_lengths_from_gff(gff_path, logger)
            else:
                # 当前版本仍默认用 GFF 估长，后续如需支持 genome.fa 可在此拓展
                seq_lengths, seq_counts = estimate_lengths_from_gff(gff_path, logger)
        except Exception as e:
            logger.error("解析 GFF 时出错（%s）：%s", species_id, str(e))
            continue

        if not seq_lengths:
            logger.warning("物种 %s 未在 GFF 中解析出任何 seqid，跳过。", species_id)
            continue

        chr_lengths_records, chr_rename_records = select_and_rename_chromosomes(
            seq_lengths, seq_counts, logger
        )

        lengths_path = OUTPUT_ROOT / f"chr_lengths_{species_id}.tsv"
        rename_path = OUTPUT_ROOT / f"chr_rename_{species_id}.tsv"

        write_tsv(
            lengths_path,
            fieldnames=["seqid_raw", "length_bp", "n_genes"],
            records=chr_lengths_records,
        )
        write_tsv(
            rename_path,
            fieldnames=[
                "species_id",
                "seqid_raw",
                "is_chromosome",
                "length_bp",
                "new_chr_name",
                "rank",
                "n_genes",
            ],
            records=[
                {"species_id": species_id, **rec}
                for rec in chr_rename_records
            ],
        )

        n_seqid_total = len(chr_lengths_records)
        n_chr_kept = sum(1 for r in chr_rename_records if r["is_chromosome"] == "yes")
        total_len_all = sum(r["length_bp"] for r in chr_lengths_records)
        total_len_kept = sum(r["length_bp"] for r in chr_rename_records if r["is_chromosome"] == "yes")
        kept_fraction = total_len_kept / total_len_all if total_len_all > 0 else 0.0

        summary_records.append(
            {
                "species_id": species_id,
                "n_seqid_total": n_seqid_total,
                "n_chromosomes_kept": n_chr_kept,
                "length_threshold_bp": CHR_LENGTH_THRESHOLD_BP,
                "total_length_kept_bp": total_len_kept,
                "total_length_all_bp": total_len_all,
                "kept_fraction": f"{kept_fraction:.4f}",
                "comment": "",
            }
        )

        logger.info(
            "物种 %s: 总 seqid=%d, 保留染色体=%d, 保留长度比例=%.2f%%",
            species_id,
            n_seqid_total,
            n_chr_kept,
            kept_fraction * 100.0,
        )

    summary_path = OUTPUT_ROOT / "chr_rename_summary.tsv"
    write_tsv(
        summary_path,
        fieldnames=[
            "species_id",
            "n_seqid_total",
            "n_chromosomes_kept",
            "length_threshold_bp",
            "total_length_kept_bp",
            "total_length_all_bp",
            "kept_fraction",
            "comment",
        ],
        records=summary_records,
    )

    logger.info("写出汇总表：%s", summary_path)
    logger.info("synteny_00 完成。")


if __name__ == "__main__":
    main()

