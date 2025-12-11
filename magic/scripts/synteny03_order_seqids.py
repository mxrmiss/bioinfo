#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny03_order_seqids.py
—— 染色体顺序推断 + seqids 生成（宏观共线性用）

当前职责（对应蓝图 Step 03）：
  1) 读取 synteny_species_meta.tsv，识别：
       - 全部 species_id
       - 唯一参考物种 ref_species（is_reference == "yes"）
  2) 读取 synteny_00_chr_rename 的 chr_rename_<species_id>.tsv，
     获取每个物种的主染色体列表及长度。
  3) 读取 synteny_02_blocks_macro/blocks_normalized 下
       <ref_short>__vs__<qry_short>.blocks.tsv，
     对每个「参考 vs 其它物种」：
       - 按 qry_chr × ref_chr 统计 blocks 总跨度（bp）以及方向；
       - 推断每条 qry_chr 的：
           dominant_ref_chr      —— 共线跨度最大的参考 Chr
           dominant_fraction     —— 该 Chr 被 dominant_ref_chr 覆盖的比例
           orientation           —— 主方向（+ / - / .）
           n_blocks              —— 共线块数量
           total_block_span_bp   —— 所有 blocks 在 qry_chr 上的总跨度
  4) 写出每个物种的染色体顺序表：
       synteny_03_chr_order/chr_order_<species_id>.tsv
  5) 根据 chr_order_* 生成：
       - seqids_species/<species_id>.seqids
       - 全局 seqids 文件
  6) 写出 chr_order_overview.tsv，概览每个物种的匹配质量。

说明：
  - 本脚本不调用 jcvi，也不再读 BLAST 文件；
    所有共线信息均来自 Step 02 的 blocks_normalized。
"""

from __future__ import annotations

import sys
import csv
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

# =========================
# 参数区（皇上可在此修改）
# =========================

# 项目根目录（默认：当前脚本所在目录的上一级 magic/）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 原始数据与上游输出目录
RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

CHR_RENAME_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"
STEP02_DIR = PROJECT_ROOT / "output" / "synteny_02_blocks_macro"
BLOCKS_NORMALIZED_DIR = STEP02_DIR / "blocks_normalized"

# 本脚本输出目录
OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_03_chr_order"
CHR_ORDER_DIR = OUTPUT_ROOT
SEQIDS_SPECIES_DIR = OUTPUT_ROOT / "seqids_species"

# 日志等级
LOG_LEVEL = "INFO"

# 判定「主参考染色体」与「混合染色体」的阈值
#   - dominant_fraction < DOMINANT_MIN_FRACTION   -> 认为没有可靠主 Chr
#   - DOMINANT_MIN_FRACTION <= fraction < DOMINANT_MIXED_CUTOFF -> mixed
DOMINANT_MIN_FRACTION = 0.30
DOMINANT_MIXED_CUTOFF = 0.70

# 物种短名映射（需与 02 脚本保持完全一致）
SPECIES_SHORT_NAME: Dict[str, str] = {
    "Sinonovacula_constricta": "Sco",
    "Sinonovacula_rivularis": "Sri",
    "Novaculina_chinensis": "Nch",
    "Panopea_generosa": "Pge",
    "Mya_arenaria": "Mar",
    "Meretrix_meretrix": "Mme",
    "Mercenaria_mercenaria": "Mmc",
    "Tegillarca_granosa": "Tgr",
    "Mytilus_edulis": "Med",
    "Mytilus_galloprovincialis": "Mga",
    "Pinctada_fucata": "Pfu",
    "Ostrea_edulis": "Oed",
    "Crassostrea_gigas": "Cgi",
    "Ctenoides_ales": "Cal",
    "Pecten_maximus": "Pma",
    "Argopecten_irradians": "Air",
}


# =========================
# 通用工具函数
# =========================

def setup_logging(log_dir: Path) -> logging.Logger:
    """初始化日志系统：同时输出到屏幕与日志文件。"""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny03_order_seqids.log"

    logger = logging.getLogger("synteny03")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    for h in list(logger.handlers):
        logger.removeHandler(h)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    ch = logging.StreamHandler(sys.stdout)

    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)

    fh.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))
    ch.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    logger.addHandler(fh)
    logger.addHandler(ch)

    logger.info("========== synteny03 — 染色体顺序 + seqids 构建 ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info("DOMINANT_MIN_FRACTION = %.2f", DOMINANT_MIN_FRACTION)
    logger.info("DOMINANT_MIXED_CUTOFF = %.2f", DOMINANT_MIXED_CUTOFF)

    return logger


def clean_output_root(output_root: Path) -> None:
    """删除旧的输出目录并重建子目录。"""
    if output_root.exists():
        shutil.rmtree(output_root)
    CHR_ORDER_DIR.mkdir(parents=True, exist_ok=True)
    SEQIDS_SPECIES_DIR.mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)


def load_species_meta(meta_file: Path, logger: logging.Logger) -> Tuple[List[str], str]:
    """读取 synteny_species_meta.tsv，返回 species_id 列表与参考物种。"""
    if not meta_file.exists():
        logger.error("物种 meta 文件不存在：%s", meta_file)
        sys.exit(1)

    species_ids: List[str] = []
    ref_species: Optional[str] = None

    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("synteny_species_meta.tsv 必须包含列：species_id")
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            species_ids.append(sid)
            is_ref = (row.get("is_reference") or "").strip().lower()
            if is_ref == "yes":
                if ref_species is not None and ref_species != sid:
                    logger.error("检测到多个 is_reference == yes 的物种：%s, %s",
                                 ref_species, sid)
                    sys.exit(1)
                ref_species = sid

    if ref_species is None:
        logger.error("synteny_species_meta.tsv 中未找到 is_reference == yes 的物种。")
        sys.exit(1)

    logger.info("总计 %d 个物种，参考物种 = %s", len(species_ids), ref_species)
    return species_ids, ref_species


def get_short_name(species_id: str) -> str:
    """获取物种短名；与 02 脚本保持一致。"""
    if species_id in SPECIES_SHORT_NAME:
        return SPECIES_SHORT_NAME[species_id]

    parts = species_id.split("_")
    if len(parts) >= 2:
        g, s = parts[0], parts[1]
        return (g[0] + s[:3]).capitalize()
    return species_id[:4]


def parse_chr_rename_for_species(
    species_id: str,
    logger: logging.Logger,
) -> List[Tuple[str, int]]:
    """
    读取 chr_rename_<species_id>.tsv，返回主染色体列表：
       [(chr_name, length_bp), ...]
    顺序按 rank 从小到大排序。
    """
    path = CHR_RENAME_DIR / f"chr_rename_{species_id}.tsv"
    if not path.exists():
        logger.error("找不到 chr_rename 文件：%s", path)
        sys.exit(1)

    rows: List[Tuple[int, str, int]] = []

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("is_chromosome") or "") != "yes":
                continue
            chr_name = (row.get("new_chr_name") or "").strip()
            if not chr_name:
                continue
            rank_str = row.get("rank") or ""
            try:
                rank = int(rank_str)
            except ValueError:
                # 若 rank 缺失，则用一个较大的值，后续按 chr 名再排序
                rank = 10_000
            length_bp = int(row.get("length_bp") or "0")
            rows.append((rank, chr_name, length_bp))

    if not rows:
        logger.warning("物种 %s 在 chr_rename 中没有标记任何主染色体。", species_id)

    rows.sort(key=lambda x: (x[0], x[1]))
    result = [(chr_name, length_bp) for _, chr_name, length_bp in rows]
    logger.info("物种 %s：chr_rename 中共解析到 %d 条主染色体。", species_id, len(result))
    return result


def parse_blocks_normalized_for_species(
    ref_species: str,
    qry_species: str,
    logger: logging.Logger,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    读取某个物种对的 blocks_normalized：
      <ref_short>__vs__<qry_short>.blocks.tsv

    返回结构：
      stats[qry_chr][ref_chr] = {
          "len_pos": float,
          "len_neg": float,
          "n_blocks": int,
      }
    其中 len_pos / len_neg 为在 qry_chr 上的跨度累积（bp）。
    """
    ref_short = get_short_name(ref_species)
    qry_short = get_short_name(qry_species)

    path = BLOCKS_NORMALIZED_DIR / f"{ref_short}__vs__{qry_short}.blocks.tsv"
    stats: Dict[str, Dict[str, Dict[str, float]]] = defaultdict(
        lambda: defaultdict(lambda: {"len_pos": 0.0, "len_neg": 0.0, "n_blocks": 0.0})
    )

    if not path.exists():
        logger.warning("blocks_normalized 缺失：%s，本物种将被视为没有共线信息。", path)
        return stats

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required_cols = {
            "ref_chr",
            "ref_start",
            "ref_end",
            "qry_chr",
            "qry_start",
            "qry_end",
            "orientation",
        }
        if not reader.fieldnames or not required_cols.issubset(reader.fieldnames):
            logger.error("blocks_normalized 表头缺失必要列：%s", path)
            return stats

        for row in reader:
            qry_chr = (row.get("qry_chr") or "").strip()
            ref_chr = (row.get("ref_chr") or "").strip()
            if not qry_chr or not ref_chr:
                continue
            try:
                qry_start = int(row.get("qry_start") or "0")
                qry_end = int(row.get("qry_end") or "0")
            except ValueError:
                continue
            if qry_end <= qry_start:
                continue
            span = float(qry_end - qry_start)
            orientation = (row.get("orientation") or "+").strip()

            cell = stats[qry_chr][ref_chr]
            if orientation == "-":
                cell["len_neg"] += span
            else:
                cell["len_pos"] += span
            cell["n_blocks"] += 1.0

    logger.info(
        "blocks_normalized(%s vs %s)：共统计到 %d 条 qry_chr。",
        ref_species,
        qry_species,
        len(stats),
    )
    return stats


def parse_chr_index_from_name(chr_name: str) -> int:
    """
    尝试从 Chr 名中解析数字部分，以便排序。
    例如：Chr1 -> 1；Chr12 -> 12；解析失败返回一个较大值。
    """
    name = chr_name
    if name.lower().startswith("chr"):
        name = name[3:]
    digits = []
    for ch in name:
        if ch.isdigit():
            digits.append(ch)
        else:
            break
    if not digits:
        return 10_000
    try:
        return int("".join(digits))
    except ValueError:
        return 10_000


# =========================
# 主流程函数
# =========================

def build_chr_order_for_species(
    species_id: str,
    ref_species: str,
    chr_list: List[Tuple[str, int]],
    blocks_stats: Optional[Dict[str, Dict[str, Dict[str, float]]]],
    logger: logging.Logger,
) -> Tuple[List[Dict[str, object]], List[str]]:
    """
    为某个物种构建 chr_order_<species_id>.tsv 所需记录。

    输入：
      - species_id    当前物种
      - ref_species   参考物种（缢蛏）
      - chr_list      来自 chr_rename，[(chr_name, length_bp), ...]
      - blocks_stats  若为参考物种则传 None，其他物种传 parse_blocks_* 的结果

    返回：
      - records  —— 写入 chr_order_<species_id>.tsv 的字典列表
      - chr_order —— 按 rank 排序后的 chr_name 列表（用于 seqids）
    """
    records: List[Dict[str, object]] = []

    # 参考物种：直接按 Chr1..N 的顺序输出
    if species_id == ref_species or blocks_stats is None:
        rank = 1
        for chr_name, length_bp in chr_list:
            rec = {
                "species_id": species_id,
                "rank": rank,
                "chr": chr_name,
                "dominant_ref_chr": chr_name if species_id == ref_species else "",
                "dominant_fraction": 1.0 if species_id == ref_species else 0.0,
                "orientation": "+",
                "n_blocks": 0,
                "total_block_span_bp": length_bp if species_id == ref_species else 0,
                "note": "reference" if species_id == ref_species else "no_synteny",
            }
            records.append(rec)
            rank += 1

        chr_order = [c for c, _ in chr_list]
        return records, chr_order

    # 其它物种：根据 blocks_stats 推断每条 Chr 的主参考 Chr
    # 先预汇总：每条 qry_chr 的总跨度与主 ref_chr
    chr_stats: Dict[str, Dict[str, object]] = {}

    for chr_name, length_bp in chr_list:
        cell_stats = blocks_stats.get(chr_name, {})
        total_span = 0.0
        best_ref_chr = ""
        best_len = 0.0
        best_len_pos = 0.0
        best_len_neg = 0.0
        total_blocks = 0.0

        for ref_chr, v in cell_stats.items():
            len_pos = v.get("len_pos", 0.0)
            len_neg = v.get("len_neg", 0.0)
            span = len_pos + len_neg
            total_span += span
            total_blocks += v.get("n_blocks", 0.0)
            if span > best_len:
                best_len = span
                best_ref_chr = ref_chr
                best_len_pos = len_pos
                best_len_neg = len_neg

        if total_span <= 0.0 or best_len <= 0.0:
            dominant_ref_chr = ""
            dominant_fraction = 0.0
            orientation = "."
        else:
            dominant_ref_chr = best_ref_chr
            dominant_fraction = best_len / total_span
            if dominant_fraction < DOMINANT_MIN_FRACTION:
                # 共线太分散，不认为有可靠主 Chr
                dominant_ref_chr = ""
                orientation = "."
            else:
                if best_len_pos >= best_len_neg:
                    orientation = "+"
                else:
                    orientation = "-"

        chr_stats[chr_name] = {
            "length_bp": length_bp,
            "dominant_ref_chr": dominant_ref_chr,
            "dominant_fraction": dominant_fraction,
            "orientation": orientation,
            "n_blocks": int(total_blocks),
            "total_block_span_bp": int(total_span),
        }

    # 按「主参考 Chr 的编号」+ 本身 Chr 编号排序，得到 rank
    def order_key(chr_name: str) -> Tuple[int, int, str]:
        d = chr_stats.get(chr_name, {})
        ref_chr = d.get("dominant_ref_chr") or ""
        ref_idx = parse_chr_index_from_name(ref_chr) if ref_chr else 10_000
        self_idx = parse_chr_index_from_name(chr_name)
        return (ref_idx, self_idx, chr_name)

    sorted_chr_names = sorted([c for c, _ in chr_list], key=order_key)

    records = []
    for rank, chr_name in enumerate(sorted_chr_names, start=1):
        d = chr_stats.get(chr_name, {})
        rec = {
            "species_id": species_id,
            "rank": rank,
            "chr": chr_name,
            "dominant_ref_chr": d.get("dominant_ref_chr", ""),
            "dominant_fraction": float(d.get("dominant_fraction", 0.0)),
            "orientation": d.get("orientation", "."),
            "n_blocks": int(d.get("n_blocks", 0)),
            "total_block_span_bp": int(d.get("total_block_span_bp", 0)),
            "note": "",
        }
        records.append(rec)

    return records, sorted_chr_names


def write_chr_order_tsv(
    species_id: str,
    records: List[Dict[str, object]],
    logger: logging.Logger,
) -> Path:
    """写出单个物种的 chr_order_<species_id>.tsv。"""
    out_path = CHR_ORDER_DIR / f"chr_order_{species_id}.tsv"
    fieldnames = [
        "species_id",
        "rank",
        "chr",
        "dominant_ref_chr",
        "dominant_fraction",
        "orientation",
        "n_blocks",
        "total_block_span_bp",
        "note",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)

    logger.info("chr_order_%s.tsv 已写出：%s", species_id, out_path)
    return out_path


def write_seqids_species(
    species_id: str,
    chr_order: List[str],
    logger: logging.Logger,
) -> Path:
    """写出 seqids_species/<species_id>.seqids。"""
    out_path = SEQIDS_SPECIES_DIR / f"{species_id}.seqids"
    with out_path.open("w", encoding="utf-8") as f:
        f.write(f">{species_id}\n")
        for chr_name in chr_order:
            f.write(f"{chr_name}\n")

    logger.info("seqids_species/%s.seqids 已写出：%s", species_id, out_path)
    return out_path


def write_global_seqids(
    species_ids: List[str],
    chr_order_map: Dict[str, List[str]],
    logger: logging.Logger,
) -> Path:
    """写出全局 seqids 文件（按照 meta 中的物种顺序）。"""
    out_path = OUTPUT_ROOT / "seqids"
    with out_path.open("w", encoding="utf-8") as f:
        for sid in species_ids:
            chrs = chr_order_map.get(sid, [])
            if not chrs:
                logger.warning("物种 %s 在 chr_order_map 中没有找到 chr 列表，跳过。", sid)
                continue
            f.write(f">{sid}\n")
            for chr_name in chrs:
                f.write(f"{chr_name}\n")

    logger.info("全局 seqids 已写出：%s", out_path)
    return out_path


def write_chr_order_overview(
    species_ids: List[str],
    ref_species: str,
    chr_order_records: Dict[str, List[Dict[str, object]]],
    logger: logging.Logger,
) -> Path:
    """写出 chr_order_overview.tsv，统计每个物种的匹配情况。"""
    out_path = OUTPUT_ROOT / "chr_order_overview.tsv"
    fieldnames = [
        "species_id",
        "n_chr",
        "n_chr_with_dominant_ref",
        "n_chr_mixed",
        "comment",
    ]

    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for sid in species_ids:
            recs = chr_order_records.get(sid, [])
            n_chr = len(recs)
            n_with_dom = 0
            n_mixed = 0

            for r in recs:
                frac = float(r.get("dominant_fraction") or 0.0)
                dom_ref = (r.get("dominant_ref_chr") or "").strip()
                if not dom_ref:
                    continue
                if frac >= DOMINANT_MIXED_CUTOFF:
                    n_with_dom += 1
                elif frac >= DOMINANT_MIN_FRACTION:
                    n_mixed += 1

            if sid == ref_species:
                comment = "reference"
            else:
                comment = ""

            writer.writerow(
                {
                    "species_id": sid,
                    "n_chr": n_chr,
                    "n_chr_with_dominant_ref": n_with_dom,
                    "n_chr_mixed": n_mixed,
                    "comment": comment,
                }
            )

    logger.info("chr_order_overview.tsv 已写出：%s", out_path)
    return out_path


# =========================
# main
# =========================

def main() -> None:
    # 1) 清空输出目录 + 初始化日志
    clean_output_root(OUTPUT_ROOT)
    logger = setup_logging(OUTPUT_ROOT / "logs")

    # 2) 读取物种列表与参考物种
    species_ids, ref_species = load_species_meta(SPECIES_META_FILE, logger)

    # 3) 为每个物种读取 chr_rename 中的主染色体列表
    chr_list_map: Dict[str, List[Tuple[str, int]]] = {}
    for sid in species_ids:
        chr_list = parse_chr_rename_for_species(sid, logger)
        chr_list_map[sid] = chr_list

    # 4) 为每个非参考物种读取 blocks_normalized 统计
    blocks_stats_map: Dict[str, Dict[str, Dict[str, Dict[str, float]]]] = {}
    for sid in species_ids:
        if sid == ref_species:
            continue
        stats = parse_blocks_normalized_for_species(ref_species, sid, logger)
        blocks_stats_map[sid] = stats

    # 5) 构建 chr_order_<species_id>.tsv 与 per-species seqids
    chr_order_map: Dict[str, List[str]] = {}
    chr_order_records: Dict[str, List[Dict[str, object]]] = {}

    for sid in species_ids:
        chr_list = chr_list_map.get(sid, [])
        if not chr_list:
            logger.warning("物种 %s 没有主染色体信息，跳过。", sid)
            continue

        if sid == ref_species:
            records, chr_order = build_chr_order_for_species(
                species_id=sid,
                ref_species=ref_species,
                chr_list=chr_list,
                blocks_stats=None,
                logger=logger,
            )
        else:
            stats = blocks_stats_map.get(sid, {})
            records, chr_order = build_chr_order_for_species(
                species_id=sid,
                ref_species=ref_species,
                chr_list=chr_list,
                blocks_stats=stats,
                logger=logger,
            )

        write_chr_order_tsv(sid, records, logger)
        write_seqids_species(sid, chr_order, logger)

        chr_order_map[sid] = chr_order
        chr_order_records[sid] = records

    # 6) 写出全局 seqids
    write_global_seqids(species_ids, chr_order_map, logger)

    # 7) 写出 chr_order_overview.tsv
    write_chr_order_overview(species_ids, ref_species, chr_order_records, logger)

    logger.info("synteny03 完成：chr_order + seqids 已全部生成。")


if __name__ == "__main__":
    main()

