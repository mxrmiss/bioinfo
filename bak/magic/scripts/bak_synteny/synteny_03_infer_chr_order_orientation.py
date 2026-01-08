#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_03_infer_chr_order_orientation.py
—— 基于 anchors_global 推断各物种染色体对应的 ALG（参考染色体）、
   染色体间对应关系、顺序和方向，并生成 JCVI 可用的 seqids 文件。

职责概述：
  1) 读取 synteny_species_meta.tsv 获取物种列表和参考物种；
  2) 读取 synteny_01_global_anchors/anchors_global.tsv：
       - 每条记录包含 og_id, species_id, chr, start, ref_chr, ref_start 等；
  3) 统计每个物种、每条 chr 对每条 ref_chr 的锚点数量与构成比例：
       - 输出 chr_alg_composition_<species_id>.tsv
  4) 对每个物种的每条 chr：
       - 找出“主导 ref_chr”（anchors 数量最多的那个 ref_chr）；
       - 计算该 chr 上主导 ref_chr 的占比（dominant_fraction）；
       - 计算该 chr 与主导 ref_chr 上 anchor 坐标的相关系数，
         判断该 chr 在绘图时应正向（+）还是反向（-）画；
       - 输出 chr_order_<species_id>.tsv
  5) 确定参考物种的染色体顺序，并据此为所有物种生成 seqids：
       - 对参考物种：按 Chr1, Chr2, Chr3... 的顺序；
       - 对其它物种：按 dominant_ref_chr 映射到参考顺序；
       - 写出 seqids/<species_id>.seqids，供 JCVI karyotype 使用；
  6) 输出 chr_order_overview.tsv 作为整体 QC。

输入：
  - raw_data/synteny_species_meta.tsv
      必须至少包含：
        species_id
        is_reference （yes / no）

  - output/synteny_01_global_anchors/anchors_global.tsv
      必须至少包含：
        og_id
        species_id
        chr
        start
        ref_chr
        ref_start

输出（每次运行前自动清空本脚本对应输出目录）：
  - output/synteny_03_chr_order/
      ├── chr_alg_composition_<species_id>.tsv
      ├── chr_order_<species_id>.tsv
      ├── seqids/<species_id>.seqids
      ├── chr_order_overview.tsv
      └── logs/synteny_03_infer_chr_order_orientation.log

核心参数：
  - MIN_DOMINANT_FRACTION   ：为了认为某条 chr 有清晰的主 ALG，需要的最小占比
  - MIN_ANCHORS_PER_CHR     ：考虑主 ALG 与方向判断所需的最小锚点数量
  - CORRELATION_THRESHOLD   ：相关系数绝对值低于该值，认为方向不可靠（orientation='?')

所有路径与参数均在脚本顶部配置，不依赖命令行参数。
"""

from __future__ import annotations

import sys
import csv
import math
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any

# =========================
# 参数区（皇上可在此修改）
# =========================

# 假定本脚本位于 magic/scripts/ 下
PROJECT_ROOT = Path(__file__).resolve().parent.parent

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

ANCHORS_GLOBAL_DIR = PROJECT_ROOT / "output" / "synteny_01_global_anchors"
ANCHORS_GLOBAL_FILE = ANCHORS_GLOBAL_DIR / "anchors_global.tsv"

OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_03_chr_order"
SEQIDS_DIR = OUTPUT_ROOT / "seqids"

LOG_LEVEL = "INFO"

# 认为某条 chr 拥有“主 ALG”的最小占比（比如 >= 0.5）
MIN_DOMINANT_FRACTION = 0.5

# 一条 chr 至少有多少 anchors 才参与判断主 ALG 和方向
MIN_ANCHORS_PER_CHR = 10

# 相关系数绝对值阈值，用于判断方向是否可靠
CORRELATION_THRESHOLD = 0.3


# =========================
# 工具函数：日志与输出目录
# =========================

def setup_logging() -> logging.Logger:
    """
    初始化日志系统：屏幕 + 文件双通道输出。
    """
    log_dir = OUTPUT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny_03_infer_chr_order_orientation.log"

    logger = logging.getLogger("synteny_chr_order")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    for h in list(logger.handlers):
        logger.removeHandler(h)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)

    logger.info("========== synteny_03 — 染色体顺序与方向推断 ==========")
    logger.info("PROJECT_ROOT           = %s", PROJECT_ROOT)
    logger.info("ANCHORS_GLOBAL_FILE    = %s", ANCHORS_GLOBAL_FILE)
    logger.info("OUTPUT_ROOT            = %s", OUTPUT_ROOT)
    logger.info("MIN_DOMINANT_FRACTION  = %.2f", MIN_DOMINANT_FRACTION)
    logger.info("MIN_ANCHORS_PER_CHR    = %d", MIN_ANCHORS_PER_CHR)
    logger.info("CORRELATION_THRESHOLD  = %.2f", CORRELATION_THRESHOLD)

    return logger


def clean_output_root(logger: logging.Logger) -> None:
    """
    每次运行前清空本脚本对应的输出目录，避免旧文件干扰。
    """
    if OUTPUT_ROOT.exists():
        logger.info("删除旧输出目录：%s", OUTPUT_ROOT)
        shutil.rmtree(OUTPUT_ROOT)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    SEQIDS_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "logs").mkdir(parents=True, exist_ok=True)


# =========================
# 工具函数：读取 meta 与 anchors
# =========================

def load_species_meta(meta_file: Path, logger: logging.Logger) -> Tuple[List[str], str]:
    """
    读取 synteny_species_meta.tsv：
      - species_id 按行顺序用于后续输出顺序；
      - is_reference=yes 的为参考物种。

    要求：
      - 至少一列：species_id
      - 一列：is_reference（yes/no），且仅有一个 yes
    """
    if not meta_file.exists():
        logger.error("meta 文件不存在：%s", meta_file)
        sys.exit(1)

    species_ids: List[str] = []
    ref_species_id: str | None = None

    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("meta 表缺少列：species_id")
            sys.exit(1)
        if "is_reference" not in reader.fieldnames:
            logger.error("meta 表缺少列：is_reference")
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            species_ids.append(sid)
            is_ref = (row.get("is_reference") or "").strip().lower()
            if is_ref == "yes":
                if ref_species_id is not None:
                    logger.error("meta 表中 is_reference=yes 超过一个：%s, %s", ref_species_id, sid)
                    sys.exit(1)
                ref_species_id = sid

    if not species_ids:
        logger.error("meta 表中没有有效的 species_id")
        sys.exit(1)
    if ref_species_id is None:
        logger.error("meta 表中未设置 is_reference=yes 的参考物种")
        sys.exit(1)

    logger.info("meta 物种数量 = %d", len(species_ids))
    logger.info("参考物种 = %s", ref_species_id)
    return species_ids, ref_species_id


def safe_int(value: Any, default: int = 0) -> int:
    """
    安全地将字符串转换为 int，失败则返回 default。
    """
    try:
        return int(value)
    except Exception:
        return default


def compute_pearson_correlation(points: List[Tuple[int, int]]) -> float:
    """
    计算 (x, y) 点集合的皮尔逊相关系数，用于判断方向：
      - r > 0 代表 x 与 y 同向排序（正向）；
      - r < 0 代表 x 与 y 反向排序（需要翻转）；
      - |r| 很小时认为方向不可靠。

    参数：
      - points: [(ref_start, species_start), ...]

    返回：
      - r: 浮点数，范围通常在 [-1, 1] 之间；
      - 若点数不足或方差为 0，则返回 0.0。
    """
    n = len(points)
    if n < 2:
        return 0.0

    xs = [float(p[0]) for p in points]
    ys = [float(p[1]) for p in points]

    mean_x = sum(xs) / n
    mean_y = sum(ys) / n

    num = 0.0
    sum_x2 = 0.0
    sum_y2 = 0.0

    for x, y in zip(xs, ys):
        dx = x - mean_x
        dy = y - mean_y
        num += dx * dy
        sum_x2 += dx * dx
        sum_y2 += dy * dy

    if sum_x2 <= 0 or sum_y2 <= 0:
        return 0.0

    denom = math.sqrt(sum_x2 * sum_y2)
    if denom == 0.0:
        return 0.0

    r = num / denom
    # 理论上 r 在 [-1, 1]，但数值误差可能略超出，做个裁剪
    if r > 1.0:
        r = 1.0
    elif r < -1.0:
        r = -1.0
    return r


def load_anchors_global(
    anchors_file: Path,
    logger: logging.Logger,
) -> Tuple[
    Dict[str, Dict[str, Dict[str, Dict[str, Any]]]],
    Dict[str, Dict[str, int]],
]:
    """
    读取 anchors_global.tsv 并构建数据结构：

    species_chr_ref_data:
      species_id -> chr -> ref_chr -> {
          "n_anchors": int,
          "coords": [(ref_start, start), ...]
      }

    species_refchr_totals:
      species_id -> ref_chr -> total_anchors（该物种该 ref_chr 的总锚点数）

    这些数据后续用于：
      - 统计 chr 与 ref_chr 的组成；
      - 判断主导 ref_chr；
      - 判断方向（相关性）。
    """
    if not anchors_file.exists():
        logger.error("未找到 anchors_global.tsv：%s", anchors_file)
        sys.exit(1)

    species_chr_ref_data: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = {}
    species_refchr_totals: Dict[str, Dict[str, int]] = {}

    n_records = 0

    logger.info("读取 anchors_global：%s", anchors_file)
    with anchors_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"species_id", "chr", "start", "ref_chr", "ref_start"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            logger.error("anchors_global.tsv 缺少必要列：%s", ",".join(sorted(missing)))
            sys.exit(1)

        for row in reader:
            n_records += 1
            sid = (row.get("species_id") or "").strip()
            chr_name = (row.get("chr") or "").strip()
            ref_chr = (row.get("ref_chr") or "").strip()
            if not sid or not chr_name or not ref_chr:
                continue

            start = safe_int(row.get("start", 0), 0)
            ref_start = safe_int(row.get("ref_start", 0), 0)

            sp_map = species_chr_ref_data.setdefault(sid, {})
            chr_map = sp_map.setdefault(chr_name, {})
            ref_entry = chr_map.setdefault(ref_chr, {"n_anchors": 0, "coords": []})

            ref_entry["n_anchors"] += 1
            ref_entry["coords"].append((ref_start, start))

            sp_ref_tot = species_refchr_totals.setdefault(sid, {})
            sp_ref_tot[ref_chr] = sp_ref_tot.get(ref_chr, 0) + 1

    logger.info("anchors_global 记录数（用于 chr/ALG 统计） = %d", n_records)
    return species_chr_ref_data, species_refchr_totals


# =========================
# 核心：计算组成、顺序和方向
# =========================

def parse_chr_numeric_index(chr_name: str) -> int:
    """
    从染色体名字中提取数字部分，用于排序：
      - 例如：Chr1 -> 1, Chr10 -> 10；
      - 若无法解析数字，则返回一个较大的值（100000），并最终按字母排序兜底。
    """
    name = chr_name.strip()
    num_part = ""
    for ch in name:
        if ch.isdigit():
            num_part += ch
    if num_part:
        try:
            return int(num_part)
        except Exception:
            return 100000
    return 100000


def compute_chr_alg_composition_and_order(
    species_ids: List[str],
    ref_species_id: str,
    species_chr_ref_data: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]],
    species_refchr_totals: Dict[str, Dict[str, int]],
    logger: logging.Logger,
) -> Tuple[
    Dict[str, List[Dict[str, Any]]],
    Dict[str, List[Dict[str, Any]]],
    Dict[str, Dict[str, Any]],
    Dict[str, int],
]:
    """
    计算：
      - 每个物种每条 chr 的 ALG 组成（chr_alg_composition_*）
      - 每个物种每条 chr 的主 ALG + 方向 + 排序信息（chr_order_*）
      - 全局概览（chr_order_overview）
      - 参考物种的 Chr 顺序映射（ref_chr_order_map）

    返回：
      chr_alg_records:
        species_id -> [ {species_id, chr, ref_chr, n_anchors, fraction_of_chr_anchors, fraction_of_alg_anchors}, ... ]

      chr_order_records:
        species_id -> [ {species_id, rank, chr, dominant_ref_chr, dominant_fraction, orientation, n_anchors_total, note}, ... ]

      overview_records_dict:
        species_id -> {species_id, n_chr, n_chr_with_clear_dominant_ref_chr, n_chr_mixed, comment}

      ref_chr_order_map:
        ref_chr -> order_index（1,2,3,...）
    """
    chr_alg_records: Dict[str, List[Dict[str, Any]]] = {}
    chr_order_records: Dict[str, List[Dict[str, Any]]] = {}
    overview_records_dict: Dict[str, Dict[str, Any]] = {}

    # 1) 先确定参考物种 ref_chr 的顺序
    ref_chr_order_map: Dict[str, int] = {}
    ref_sp_chrs = species_chr_ref_data.get(ref_species_id, {})
    if not ref_sp_chrs:
        logger.error("参考物种 %s 在 anchors_global 中没有任何 chr 记录", ref_species_id)
        sys.exit(1)

    # 参考物种的 Chr 列表
    ref_chr_list = sorted(ref_sp_chrs.keys(), key=lambda c: (parse_chr_numeric_index(c), c))
    for idx, chr_name in enumerate(ref_chr_list, start=1):
        ref_chr_order_map[chr_name] = idx

    logger.info("参考物种 %s 的 Chr 顺序：%s", ref_species_id, ", ".join(ref_chr_list))

    # 2) 逐物种计算
    for sid in species_ids:
        sp_chr_map = species_chr_ref_data.get(sid, {})
        sp_ref_tot = species_refchr_totals.get(sid, {})

        chr_list = sorted(sp_chr_map.keys(), key=lambda c: (parse_chr_numeric_index(c), c))
        sp_chr_alg_list: List[Dict[str, Any]] = []
        sp_chr_order_list: List[Dict[str, Any]] = []

        n_chr_total = len(chr_list)
        n_chr_clear = 0
        n_chr_mixed = 0

        logger.info("------ 物种 %s：chr 数量 = %d ------", sid, n_chr_total)

        # 先计算每个 chr 的组成与方向信息（但 rank 要等知道 ref_chr_order_map 后统一排序）
        temp_chr_info: Dict[str, Dict[str, Any]] = {}

        for chr_name in chr_list:
            ref_map = sp_chr_map.get(chr_name, {})
            if not ref_map:
                # 理论上不会出现，因为没有 anchors 不会出现在 sp_chr_map
                continue

            # 这一条 chr 的总 anchors 数
            total_anchors_chr = sum(entry["n_anchors"] for entry in ref_map.values())
            if total_anchors_chr <= 0:
                continue

            # 统计 ALG 组成
            for ref_chr, entry in ref_map.items():
                n_anch = entry["n_anchors"]
                if n_anch <= 0:
                    continue
                total_for_alg = sp_ref_tot.get(ref_chr, 0)
                frac_chr = n_anch / total_anchors_chr if total_anchors_chr > 0 else 0.0
                frac_alg = n_anch / total_for_alg if total_for_alg > 0 else 0.0

                sp_chr_alg_list.append(
                    {
                        "species_id": sid,
                        "chr": chr_name,
                        "ref_chr": ref_chr,
                        "n_anchors": n_anch,
                        "fraction_of_chr_anchors": f"{frac_chr:.4f}",
                        "fraction_of_alg_anchors": f"{frac_alg:.44f}",
                    }
                )

            # 找主导 ref_chr
            dominant_ref_chr = None
            dominant_n = 0
            for ref_chr, entry in ref_map.items():
                n_anch = entry["n_anchors"]
                if n_anch > dominant_n:
                    dominant_n = n_anch
                    dominant_ref_chr = ref_chr

            dominant_fraction = (
                dominant_n / total_anchors_chr if total_anchors_chr > 0 else 0.0
            )

            # 判断方向
            orientation = "?"
            note_parts: List[str] = []

            if total_anchors_chr < MIN_ANCHORS_PER_CHR:
                note_parts.append(f"few_anchors(<{MIN_ANCHORS_PER_CHR})")

            if dominant_ref_chr is None or dominant_n <= 0:
                note_parts.append("no_dominant_ref_chr")
            else:
                if dominant_fraction < MIN_DOMINANT_FRACTION:
                    note_parts.append(f"mixed_ref_chr(<{MIN_DOMINANT_FRACTION:.2f})")
                # 只有 anchor 足够且有主导 ALG 时才尝试算相关性
                if total_anchors_chr >= MIN_ANCHORS_PER_CHR and dominant_fraction >= 0.0:
                    coords = ref_map[dominant_ref_chr]["coords"]
                    r = compute_pearson_correlation(coords)
                    if abs(r) < CORRELATION_THRESHOLD:
                        orientation = "?"
                        note_parts.append(f"weak_corr(|r|<{CORRELATION_THRESHOLD:.2f})")
                    else:
                        orientation = "+" if r >= 0 else "-"
                        note_parts.append(f"corr={r:.3f}")

            # 统计“清晰”的 chr
            is_clear = (
                dominant_ref_chr is not None
                and dominant_fraction >= MIN_DOMINANT_FRACTION
                and total_anchors_chr >= MIN_ANCHORS_PER_CHR
            )
            if is_clear:
                n_chr_clear += 1
            else:
                n_chr_mixed += 1

            temp_chr_info[chr_name] = {
                "species_id": sid,
                "chr": chr_name,
                "dominant_ref_chr": dominant_ref_chr or "NA",
                "dominant_fraction": dominant_fraction,
                "orientation": orientation,
                "n_anchors_total": total_anchors_chr,
                "note": ";".join(note_parts) if note_parts else "",
            }

        # 根据参考物种的 ref_chr 顺序，确定本物种 chr 的绘图顺序（rank）
        # 逻辑：
        #   1) 优先：dominant_ref_chr 在 ref_chr_order_map 中，且满足“清晰”条件；
        #   2) 剩下的 chr 排在后面，按 anchor 数量降序 + chr 名排序。
        good_chrs: List[Tuple[Tuple[int, int, str], Dict[str, Any]]] = []
        bad_chrs: List[Tuple[Tuple[int, int, str], Dict[str, Any]]] = []

        for chr_name in chr_list:
            info = temp_chr_info.get(chr_name)
            if info is None:
                continue
            dom_ref = info["dominant_ref_chr"]
            total_anch = info["n_anchors_total"]
            dom_frac = info["dominant_fraction"]

            if (
                dom_ref in ref_chr_order_map
                and dom_ref != "NA"
                and dom_frac >= MIN_DOMINANT_FRACTION
                and total_anch >= MIN_ANCHORS_PER_CHR
            ):
                key = (ref_chr_order_map[dom_ref], -total_anch, chr_name)
                good_chrs.append((key, info))
            else:
                key = (999999, -total_anch, chr_name)
                bad_chrs.append((key, info))

        good_chrs.sort(key=lambda x: x[0])
        bad_chrs.sort(key=lambda x: x[0])

        ranked_infos: List[Dict[str, Any]] = []
        rank_counter = 0

        for key, info in good_chrs + bad_chrs:
            rank_counter += 1
            record = {
                "species_id": sid,
                "rank": rank_counter,
                "chr": info["chr"],
                "dominant_ref_chr": info["dominant_ref_chr"],
                "dominant_fraction": f"{info['dominant_fraction']:.4f}",
                "orientation": info["orientation"],
                "n_anchors_total": info["n_anchors_total"],
                "note": info["note"],
            }
            ranked_infos.append(record)

        chr_alg_records[sid] = sp_chr_alg_list
        chr_order_records[sid] = ranked_infos

        comment = ""
        if n_chr_clear == 0:
            comment = "no_clear_chr"
        elif n_chr_mixed > 0:
            comment = "mixed_exists"
        else:
            comment = "all_clear"

        overview_records_dict[sid] = {
            "species_id": sid,
            "n_chr": n_chr_total,
            "n_chr_with_clear_dominant_ref_chr": n_chr_clear,
            "n_chr_mixed": n_chr_mixed,
            "comment": comment,
        }

        logger.info(
            "物种 %s：n_chr = %d, 清晰 chr = %d, 其它 chr = %d, comment = %s",
            sid,
            n_chr_total,
            n_chr_clear,
            n_chr_mixed,
            comment,
        )

    return chr_alg_records, chr_order_records, overview_records_dict, ref_chr_order_map


# =========================
# 写出结果文件
# =========================

def write_chr_alg_composition(
    chr_alg_records: Dict[str, List[Dict[str, Any]]],
    logger: logging.Logger,
) -> None:
    """
    写出每个物种的 chr_alg_composition_<species_id>.tsv。
    """
    for sid, records in chr_alg_records.items():
        out_path = OUTPUT_ROOT / f"chr_alg_composition_{sid}.tsv"
        fieldnames = [
            "species_id",
            "chr",
            "ref_chr",
            "n_anchors",
            "fraction_of_chr_anchors",
            "fraction_of_alg_anchors",
        ]
        with out_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for rec in records:
                writer.writerow(rec)
        logger.info("写出 ALG 组成表：%s （记录数 = %d）", out_path, len(records))


def write_chr_order_and_seqids(
    chr_order_records: Dict[str, List[Dict[str, Any]]],
    logger: logging.Logger,
) -> None:
    """
    写出每个物种的 chr_order_<sid>.tsv 与 seqids/<sid>.seqids。
    """
    for sid, records in chr_order_records.items():
        # 1) 写 chr_order_<sid>.tsv
        order_path = OUTPUT_ROOT / f"chr_order_{sid}.tsv"
        fieldnames = [
            "species_id",
            "rank",
            "chr",
            "dominant_ref_chr",
            "dominant_fraction",
            "orientation",
            "n_anchors_total",
            "note",
        ]
        with order_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for rec in records:
                writer.writerow(rec)
        logger.info("写出 chr 顺序表：%s （记录数 = %d）", order_path, len(records))

        # 2) 写 seqids/<sid>.seqids
        seqids_path = SEQIDS_DIR / f"{sid}.seqids"
        with seqids_path.open("w", encoding="utf-8", newline="") as f:
            for rec in records:
                chr_name = rec["chr"]
                orient = (rec["orientation"] or "").strip()
                if orient not in {"+", "-"}:
                    # 若方向不确定，则不加符号，由后续绘图或手动调整
                    line = chr_name
                else:
                    line = f"{chr_name}{orient}"
                f.write(line + "\n")
        logger.info("写出 seqids 文件：%s", seqids_path)


def write_overview(
    overview_records_dict: Dict[str, Dict[str, Any]],
    logger: logging.Logger,
) -> None:
    """
    写出 chr_order_overview.tsv。
    """
    out_path = OUTPUT_ROOT / "chr_order_overview.tsv"
    fieldnames = [
        "species_id",
        "n_chr",
        "n_chr_with_clear_dominant_ref_chr",
        "n_chr_mixed",
        "comment",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sid in sorted(overview_records_dict.keys()):
            writer.writerow(overview_records_dict[sid])
    logger.info("写出 chr_order_overview.tsv：%s", out_path)


# =========================
# 主流程
# =========================

def main() -> None:
    logger = setup_logging()
    clean_output_root(logger)

    # 1) 读取物种顺序和参考物种
    species_ids, ref_species_id = load_species_meta(SPECIES_META_FILE, logger)

    # 2) 读取 anchors_global，构建 species_chr_ref_data & species_refchr_totals
    species_chr_ref_data, species_refchr_totals = load_anchors_global(
        ANCHORS_GLOBAL_FILE, logger
    )

    # 3) 计算 ALG 组成、chr 顺序和方向、概览、以及 ref_chr 顺序映射
    (
        chr_alg_records,
        chr_order_records,
        overview_records_dict,
        ref_chr_order_map,
    ) = compute_chr_alg_composition_and_order(
        species_ids,
        ref_species_id,
        species_chr_ref_data,
        species_refchr_totals,
        logger,
    )

    logger.info("参考物种 %s 的 ref_chr 顺序映射：%s", ref_species_id, ref_chr_order_map)

    # 4) 写出各类结果文件
    write_chr_alg_composition(chr_alg_records, logger)
    write_chr_order_and_seqids(chr_order_records, logger)
    write_overview(overview_records_dict, logger)

    logger.info("synteny_03 完成。")


if __name__ == "__main__":
    main()

