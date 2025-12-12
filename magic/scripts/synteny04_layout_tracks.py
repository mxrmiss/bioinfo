#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny04_layout_tracks.py
—— 轨道 layout（Step 04）

当前职责（对应蓝图 Step 04）：
  1) 基于 synteny_species_meta.tsv 确定物种顺序与分组信息；
  2) 为每个物种分配纵向轨道：
       - order_index（从上到下）
       - y_center（0~1 之间）
       - track_height（轨道高度）
       - short_label / plot_label / include_in_synteny
  3) 基于 group（digger / non_digger）生成纵向背景区间：
       - layout_background_blocks.tsv（仅包含 group 与 y 范围，不含颜色信息）
  4) 汇总 layout 信息写出：
       - layout_species_tracks.tsv
       - layout_base_tracks.tsv

说明：
  - 本脚本不调用 jcvi，也不读写 anchors/simple；
  - 只依赖 Step 00 + Step 03 的结果，用于驱动 Step 05 生成 layout 文件，
    所有具体颜色由后续出图脚本统一分配。
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

# 原始 meta 表、上游输出目录
RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

CHR_RENAME_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"
CHR_ORDER_DIR = PROJECT_ROOT / "output" / "synteny_03_chr_order"

# 本脚本输出目录
OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_04_layout"

# 日志等级
LOG_LEVEL = "INFO"

# 纵向布局参数（0~1 范围）
# 整个图像纵向可用空间为 [VERTICAL_MARGIN, 1 - VERTICAL_MARGIN]
VERTICAL_MARGIN = 0.06

# 轨道高度控制：
#   若 FIXED_TRACK_HEIGHT 为 None，则自动根据物种数均匀分配高度：
#      track_height = step * AUTO_TRACK_HEIGHT_FRACTION
#   若 FIXED_TRACK_HEIGHT 为浮点数（例如 0.04），则所有物种使用同一高度。
FIXED_TRACK_HEIGHT: Optional[float] = None
AUTO_TRACK_HEIGHT_FRACTION = 0.6

# 背景块纵向扩展（在最上/最下物种轨道基础上略微扩展）
BACKGROUND_PADDING = 0.01

# 物种短名映射表（运行时自动生成，不再手动维护）
SHORT_NAME_MAP: Dict[str, str] = {}


# =========================
# 通用工具函数
# =========================

def setup_logging(log_dir: Path) -> logging.Logger:
    """
    初始化日志系统：同时输出到屏幕与日志文件。
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny04_layout_tracks.log"

    logger = logging.getLogger("synteny04")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清空旧 handler，避免重复
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

    logger.info("========== synteny04 — layout ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info("VERTICAL_MARGIN = %.3f", VERTICAL_MARGIN)
    logger.info("FIXED_TRACK_HEIGHT = %s", str(FIXED_TRACK_HEIGHT))
    logger.info("AUTO_TRACK_HEIGHT_FRACTION = %.3f", AUTO_TRACK_HEIGHT_FRACTION)

    return logger


def clean_output_root(output_root: Path) -> None:
    """
    删除旧的输出目录并重建子目录。
    """
    if output_root.exists():
        shutil.rmtree(output_root)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)


def get_short_name(species_id: str) -> str:
    """
    获取物种短名；优先使用自动生成的 SHORT_NAME_MAP。
    若未找到，则退化为首字母组合。
    """
    if species_id in SHORT_NAME_MAP:
        return SHORT_NAME_MAP[species_id]

    parts = species_id.split("_")
    if len(parts) >= 2:
        g, s = parts[0], parts[1]
        return (g[0] + s[:3]).capitalize()
    return species_id[:4]


# =========================
# 读取 meta 与 chr_order / chr_rename
# =========================

def load_species_meta(meta_file: Path, logger: logging.Logger) -> Tuple[List[Dict[str, str]], str]:
    """
    读取 synteny_species_meta.tsv，保持原始行顺序。

    返回：
      meta_rows: 每行为一个 dict，至少包含：
         - species_id
         - group
         - is_reference
      ref_species: 唯一参考物种 ID
    """
    if not meta_file.exists():
        logger.error("物种 meta 文件不存在：%s", meta_file)
        sys.exit(1)

    meta_rows: List[Dict[str, str]] = []
    ref_species: Optional[str] = None

    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("synteny_species_meta.tsv 必须包含列：species_id")
            sys.exit(1)

        if "group" not in reader.fieldnames:
            logger.error("synteny_species_meta.tsv 必须包含列：group")
            sys.exit(1)

        if "is_reference" not in reader.fieldnames:
            logger.error("synteny_species_meta.tsv 必须包含列：is_reference")
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            group = (row.get("group") or "").strip()
            is_ref = (row.get("is_reference") or "").strip().lower()

            meta_rows.append(
                {
                    "species_id": sid,
                    "group": group,
                    "is_reference": is_ref,
                }
            )

            if is_ref == "yes":
                if ref_species is not None and ref_species != sid:
                    logger.error("检测到多个 is_reference == yes 的物种：%s, %s",
                                 ref_species, sid)
                    sys.exit(1)
                ref_species = sid

    if ref_species is None:
        logger.error("synteny_species_meta.tsv 中未找到 is_reference == yes 的物种。")
        sys.exit(1)

    logger.info("物种 meta 读取完成：总计 %d 个物种，参考物种 = %s",
                len(meta_rows), ref_species)
    return meta_rows, ref_species


def build_short_name_map(species_ids: List[str], logger: logging.Logger) -> Dict[str, str]:
    """
    根据 species_id 自动生成短名：
      - 初始规则：Genus[0] + species[:2]，例如：
          Sinonovacula_constricta -> Sco
      - 若不同物种间出现重名，则自动增加更多 species 字母，
        直到在当前物种集合中唯一。

    生成的映射会写入日志，以便检查。
    """
    mapping: Dict[str, str] = {}
    used: Dict[str, str] = {}

    for sid in species_ids:
        parts = sid.split("_")
        if len(parts) >= 2:
            genus, species = parts[0], parts[1]
        else:
            genus, species = sid, ""

        if species:
            base = genus[0] + species[:2]
        else:
            base = (genus[:3] or sid[:3])

        cand = base.capitalize()
        extra = 2
        while cand in used and used[cand] != sid:
            extra += 1
            if species and extra < len(species):
                cand = (genus[0] + species[:extra]).capitalize()
            else:
                cand = (genus[:min(len(genus), extra)] + species[:1]).capitalize()
                if cand in used and used[cand] != sid:
                    cand = (base + str(extra)).capitalize()

        mapping[sid] = cand
        used[cand] = sid

    logger.info("物种短名映射：%s",
                ", ".join(f"{k} -> {v}" for k, v in mapping.items()))
    return mapping


def check_chr_order_files(meta_rows: List[Dict[str, str]], logger: logging.Logger) -> None:
    """
    简单检查 Step 03 输出的 chr_order_<species_id>.tsv 是否存在。
    不深入解析内容，只做存在性检查，避免后续 Step 05 报错才发现。
    """
    missing: List[str] = []
    for row in meta_rows:
        sid = row["species_id"]
        path = CHR_ORDER_DIR / f"chr_order_{sid}.tsv"
        if not path.exists():
            missing.append(str(path))

    if missing:
        logger.error("以下 chr_order_<species_id>.tsv 文件缺失，请先完成 Step 03：")
        for p in missing:
            logger.error("  - %s", p)
        sys.exit(1)
    else:
        logger.info("Step 03 的 chr_order_<species_id>.tsv 已全部就绪。")


# =========================
# 轨道 layout 计算
# =========================

def compute_track_layout(meta_rows: List[Dict[str, str]], logger: logging.Logger
                         ) -> List[Dict[str, object]]:
    """
    根据物种 meta 行顺序，计算每个物种的轨道参数：

      - order_index（从 1 开始，meta 第一行在最上）
      - y_center（0~1）
      - track_height
      - short_label / plot_label / group / include_in_synteny

    返回 layout_rows 列表，每行为一个 dict。
    """
    n = len(meta_rows)
    if n == 0:
        logger.error("物种数量为 0，无法生成轨道 layout。")
        sys.exit(1)

    if VERTICAL_MARGIN < 0 or VERTICAL_MARGIN >= 0.5:
        logger.error("VERTICAL_MARGIN 参数不合理（应在 [0, 0.5) 之间）：%.3f", VERTICAL_MARGIN)
        sys.exit(1)

    available_height = 1.0 - 2 * VERTICAL_MARGIN
    if available_height <= 0:
        logger.error("VERTICAL_MARGIN 过大，导致可用纵向高度 <= 0。")
        sys.exit(1)

    step = available_height / float(n)
    logger.info("纵向布局：n_species = %d, available_height = %.3f, step = %.3f",
                n, available_height, step)

    if FIXED_TRACK_HEIGHT is not None:
        track_height = float(FIXED_TRACK_HEIGHT)
    else:
        track_height = step * AUTO_TRACK_HEIGHT_FRACTION

    if track_height <= 0 or track_height >= 1.0:
        logger.error("track_height 参数不合理（应在 (0,1) 之间）：%.3f", track_height)
        sys.exit(1)

    logger.info("轨道高度设置：track_height = %.3f", track_height)

    layout_rows: List[Dict[str, object]] = []

    # meta 行的顺序 = 从上到下的物种顺序
    for idx, row in enumerate(meta_rows, start=1):
        sid = row["species_id"]
        group = row["group"]
        is_ref = row["is_reference"]

        # y_center：越靠上的物种 y 越大（图上从上到下）
        # order_index = 1 -> 最上
        y_center = 1.0 - VERTICAL_MARGIN - step * (idx - 0.5)

        short_label = get_short_name(sid)
        plot_label = sid.replace("_", " ")

        include_in_synteny = "yes"  # 当前版本默认全部纳入共线性图

        layout_rows.append(
            {
                "species_id": sid,
                "short_label": short_label,
                "plot_label": plot_label,
                "group": group,
                "order_index": idx,
                "y_center": y_center,
                "track_height": track_height,
                "is_reference": is_ref,
                "include_in_synteny": include_in_synteny,
            }
        )

    logger.info("轨道 layout 已为 %d 个物种生成。", len(layout_rows))
    return layout_rows


# =========================
# 背景块生成（按 group）
# =========================

def build_background_blocks(
    layout_rows: List[Dict[str, object]],
    logger: logging.Logger,
) -> List[Dict[str, object]]:
    """
    根据物种的 group 与轨道位置，生成背景区间（不含颜色）：

      - 对每个 group（如 digger / non_digger）：
         - 找到该组所有物种的 y_center ± track_height/2 的范围；
         - 在该范围基础上向外扩展 BACKGROUND_PADDING；
         - 记录 group 与纵向范围，颜色由后续出图脚本决定。
    """
    by_group: Dict[str, List[Tuple[float, float]]] = defaultdict(list)
    for row in layout_rows:
        group = str(row["group"])
        y_center = float(row["y_center"])
        th = float(row["track_height"])
        y_min = y_center - th / 2.0
        y_max = y_center + th / 2.0
        by_group[group].append((y_min, y_max))

    blocks: List[Dict[str, object]] = []
    block_id_counter = 1

    for group, intervals in by_group.items():
        if not intervals:
            continue
        ys_min = min(v[0] for v in intervals)
        ys_max = max(v[1] for v in intervals)

        y_min = max(0.0, ys_min - BACKGROUND_PADDING)
        y_max = min(1.0, ys_max + BACKGROUND_PADDING)

        block_id = f"bg_{block_id_counter}"
        block_id_counter += 1

        blocks.append(
            {
                "block_id": block_id,
                "group": group,
                "y_min": y_min,
                "y_max": y_max,
                "note": f"group={group}",
            }
        )

    logger.info("背景块生成完成：共 %d 个 group 背景块。", len(blocks))
    return blocks


# =========================
# 写出 TSV 文件
# =========================

def write_layout_species_tracks(
    layout_rows: List[Dict[str, object]],
    out_path: Path,
    logger: logging.Logger,
) -> None:
    """
    写出 layout_species_tracks.tsv：
      - species_id
      - short_label
      - plot_label
      - group
      - order_index
      - y_center
      - track_height
      - include_in_synteny
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "species_id",
        "short_label",
        "plot_label",
        "group",
        "order_index",
        "y_center",
        "track_height",
        "include_in_synteny",
        "is_reference",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in layout_rows:
            writer.writerow(
                {
                    "species_id": row["species_id"],
                    "short_label": row["short_label"],
                    "plot_label": row["plot_label"],
                    "group": row["group"],
                    "order_index": int(row["order_index"]),
                    "y_center": float(row["y_center"]),
                    "track_height": float(row["track_height"]),
                    "include_in_synteny": row["include_in_synteny"],
                    "is_reference": row["is_reference"],
                }
            )
    logger.info("layout_species_tracks.tsv 已写出：%s", out_path)


def write_layout_base_tracks(
    layout_rows: List[Dict[str, object]],
    out_path: Path,
    logger: logging.Logger,
) -> None:
    """
    写出 layout_base_tracks.tsv：

      - track_index（与 order_index 保持一致）
      - species_id
      - short_label
      - plot_label
      - group
      - order_index
      - y_center
      - track_height
      - include_in_synteny
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "track_index",
        "species_id",
        "short_label",
        "plot_label",
        "group",
        "order_index",
        "y_center",
        "track_height",
        "include_in_synteny",
        "is_reference",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in layout_rows:
            order_index = int(row["order_index"])
            writer.writerow(
                {
                    "track_index": order_index,
                    "species_id": row["species_id"],
                    "short_label": row["short_label"],
                    "plot_label": row["plot_label"],
                    "group": row["group"],
                    "order_index": order_index,
                    "y_center": float(row["y_center"]),
                    "track_height": float(row["track_height"]),
                    "include_in_synteny": row["include_in_synteny"],
                    "is_reference": row["is_reference"],
                }
            )
    logger.info("layout_base_tracks.tsv 已写出：%s", out_path)


def write_background_blocks(
    bg_blocks: List[Dict[str, object]],
    out_path: Path,
    logger: logging.Logger,
) -> None:
    """
    写出 layout_background_blocks.tsv：
      - block_id
      - group
      - y_min
      - y_max
      - note
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "block_id",
        "group",
        "y_min",
        "y_max",
        "note",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in bg_blocks:
            writer.writerow(
                {
                    "block_id": row["block_id"],
                    "group": row["group"],
                    "y_min": float(row["y_min"]),
                    "y_max": float(row["y_max"]),
                    "note": row["note"],
                }
            )
    logger.info("layout_background_blocks.tsv 已写出：%s", out_path)


# =========================
# 主流程
# =========================

def main() -> None:
    # 1) 清空输出目录 + 初始化日志
    clean_output_root(OUTPUT_ROOT)
    logger = setup_logging(OUTPUT_ROOT / "logs")

    # 2) 读取物种 meta，并确认参考物种
    meta_rows, ref_species = load_species_meta(SPECIES_META_FILE, logger)
    logger.info("参考物种（仅用于后续步骤配色等）：%s", ref_species)

    # 2b) 自动构建物种短名映射
    species_ids = [row["species_id"] for row in meta_rows]
    global SHORT_NAME_MAP
    SHORT_NAME_MAP = build_short_name_map(species_ids, logger)

    # 3) 检查 Step 03 的 chr_order_<species_id>.tsv 是否存在
    check_chr_order_files(meta_rows, logger)

    # 4) 计算每个物种的轨道布局（纵向 y_center + track_height）
    layout_rows = compute_track_layout(meta_rows, logger)

    # 5) 基于 group 生成背景区间（不含颜色）
    bg_blocks = build_background_blocks(layout_rows, logger)

    # 6) 写出各类 TSV 文件
    layout_species_path = OUTPUT_ROOT / "layout_species_tracks.tsv"
    layout_base_path = OUTPUT_ROOT / "layout_base_tracks.tsv"
    bg_blocks_path = OUTPUT_ROOT / "layout_background_blocks.tsv"

    write_layout_species_tracks(layout_rows, layout_species_path, logger)
    write_layout_base_tracks(layout_rows, layout_base_path, logger)
    write_background_blocks(bg_blocks, bg_blocks_path, logger)

    logger.info("synteny04 完成：layout_species_tracks + background_blocks + layout_base_tracks 均已生成。")


if __name__ == "__main__":
    main()

