#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_04_build_layout_and_colors.py
—— 为 16 物种线性染色体共线性堆叠图构建布局与颜色表

职责概要：
  1) 读取物种 meta（synteny_species_meta.tsv）与 chr_order_<species>.tsv：
       - 确定物种的垂直顺序（tracks）
       - 确定参考物种的染色体列表（ref_chr）
  2) 自动为每个物种计算：
       - y_center：在整张图中的垂直中心位置（0~FIGURE_HEIGHT）
       - track_height：该物种染色体“带”的高度
       - short_label：短标签（如 Sco, Sri）
       - plot_label：长标签（如 S.constricta）
  3) 为 digger / non_digger 两类物种生成背景 block：
       - y_min / y_max 主要用于后续 Inkscape / AI 调整
  4) 为参考物种的每条 ref_chr 分配颜色：
       - 使用皇上的“御用配色”循环分配给 Chr1..ChrN
       - 生成 color_table_ref_chr.tsv

输入：
  - raw_data/synteny_species_meta.tsv
      当前至少包含：
        species_id
        group         (digger / non_digger / 其它)
        is_reference  (yes / no)
  - output/synteny_03_chr_order/chr_order_<species_id>.tsv
      至少包含字段：
        species_id
        rank
        chr
        dominant_ref_chr
        dominant_fraction
        orientation
        n_anchors_total
        note

输出（每次运行前清空本脚本对应输出目录）：
  - output/synteny_04_layout/
      ├── layout_species_tracks.tsv
      ├── layout_background_blocks.tsv
      ├── color_table_ref_chr.tsv
      └── logs/synteny_04_build_layout_and_colors.log

说明：
  - 本脚本不直接调用 JCVI，只生成后续绘图需要的布局与颜色元数据；
  - 不写死物种名称，完全由 synteny_species_meta.tsv 驱动；
  - y 坐标为相对坐标（0 ~ FIGURE_HEIGHT），方便后续在任意画布中映射。
"""

from __future__ import annotations

import sys
import csv
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any

# =========================
# 参数区（皇上可在此修改）
# =========================

# 项目根目录（假定本脚本位于 magic/scripts/ 下）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

CHR_ORDER_DIR = PROJECT_ROOT / "output" / "synteny_03_chr_order"

OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_04_layout"

LOG_LEVEL = "INFO"

# 图像总高度（相对单位），所有 y 坐标都会被压缩到 [0, FIGURE_HEIGHT]
FIGURE_HEIGHT = 1.0

# 期望的单个 track 高度 & track 间距（会按物种数量自动缩放）
TRACK_HEIGHT_BASE = 0.03
TRACK_VERTICAL_SPACING_BASE = 0.02

# 背景块在上下方向额外扩展的 margin（相对 FIGURE_HEIGHT 的比例）
BACKGROUND_MARGIN_FRACTION = 0.01

# 皇上御用配色（用于参考物种的 ref_chr 颜色循环）
ROYAL_COLORS = [
    "#E64B35",  # Maroon
    "#4DBBD5",  # SkyBlue
    "#00A087",  # Teal
    "#3C5488",  # Navy
    "#F39B7F",  # Light Orange
    "#8491B4",  # Slate Blue
    "#808180",  # Dark Gray
]

# 挖掘组 / 非挖掘组背景颜色（稍微偏淡，用于打底）
DIGGER_BG_COLOR = "#E6FFF0"
NON_DIGGER_BG_COLOR = "#FFE6F0"
OTHER_BG_COLOR = "#F0F0F0"


# =========================
# 基础工具：日志 & 输出目录
# =========================

def setup_logging() -> logging.Logger:
    """
    初始化日志系统：屏幕 + 文件双通道输出。
    """
    log_dir = OUTPUT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny_04_build_layout_and_colors.log"

    logger = logging.getLogger("synteny_layout_colors")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清理旧 handler
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

    logger.info("========== synteny_04 — layout & color table ==========")
    logger.info("PROJECT_ROOT      = %s", PROJECT_ROOT)
    logger.info("CHR_ORDER_DIR     = %s", CHR_ORDER_DIR)
    logger.info("OUTPUT_ROOT       = %s", OUTPUT_ROOT)
    logger.info("FIGURE_HEIGHT     = %.3f", FIGURE_HEIGHT)

    return logger


def clean_output_root(logger: logging.Logger) -> None:
    """
    按皇上要求，每次运行前清空本脚本对应输出目录，再重新创建。
    """
    if OUTPUT_ROOT.exists():
        logger.info("删除旧输出目录：%s", OUTPUT_ROOT)
        shutil.rmtree(OUTPUT_ROOT)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "logs").mkdir(parents=True, exist_ok=True)


# =========================
# 工具：解析 meta & chr_order
# =========================

def load_species_meta(meta_file: Path, logger: logging.Logger) -> Tuple[List[Dict[str, str]], str]:
    """
    读取 synteny_species_meta.tsv：
      - 保留原始行顺序，作为物种的垂直堆叠顺序；
      - 找到参考物种 ref_species_id。

    当前已知列：
      - species_id
      - group        (digger / non_digger / 其它)
      - is_reference (yes / no)

    返回：
      - species_info_list: 每个元素是 {"species_id", "group", "is_reference"}；
      - ref_species_id: 参考物种 ID。
    """
    if not meta_file.exists():
        logger.error("meta 文件不存在：%s", meta_file)
        sys.exit(1)

    species_info_list: List[Dict[str, str]] = []
    ref_species_id: str | None = None

    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("meta 表缺少列：species_id")
            sys.exit(1)
        if "group" not in reader.fieldnames:
            logger.error("meta 表缺少列：group")
            sys.exit(1)
        if "is_reference" not in reader.fieldnames:
            logger.error("meta 表缺少列：is_reference")
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            group = (row.get("group") or "").strip()
            is_ref = (row.get("is_reference") or "").strip().lower()
            species_info_list.append(
                {
                    "species_id": sid,
                    "group": group,
                    "is_reference": is_ref,
                }
            )
            if is_ref == "yes":
                if ref_species_id is not None:
                    logger.error("meta 表中 is_reference=yes 超过一个：%s, %s", ref_species_id, sid)
                    sys.exit(1)
                ref_species_id = sid

    if not species_info_list:
        logger.error("meta 表中没有有效的 species_id")
        sys.exit(1)
    if ref_species_id is None:
        logger.error("未在 meta 表中找到 is_reference=yes 的参考物种")
        sys.exit(1)

    logger.info("物种总数（用于 synteny） = %d", len(species_info_list))
    logger.info("参考物种 = %s", ref_species_id)
    return species_info_list, ref_species_id


def parse_chr_numeric_index(chr_name: str) -> int:
    """
    从染色体名中提取数字部分用于排序：
      - Chr1 -> 1
      - Chr10 -> 10
      - 若无法提取数字，则返回一个较大的值以放在后面。
    """
    name = chr_name.strip()
    digits = "".join(ch for ch in name if ch.isdigit())
    if digits:
        try:
            return int(digits)
        except Exception:
            return 100000
    return 100000


def auto_short_label(species_id: str) -> str:
    """
    自动生成 short_label，例如：
      Sinonovacula_constricta -> Sco
      Argopecten_irradians    -> Ari
      Mytilus_edulis          -> Med
    简单规则：
      - 拆分为属名 + 种名；
      - 属名首字母 + 种名前两字母；
      - 若格式不规范，则退化为 species_id 的前三个非下划线字符。
    """
    parts = species_id.split("_")
    if len(parts) >= 2:
        genus = parts[0]
        species = parts[1]
        g = genus[0].upper() if genus else "X"
        s1 = species[0].lower() if len(species) >= 1 else "x"
        s2 = species[1].lower() if len(species) >= 2 else ""
        return f"{g}{s1}{s2}"
    # fallback
    compact = "".join(ch for ch in species_id if ch.isalpha())
    if len(compact) >= 3:
        return compact[:3]
    return species_id[:3]


def auto_plot_label(species_id: str) -> str:
    """
    自动生成 plot_label，例如：
      Sinonovacula_constricta -> S.constricta
      Argopecten_irradians    -> A.irradians
    """
    parts = species_id.split("_")
    if len(parts) >= 2:
        genus = parts[0]
        species = parts[1]
        initial = genus[0].upper() if genus else "X"
        return f"{initial}.{species}"
    return species_id


def load_chr_order_for_species(
    species_id: str,
    logger: logging.Logger,
) -> List[Dict[str, Any]]:
    """
    读取 chr_order_<species_id>.tsv，返回按 rank 排序的记录列表。
    """
    path = CHR_ORDER_DIR / f"chr_order_{species_id}.tsv"
    if not path.exists():
        logger.error("未找到 chr_order 文件：%s", path)
        sys.exit(1)

    records: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"species_id", "rank", "chr", "dominant_ref_chr", "orientation"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            logger.error("chr_order_%s.tsv 缺少必要列：%s", species_id, ",".join(sorted(missing)))
            sys.exit(1)

        for row in reader:
            try:
                row["rank"] = int(row.get("rank", "0"))
            except Exception:
                row["rank"] = 0
            records.append(row)

    # 按 rank 正序排序
    records.sort(key=lambda r: r["rank"])
    logger.info("读取 chr_order_%s.tsv：染色体数 = %d", species_id, len(records))
    return records


# =========================
# 计算 tracks 布局 & 背景块 & 颜色表
# =========================

def compute_track_positions(
    species_info_list: List[Dict[str, str]],
    logger: logging.Logger,
) -> Dict[str, Dict[str, Any]]:
    """
    为每个物种计算 track 的 y_center 与 track_height。

    逻辑：
      - 将物种视为一条垂直堆叠链，按 meta 中出现顺序从上到下；
      - 使用 TRACK_HEIGHT_BASE 和 TRACK_VERTICAL_SPACING_BASE 计算“期望总高度”；
      - 然后根据 FIGURE_HEIGHT 做等比例缩放；
      - 确保所有 tracks 都被压缩到 [0, FIGURE_HEIGHT] 范围内。
    """
    n = len(species_info_list)
    if n <= 0:
        logger.error("没有物种，无法计算 tracks 布局")
        sys.exit(1)

    # 期望的总高度（含上下间隔）
    total_height_nominal = (
        n * (TRACK_HEIGHT_BASE + TRACK_VERTICAL_SPACING_BASE)
        + TRACK_VERTICAL_SPACING_BASE
    )
    if total_height_nominal <= 0:
        total_height_nominal = 1.0

    scale = FIGURE_HEIGHT / total_height_nominal

    track_height = TRACK_HEIGHT_BASE * scale
    step = (TRACK_HEIGHT_BASE + TRACK_VERTICAL_SPACING_BASE) * scale
    top_margin = TRACK_VERTICAL_SPACING_BASE * scale

    logger.info("tracks 布局缩放比例 scale = %.4f", scale)
    logger.info("实际 track_height = %.4f, step = %.4f", track_height, step)

    track_info: Dict[str, Dict[str, Any]] = {}

    # i=0 对应 meta 的第一行物种（最顶部）
    for i, sp in enumerate(species_info_list):
        sid = sp["species_id"]
        # 从上往下：第 i 条 track 的中心
        y_center = FIGURE_HEIGHT - (top_margin + track_height / 2.0 + i * step)

        info = {
            "species_id": sid,
            "group": sp["group"],
            "order_index": i + 1,
            "y_center": y_center,
            "track_height": track_height,
            "short_label": auto_short_label(sid),
            "plot_label": auto_plot_label(sid),
        }
        track_info[sid] = info
        logger.info(
            "track: %s  order_index=%d  y_center=%.4f  track_height=%.4f  group=%s",
            sid,
            info["order_index"],
            info["y_center"],
            info["track_height"],
            info["group"],
        )

    return track_info


def build_background_blocks(
    track_info: Dict[str, Dict[str, Any]],
    logger: logging.Logger,
) -> List[Dict[str, Any]]:
    """
    根据 group（digger / non_digger / 其它），构建背景块：

    block 粒度：
      - digger 物种整体一个块；
      - non_digger 物种整体一个块；
      - 其它分组可选（暂不强制分块）。

    每个 block 输出：
      - block_id
      - group
      - y_min
      - y_max
      - fill_color
      - note
    """
    # 按 group 聚合 y 范围
    group_to_ys: Dict[str, List[Tuple[float, float]]] = {}
    for sid, info in track_info.items():
        group = (info.get("group") or "").strip() or "unknown"
        yc = float(info["y_center"])
        h = float(info["track_height"])
        y_min = yc - h / 2.0
        y_max = yc + h / 2.0
        group_to_ys.setdefault(group, []).append((y_min, y_max))

    blocks: List[Dict[str, Any]] = []
    margin = FIGURE_HEIGHT * BACKGROUND_MARGIN_FRACTION

    def pick_color(group: str) -> str:
        g = group.lower()
        if g == "digger":
            return DIGGER_BG_COLOR
        if g == "non_digger":
            return NON_DIGGER_BG_COLOR
        return OTHER_BG_COLOR

    block_id_counter = 0
    for group, ymin_ymax_list in group_to_ys.items():
        if not ymin_ymax_list:
            continue
        ymin = min(y[0] for y in ymin_ymax_list) - margin
        ymax = max(y[1] for y in ymin_ymax_list) + margin

        # 保证仍然在 [0, FIGURE_HEIGHT] 范围内
        ymin = max(0.0, ymin)
        ymax = min(FIGURE_HEIGHT, ymax)

        block_id_counter += 1
        bid = f"block_{block_id_counter:02d}"

        blocks.append(
            {
                "block_id": bid,
                "group": group,
                "y_min": f"{ymin:.4f}",
                "y_max": f"{ymax:.4f}",
                "fill_color": pick_color(group),
                "note": "",
            }
        )

        logger.info(
            "背景块：%s  group=%s  y_min=%.4f  y_max=%.4f",
            bid,
            group,
            ymin,
            ymax,
        )

    return blocks


def build_ref_chr_color_table(
    ref_species_id: str,
    chr_order_records: Dict[str, List[Dict[str, Any]]],
    logger: logging.Logger,
) -> List[Dict[str, Any]]:
    """
    为参考物种的每条染色体（ref_chr）分配颜色。

    逻辑：
      - 读取 chr_order_<ref>.tsv 中所有 chr 名；
      - 按 Chr 数字排序；
      - 按顺序循环使用皇上的御用配色；
      - 输出 color_table_ref_chr.tsv：
          ref_chr, color_hex, ALG_id, note
    """
    if ref_species_id not in chr_order_records:
        logger.error("参考物种 %s 未在 chr_order_records 中找到", ref_species_id)
        sys.exit(1)

    records = chr_order_records[ref_species_id]
    # 提取所有 chr（避免重复）
    chr_set = {r["chr"] for r in records if (r.get("chr") or "").strip()}
    chr_list = sorted(chr_set, key=lambda c: (parse_chr_numeric_index(c), c))

    logger.info("参考物种 %s 染色体列表用于配色：%s", ref_species_id, ", ".join(chr_list))

    color_table: List[Dict[str, Any]] = []
    n_colors = len(ROYAL_COLORS)
    if n_colors == 0:
        logger.error("ROYAL_COLORS 为空，无法分配颜色")
        sys.exit(1)

    for i, chr_name in enumerate(chr_list):
        color = ROYAL_COLORS[i % n_colors]
        color_table.append(
            {
                "ref_chr": chr_name,
                "color_hex": color,
                "ALG_id": chr_name,
                "note": "",
            }
        )
        logger.info("ref_chr=%s  color=%s", chr_name, color)

    return color_table


# =========================
# 写出结果文件
# =========================

def write_layout_species_tracks(
    track_info: Dict[str, Dict[str, Any]],
    species_info_list: List[Dict[str, str]],
    logger: logging.Logger,
) -> None:
    """
    写出 layout_species_tracks.tsv：
      - species_id
      - short_label
      - plot_label
      - order_index
      - group
      - y_center
      - track_height
      - include_in_synteny
    """
    out_path = OUTPUT_ROOT / "layout_species_tracks.tsv"
    fieldnames = [
        "species_id",
        "short_label",
        "plot_label",
        "order_index",
        "group",
        "y_center",
        "track_height",
        "include_in_synteny",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        # 按 meta 中顺序输出
        for sp in species_info_list:
            sid = sp["species_id"]
            info = track_info[sid]
            writer.writerow(
                {
                    "species_id": sid,
                    "short_label": info["short_label"],
                    "plot_label": info["plot_label"],
                    "order_index": info["order_index"],
                    "group": info["group"],
                    "y_center": f"{info['y_center']:.4f}",
                    "track_height": f"{info['track_height']:.4f}",
                    "include_in_synteny": "yes",
                }
            )
    logger.info("写出 layout_species_tracks.tsv：%s", out_path)


def write_layout_background_blocks(
    blocks: List[Dict[str, Any]],
    logger: logging.Logger,
) -> None:
    """
    写出 layout_background_blocks.tsv：
      - block_id
      - group
      - y_min
      - y_max
      - fill_color
      - note
    """
    out_path = OUTPUT_ROOT / "layout_background_blocks.tsv"
    fieldnames = ["block_id", "group", "y_min", "y_max", "fill_color", "note"]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in blocks:
            writer.writerow(rec)
    logger.info("写出 layout_background_blocks.tsv：%s", out_path)


def write_color_table_ref_chr(
    color_table: List[Dict[str, Any]],
    logger: logging.Logger,
) -> None:
    """
    写出 color_table_ref_chr.tsv：
      - ref_chr
      - color_hex
      - ALG_id
      - note
    """
    out_path = OUTPUT_ROOT / "color_table_ref_chr.tsv"
    fieldnames = ["ref_chr", "color_hex", "ALG_id", "note"]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in color_table:
            writer.writerow(rec)
    logger.info("写出 color_table_ref_chr.tsv：%s", out_path)


# =========================
# 主流程
# =========================

def main() -> None:
    logger = setup_logging()
    clean_output_root(logger)

    # 1) 读取物种 meta + 找出参考物种
    species_info_list, ref_species_id = load_species_meta(SPECIES_META_FILE, logger)

    # 2) 为每个物种读取 chr_order 文件（用于配色 & 校对）
    chr_order_records: Dict[str, List[Dict[str, Any]]] = {}
    for sp in species_info_list:
        sid = sp["species_id"]
        chr_order_records[sid] = load_chr_order_for_species(sid, logger)

    # 3) 计算 tracks 垂直布局（y_center & track_height）
    track_info = compute_track_positions(species_info_list, logger)

    # 4) 生成 digger / non_digger 等背景块
    blocks = build_background_blocks(track_info, logger)

    # 5) 为参考物种的染色体生成颜色表
    color_table = build_ref_chr_color_table(ref_species_id, chr_order_records, logger)

    # 6) 写出各种 layout / color 表
    write_layout_species_tracks(track_info, species_info_list, logger)
    write_layout_background_blocks(blocks, logger)
    write_color_table_ref_chr(color_table, logger)

    logger.info("synteny_04 完成。")


if __name__ == "__main__":
    main()

