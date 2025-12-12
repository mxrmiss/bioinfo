#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny05_plot_karyotype.py
—— JCVI karyotype 多物种宏观共线性最终绘图脚本（05 终版 · 颜色完全由本脚本决定）

职责概述：
  1) 读取 layout_species_tracks.tsv，确定物种顺序、参考物种、轨道位置；
  2) 读取 synteny_03 产生的各物种 seqids，合并生成全局 seqids 文件；
  3) 读取 synteny_02/anchors_for_links 中的 *.anchors.simple.filtered，
     保持 6 列基因锚点格式，同时将第 5 列 score 统一转换为整数；
  4) 为每个物种生成 JCVI 所需 BED：
       - 先写入 geneorder 中的所有基因条目（供锚点使用）；
       - 再写入整条染色体的灰色背景条（仅参考物种显示编号，其余物种不显示编号）；
       - 再写入彩色宏观色带（颜色完全由本脚本按 Chr1..ChrN 指定）；
  5) 构建 layout 文件：一行一个物种轨道，最后追加多组 e 边信息，对接 simple 文件；
  6) 调用 `python -m jcvi.graphics.karyotype` 生成最终 PDF。

重要约定：
  - 所有颜色由本脚本内部的 `CHR_COLOR_PALETTE` 控制，不再依赖 02 / 04 的 color 列；
  - 参考物种染色体编号使用纯数字（1..N），其他物种染色体编号一律隐藏（"."）；
  - 左侧物种名到染色体条带的水平距离通过 TRACK_XSTART / TRACK_XEND 手动调节；
  - simple 文件保持 6 列 (geneA, geneB, chrA, chrB, score, orientation)，score 强制为整数。
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# =========================
# 参数自定义区（皇上可在此手动修改）
# =========================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 上游目录
STEP00_RENAME_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"
STEP01_GENEORDER_DIR = PROJECT_ROOT / "output" / "synteny_01_mcscan_catalog" / "geneorder"
STEP02_DIR = PROJECT_ROOT / "output" / "synteny_02_blocks_macro"
PAINT_SEGMENTS_DIR = STEP02_DIR / "paint_segments"
ANCHORS_FOR_LINKS_DIR = STEP02_DIR / "anchors_for_links"
STEP03_DIR = PROJECT_ROOT / "output" / "synteny_03_chr_order"
SEQIDS_SPECIES_DIR = STEP03_DIR / "seqids_species"
STEP04_DIR = PROJECT_ROOT / "output" / "synteny_04_layout"
LAYOUT_TRACKS_FILE = STEP04_DIR / "layout_species_tracks.tsv"

# 输出目录
OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_05_plot"
BED_DIR = OUTPUT_ROOT / "bed_noheader"
SIMPLE_DIR = OUTPUT_ROOT / "simple"
LOG_DIR = OUTPUT_ROOT / "logs"

# JCVI 调用参数
JCVI_PYTHON = "python"
OUTPUT_PDF_NAME = "synteny_linear.pdf"
FIG_WIDTH = 10
FIG_HEIGHT = 7

# 字体控制：True 则尝试使用 FONT_FAMILY；False 则交给 matplotlib 默认
USE_CUSTOM_FONT = True
FONT_FAMILY = "Arial"

# 物种名与染色体条带的水平距离控制
# 数值越小，条带越靠左（离物种名更近）；越大则条带整体往右移动
TRACK_XSTART = 0.080
TRACK_XEND = 0.950

# 染色体颜色调色板（按 Chr1..ChrN 顺序循环）
# 使用皇上的御用配色（顺序可自行调整）
CHR_COLOR_PALETTE = [
    "#E64B35",  # Maroon
    "#4DBBD5",  # SkyBlue
    "#00A087",  # Teal
    "#3C5488",  # Navy
    "#F39B7F",  # Light Orange
    "#8491B4",  # Slate Blue
    "#808180",  # Dark Gray
]

# 日志级别
LOG_LEVEL = "INFO"


# =========================
# 工具函数
# =========================

def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / "synteny05_plot_karyotype.log"

    logger = logging.getLogger("synteny05")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清空旧 handler，避免重复输出
    for h in list(logger.handlers):
        logger.removeHandler(h)

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def clean_output_root(logger: logging.Logger) -> None:
    if OUTPUT_ROOT.exists():
        shutil.rmtree(OUTPUT_ROOT)
    BED_DIR.mkdir(parents=True, exist_ok=True)
    SIMPLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("OUTPUT_ROOT 已清空并重建：%s", OUTPUT_ROOT)


def load_layout_tracks(logger: logging.Logger) -> Tuple[List[Dict[str, str]], str]:
    """
    读取 layout_species_tracks.tsv，得到物种顺序和参考物种 ID。
    """
    if not LAYOUT_TRACKS_FILE.exists():
        logger.error("layout_tracks 文件不存在：%s", LAYOUT_TRACKS_FILE)
        sys.exit(1)

    rows: List[Dict[str, str]] = []
    ref_species: Optional[str] = None

    with LAYOUT_TRACKS_FILE.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(r)
            if r.get("is_reference", "").lower() == "yes":
                ref_species = r["species_id"]

    if not rows:
        logger.error("layout_tracks 文件为空：%s", LAYOUT_TRACKS_FILE)
        sys.exit(1)
    if ref_species is None:
        logger.error("layout_tracks 中未找到参考物种（is_reference==yes）")
        sys.exit(1)

    # 按 order_index 排序
    rows.sort(key=lambda x: int(x.get("order_index", "0")))
    logger.info("共读取物种轨道 %d 条，参考物种 = %s", len(rows), ref_species)
    return rows, ref_species


def get_short_label(species_id: str, track_rows: List[Dict[str, str]]) -> str:
    for r in track_rows:
        if r["species_id"] == species_id:
            sl = r.get("short_label", "").strip()
            if sl:
                return sl
    # 兜底：物种名首字母 + 两个字符
    parts = species_id.split("_")
    if len(parts) >= 2:
        return f"{parts[0][0]}{parts[1][:2]}"
    return species_id[:3]


def build_seqids_file(track_rows: List[Dict[str, str]], logger: logging.Logger) -> None:
    """
    将 synteny_03 中各物种的 *.seqids 合并成一个总的 seqids 文件，
    一行一个物种，对应 layout 的轨道顺序。
    """
    out_path = OUTPUT_ROOT / "seqids"
    with out_path.open("w", encoding="utf-8") as fout:
        for row in track_rows:
            sid = row["species_id"]
            seqids_path = SEQIDS_SPECIES_DIR / f"{sid}.seqids"
            chroms: List[str] = []
            if seqids_path.exists():
                with seqids_path.open("r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith(">"):
                            continue
                        chroms.append(line)
            logger.info("  seqids[%s]：%d 条 Chr。", sid, len(chroms))
            fout.write(",".join(chroms) + "\n")
    logger.info("  seqids 写出：%s", out_path)


def collect_ref_chr_and_colors(
    logger: logging.Logger,
) -> Dict[str, str]:
    """
    从所有 paint_segments_*.tsv 中收集 ref_chr，并按 Chr1..ChrN 的顺序
    为每条参考染色体分配颜色（完全由本脚本决定，不依赖上游颜色）。
    """
    ref_chrs: set[str] = set()

    if not PAINT_SEGMENTS_DIR.exists():
        logger.error("paint_segments 目录不存在：%s", PAINT_SEGMENTS_DIR)
        sys.exit(1)

    for tsv in sorted(PAINT_SEGMENTS_DIR.glob("paint_segments_*.tsv")):
        with tsv.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rc = row.get("ref_chr", "").strip()
                if rc:
                    ref_chrs.add(rc)

    if not ref_chrs:
        logger.error("在 paint_segments 中未找到任何 ref_chr，请检查上游 02 步。")
        sys.exit(1)

    # 尝试按 Chr 前缀和数字排序
    def chr_sort_key(x: str) -> Tuple[int, str]:
        s = x.replace("chr", "").replace("Chr", "")
        try:
            return (0, int(s))
        except ValueError:
            return (1, x)

    sorted_ref = sorted(ref_chrs, key=chr_sort_key)
    logger.info(
        "从 paint_segments 中共收集到 %d 个 ref_chr：%s",
        len(sorted_ref),
        ", ".join(sorted_ref),
    )

    color_map: Dict[str, str] = {}
    for idx, chr_name in enumerate(sorted_ref):
        color = CHR_COLOR_PALETTE[idx % len(CHR_COLOR_PALETTE)]
        color_map[chr_name] = color

    logger.info(
        "参考染色体颜色映射（完全由 05 决定）：%s",
        ", ".join(f"{c}→{h}" for c, h in color_map.items()),
    )
    return color_map


def copy_simple_files(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    logger: logging.Logger,
) -> None:
    """
    从 02 步的 anchors_for_links 中读取 *.anchors.simple.filtered，
    保持 6 列结构，只是把第 5 列 score 强制转换为整数，写入 05/simple 中。
    """
    ref_short = get_short_label(ref_species, track_rows)
    total_kept = 0

    logger.info("Step 1: 复制 simple 文件（保持 6 列 gene 结构）...")

    for row in track_rows:
        sid = row["species_id"]
        if sid == ref_species:
            continue

        qry_short = get_short_label(sid, track_rows)
        src = ANCHORS_FOR_LINKS_DIR / f"{ref_short}__vs__{qry_short}.anchors.simple.filtered"
        dst = SIMPLE_DIR / f"{ref_short}.{qry_short}.simple"

        if not src.exists():
            logger.warning("  simple 源文件缺失，跳过：%s", src)
            continue

        kept = 0
        skipped = 0
        with src.open("r", encoding="utf-8") as fin, dst.open(
            "w", encoding="utf-8"
        ) as fout:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    skipped += 1
                    continue

                gene_a, gene_b, chr_a, chr_b, score_raw, orient = parts[:6]

                # 将 score 强制转换为整数
                try:
                    score_int = int(float(score_raw))
                except ValueError:
                    logger.warning(
                        "  simple 行 score 无法转换为整数（%s），使用 1 代替：%s",
                        score_raw,
                        line,
                    )
                    score_int = 1

                fout.write(
                    "\t".join(
                        [gene_a, gene_b, chr_a, chr_b, str(score_int), orient]
                    )
                    + "\n"
                )
                kept += 1

        total_kept += kept
        logger.info(
            "  simple 写出：%s（保留 %d 条锚点，跳过 %d 条格式不合规行）。",
            dst,
            kept,
            skipped,
        )

    logger.info("simple 文件合计保留锚点 %d 条。", total_kept)


def load_chr_lengths(species_id: str) -> Dict[str, int]:
    """
    从 synteny_00 的 chr_rename_{species}.tsv 中读取染色体长度。
    只保留 is_chromosome == "yes" 的 Chr。
    """
    lengths: Dict[str, int] = {}
    tsv = STEP00_RENAME_DIR / f"chr_rename_{species_id}.tsv"
    if not tsv.exists():
        return lengths

    with tsv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("is_chromosome", "").lower() != "yes":
                continue
            name = row.get("new_chr_name") or row.get("chr_name") or row.get("seqid")
            if not name:
                continue
            try:
                length = int(row["length_bp"])
            except (KeyError, ValueError):
                continue
            lengths[name] = length
    return lengths


def load_geneorder(species_id: str) -> List[Tuple[str, int, int, str]]:
    """
    读取 synteny_01 产生的 geneorder BED：
      chr, start, end, gene_id, ...
    返回列表，用于写入 BED。
    """
    bed_path = STEP01_GENEORDER_DIR / f"{species_id}.bed"
    records: List[Tuple[str, int, int, str]] = []
    if not bed_path.exists():
        return records

    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            gene_id = parts[3]
            records.append((chrom, start, end, gene_id))
    return records


def load_paint_segments(
    species_id: str,
    logger: logging.Logger,
) -> List[Tuple[str, int, int, str]]:
    """
    读取 paint_segments_{species}.tsv：
      只依赖列：chr, start, end, ref_chr
    颜色在 05 中根据 ref_chr → color_map 计算。
    返回列表 (chr, start, end, ref_chr)。
    """
    tsv = PAINT_SEGMENTS_DIR / f"paint_segments_{species_id}.tsv"
    segs: List[Tuple[str, int, int, str]] = []
    if not tsv.exists():
        logger.warning("paint_segments 文件缺失：%s", tsv)
        return segs

    with tsv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("chr", "").strip()
            ref_chr = row.get("ref_chr", "").strip()
            if not chrom or not ref_chr:
                continue
            try:
                start = int(row.get("start", "0"))
                end = int(row.get("end", "0"))
            except ValueError:
                continue
            segs.append((chrom, start, end, ref_chr))

    logger.info("  paint_segments 加载：%s → %d 条宏观色带。", species_id, len(segs))
    return segs


def generate_bed_files(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    ref_chr_color_map: Dict[str, str],
    logger: logging.Logger,
) -> None:
    """
    为每个物种生成 BED：
      1) 所有基因（来自 synteny_01 geneorder）；
      2) 染色体灰色背景条；
      3) 彩色宏观色带（使用 ref_chr → color 映射）。
    """
    logger.info("Step 2: 生成 BED 文件（基因 + 灰底染色体 + 彩色宏观色带）...")

    for row in track_rows:
        sid = row["species_id"]
        out_bed = BED_DIR / f"{sid}.bed"

        gene_records = load_geneorder(sid)
        chr_lengths = load_chr_lengths(sid)
        segs = load_paint_segments(sid, logger)

        is_ref = sid == ref_species

        with out_bed.open("w", encoding="utf-8") as fout:
            # 1) 写入所有基因（供锚点基因查找使用）
            for chrom, start, end, gene_id in gene_records:
                fout.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

            # 2) 整条染色体灰色背景条
            for chrom, length in chr_lengths.items():
                if is_ref:
                    # 参考物种：Chr1 -> 1
                    label = chrom.replace("Chr", "").replace("chr", "")
                else:
                    # 非参考物种：不显示编号
                    label = "."
                fout.write(f"{chrom}\t0\t{length}\t{label}\twhitesmoke\n")

            # 3) 宏观色带，颜色按 ref_chr 映射
            for chrom, start, end, ref_chr in segs:
                color = ref_chr_color_map.get(ref_chr, "#D3D3D3")
                fout.write(f"{chrom}\t{start}\t{end}\t.\t{color}\n")

        logger.info(
            "  写出：%s（基因 = %d，背景 Chr = %d，宏观色带 = %d）。",
            out_bed,
            len(gene_records),
            len(chr_lengths),
            len(segs),
        )


def format_species_label(species_id: str) -> str:
    """
    把 Sinonovacula_constricta 格式化为 S. constricta。
    """
    parts = species_id.split("_")
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return species_id


def build_layout_file(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    logger: logging.Logger,
) -> None:
    """
    生成 layout 文件：
      - 每个物种一行轨道；
      - 最后追加 edges，按 ref ↔ 其它物种连线。
    """
    layout_path = OUTPUT_ROOT / "layout"
    ref_short = get_short_label(ref_species, track_rows)

    # 记录 ref 轨道索引
    ref_index = None
    for i, r in enumerate(track_rows):
        if r["species_id"] == ref_species:
            ref_index = i
            break
    if ref_index is None:
        logger.error("未在 track_rows 中找到参考物种：%s", ref_species)
        sys.exit(1)

    with layout_path.open("w", encoding="utf-8") as f:
        f.write("# y, xstart, xend, rotation, color, label, va, bed\n")

        for idx, row in enumerate(track_rows):
            sid = row["species_id"]
            y_center = float(row.get("y_center", "0.5") or "0.5")
            label = format_species_label(sid)

            # 物种名颜色：默认黑色，参考物种用红色（E64B35）
            if sid == ref_species:
                color = "#E64B35"
            else:
                color = "black"

            # 物种名垂直对齐统一用 top
            f.write(
                f"{y_center:.3f}, {TRACK_XSTART:.3f}, {TRACK_XEND:.3f}, "
                f"0, {color}, {label}, top, bed_noheader/{sid}.bed\n"
            )

        # edges
        f.write("# edges\n")
        for idx, row in enumerate(track_rows):
            sid = row["species_id"]
            if sid == ref_species:
                continue
            qry_short = get_short_label(sid, track_rows)
            simple_path = SIMPLE_DIR / f"{ref_short}.{qry_short}.simple"
            if simple_path.exists():
                f.write(f"e, {ref_index}, {idx}, simple/{ref_short}.{qry_short}.simple\n")

    logger.info("  layout 写出：%s", layout_path)


def call_jcvi(logger: logging.Logger) -> None:
    """
    调用 jcvi.graphics.karyotype 生成 PDF。
    """
    cmd: List[str] = [
        JCVI_PYTHON,
        "-m",
        "jcvi.graphics.karyotype",
        "seqids",
        "layout",
        "--format",
        "pdf",
        f"--figsize={FIG_WIDTH}x{FIG_HEIGHT}",
        "--notex",
    ]
    if USE_CUSTOM_FONT and FONT_FAMILY:
        cmd.extend(["--font", FONT_FAMILY])

    logger.info("Step 4: 调用 JCVI 绘图生成 PDF...")
    logger.info("运行命令：%s", " ".join(cmd))

    proc = subprocess.run(
        cmd,
        cwd=str(OUTPUT_ROOT),
        text=True,
        capture_output=True,
    )
    if proc.stdout:
        logger.info("[STDOUT]\n%s", proc.stdout.strip())
    if proc.stderr:
        logger.info("[STDERR]\n%s", proc.stderr.strip())

    if proc.returncode != 0:
        logger.error("synteny05 — JCVI 绘图失败，请根据日志排查。")
    else:
        logger.info("synteny05 — 运行完成。最终图像：%s", OUTPUT_ROOT / OUTPUT_PDF_NAME)


# =========================
# Main
# =========================

def main() -> None:
    logger = setup_logging()
    logger.info("========== synteny05 — JCVI karyotype 最终绘图 ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info(
        "USE_CUSTOM_FONT = %s (font=%s)", "True" if USE_CUSTOM_FONT else "False", FONT_FAMILY
    )
    logger.info(
        "TRACK_XSTART = %.3f, TRACK_XEND = %.3f", TRACK_XSTART, TRACK_XEND
    )

    clean_output_root(logger)

    # 1) 读取物种轨道信息
    track_rows, ref_species = load_layout_tracks(logger)

    # 2) 生成 seqids
    logger.info("Step 0: 生成全局 seqids 文件...")
    build_seqids_file(track_rows, logger)

    # 3) 构建参考染色体颜色表
    ref_chr_color_map = collect_ref_chr_and_colors(logger)

    # 4) 复制 simple 文件（score 转为整数，保持 6 列）
    copy_simple_files(track_rows, ref_species, logger)

    # 5) 生成每个物种的 BED（基因 + 灰底 Chr + 彩色宏观色带）
    generate_bed_files(track_rows, ref_species, ref_chr_color_map, logger)

    # 6) 构建 layout
    logger.info("Step 3: 生成 layout 文件...")
    build_layout_file(track_rows, ref_species, logger)

    # 7) 调用 JCVI 绘图
    call_jcvi(logger)


if __name__ == "__main__":
    main()

