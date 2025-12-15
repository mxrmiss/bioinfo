#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny05_plot_karyotype.py
—— JCVI karyotype 多物种宏观共线性最终绘图脚本（05 终版 · 颜色由 02 的规则复用 + 05 负责落地出图）

目标（皇上钦定）：
  - ribbons 继续彩色（复用 02 的 ref_chr → 颜色规则，且 simple 中的 #HEX*gene 颜色前缀必须被保留）；
  - 去掉“空壳两端”：采用方案 A —— 不再绘制“整条染色体背景条”，只保留中间的实心色带；
  - 圆圈与圆圈数字：要么全有，要么全无（脚本内只保留一个 circles 开关，不再拆分成多个开关）。

职责概述：
  1) 读取 layout_species_tracks.tsv，确定物种顺序、参考物种、轨道位置；
  2) 读取 synteny_03 产生的各物种 seqids，合并生成全局 seqids 文件；
  3) 读取 synteny_02/anchors_for_links 中的 *.anchors.simple.filtered，
     生成 JCVI 需要的 *.simple：
       - 不把 “#HEX*gene” 当注释；
       - 仅保留 len(parts) >= 6 且第 5 列 score 可转数字的行；
       - 将 score 统一写为整数；
  4) 为每个物种生成 JCVI 所需 BED：
       - geneorder 中所有基因条目（供锚点定位）；
       - 彩色宏观色带（paint_segments，颜色复用 02 的 ref_chr 映射）；
     重要稳健性（颜色注释问题）：
       - 若 simple 的参考基因带 '#HEX*' 前缀，参考物种 BED 的 gene 行最稳方式：
         “同时写 plain 和 #HEX* 两份”，避免上游混前缀导致找不到基因。
  5) 构建 layout 文件并追加 edges；
  6) 调用 `python -m jcvi.graphics.karyotype` 生成最终 PDF：
       - circles：只用一个开关控制（要么都有，要么都没有）。
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
FIG_WIDTH = 16
FIG_HEIGHT = 7

# 字体控制：True 则尝试使用 FONT_FAMILY；False 则交给 matplotlib 默认
USE_CUSTOM_FONT = True
FONT_FAMILY = "Arial"

# circles（圆圈+圆圈数字）：True=显示；False=全部去掉
SHOW_CIRCLES = False

# ---------- 物种名左侧空间 & 名称到染色体的水平距离（皇上可调） ----------
TRACK_XSTART_BASE = 0.120
NAME_TO_TRACK_GAP = 0.001
TRACK_XEND = 0.950
TRACK_XSTART = TRACK_XSTART_BASE + NAME_TO_TRACK_GAP

# 染色体颜色调色板（按 Chr1..ChrN 顺序循环）
CHR_COLOR_PALETTE = [
    "#F79D93",  # A - 干燥玫瑰 (Dusty Rose) - 比原红更柔和
    "#F6CD96",  # B - 杏仁奶油 (Apricot) - 去掉了刺眼的橙光
    "#9FD5CB",  # D - 灰豆绿 (Sage Green) - 很有质感的植物绿
    "#95C8F2",  # F - 雾霾蓝 (Hazy Blue) - 非常经典的科研蓝
    "#C8C0F0",  # G - 丁香紫 (Lilac) - 淡雅的紫色
    "#F6C6E7"   # H - 柔粉色 (Soft Pink) - 少女感但稳重
]

DEFAULT_GRAY = "#D3D3D3"
REF_LABEL_COLOR = "#E64B35"

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


def _norm_key(k: str) -> str:
    return str(k).strip().strip("\ufeff").replace("\r", "")


def _norm_val(v: str) -> str:
    return str(v).strip().replace("\r", "")


def _normalize_row(row: Dict[str, str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for k, v in row.items():
        out[_norm_key(k)] = _norm_val(v)
    return out


def load_layout_tracks(logger: logging.Logger) -> Tuple[List[Dict[str, str]], str]:
    """
    读取 layout_species_tracks.tsv，得到物种顺序和参考物种 ID。
    兼容 Windows CRLF（字段名/字段值可能带 \\r）。
    """
    if not LAYOUT_TRACKS_FILE.exists():
        logger.error("layout_tracks 文件不存在：%s", LAYOUT_TRACKS_FILE)
        sys.exit(1)

    rows: List[Dict[str, str]] = []
    ref_species: Optional[str] = None

    with LAYOUT_TRACKS_FILE.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rr = _normalize_row(r)
            rows.append(rr)
            if (rr.get("is_reference", "") or "").lower() == "yes":
                ref_species = rr.get("species_id", "").strip()

    if not rows:
        logger.error("layout_tracks 文件为空：%s", LAYOUT_TRACKS_FILE)
        sys.exit(1)
    if ref_species is None or not ref_species.strip():
        logger.error("layout_tracks 中未找到参考物种（is_reference==yes）")
        sys.exit(1)

    def _safe_int(x: str) -> int:
        try:
            return int(str(x).strip())
        except ValueError:
            return 0

    rows.sort(key=lambda x: _safe_int(x.get("order_index", "0")))
    logger.info("共读取物种轨道 %d 条，参考物种 = %s", len(rows), ref_species)
    return rows, ref_species


def get_short_label(species_id: str, track_rows: List[Dict[str, str]]) -> str:
    for r in track_rows:
        if (r.get("species_id") or "").strip() == species_id:
            sl = (r.get("short_label") or "").strip()
            if sl:
                return sl
    parts = species_id.split("_")
    if len(parts) >= 2:
        return f"{parts[0][0]}{parts[1][:2]}"
    return species_id[:3]


def build_seqids_file(track_rows: List[Dict[str, str]], logger: logging.Logger) -> None:
    out_path = OUTPUT_ROOT / "seqids"
    with out_path.open("w", encoding="utf-8") as fout:
        for row in track_rows:
            sid = (row.get("species_id") or "").strip()
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


def collect_ref_chr_and_colors(logger: logging.Logger) -> Dict[str, str]:
    ref_chrs: set[str] = set()

    if not PAINT_SEGMENTS_DIR.exists():
        logger.error("paint_segments 目录不存在：%s", PAINT_SEGMENTS_DIR)
        sys.exit(1)

    for tsv in sorted(PAINT_SEGMENTS_DIR.glob("paint_segments_*.tsv")):
        with tsv.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                row = _normalize_row(row)
                rc = (row.get("ref_chr") or "").strip()
                if rc:
                    ref_chrs.add(rc)

    if not ref_chrs:
        logger.error("在 paint_segments 中未找到任何 ref_chr，请检查上游 02 步。")
        sys.exit(1)

    def chr_sort_key(x: str) -> Tuple[int, int, str]:
        s = x.replace("chr", "").replace("Chr", "")
        try:
            return (0, int(s), x)
        except ValueError:
            return (1, 10**9, x)

    sorted_ref = sorted(ref_chrs, key=chr_sort_key)
    logger.info(
        "从 paint_segments 中共收集到 %d 个 ref_chr：%s",
        len(sorted_ref),
        ", ".join(sorted_ref),
    )

    color_map: Dict[str, str] = {}
    for idx, chr_name in enumerate(sorted_ref):
        color_map[chr_name] = CHR_COLOR_PALETTE[idx % len(CHR_COLOR_PALETTE)]

    logger.info(
        "参考染色体颜色映射（复用 synteny02 规则）：%s",
        ", ".join(f"{c}#{h.lstrip('#')}" for c, h in color_map.items()),
    )
    return color_map


def _first_valid_record_from_file(path: Path) -> Optional[List[str]]:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 6:
                return parts[:6]
    return None


def detect_ref_gene_prefix_needed(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    logger: logging.Logger,
) -> bool:
    ref_short = get_short_label(ref_species, track_rows)

    for row in track_rows:
        sid = (row.get("species_id") or "").strip()
        if sid == ref_species:
            continue
        qry_short = get_short_label(sid, track_rows)
        src = ANCHORS_FOR_LINKS_DIR / f"{ref_short}__vs__{qry_short}.anchors.simple.filtered"
        parts = _first_valid_record_from_file(src)
        if not parts:
            continue
        gene_a = parts[0]
        if gene_a.startswith("#") and "*" in gene_a and len(gene_a.split("*", 1)[0]) == 7:
            logger.info("检测到 simple 参考基因带颜色前缀（#HEX*）：True（示例：%s）", gene_a)
            return True
        logger.info("检测到 simple 参考基因带颜色前缀（#HEX*）：False（示例：%s）", gene_a)
        return False

    logger.warning("未能从任何 filtered 文件中抽到有效记录，默认不使用 #HEX* 前缀。")
    return False


def copy_simple_files(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    logger: logging.Logger,
) -> None:
    """
    关键：不把行首 '#' 当注释（因为 #HEX*gene 是数据）。
    仅收：len(parts)>=6 且第 5 列 score 可转数字。
    """
    ref_short = get_short_label(ref_species, track_rows)
    total_kept = 0
    total_skipped = 0

    logger.info("Step 1: 复制 simple 文件（列数>=6 且 score 可转数字才收）...")

    for row in track_rows:
        sid = (row.get("species_id") or "").strip()
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

        with src.open("r", encoding="utf-8") as fin, dst.open("w", encoding="utf-8") as fout:
            for line in fin:
                line = line.strip()
                if not line:
                    continue

                parts = line.split()
                if len(parts) < 6:
                    skipped += 1
                    continue

                gene_a, gene_b, chr_a, chr_b, score_raw, orient = parts[:6]

                try:
                    score_float = float(score_raw)
                except ValueError:
                    skipped += 1
                    continue

                score_int = int(score_float)
                fout.write("\t".join([gene_a, gene_b, chr_a, chr_b, str(score_int), orient]) + "\n")
                kept += 1

        total_kept += kept
        total_skipped += skipped
        logger.info("  simple 写出：%s（保留 %d，跳过 %d）。", dst, kept, skipped)

    logger.info("simple 文件合计保留锚点 %d 条，跳过 %d 条。", total_kept, total_skipped)


def load_geneorder(species_id: str) -> List[Tuple[str, int, int, str]]:
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
            chrom = parts[0].strip()
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            gene_id = parts[3].strip().replace("\r", "")
            records.append((chrom, start, end, gene_id))
    return records


def load_paint_segments(species_id: str, logger: logging.Logger) -> List[Tuple[str, int, int, str]]:
    tsv = PAINT_SEGMENTS_DIR / f"paint_segments_{species_id}.tsv"
    segs: List[Tuple[str, int, int, str]] = []
    if not tsv.exists():
        logger.warning("paint_segments 文件缺失：%s", tsv)
        return segs

    with tsv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = _normalize_row(row)
            chrom = (row.get("chr") or "").strip()
            ref_chr = (row.get("ref_chr") or "").strip()
            if not chrom or not ref_chr:
                continue
            try:
                start = int(str(row.get("start", "0")).strip())
                end = int(str(row.get("end", "0")).strip())
            except ValueError:
                continue
            segs.append((chrom, start, end, ref_chr))

    logger.info("  paint_segments 加载：%s → %d 条宏观色带。", species_id, len(segs))
    return segs


def generate_bed_files(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    ref_chr_color_map: Dict[str, str],
    ref_gene_need_prefix: bool,
    logger: logging.Logger,
) -> None:
    """
    方案 A（皇上选择）：不再绘制“整条染色体背景条”，只保留中间实心色带（paint segments）。

    颜色注释稳健性：
      - 参考物种 gene 行：如果需要前缀，则“同时写两份”（plain + #HEX*），避免 simple 混前缀导致找不到。
    """
    logger.info("Step 2: 生成 BED 文件（基因 + 彩色宏观色带）...")
    logger.info("参考物种基因是否需要 '#HEX*' 前缀：%s", "True" if ref_gene_need_prefix else "False")

    for row in track_rows:
        sid = (row.get("species_id") or "").strip()
        out_bed = BED_DIR / f"{sid}.bed"

        gene_records = load_geneorder(sid)
        segs = load_paint_segments(sid, logger)

        is_ref = sid == ref_species

        n_prefixed = 0
        n_missing_color = 0
        n_plain = 0

        with out_bed.open("w", encoding="utf-8") as fout:
            # 1) gene 行（供锚点定位使用）
            for chrom, start, end, gene_id in gene_records:
                if is_ref and ref_gene_need_prefix:
                    # 先写 plain（兜底）
                    fout.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")
                    n_plain += 1

                    # 再写 prefixed（给 simple 中 #HEX*gene 找到用）
                    color = ref_chr_color_map.get(chrom)
                    if not color:
                        color = DEFAULT_GRAY
                        n_missing_color += 1
                    fout.write(f"{chrom}\t{start}\t{end}\t{color}*{gene_id}\n")
                    n_prefixed += 1
                else:
                    fout.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

            chrom_max_end: Dict[str, int] = {}
            for chrom, start, end, ref_chr in segs:
                prev = chrom_max_end.get(chrom, 0)
                if end > prev:
                    chrom_max_end[chrom] = end
            if not chrom_max_end:
                for chrom, start, end, gene_id in gene_records:
                    prev = chrom_max_end.get(chrom, 0)
                    if end > prev:
                        chrom_max_end[chrom] = end
            for chrom, max_end in chrom_max_end.items():
                if max_end > 0:
                    fout.write(f"{chrom}\t0\t{max_end}\t.\t#D9D9D9\n")

            # 2) 宏观色带（paint segments）：只画实心的中间部分
            for chrom, start, end, ref_chr in segs:
                color = ref_chr_color_map.get(ref_chr, DEFAULT_GRAY)
                fout.write(f"{chrom}\t{start}\t{end}\t.\t{color}\n")

        if is_ref and ref_gene_need_prefix:
            logger.info(
                "  写出：%s（基因原始=%d；plain=%d；prefixed=%d；宏观色带=%d；找不到Chr颜色=%d）。",
                out_bed, len(gene_records), n_plain, n_prefixed, len(segs), n_missing_color
            )
        else:
            logger.info(
                "  写出：%s（基因=%d；宏观色带=%d）。",
                out_bed, len(gene_records), len(segs)
            )


def format_species_label(species_id: str) -> str:
    parts = species_id.split("_")
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return species_id


def build_layout_file(
    track_rows: List[Dict[str, str]],
    ref_species: str,
    logger: logging.Logger,
) -> None:
    layout_path = OUTPUT_ROOT / "layout"
    ref_short = get_short_label(ref_species, track_rows)

    ref_index = None
    for i, r in enumerate(track_rows):
        if (r.get("species_id") or "").strip() == ref_species:
            ref_index = i
            break
    if ref_index is None:
        logger.error("未在 track_rows 中找到参考物种：%s", ref_species)
        sys.exit(1)

    with layout_path.open("w", encoding="utf-8") as f:
        f.write("# y, xstart, xend, rotation, color, label, va, bed\n")

        for idx, row in enumerate(track_rows):
            sid = (row.get("species_id") or "").strip()
            y_center = float((row.get("y_center") or "0.5").strip() or "0.5")
            label = format_species_label(sid)

            if sid == ref_species:
                color = REF_LABEL_COLOR
            else:
                color = "black"

            f.write(
                f"{y_center:.3f}, {TRACK_XSTART:.3f}, {TRACK_XEND:.3f}, "
                f"0, {color}, {label}, top, bed_noheader/{sid}.bed\n"
            )

        f.write("# edges\n")
        for idx, row in enumerate(track_rows):
            sid = (row.get("species_id") or "").strip()
            if sid == ref_species:
                continue
            qry_short = get_short_label(sid, track_rows)
            simple_path = SIMPLE_DIR / f"{ref_short}.{qry_short}.simple"
            if simple_path.exists():
                f.write(f"e, {ref_index}, {idx}, simple/{ref_short}.{qry_short}.simple\n")

    logger.info("  layout 写出：%s", layout_path)


def call_jcvi(logger: logging.Logger) -> None:
    """
    circles（圆圈+数字）只用一个开关控制：
      - SHOW_CIRCLES=False：固定加 --nocircles（全无）
      - SHOW_CIRCLES=True ：不加 --nocircles（全有）
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
#        "--chrstyle",
#        "rect",
    ]

    if not SHOW_CIRCLES:
        cmd.append("--nocircles")

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
        logger.info("synteny05 — 运行完成。请在输出目录查看生成的 PDF：%s", OUTPUT_ROOT)


# =========================
# Main
# =========================

def main() -> None:
    logger = setup_logging()
    logger.info("========== synteny05 — JCVI karyotype 最终绘图 ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info("USE_CUSTOM_FONT = %s (font=%s)", "True" if USE_CUSTOM_FONT else "False", FONT_FAMILY)
    logger.info("TRACK_XSTART_BASE = %.3f", TRACK_XSTART_BASE)
    logger.info("NAME_TO_TRACK_GAP = %.3f", NAME_TO_TRACK_GAP)
    logger.info("TRACK_XSTART = %.3f, TRACK_XEND = %.3f", TRACK_XSTART, TRACK_XEND)
    logger.info("SHOW_CIRCLES = %s", "True" if SHOW_CIRCLES else "False")

    clean_output_root(logger)

    # 1) 读取物种轨道信息
    track_rows, ref_species = load_layout_tracks(logger)

    # 2) 生成 seqids
    logger.info("Step 0: 生成全局 seqids 文件...")
    build_seqids_file(track_rows, logger)

    # 3) 构建参考染色体颜色表（复用 02 规则）
    ref_chr_color_map = collect_ref_chr_and_colors(logger)

    # 4) 检测参考基因是否使用 #HEX* 前缀
    ref_gene_need_prefix = detect_ref_gene_prefix_needed(track_rows, ref_species, logger)

    # 5) 复制 simple 文件（列数>=6 且 score 可转数字才收）
    copy_simple_files(track_rows, ref_species, logger)

    # 6) 生成每个物种的 BED（基因 + 彩色宏观色带）
    generate_bed_files(track_rows, ref_species, ref_chr_color_map, ref_gene_need_prefix, logger)

    # 7) 构建 layout
    logger.info("Step 3: 生成 layout 文件...")
    build_layout_file(track_rows, ref_species, logger)

    # 8) 调用 JCVI 绘图
    call_jcvi(logger)


if __name__ == "__main__":
    main()

