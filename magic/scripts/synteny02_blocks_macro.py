#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny02_blocks_macro.py
—— JCVI MCscan 结果整理：blocks 标准化 + 染色体 segments + anchors.simple 过滤

当前职责（对应蓝图的 Step 02）：
  1) 02a：统一调用 jcvi 为“参考物种 vs 其它物种”生成 anchors.simple：
       - 在 raw_jcvi/ 下为每个物种准备 <short>.bed + <short>.pep 软链接；
       - 调用：
           python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names ref_short qry_short
           python -m jcvi.compara.synteny screen --minspan=MINSPAN --simple ref_short.qry_short.anchors ref_short.qry_short.anchors.new
       - 得到：raw_jcvi/<ref_short>.<qry_short>.anchors.simple
  2) 02b：基于 anchors.simple + geneorder 生成 blocks_normalized：
       - blocks_normalized/<ref_short>__vs__<qry_short>.blocks.tsv
       - pairwise_block_summary.tsv
       - paint_segments/paint_segments_<species_id>.tsv
  3) 02c：从 anchors.simple 中筛选“画细线”的简化锚点：
       - anchors_for_links/<ref_short>__vs__<qry_short>.anchors.simple.filtered

说明：
  - 本脚本使用的共线信息全部来自 jcvi 输出的 anchors.simple，
    不依赖 synteny_01 的 BLAST 文件（01 的 BLAST 可用于后续其它分析）。
  - paint_segments 仅负责各物种染色体上“对应参考 Chr 的区段”坐标信息，
    不再包含颜色信息；颜色由后续出图脚本统一分配。
"""

from __future__ import annotations

import sys
import csv
import shutil
import logging
import subprocess
import os
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
PROTEIN_DIR = RAW_DATA_DIR / "protein"

CHR_RENAME_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"
STEP01_DIR = PROJECT_ROOT / "output" / "synteny_01_mcscan_catalog"
STEP01_GENEORDER_DIR = STEP01_DIR / "geneorder"

# 本脚本输出目录
OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_02_blocks_macro"
RAW_JCVI_DIR = OUTPUT_ROOT / "raw_jcvi"
BLOCKS_NORMALIZED_DIR = OUTPUT_ROOT / "blocks_normalized"
PAINT_SEGMENTS_DIR = OUTPUT_ROOT / "paint_segments"
ANCHORS_FOR_LINKS_DIR = OUTPUT_ROOT / "anchors_for_links"

# 是否在脚本内调用 jcvi 生成 anchors.simple
# 当前设为 True：00 → 01 → 02 一条龙，无需额外手动步骤。
RUN_JCVI_PIPELINE = True

# jcvi 相关参数
JCVI_PYTHON = "python"   # jcvi 所在的 python 命令
JCVI_CSCORE = 0.7        # compara.catalog ortholog --cscore
JCVI_MINSPAN = 30        # synteny.screen --minspan

# “画细线”时每条参考染色体最多保留多少个共线区块
MAX_LINK_BLOCKS_PER_REF_CHR = 300

# paint_segments 中保留的最小区段长度（bp）
MIN_SEGMENT_LENGTH_BP = 100_000

# 日志等级
LOG_LEVEL = "INFO"

# jcvi dotplot 字体设置：
#   - DOTPLOT_SANS_SERIF_FONTS：写入 matplotlibrc 的 sans-serif 字体列表（按顺序偏好）。
#   - ENABLE_DOTPLOT_FONT_TWEAK：
#       True  时：为 jcvi 子进程设置 MPLCONFIGDIR，使 dotplot 使用这些字体。
#       False 时：完全不干预，使用 jcvi / matplotlib 默认字体。
DOTPLOT_SANS_SERIF_FONTS: List[str] = ["Arial", "DejaVu Sans"]
ENABLE_DOTPLOT_FONT_TWEAK = True

# 供 run_cmd 使用的 MPLCONFIGDIR（若未成功创建则为 None）
MPLCONFIG_DIR_FOR_JCVI: Optional[Path] = None

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
    log_file = log_dir / "synteny02_blocks_macro.log"

    logger = logging.getLogger("synteny02")
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

    logger.info("========== synteny02 — JCVI blocks 整理 ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info("RUN_JCVI_PIPELINE                = %s", RUN_JCVI_PIPELINE)
    logger.info("MAX_LINK_BLOCKS_PER_REF_CHR      = %d", MAX_LINK_BLOCKS_PER_REF_CHR)
    logger.info("MIN_SEGMENT_LENGTH_BP            = %d", MIN_SEGMENT_LENGTH_BP)
    logger.info("JCVI_MINSPAN (screen --minspan)  = %d", JCVI_MINSPAN)
    logger.info("JCVI_CSCORE  (ortholog --cscore) = %.2f", JCVI_CSCORE)
    logger.info("DOTPLOT_SANS_SERIF_FONTS         = %s", ", ".join(DOTPLOT_SANS_SERIF_FONTS))
    logger.info("ENABLE_DOTPLOT_FONT_TWEAK        = %s", ENABLE_DOTPLOT_FONT_TWEAK)

    return logger


def clean_output_root(output_root: Path) -> None:
    """
    删除旧的输出目录并重建子目录。
    """
    if output_root.exists():
        shutil.rmtree(output_root)
    RAW_JCVI_DIR.mkdir(parents=True, exist_ok=True)
    BLOCKS_NORMALIZED_DIR.mkdir(parents=True, exist_ok=True)
    PAINT_SEGMENTS_DIR.mkdir(parents=True, exist_ok=True)
    ANCHORS_FOR_LINKS_DIR.mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)


def load_species_meta(meta_file: Path, logger: logging.Logger) -> Tuple[List[str], str]:
    """
    读取 synteny_species_meta.tsv，返回 species_id 列表和唯一参考物种。
    """
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

        # 基础前缀：首字母 + 种名前两位（若不足则尽量多）
        if species:
            base = genus[0] + species[:2]
        else:
            base = (genus[:3] or sid[:3])

        cand = base.capitalize()
        # 若已被占用，则逐步延长种名部分，直至唯一
        extra = 2
        while cand in used and used[cand] != sid:
            extra += 1
            if species and extra < len(species):
                cand = (genus[0] + species[:extra]).capitalize()
            else:
                # 退而求其次，拼接 genus 更多字母或序号
                cand = (genus[:min(len(genus), extra)] + species[:1]).capitalize()
                if cand in used and used[cand] != sid:
                    cand = (base + str(extra)).capitalize()

        mapping[sid] = cand
        used[cand] = sid

    logger.info("物种短名映射：%s",
                ", ".join(f"{k} -> {v}" for k, v in mapping.items()))
    return mapping


def get_short_name(species_id: str) -> str:
    """
    获取物种短名；优先使用自动生成的 SHORT_NAME_MAP。
    """
    if species_id in SHORT_NAME_MAP:
        return SHORT_NAME_MAP[species_id]

    parts = species_id.split("_")
    if len(parts) >= 2:
        g, s = parts[0], parts[1]
        return (g[0] + s[:2]).capitalize()
    return species_id[:4]


def prepare_mpl_config_for_jcvi(logger: logging.Logger) -> Optional[Path]:
    """
    如有需要，为 jcvi 子进程创建专用 MPLCONFIGDIR：
      - 在 matplotlibrc 中设置：
            font.family    : sans-serif
            font.sans-serif: <DOTPLOT_SANS_SERIF_FONTS 列表>
      - 通过环境变量 MPLCONFIGDIR 传递给 jcvi 使用。

    不再检测字体是否真实存在，若缺失则由 matplotlib 自己回退并可能打印警告。
    """
    if not ENABLE_DOTPLOT_FONT_TWEAK:
        logger.info("不调整 jcvi dotplot 字体（ENABLE_DOTPLOT_FONT_TWEAK = False）。")
        return None

    fonts = [f.strip() for f in DOTPLOT_SANS_SERIF_FONTS if f.strip()]
    if not fonts:
        logger.info("DOTPLOT_SANS_SERIF_FONTS 为空，跳过 dotplot 字体设置。")
        return None

    fonts_str = ", ".join(fonts)

    mpl_dir = OUTPUT_ROOT / "mplconfig_for_jcvi"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    rc_file = mpl_dir / "matplotlibrc"
    with rc_file.open("w", encoding="utf-8") as f:
        f.write("font.family: sans-serif\n")
        f.write(f"font.sans-serif: {fonts_str}\n")

    logger.info(
        "已为 jcvi dotplot 创建 MPLCONFIGDIR：%s（font.sans-serif = [%s]）",
        mpl_dir,
        fonts_str,
    )
    return mpl_dir


def run_cmd(cmd: List[str], logger: logging.Logger, cwd: Optional[Path] = None) -> int:
    """
    统一封装外部命令调用。
    如 MPLCONFIG_DIR_FOR_JCVI 非空，则为子进程设置 MPLCONFIGDIR（用于控制 dotplot 字体）。
    """
    cmd_str = " ".join(str(x) for x in cmd)
    logger.info("运行命令：%s", cmd_str)

    env = None
    if MPLCONFIG_DIR_FOR_JCVI is not None:
        env = dict(os.environ)
        env["MPLCONFIGDIR"] = str(MPLCONFIG_DIR_FOR_JCVI)

    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError as e:
        logger.error("命令不存在：%s", e)
        return 1

    if result.stdout:
        logger.info("[STDOUT]\n%s\n", result.stdout.strip())
    if result.stderr:
        logger.info("[STDERR]\n%s\n", result.stderr.strip())

    if result.returncode != 0:
        logger.error("命令执行失败，退出码：%d", result.returncode)
    else:
        logger.info("命令执行成功。")
    return result.returncode


# =========================
# 02a — jcvi 调用生成 anchors.simple
# =========================

def prepare_jcvi_inputs_for_species(
    species_id: str,
    short_name: str,
    logger: logging.Logger,
) -> None:
    """
    为 jcvi 准备 <short>.bed + <short>.pep 软链接：
      - bed 直接使用 synteny_01 的 geneorder bed；
      - pep 使用 raw_data/protein/<species_id>.faa。
    """
    # bed
    src_bed = STEP01_GENEORDER_DIR / f"{species_id}.bed"
    if not src_bed.exists():
        logger.error("找不到 geneorder bed：%s", src_bed)
        sys.exit(1)

    dst_bed = RAW_JCVI_DIR / f"{short_name}.bed"
    if dst_bed.exists():
        dst_bed.unlink()
    dst_bed.symlink_to(src_bed)

    # pep
    src_faa = PROTEIN_DIR / f"{species_id}.faa"
    if not src_faa.exists():
        logger.error("找不到蛋白文件：%s", src_faa)
        sys.exit(1)

    dst_pep = RAW_JCVI_DIR / f"{short_name}.pep"
    if dst_pep.exists():
        dst_pep.unlink()
    dst_pep.symlink_to(src_faa)

    logger.info("为物种 %s 准备 jcvi 输入：%s, %s", species_id, dst_bed, dst_pep)


def run_jcvi_for_pair(
    ref_species: str,
    qry_species: str,
    logger: logging.Logger,
) -> None:
    """
    为 ref vs qry 生成 anchors.simple：
      1) 准备 <ref_short>.bed/.pep 和 <qry_short>.bed/.pep；
      2) jcvi.compara.catalog ortholog（生成 ref.qry.anchors 等）；
      3) jcvi.compara.synteny screen --simple（生成 ref.qry.anchors.simple）。
    """
    ref_short = get_short_name(ref_species)
    qry_short = get_short_name(qry_species)

    anchors_simple = RAW_JCVI_DIR / f"{ref_short}.{qry_short}.anchors.simple"
    if anchors_simple.exists():
        logger.info("检测到已有 anchors.simple：%s，本次不再重复运行 jcvi。", anchors_simple)
        return

    logger.info("====== JCVI pipeline：%s (%s) vs %s (%s) ======",
                ref_species, ref_short, qry_species, qry_short)

    prepare_jcvi_inputs_for_species(ref_species, ref_short, logger)
    prepare_jcvi_inputs_for_species(qry_species, qry_short, logger)

    # 1) compara.catalog ortholog
    cmd_catalog = [
        JCVI_PYTHON, "-m", "jcvi.compara.catalog", "ortholog",
        "--dbtype", "prot",
        "--no_strip_names",
        f"--cscore={JCVI_CSCORE}",
        ref_short,
        qry_short,
    ]
    rc = run_cmd(cmd_catalog, logger, cwd=RAW_JCVI_DIR)
    if rc != 0:
        logger.error("jcvi.compara.catalog ortholog 运行失败（%s vs %s），终止该物种对。",
                     ref_short, qry_short)
        return

    anchors_prefix = RAW_JCVI_DIR / f"{ref_short}.{qry_short}.anchors"
    anchors_new = RAW_JCVI_DIR / f"{ref_short}.{qry_short}.anchors.new"

    # 2) synteny.screen --simple
    cmd_screen = [
        JCVI_PYTHON, "-m", "jcvi.compara.synteny", "screen",
        f"--minspan={JCVI_MINSPAN}",
        "--simple",
        str(anchors_prefix),
        str(anchors_new),
    ]
    rc = run_cmd(cmd_screen, logger, cwd=RAW_JCVI_DIR)
    if rc != 0:
        logger.error("jcvi.compara.synteny screen 运行失败（%s vs %s）。", ref_short, qry_short)
        return

    if anchors_simple.exists():
        logger.info("anchors.simple 生成成功：%s", anchors_simple)
    else:
        logger.warning("jcvi 运行结束但未找到 anchors.simple：%s", anchors_simple)


# =========================
# 读取 geneorder / anchors.simple
# =========================

def load_geneorder_bed(
    species_id: str,
    logger: logging.Logger,
) -> Tuple[Dict[str, Tuple[str, int, int]], Dict[str, List[Tuple[int, int, str]]]]:
    """
    读取 synteny_01_mcscan_catalog/geneorder/<species_id>.bed

    返回：
      gene_info: gene_id -> (chr, start, end)
      chr_intervals: chr -> [(start, end, gene_id), ...]，按 start 排序
    """
    bed_path = STEP01_GENEORDER_DIR / f"{species_id}.bed"
    if not bed_path.exists():
        logger.error("找不到 geneorder bed：%s", bed_path)
        sys.exit(1)

    gene_info: Dict[str, Tuple[str, int, int]] = {}
    chr_intervals: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)

    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chr_name, start_str, end_str, gene_id = parts[0], parts[1], parts[2], parts[3]
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue
            gene_info[gene_id] = (chr_name, start, end)
            chr_intervals[chr_name].append((start, end, gene_id))

    for chr_name in chr_intervals:
        chr_intervals[chr_name].sort(key=lambda x: x[0])

    logger.info("geneorder(%s): 共解析到 %d 个基因，%d 条染色体。",
                species_id, len(gene_info), len(chr_intervals))
    return gene_info, chr_intervals


def parse_anchors_simple(
    anchors_simple_path: Path,
    logger: logging.Logger,
) -> List[Tuple[str, str, str, str, float, str]]:
    """
    解析 jcvi 的 *.anchors.simple 文件。

    约定每行 6 列：
      ref_start_gene  ref_end_gene  qry_start_gene  qry_end_gene  score  orientation
    """
    records: List[Tuple[str, str, str, str, float, str]] = []

    if not anchors_simple_path.exists():
        logger.error("anchors.simple 文件不存在：%s", anchors_simple_path)
        return records

    with anchors_simple_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            ref_start_gene, ref_end_gene, qry_start_gene, qry_end_gene = parts[0], parts[1], parts[2], parts[3]
            try:
                score = float(parts[4])
            except ValueError:
                score = 0.0
            orientation = parts[5]
            records.append((ref_start_gene, ref_end_gene, qry_start_gene, qry_end_gene, score, orientation))

    logger.info("读取 anchors.simple：%s，共 %d 个共线区块。", anchors_simple_path.name, len(records))
    return records


# =========================
# 02b — blocks_normalized + paint_segments
# =========================

def build_blocks_normalized_and_segments_for_pair(
    ref_species: str,
    qry_species: str,
    ref_gene_info: Dict[str, Tuple[str, int, int]],
    qry_gene_info: Dict[str, Tuple[str, int, int]],
    anchors_records: List[Tuple[str, str, str, str, float, str]],
    blocks_path: Path,
) -> Tuple[Dict[str, List[Tuple[int, int, str, int]]], Dict[str, object]]:
    """
    把某一物种对的 anchors.simple 转换成：
      - blocks_normalized/<ref>__vs__<qry>.blocks.tsv
      - qry 物种上的原始 segments_by_chr（后续再合并）
      - 该物种对的 summary 记录
    """
    segments_by_chr: Dict[str, List[Tuple[int, int, str, int]]] = defaultdict(list)

    blocks_fieldnames = [
        "ref_species",
        "qry_species",
        "ref_chr",
        "ref_start",
        "ref_end",
        "qry_chr",
        "qry_start",
        "qry_end",
        "orientation",
        "n_genes",
        "score",
    ]
    blocks_path.parent.mkdir(parents=True, exist_ok=True)
    f_blocks = blocks_path.open("w", encoding="utf-8", newline="")
    writer = csv.DictWriter(f_blocks, fieldnames=blocks_fieldnames, delimiter="\t")
    writer.writeheader()

    n_blocks_raw = len(anchors_records)
    n_blocks_kept = 0
    span_ref_list: List[int] = []
    span_qry_list: List[int] = []

    for ref_start_gene, ref_end_gene, qry_start_gene, qry_end_gene, score, orientation in anchors_records:
        if ref_start_gene not in ref_gene_info or ref_end_gene not in ref_gene_info:
            continue
        if qry_start_gene not in qry_gene_info or qry_end_gene not in qry_gene_info:
            continue

        ref_chr1, ref_s1, ref_e1 = ref_gene_info[ref_start_gene]
        ref_chr2, ref_s2, ref_e2 = ref_gene_info[ref_end_gene]
        qry_chr1, qry_s1, qry_e1 = qry_gene_info[qry_start_gene]
        qry_chr2, qry_s2, qry_e2 = qry_gene_info[qry_end_gene]

        # 起止基因不在同一 Chr 上就跳过
        if ref_chr1 != ref_chr2 or qry_chr1 != qry_chr2:
            continue

        ref_chr = ref_chr1
        qry_chr = qry_chr1

        ref_start = min(ref_s1, ref_s2, ref_e1, ref_e2)
        ref_end = max(ref_s1, ref_s2, ref_e1, ref_e2)
        qry_start = min(qry_s1, qry_s2, qry_e1, qry_e2)
        qry_end = max(qry_s1, qry_s2, qry_e1, qry_e2)

        if ref_end <= ref_start or qry_end <= qry_start:
            continue

        n_genes = 2  # 这里用 2 占位即可，paint_segments 主要看长度

        writer.writerow(
            {
                "ref_species": ref_species,
                "qry_species": qry_species,
                "ref_chr": ref_chr,
                "ref_start": ref_start,
                "ref_end": ref_end,
                "qry_chr": qry_chr,
                "qry_start": qry_start,
                "qry_end": qry_end,
                "orientation": orientation,
                "n_genes": n_genes,
                "score": score,
            }
        )
        n_blocks_kept += 1
        span_ref_list.append(ref_end - ref_start)
        span_qry_list.append(qry_end - qry_start)

        segments_by_chr[qry_chr].append((qry_start, qry_end, ref_chr, n_genes))

    f_blocks.close()

    def _median(lst: List[int]) -> float:
        if not lst:
            return 0.0
        lst_sorted = sorted(lst)
        n = len(lst_sorted)
        mid = n // 2
        if n % 2 == 1:
            return float(lst_sorted[mid])
        return (lst_sorted[mid - 1] + lst_sorted[mid]) / 2.0

    summary_rec = {
        "ref_species": ref_species,
        "qry_species": qry_species,
        "n_blocks_raw": n_blocks_raw,
        "n_blocks_kept": n_blocks_kept,
        "mean_block_span_ref": sum(span_ref_list) / len(span_ref_list) if span_ref_list else 0.0,
        "median_block_span_ref": _median(span_ref_list),
        "mean_block_span_qry": sum(span_qry_list) / len(span_qry_list) if span_qry_list else 0.0,
        "median_block_span_qry": _median(span_qry_list),
        "comment": "",
    }

    return segments_by_chr, summary_rec


def merge_segments_for_species(
    species_id: str,
    segments_by_chr: Dict[str, List[Tuple[int, int, str, int]]],
    logger: logging.Logger,
) -> List[Dict[str, object]]:
    """
    将某物种在各 Chr 上的原始 segments 合并：
      - 同一 Chr 内按 start 排序；
      - ref_chr 相同且重叠/相邻的区段合并；
      - 丢弃长度 < MIN_SEGMENT_LENGTH_BP 的碎片。
    """
    records: List[Dict[str, object]] = []

    for chr_name, segs in segments_by_chr.items():
        if not segs:
            continue
        segs_sorted = sorted(segs, key=lambda x: x[0])

        merged: List[Tuple[int, int, str, int]] = []
        cur_start, cur_end, cur_ref_chr, cur_genes = segs_sorted[0]

        for s, e, ref_chr, n_genes in segs_sorted[1:]:
            if ref_chr == cur_ref_chr and s <= cur_end:
                cur_end = max(cur_end, e)
                cur_genes += n_genes
            else:
                merged.append((cur_start, cur_end, cur_ref_chr, cur_genes))
                cur_start, cur_end, cur_ref_chr, cur_genes = s, e, ref_chr, n_genes
        merged.append((cur_start, cur_end, cur_ref_chr, cur_genes))

        for s, e, ref_chr, n_genes in merged:
            if e - s < MIN_SEGMENT_LENGTH_BP:
                continue
            records.append(
                {
                    "species_id": species_id,
                    "chr": chr_name,
                    "start": int(s),
                    "end": int(e),
                    "ref_chr": ref_chr,
                    "n_genes_in_block": int(n_genes),
                }
            )

    logger.info("物种 %s：paint_segments 共 %d 条 segments。", species_id, len(records))
    return records


def build_reference_segments_from_chr_rename(
    ref_species: str,
    logger: logging.Logger,
) -> List[Dict[str, object]]:
    """
    参考物种的染色体 segments：
      - 每条 Chr 作为一个完整区段（0 ~ length_bp）；
      - ref_chr = Chr 名本身。
    """
    chr_rename_file = CHR_RENAME_DIR / f"chr_rename_{ref_species}.tsv"
    if not chr_rename_file.exists():
        logger.error("找不到参考物种 chr_rename 文件：%s", chr_rename_file)
        sys.exit(1)

    records: List[Dict[str, object]] = []

    with chr_rename_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("is_chromosome") or "") != "yes":
                continue
            chr_name = (row.get("new_chr_name") or "").strip()
            if not chr_name:
                continue
            length_bp = int(row.get("length_bp") or "0")
            if length_bp <= 0:
                continue
            records.append(
                {
                    "species_id": ref_species,
                    "chr": chr_name,
                    "start": 0,
                    "end": length_bp,
                    "ref_chr": chr_name,
                    "n_genes_in_block": int(row.get("n_genes") or 0),
                }
            )

    logger.info("参考物种 %s：paint_segments 共 %d 条 Chr 段。", ref_species, len(records))
    return records


# =========================
# 02c — anchors.simple.filtered（画细线）
# =========================

def filter_anchors_for_links(
    anchors_records: List[Tuple[str, str, str, str, float, str]],
    ref_gene_info: Dict[str, Tuple[str, int, int]],
) -> List[Tuple[str, str, str, str, float, str]]:
    """
    按“每条参考 Chr 最多保留 N 个 block”的规则，从 anchors.simple 中挑选
    适合画细线的记录。
    """
    items_with_chr = []
    for rec in anchors_records:
        ref_start_gene = rec[0]
        if ref_start_gene not in ref_gene_info:
            continue
        ref_chr = ref_gene_info[ref_start_gene][0]
        items_with_chr.append((ref_chr, rec))

    # 按 score 降序
    items_with_chr.sort(key=lambda x: x[1][4], reverse=True)

    kept: List[Tuple[str, str, str, str, float, str]] = []
    chr_counts: Dict[str, int] = defaultdict(int)

    for ref_chr, rec in items_with_chr:
        if chr_counts[ref_chr] >= MAX_LINK_BLOCKS_PER_REF_CHR:
            continue
        chr_counts[ref_chr] += 1
        kept.append(rec)

    return kept


def write_filtered_anchors_simple(
    out_path: Path,
    records: List[Tuple[str, str, str, str, float, str]],
) -> None:
    """
    将筛选后的 anchors.simple 记录写出为 6 列文件。
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for ref_start_gene, ref_end_gene, qry_start_gene, qry_end_gene, score, orientation in records:
            f.write(
                f"{ref_start_gene}\t{ref_end_gene}\t"
                f"{qry_start_gene}\t{qry_end_gene}\t"
                f"{score}\t{orientation}\n"
            )


# =========================
# 主流程
# =========================

def main() -> None:
    # 1) 清空输出目录 + 初始化日志
    clean_output_root(OUTPUT_ROOT)
    logger = setup_logging(OUTPUT_ROOT / "logs")

    # 2) 读取物种列表与参考物种
    species_ids, ref_species = load_species_meta(SPECIES_META_FILE, logger)

    # 2b) 自动构建物种短名映射
    global SHORT_NAME_MAP
    SHORT_NAME_MAP = build_short_name_map(species_ids, logger)
    ref_short = get_short_name(ref_species)

    # 2c) 准备 jcvi dotplot 字体配置（若启用）
    global MPLCONFIG_DIR_FOR_JCVI
    MPLCONFIG_DIR_FOR_JCVI = prepare_mpl_config_for_jcvi(logger)

    # 3) geneorder 预载
    gene_info_cache: Dict[str, Dict[str, Tuple[str, int, int]]] = {}
    chr_intervals_cache: Dict[str, Dict[str, List[Tuple[int, int, str]]]] = {}
    for sid in species_ids:
        gi, ci = load_geneorder_bed(sid, logger)
        gene_info_cache[sid] = gi
        chr_intervals_cache[sid] = ci

    # 4) 根据需要调用 jcvi 生成 anchors.simple
    if RUN_JCVI_PIPELINE:
        logger.info("RUN_JCVI_PIPELINE = True，将在本脚本中统一调用 jcvi 生成 anchors.simple。")
        for sid in species_ids:
            if sid == ref_species:
                continue
            run_jcvi_for_pair(ref_species, sid, logger)
    else:
        logger.info("RUN_JCVI_PIPELINE = False，假定 raw_jcvi/ 下已有 anchors.simple。")

    # 5) 针对每个物种对构建 blocks_normalized + anchors.simple.filtered
    pair_summary: List[Dict[str, object]] = []
    all_segments_per_species: Dict[str, Dict[str, List[Tuple[int, int, str, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for sid in species_ids:
        if sid == ref_species:
            continue
        qry_species = sid
        qry_short = get_short_name(qry_species)

        anchors_simple_path = RAW_JCVI_DIR / f"{ref_short}.{qry_short}.anchors.simple"
        if not anchors_simple_path.exists():
            logger.error("缺少 anchors.simple：%s，跳过该物种对。", anchors_simple_path)
            continue

        anchors_records = parse_anchors_simple(anchors_simple_path, logger)
        if not anchors_records:
            logger.warning("anchors.simple 为空：%s，跳过该物种对。", anchors_simple_path)
            continue

        ref_gene_info = gene_info_cache[ref_species]
        qry_gene_info = gene_info_cache[qry_species]

        blocks_path = BLOCKS_NORMALIZED_DIR / f"{ref_short}__vs__{qry_short}.blocks.tsv"
        segments_by_chr, summary_rec = build_blocks_normalized_and_segments_for_pair(
            ref_species=ref_species,
            qry_species=qry_species,
            ref_gene_info=ref_gene_info,
            qry_gene_info=qry_gene_info,
            anchors_records=anchors_records,
            blocks_path=blocks_path,
        )
        pair_summary.append(summary_rec)
        all_segments_per_species[qry_species] = segments_by_chr
        logger.info(
            "blocks_normalized 已写出：%s（kept=%d / raw=%d）",
            blocks_path,
            summary_rec["n_blocks_kept"],
            summary_rec["n_blocks_raw"],
        )

        # 02c：生成 anchors.simple.filtered
        filtered_records = filter_anchors_for_links(anchors_records, ref_gene_info)
        out_filtered = ANCHORS_FOR_LINKS_DIR / f"{ref_short}__vs__{qry_short}.anchors.simple.filtered"
        write_filtered_anchors_simple(out_filtered, filtered_records)
        logger.info(
            "anchors.simple.filtered 已写出：%s（保留 %d / %d 个 blocks）",
            out_filtered,
            len(filtered_records),
            len(anchors_records),
        )

    # 6) pairwise_block_summary.tsv
    summary_path = OUTPUT_ROOT / "pairwise_block_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as f_sum:
        fieldnames = [
            "ref_species",
            "qry_species",
            "n_blocks_raw",
            "n_blocks_kept",
            "mean_block_span_ref",
            "median_block_span_ref",
            "mean_block_span_qry",
            "median_block_span_qry",
            "comment",
        ]
        writer = csv.DictWriter(f_sum, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in pair_summary:
            writer.writerow(rec)
    logger.info("pairwise_block_summary.tsv 已写出：%s", summary_path)

    # 7) paint_segments：参考物种 + 其它物种
    paint_fieldnames = [
        "species_id",
        "chr",
        "start",
        "end",
        "ref_chr",
        "n_genes_in_block",
    ]

    # 7.1 参考物种
    ref_segments_records = build_reference_segments_from_chr_rename(
        ref_species, logger
    )
    paint_ref_path = PAINT_SEGMENTS_DIR / f"paint_segments_{ref_species}.tsv"
    with paint_ref_path.open("w", encoding="utf-8", newline="") as f_ref:
        writer = csv.DictWriter(f_ref, fieldnames=paint_fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in ref_segments_records:
            writer.writerow(rec)
    logger.info("参考物种 paint_segments 文件写出：%s", paint_ref_path)

    # 7.2 其它物种
    for sid in species_ids:
        if sid == ref_species:
            continue
        segs_by_chr = all_segments_per_species.get(sid, {})
        records = merge_segments_for_species(
            species_id=sid,
            segments_by_chr=segs_by_chr,
            logger=logger,
        )
        out_path = PAINT_SEGMENTS_DIR / f"paint_segments_{sid}.tsv"
        with out_path.open("w", encoding="utf-8", newline="") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=paint_fieldnames, delimiter="\t")
            writer.writeheader()
            for rec in records:
                writer.writerow(rec)
        logger.info("物种 %s 的 paint_segments 文件写出：%s", sid, out_path)

    logger.info("synteny02 完成：blocks_normalized + paint_segments + anchors.simple.filtered 均已生成。")


if __name__ == "__main__":
    main()

