#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_02_build_pairwise_and_bed.py
—— 从 anchors_global.tsv 构建：
     1) 链式 pairwise anchors.simple（2 列基因配对）
     2) 每物种 bed 文件
     3) 每物种对的“区段级” simple_blocks（6 列 block simple，使用 gene_id）

职责（v3·修正版）：
  1) 读取 synteny_species_meta.tsv，确定参与 synteny 的物种列表及顺序：
       - meta 中的 species_id（按行顺序作为绘图顺序）
       - 将 meta 中 16 个贝类按出现顺序组成一条“物种链”
         例如：[Ari, Pma, Cal, ... , Sri]，再生成相邻物种 pair：
           (Ari, Pma), (Pma, Cal), ..., (Sco, Sri)
  2) 读取 synteny_01 的 anchors_global.tsv：
       - 这是长表，每行是一个 (OG, 物种, 基因) 的记录；
       - 包含：og_id, species_id, id_core, gene_id, chr, start, end, ref_chr 等。
       - 我们只使用 synteny_01 已经筛好的锚点，不再做额外过滤。
  3) 为每个物种写出 bed 文件：
       - chr, start, end, gene_id, full_id, og_id, ref_chr, species_id
       - gene_id 使用 id_core（与 OrthoFinder / anchors.simple 对齐），
         full_id 保留 GFF 中解析出的原始 ID。
  4) 为链式物种对写出 2 列 anchors.simple：
       - 每个 OG 在 (A, B) 中都存在时，写出：
             id_core_A [TAB] id_core_B
       - anchors.simple 中的 ID 与 bed 的 gene_id 完全一致。
  5) 基于 bed + anchors.simple，为每个物种对自动聚类“共线区段”，写出：
       - simple_blocks/<A>__vs__<B>.anchors.simple
       - 这是 JCVI block 模式可直接使用的 6 列 simple 文件，每行代表一个 synteny block：
             geneA_start   geneA_end   geneB_start   geneB_end   score   orientation
         其中 geneA_start/geneA_end/geneB_start/geneB_end 都是 gene_id，
         JCVI 会用这些基因 ID 去 bed 文件中找坐标。
  6) 输出 pairwise_summary.tsv，用于快速查看每个物种对的锚点与 block 数量。

输入：
  - raw_data/synteny_species_meta.tsv
  - output/synteny_01_global_anchors/anchors_global.tsv

输出（每次运行前自动删除整个输出目录，重新创建）：
  - output/synteny_02_pairwise_bed/
      ├── bed/
      │     ├── <species_id>.bed
      │     └── ...
      ├── anchors/
      │     ├── <A>__vs__<B>.anchors.simple   （2 列基因配对）
      │     └── ...
      ├── simple_blocks/
      │     ├── <A>__vs__<B>.anchors.simple   （6 列 block simple，gene_id）
      │     └── ...
      ├── pairwise_summary.tsv
      └── logs/synteny_02_build_pairwise_and_bed.log

说明：
  - 物种顺序完全按照 synteny_species_meta.tsv 中的行顺序；
  - 只使用 synteny_01 中已经筛选好的单拷贝 / 高质量锚点，不再重复过滤；
  - 代码不依赖命令行参数，所有路径与参数在顶部统一配置，方便皇上修改。
"""

from __future__ import annotations

import sys
import csv
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any
import statistics

# =========================
# 参数区（皇上可在此修改）
# =========================

# 项目根目录（假定本脚本位于 magic/scripts/ 下）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"

SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

ANCHORS_GLOBAL_DIR = PROJECT_ROOT / "output" / "synteny_01_global_anchors"
ANCHORS_GLOBAL_FILE = ANCHORS_GLOBAL_DIR / "anchors_global.tsv"

OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_02_pairwise_bed"
BED_DIR = OUTPUT_ROOT / "bed"
ANCHORS_DIR = OUTPUT_ROOT / "anchors"
SIMPLE_BLOCK_DIR = OUTPUT_ROOT / "simple_blocks"

LOG_LEVEL = "INFO"

# 最小锚点数量阈值：如果某个物种对的锚点少于该值，在 summary 里提示（不会过滤）
MIN_ANCHORS_WARN_THRESHOLD = 50

# —— 区段级 simple_blocks 聚类参数 —— #
# 最少多少个 anchor 才认为是一个 block
MIN_ANCHORS_PER_BLOCK = 5

# 在参考物种 A 上相邻 anchor 允许的最大距离（bp）
MAX_GAP_A_BP = 5_000_000

# 在物种 B 上相邻 anchor 允许的最大距离（bp）
MAX_GAP_B_BP = 5_000_000

# block 在任一物种上的最小物理跨度（bp），小于该值的 block 丢弃；设为 0 则不限制
MIN_BLOCK_SPAN_BP = 0


# =========================
# 基础工具函数
# =========================


def get_logger() -> logging.Logger:
    """
    配置日志系统：
      - 屏幕 + 日志文件双通道输出；
      - 统一时间格式。
    """
    log_dir = OUTPUT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny_02_build_pairwise_and_bed.log"

    logger = logging.getLogger("synteny_pairwise_bed")
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

    logger.info("========== synteny_02 — pairwise anchors + bed + simple_blocks ==========")
    logger.info("PROJECT_ROOT      = %s", PROJECT_ROOT)
    logger.info("RAW_DATA_DIR      = %s", RAW_DATA_DIR)
    logger.info("ANCHORS_GLOBAL    = %s", ANCHORS_GLOBAL_FILE)
    logger.info("OUTPUT_ROOT       = %s", OUTPUT_ROOT)

    return logger


def clean_output_root(logger: logging.Logger) -> None:
    """
    在每次运行脚本前，删除本脚本对应的输出目录，
    再重新创建，确保不会被旧结果污染。
    """
    if OUTPUT_ROOT.exists():
        logger.info("删除旧输出目录：%s", OUTPUT_ROOT)
        shutil.rmtree(OUTPUT_ROOT)
    BED_DIR.mkdir(parents=True, exist_ok=True)
    ANCHORS_DIR.mkdir(parents=True, exist_ok=True)
    SIMPLE_BLOCK_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "logs").mkdir(parents=True, exist_ok=True)


# =========================
# 数据读取与构建
# =========================


def load_species_order(meta_file: Path, logger: logging.Logger) -> List[str]:
    """
    读取 synteny_species_meta.tsv，获取 synteny 物种顺序。

    当前简化版 meta 至少包含：
      - species_id

    行顺序即物种的垂直顺序（即图中的堆叠顺序，也用于构建链式物种对）。
    """
    if not meta_file.exists():
        logger.error("species meta 文件不存在：%s", meta_file)
        sys.exit(1)

    species_ids: List[str] = []

    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "species_id" not in (reader.fieldnames or []):
            logger.error("species meta 文件缺少 species_id 列：%s", meta_file)
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            species_ids.append(sid)

    if len(species_ids) < 2:
        logger.error("参与 synteny 的物种数量 < 2，无法构建物种对。")
        sys.exit(1)

    logger.info("synteny 物种顺序（共 %d 个）：%s", len(species_ids), ", ".join(species_ids))
    return species_ids


def load_anchors_global(
    anchors_file: Path,
    logger: logging.Logger,
) -> Tuple[Dict[str, Dict[str, dict]], Dict[str, List[dict]]]:
    """
    读取 anchors_global.tsv，返回两个结构：

    1) og_to_species:
         og_id -> { species_id -> record }
       用于构建 pairwise anchors.simple。

    2) species_to_records:
         species_id -> [record, record, ...]
       用于构建 per-species bed 文件。

    每个 record 至少包含字段：
      - og_id
      - species_id
      - id_core
      - gene_id
      - chr
      - start
      - end
      - ref_chr
    """
    if not anchors_file.exists():
        logger.error("anchors_global.tsv 文件不存在：%s", anchors_file)
        sys.exit(1)

    og_to_species: Dict[str, Dict[str, dict]] = {}
    species_to_records: Dict[str, List[dict]] = {}

    required = {
        "og_id",
        "species_id",
        "id_core",
        "gene_id",
        "chr",
        "start",
        "end",
        "ref_chr",
    }

    n_lines = 0
    with anchors_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            logger.error("anchors_global.tsv 缺少必要列：%s", ",".join(sorted(missing)))
            sys.exit(1)

        for row in reader:
            n_lines += 1
            og_id = (row.get("og_id") or "").strip()
            sid = (row.get("species_id") or "").strip()
            if not og_id or not sid:
                continue

            og_map = og_to_species.setdefault(og_id, {})
            og_map[sid] = row

            species_to_records.setdefault(sid, []).append(row)

    logger.info("anchors_global 总记录数 = %d", n_lines)
    logger.info("anchors_global 不同 OG 数 = %d", len(og_to_species))
    return og_to_species, species_to_records


# =========================
# 核心构建：bed
# =========================


def build_bed_files(
    species_ids: List[str],
    species_to_records: Dict[str, List[dict]],
    logger: logging.Logger,
) -> None:
    """
    为每个物种生成 bed 文件：

    bed 文件格式（带 header，制表符分隔）：
      - chr        重命名后的染色体名（Chr1, Chr2, ...）
      - start      起始坐标（从 anchors_global 直接继承）
      - end        终止坐标
      - gene_id    用于 JCVI 的基因 ID，这里统一使用 id_core（即最后一个“|”之后的 ID）
      - full_id    原始 GFF 中解析的 gene_id（transcript_id / orig_transcript_id / ID）
      - og_id      Orthogroup ID
      - ref_chr    该 OG 在参考物种中的 ALG 染色体（用于后续分析）
      - species_id 物种 ID（方便调试）

    注意：
      - anchors_simple 中会使用 gene_id（即 id_core），因此 bed 的 gene_id 必须与 id_core 一致；
      - full_id 仅作信息保留，不参与 JCVI 锚点匹配。
    """
    for sid in species_ids:
        records = species_to_records.get(sid, [])
        if not records:
            logger.warning("物种 %s 在 anchors_global 中没有记录，将不会生成 bed。", sid)
            continue

        out_path = BED_DIR / f"{sid}.bed"
        n_written = 0
        with out_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(
                [
                    "chr",
                    "start",
                    "end",
                    "gene_id",
                    "full_id",
                    "og_id",
                    "ref_chr",
                    "species_id",
                ]
            )
            for rec in records:
                chr_name = (rec.get("chr") or "").strip()
                start = rec.get("start")
                end = rec.get("end")
                id_core = (rec.get("id_core") or "").strip()
                full_id = (rec.get("gene_id") or "").strip()
                og_id = (rec.get("og_id") or "").strip()
                ref_chr = (rec.get("ref_chr") or "").strip()

                if not chr_name or not id_core:
                    continue

                writer.writerow(
                    [
                        chr_name,
                        start,
                        end,
                        id_core,  # gene_id = id_core
                        full_id,  # full_id = GFF 中解析的 ID
                        og_id,
                        ref_chr,
                        sid,
                    ]
                )
                n_written += 1

        logger.info("写出 bed 文件：%s （记录数 = %d）", out_path, n_written)


def build_gene_position_map(
    species_ids: List[str],
    logger: logging.Logger,
) -> Dict[str, Dict[str, Tuple[str, int]]]:
    """
    从 bed 文件中构建基因 -> (chr, 中心坐标) 的映射，用于后续 block 聚类。

    返回：
      gene_pos_map:
        species_id -> { gene_id -> (chr, mid_pos) }
    """
    gene_pos_map: Dict[str, Dict[str, Tuple[str, int]]] = {}

    for sid in species_ids:
        bed_file = BED_DIR / f"{sid}.bed"
        if not bed_file.exists():
            logger.warning("找不到 bed 文件（跳过构建位置信息）：%s", bed_file)
            continue

        species_map: Dict[str, Tuple[str, int]] = {}
        with bed_file.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                chr_name = (row.get("chr") or "").strip()
                gene_id = (row.get("gene_id") or "").strip()
                if not chr_name or not gene_id:
                    continue
                try:
                    start = int(float(row.get("start", "0")))
                    end = int(float(row.get("end", "0")))
                except ValueError:
                    continue
                mid = (start + end) // 2
                species_map[gene_id] = (chr_name, mid)

        gene_pos_map[sid] = species_map
        logger.info(
            "物种 %s：从 bed 解析到 %d 个基因的位置信息。",
            sid,
            len(species_map),
        )

    return gene_pos_map


# =========================
# 核心构建：pairwise anchors + blocks
# =========================


def build_blocks_for_pair(
    species_a: str,
    species_b: str,
    gene_pos_map: Dict[str, Dict[str, Tuple[str, int]]],
    anchors_file: Path,
    out_block_file: Path,
    logger: logging.Logger,
) -> Tuple[int, List[int]]:
    """
    基于 2 列 anchors.simple + bed 中的坐标，为一个物种对构建“区段级” simple_blocks。

    简单聚类策略：
      1) 从 anchors.simple 读取 (geneA, geneB)；
      2) 利用 bed 的 gene_pos_map，得到：
           (chrA, posA, chrB, posB)
         若任一基因缺少位置信息，则丢弃该锚点；
      3) 按 (chrA, chrB) 分组，每组内按 posA 升序排序；
      4) 顺序扫描，相邻 anchor 若满足：
           |posA_i - posA_{i-1}| <= MAX_GAP_A_BP
           且 |posB_i - posB_{i-1}| <= MAX_GAP_B_BP
         则视为同一 block，反之则开启新 block；
      5) 对每个 block：
           - 若锚点数量 < MIN_ANCHORS_PER_BLOCK 丢弃；
           - 若跨度 < MIN_BLOCK_SPAN_BP（可选）丢弃；
         否则写出 6 列 simple（全部为 gene_id）：
           geneA_start   geneA_end   geneB_start   geneB_end   score   orientation
         其中：
           - geneA_start / geneA_end：按 posA 排序后的首、尾锚点的 geneA；
           - geneB_start / geneB_end：
                若 orientation = "+"，取首、尾锚点的 geneB；
                若 orientation = "-"，取尾、首锚点的 geneB；
           - score = block 中锚点数量；
           - orientation = "+" 或 "-"，由 posB 整体趋势决定。
    """
    if not anchors_file.exists():
        logger.warning(
            "物种对 %s vs %s 的 anchors.simple 不存在，无法构建 blocks：%s",
            species_a,
            species_b,
            anchors_file,
        )
        out_block_file.touch()
        return 0, []

    pos_a = gene_pos_map.get(species_a, {})
    pos_b = gene_pos_map.get(species_b, {})

    # 先按 (chrA, chrB) 分组收集 anchor 位置信息
    grouped: Dict[Tuple[str, str], List[Tuple[str, int, str, int]]] = {}
    n_total_anchors = 0
    n_missing_pos = 0

    with anchors_file.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            g_a, g_b = parts[0].strip(), parts[1].strip()
            if not g_a or not g_b:
                continue

            info_a = pos_a.get(g_a)
            info_b = pos_b.get(g_b)
            if info_a is None or info_b is None:
                n_missing_pos += 1
                continue

            chr_a, x = info_a
            chr_b, y = info_b
            grouped.setdefault((chr_a, chr_b), []).append((g_a, x, g_b, y))
            n_total_anchors += 1

    if n_total_anchors == 0:
        logger.warning(
            "物种对 %s vs %s：在 bed 中找不到任何有效坐标的 anchors（anchors=%d, missing_pos=%d），不生成 blocks。",
            species_a,
            species_b,
            n_total_anchors,
            n_missing_pos,
        )
        out_block_file.touch()
        return 0, []

    out_block_file.parent.mkdir(parents=True, exist_ok=True)
    n_blocks = 0
    block_sizes: List[int] = []

    with out_block_file.open("w", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")

        for (_chr_a, _chr_b), entries in grouped.items():
            # entries: List[(geneA, posA, geneB, posB)]
            if len(entries) < MIN_ANCHORS_PER_BLOCK:
                continue

            entries.sort(key=lambda t: t[1])  # 按 posA 排序
            current_block: List[Tuple[str, int, str, int]] = [entries[0]]
            prev_x, prev_y = entries[0][1], entries[0][3]

            def flush_block(block: List[Tuple[str, int, str, int]]) -> None:
                nonlocal n_blocks
                if len(block) < MIN_ANCHORS_PER_BLOCK:
                    return
                xs = [p[1] for p in block]
                ys = [p[3] for p in block]
                span_a = max(xs) - min(xs)
                span_b = max(ys) - min(ys)
                if MIN_BLOCK_SPAN_BP > 0 and (span_a < MIN_BLOCK_SPAN_BP and span_b < MIN_BLOCK_SPAN_BP):
                    return

                # orientation：看 B 物种的整体趋势
                orient = "+" if block[-1][3] >= block[0][3] else "-"

                # A 物种的起止基因：按 posA 排序的首尾
                geneA_start = block[0][0]
                geneA_end = block[-1][0]

                # B 物种的起止基因：根据方向决定
                if orient == "+":
                    geneB_start = block[0][2]
                    geneB_end = block[-1][2]
                else:
                    geneB_start = block[-1][2]
                    geneB_end = block[0][2]

                score = len(block)

                writer.writerow(
                    [
                        geneA_start,
                        geneA_end,
                        geneB_start,
                        geneB_end,
                        score,
                        orient,
                    ]
                )
                n_blocks += 1
                block_sizes.append(score)

            for (g_a, x, g_b, y) in entries[1:]:
                if abs(x - prev_x) <= MAX_GAP_A_BP and abs(y - prev_y) <= MAX_GAP_B_BP:
                    current_block.append((g_a, x, g_b, y))
                else:
                    flush_block(current_block)
                    current_block = [(g_a, x, g_b, y)]
                prev_x, prev_y = x, y

            flush_block(current_block)

    logger.info(
        "物种对 %s vs %s：anchors 有效=%d，缺少坐标=%d，生成 blocks=%d",
        species_a,
        species_b,
        n_total_anchors,
        n_missing_pos,
        n_blocks,
    )
    return n_blocks, block_sizes


def build_pairwise_anchors(
    species_ids: List[str],
    og_to_species: Dict[str, Dict[str, dict]],
    gene_pos_map: Dict[str, Dict[str, Tuple[str, int]]],
    logger: logging.Logger,
) -> List[Dict[str, Any]]:
    """
    构建链式物种对的 anchors.simple 文件 + 区段级 simple_blocks：

    物种链：
      [S1, S2, S3, ..., Sn]
    生成物种对：
      (S1, S2), (S2, S3), ..., (S_{n-1}, S_n)

    对每个 pair (A, B)：
      1) 遍历所有 og_id：
           若该 OG 同时有物种 A 和 B 的记录，则生成一条锚点：
             anchors.simple 行： id_core_A  TAB  id_core_B
      2) 用 anchors.simple + bed 位置信息，聚类区段，生成 simple_blocks。
    """
    pairs: List[Tuple[str, str]] = []
    for i in range(len(species_ids) - 1):
        a = species_ids[i]
        b = species_ids[i + 1]
        pairs.append((a, b))

    logger.info("链式物种对（共 %d 对）：", len(pairs))
    for a, b in pairs:
        logger.info("  %s  <==>  %s", a, b)

    summary_records: List[Dict[str, Any]] = []

    for a, b in pairs:
        pair_id = f"{a}__vs__{b}"
        anchors_path = ANCHORS_DIR / f"{pair_id}.anchors.simple"

        n_anchors = 0
        n_ogs_pair = 0

        with anchors_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")

            for og_id, species_map in og_to_species.items():
                rec_a = species_map.get(a)
                rec_b = species_map.get(b)
                if rec_a is None or rec_b is None:
                    continue

                id_a = (rec_a.get("id_core") or "").strip()
                id_b = (rec_b.get("id_core") or "").strip()
                if not id_a or not id_b:
                    continue

                writer.writerow([id_a, id_b])
                n_anchors += 1
                n_ogs_pair += 1

        warn_flag = ""
        if n_anchors < MIN_ANCHORS_WARN_THRESHOLD:
            warn_flag = f"<{MIN_ANCHORS_WARN_THRESHOLD}_anchors"

        logger.info(
            "写出 anchors.simple：%s （anchors = %d, ogs = %d）%s",
            anchors_path,
            n_anchors,
            n_ogs_pair,
            f"  [WARNING: {warn_flag}]" if warn_flag else "",
        )

        block_file = SIMPLE_BLOCK_DIR / f"{pair_id}.anchors.simple"
        n_blocks, block_sizes = build_blocks_for_pair(
            a,
            b,
            gene_pos_map,
            anchors_path,
            block_file,
            logger,
        )
        median_block_size = statistics.median(block_sizes) if block_sizes else 0
        max_block_size = max(block_sizes) if block_sizes else 0

        summary_records.append(
            {
                "pair_id": pair_id,
                "species_a": a,
                "species_b": b,
                "n_anchors": n_anchors,
                "n_ogs": n_ogs_pair,
                "n_blocks": n_blocks,
                "median_block_size": median_block_size,
                "max_block_size": max_block_size,
                "comment": warn_flag,
            }
        )

    return summary_records


def write_pairwise_summary(
    summary_records: List[Dict[str, Any]],
    logger: logging.Logger,
) -> None:
    """
    将 pairwise_summary.tsv 写出，便于皇上快速查看每对物种的锚点与 block 概况。
    """
    out_path = OUTPUT_ROOT / "pairwise_summary.tsv"
    fieldnames = [
        "pair_id",
        "species_a",
        "species_b",
        "n_anchors",
        "n_ogs",
        "n_blocks",
        "median_block_size",
        "max_block_size",
        "comment",
    ]

    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in summary_records:
            writer.writerow(rec)

    logger.info("写出 pairwise_summary.tsv：%s", out_path)


# =========================
# 主流程
# =========================


def main() -> None:
    logger = get_logger()

    clean_output_root(logger)

    species_ids = load_species_order(SPECIES_META_FILE, logger)

    og_to_species, species_to_records = load_anchors_global(ANCHORS_GLOBAL_FILE, logger)

    build_bed_files(species_ids, species_to_records, logger)

    gene_pos_map = build_gene_position_map(species_ids, logger)

    summary_records = build_pairwise_anchors(
        species_ids,
        og_to_species,
        gene_pos_map,
        logger,
    )

    write_pairwise_summary(summary_records, logger)

    logger.info("synteny_02 完成。")


if __name__ == "__main__":
    main()

