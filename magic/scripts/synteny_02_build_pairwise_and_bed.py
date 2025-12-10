#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_02_build_pairwise_and_bed.py
—— 从 anchors_global.tsv 构建 pairwise anchors.simple + per-species bed 文件（链式共线性）

职责：
  1) 读取 synteny_species_meta.tsv，确定参与 synteny 的物种列表及顺序：
       - 当前简化 meta：只需要 species_id（按行顺序作为绘图顺序）
       - 将 meta 中 16 个贝类按出现顺序组成一条“物种链”
         例如：[Ari, Pma, Cal, ... , Sri]，再生成相邻物种 pair：
           (Ari, Pma), (Pma, Cal), ..., (Sco, Sri)
  2) 读取 synteny_01 的 anchors_global.tsv：
       - 这是长表，每行是一个 (OG, 物种, 基因) 的记录；
       - 包含：og_id, species_id, id_core, gene_id, chr, start, end, ref_chr 等。
  3) 为每个物种生成 bed 文件（供 JCVI 使用）：
       - 只保留 anchors_global 中出现过的基因（主染色体上，已在 01 号脚本过滤）
       - 使用 id_core 作为 bed 的 gene_id（供 anchors.simple 对接）
  4) 为每对相邻物种生成 anchors.simple：
       - 对每个 OG：
           若该 OG 同时在物种 A 与 B 中都有成员，则生成一条锚点：
             id_core_A [TAB] id_core_B
       - anchors.simple 中的 ID 与 bed 的 gene_id 完全一致。
  5) 输出 pairwise_summary.tsv，用于快速查看每个物种对的锚点数量。

输入：
  - raw_data/synteny_species_meta.tsv
  - output/synteny_01_global_anchors/anchors_global.tsv

输出（每次运行前自动删除整个输出目录，重新创建）：
  - output/synteny_02_pairwise_bed/
      ├── bed/
      │     ├── <species_id>.bed
      │     └── ...
      ├── anchors/
      │     ├── <A>__vs__<B>.anchors.simple
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
from typing import Dict, List, Tuple

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

LOG_LEVEL = "INFO"

# 最小锚点数量阈值：如果某个物种对的锚点少于该值，可以在 summary 里提示（不会过滤）
MIN_ANCHORS_WARN_THRESHOLD = 50


# =========================
# 基础工具函数
# =========================

def setup_logging() -> logging.Logger:
    """
    初始化日志系统：
      - 屏幕 + 文件双通道输出；
      - 日志文件放在 output/synteny_02_pairwise_bed/logs/ 下。
    """
    log_dir = OUTPUT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny_02_build_pairwise_and_bed.log"

    logger = logging.getLogger("synteny_pairwise_bed")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清理旧 handler，避免重复输出
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

    logger.info("========== synteny_02 — pairwise anchors + bed ==========")
    logger.info("PROJECT_ROOT      = %s", PROJECT_ROOT)
    logger.info("RAW_DATA_DIR      = %s", RAW_DATA_DIR)
    logger.info("ANCHORS_GLOBAL    = %s", ANCHORS_GLOBAL_FILE)
    logger.info("OUTPUT_ROOT       = %s", OUTPUT_ROOT)

    return logger


def clean_output_root(logger: logging.Logger) -> None:
    """
    按皇上要求：
      在每次运行脚本前，删除本脚本对应的输出目录，
      再重新创建，确保不会被旧结果污染。
    """
    if OUTPUT_ROOT.exists():
        logger.info("删除旧输出目录：%s", OUTPUT_ROOT)
        shutil.rmtree(OUTPUT_ROOT)
    BED_DIR.mkdir(parents=True, exist_ok=True)
    ANCHORS_DIR.mkdir(parents=True, exist_ok=True)
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
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("meta 表必须包含列：species_id")
            sys.exit(1)

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            species_ids.append(sid)

    if not species_ids:
        logger.error("meta 表中没有有效 species_id")
        sys.exit(1)

    logger.info("读取 meta 物种顺序（共 %d 个）：%s", len(species_ids), ", ".join(species_ids))
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
        logger.error("未找到 anchors_global.tsv：%s", anchors_file)
        sys.exit(1)

    og_to_species: Dict[str, Dict[str, dict]] = {}
    species_to_records: Dict[str, List[dict]] = {}

    n_lines = 0

    logger.info("读取 anchors_global：%s", anchors_file)
    with anchors_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
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

            # 记录到 og_to_species
            og_map = og_to_species.setdefault(og_id, {})
            # 单拷贝 OG 情况下，每个 (og_id, species_id) 应该只有一个记录
            # 如出现多个，这里后写会覆盖先写，问题不大
            og_map[sid] = row

            # 记录到 species_to_records
            species_to_records.setdefault(sid, []).append(row)

    logger.info("anchors_global 总记录数 = %d", n_lines)
    logger.info("anchors_global 不同 OG 数 = %d", len(og_to_species))
    return og_to_species, species_to_records


# =========================
# 核心构建：bed + anchors
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
            logger.warning("物种 %s 在 anchors_global 中没有记录，将不会生成 bed 文件。", sid)
            continue

        out_path = BED_DIR / f"{sid}.bed"
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
                        id_core,   # gene_id = id_core
                        full_id,   # full_id = GFF 中解析的 ID
                        og_id,
                        ref_chr,
                        sid,
                    ]
                )

        logger.info("写出 bed 文件：%s （记录数 = %d）", out_path, len(records))


def build_pairwise_anchors(
    species_ids: List[str],
    og_to_species: Dict[str, Dict[str, dict]],
    logger: logging.Logger,
) -> List[Dict[str, object]]:
    """
    构建链式物种对的 anchors.simple 文件：

    物种链：
      [S1, S2, S3, ..., Sn]
    生成物种对：
      (S1, S2), (S2, S3), ..., (S_{n-1}, S_n)

    对每个 pair (A, B)：
      遍历所有 og_id：
        若该 OG 同时有物种 A 和 B 的记录，则生成一条锚点：
          anchors.simple 行： id_core_A  TAB  id_core_B

    返回：
      pairwise_summary 记录列表，用于写 pairwise_summary.tsv。
    """
    # 构建链式物种对
    pairs: List[Tuple[str, str]] = []
    for i in range(len(species_ids) - 1):
        a = species_ids[i]
        b = species_ids[i + 1]
        pairs.append((a, b))

    logger.info("链式物种对（共 %d 对）：", len(pairs))
    for a, b in pairs:
        logger.info("  %s  <==>  %s", a, b)

    # 为每对物种创建 anchors.simple
    summary_records: List[Dict[str, object]] = []

    for a, b in pairs:
        pair_id = f"{a}__vs__{b}"
        anchors_path = ANCHORS_DIR / f"{pair_id}.anchors.simple"

        n_anchors = 0
        n_ogs_pair = 0

        with anchors_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")

            # 遍历所有 OG
            for og_id, species_map in og_to_species.items():
                rec_a = species_map.get(a)
                rec_b = species_map.get(b)
                if rec_a is None or rec_b is None:
                    continue

                # 提取 id_core（必须与 bed 中的 gene_id 一致）
                id_a = (rec_a.get("id_core") or "").strip()
                id_b = (rec_b.get("id_core") or "").strip()
                if not id_a or not id_b:
                    continue

                writer.writerow([id_a, id_b])
                n_anchors += 1
                n_ogs_pair += 1  # 在单拷贝 OG 情况下，一对物种 per OG 只会有一条锚点

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

        summary_records.append(
            {
                "pair_id": pair_id,
                "species_a": a,
                "species_b": b,
                "n_anchors": n_anchors,
                "n_ogs": n_ogs_pair,
                "comment": warn_flag,
            }
        )

    return summary_records


def write_pairwise_summary(
    summary_records: List[Dict[str, object]],
    logger: logging.Logger,
) -> None:
    """
    写出 pairwise_summary.tsv：
      - pair_id
      - species_a
      - species_b
      - n_anchors
      - n_ogs
      - comment
    """
    out_path = OUTPUT_ROOT / "pairwise_summary.tsv"
    fieldnames = ["pair_id", "species_a", "species_b", "n_anchors", "n_ogs", "comment"]
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
    logger = setup_logging()
    clean_output_root(logger)

    # 1) 读物种顺序（按 meta 行顺序）
    species_ids = load_species_order(SPECIES_META_FILE, logger)

    # 2) 读取 anchors_global，构建 og_to_species + species_to_records
    og_to_species, species_to_records = load_anchors_global(ANCHORS_GLOBAL_FILE, logger)

    # 3) 为每个物种生成 bed 文件
    build_bed_files(species_ids, species_to_records, logger)

    # 4) 构建链式 pairwise anchors.simple
    summary_records = build_pairwise_anchors(species_ids, og_to_species, logger)

    # 5) 写出 pairwise_summary.tsv
    write_pairwise_summary(summary_records, logger)

    logger.info("synteny_02 完成。")


if __name__ == "__main__":
    main()

