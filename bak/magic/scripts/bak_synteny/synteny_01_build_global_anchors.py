#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_01_build_global_anchors.py
—— 从 OrthoFinder + GFF + chr_rename 构建全局锚点表（v5 · transcript_id 家族优先版）

核心职责：
  1) 读取 synteny_species_meta.tsv，锁定物种列表和参考物种；
  2) 读取 synteny_00 的 chr_rename_<species_id>.tsv，确定主染色体并重命名；
  3) 为每个物种：
       - 自动寻找 GFF 文件；
       - 解析 GFF 中 mRNA/transcript 的坐标；
       - ★ 优先使用各种 *transcript_id 字段（transcript_id、orig_transcript_id、以及未来的 xxxtranscript_id）作为转录本 ID；
       - 若所有 transcript_id 家族都不存在，退回使用 ID= 字段；
       - 然后对该 ID 应用“最后一个 | 之后”为 id_core 的规则；
       - 只保留位于主染色体上的转录本；
       - 建立 id_core → 基因坐标信息映射。
  4) 读取 OrthoFinder Orthogroups.tsv + Orthogroups_SingleCopyOrthologues.txt：
       - 从每个蛋白 ID 提取 id_core（最后一个“|”之后的字符串）；
       - 用 id_core 关联到各物种 GFF 转录本；
       - 以参考物种成员定义 ref_chr、ref_start、ref_end；
       - 生成 anchors_global.tsv（基因级锚点表）。
  5) 输出：
       - anchors_global.tsv
       - id_mapping_qc.tsv
       - anchors_global_summary.tsv
       - 完整日志。

重要约定：
  - OrthoFinder 蛋白 ID：最后一个“|”之后为 id_core；
  - GFF：优先 transcript_id 家族，再退回 ID=；
  - 然后同样用“最后一个 | 之后”作为 id_core；
  - 只有主染色体上的转录本参与锚点构建；
  - 仅保留物种数 ≥ MIN_SPECIES_PER_OG 且包含参考物种的 OG。
"""

from __future__ import annotations

import os
import sys
import csv
import gzip
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set

# =========================
# 参数区（皇上可在此修改）
# =========================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
GFF_DIR = RAW_DATA_DIR / "gff"

# ★ 确保这里指向 Orthofinder Orthogroups 目录
ORTHOFINDER_DIR = RAW_DATA_DIR / "orthofinder"

SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"
CHR_RENAME_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"
OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_01_global_anchors"

ORTHOGROUPS_FILE = ORTHOFINDER_DIR / "Orthogroups.tsv"
ORTHOGROUPS_SC_FILE = ORTHOFINDER_DIR / "Orthogroups_SingleCopyOrthologues.txt"

USE_SINGLE_COPY_ONLY = True
MIN_SPECIES_PER_OG = 10
LOG_LEVEL = "INFO"

# =========================
# 工具函数
# =========================

def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny_01_build_global_anchors.log"

    logger = logging.getLogger("synteny_global_anchors")
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

    logger.info("========== synteny_01 — 全局锚点构建 ==========")
    logger.info("PROJECT_ROOT            = %s", PROJECT_ROOT)
    logger.info("RAW_DATA_DIR            = %s", RAW_DATA_DIR)
    logger.info("ORTHOFINDER_DIR         = %s", ORTHOFINDER_DIR)
    logger.info("CHR_RENAME_DIR          = %s", CHR_RENAME_DIR)
    logger.info("USE_SINGLE_COPY_ONLY    = %s", USE_SINGLE_COPY_ONLY)
    logger.info("MIN_SPECIES_PER_OG      = %d", MIN_SPECIES_PER_OG)

    return logger


def clean_output_root(output_root: Path, logger: logging.Logger) -> None:
    if output_root.exists():
        logger.info("删除旧输出目录：%s", output_root)
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)


def load_species_meta(meta_file: Path, logger: logging.Logger) -> tuple[list[str], str]:
    if not meta_file.exists():
        logger.error("species meta 文件不存在：%s", meta_file)
        sys.exit(1)

    species_ids: list[str] = []
    ref_species_id: str | None = None

    logger.info("读取物种 meta：%s", meta_file)
    with meta_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            logger.error("meta 表必须包含列：species_id")
            sys.exit(1)
        has_ref = "is_reference" in reader.fieldnames

        for row in reader:
            sid = (row.get("species_id") or "").strip()
            if not sid:
                continue
            species_ids.append(sid)
            if has_ref:
                flag = (row.get("is_reference") or "").strip().lower()
                if flag == "yes":
                    if ref_species_id is not None:
                        logger.error(
                            "meta 表中 is_reference=yes 超过一个：%s, %s",
                            ref_species_id,
                            sid,
                        )
                        sys.exit(1)
                    ref_species_id = sid

    if not species_ids:
        logger.error("meta 表中没有有效 species_id")
        sys.exit(1)

    if ref_species_id is None:
        logger.error("meta 表未指定 is_reference=yes 的参考物种")
        sys.exit(1)

    logger.info("synteny 物种数量 = %d", len(species_ids))
    logger.info("参考物种 = %s", ref_species_id)
    return species_ids, ref_species_id


def find_gff_for_species(species_id: str, logger: logging.Logger) -> Path | None:
    candidates = [
        GFF_DIR / f"{species_id}.gff3",
        GFF_DIR / f"{species_id}.gff",
        GFF_DIR / f"{species_id}.gff3.gz",
        GFF_DIR / f"{species_id}.gff.gz",
    ]
    for p in candidates:
        if p.exists():
            logger.info("物种 %s 使用 GFF：%s", species_id, p)
            return p

    try:
        for fname in os.listdir(GFF_DIR):
            for suffix in (".gff3", ".gff", ".gff3.gz", ".gff.gz"):
                if fname == f"{species_id}{suffix}":
                    path = GFF_DIR / fname
                    logger.info("物种 %s 使用 GFF（宽松匹配）：%s", species_id, path)
                    return path
    except FileNotFoundError:
        logger.error("GFF 目录不存在：%s", GFF_DIR)
        return None

    logger.error("未找到物种 %s 对应的 GFF 文件", species_id)
    return None


def parse_attributes_for_best_id(attr_str: str) -> str | None:
    """
    从 GFF 属性列选择“最适合作为转录本 ID”的字段。

    优先级：
      1) 所有 key.lower().endswith("transcript_id") 的字段（transcript_id 家族）：
         - 优先 transcript_id
         - 其次 orig_transcript_id
         - 再其次其它 *transcript_id（按 key 名排序保证稳定）
      2) 若无 transcript_id 家族，退回 ID=

    示例 Teg：
      ID=rna-Teg00024835-RA;orig_transcript_id=gnl|WGS:JARBDR|Teg00024835-RA;...
      -> 先命中 orig_transcript_id
      -> 返回 "gnl|WGS:JARBDR|Teg00024835-RA"
      -> 再经 extract_id_core() => "Teg00024835-RA"
    """
    if not attr_str:
        return None
    raw_fields = attr_str.split(";")
    kv: Dict[str, str] = {}
    for field in raw_fields:
        field = field.strip()
        if not field:
            continue
        if "=" not in field:
            continue
        k, v = field.split("=", 1)
        kv[k.strip()] = v.strip()

    # 找出所有 *transcript_id 键
    transcript_like_keys = [
        k for k in kv.keys()
        if k.lower().endswith("transcript_id")
    ]

    if transcript_like_keys:
        # 优先 transcript_id
        for key in ("transcript_id", "Transcript_id", "TRANSCRIPT_ID"):
            if key in kv and kv[key]:
                return kv[key]
        # 其次 orig_transcript_id
        for key in ("orig_transcript_id", "Orig_transcript_id", "ORIG_TRANSCRIPT_ID"):
            if key in kv and kv[key]:
                return kv[key]
        # 其它 *transcript_id，按 key 名排序保证稳定
        for key in sorted(transcript_like_keys):
            if kv.get(key):
                return kv[key]

    # 若没有任何 transcript_id 家族，则退回 ID
    if "ID" in kv and kv["ID"]:
        return kv["ID"]

    return None


def extract_id_core(full_id: str) -> str:
    """
    从完整 ID 中提取 id_core：
      - 若包含竖线 '|'，取最后一个 '|' 之后的部分；
      - 否则直接返回原始 ID。
    """
    if "|" in full_id:
        return full_id.split("|")[-1]
    return full_id


def load_chr_rename_for_species(
    species_id: str,
    logger: logging.Logger,
) -> Dict[str, Tuple[str | None, bool]]:
    path = CHR_RENAME_DIR / f"chr_rename_{species_id}.tsv"
    mapping: Dict[str, Tuple[str | None, bool]] = {}
    if not path.exists():
        logger.error("未找到 chr_rename 文件：%s", path)
        return mapping

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"seqid_raw", "is_chromosome", "new_chr_name"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            logger.error("chr_rename_%s.tsv 缺少列：%s", species_id, ",".join(sorted(missing)))
            return mapping
        for row in reader:
            seqid_raw = (row.get("seqid_raw") or "").strip()
            if not seqid_raw:
                continue
            is_chr = (row.get("is_chromosome") or "").strip().lower() == "yes"
            new_chr = (row.get("new_chr_name") or "").strip() or None
            mapping[seqid_raw] = (new_chr, is_chr)

    n_chr = sum(1 for v in mapping.values() if v[1] and v[0])
    logger.info("物种 %s：主染色体数量 = %d", species_id, n_chr)
    return mapping


def build_gene_mapping_for_species(
    species_id: str,
    chr_map: Dict[str, Tuple[str | None, bool]],
    logger: logging.Logger,
) -> Dict[str, Dict[str, object]]:
    gff_path = find_gff_for_species(species_id, logger)
    if gff_path is None:
        return {}

    is_gz = gff_path.suffix == ".gz"
    open_func = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    target_types = {"mrna", "transcript"}

    gene_map: Dict[str, Dict[str, object]] = {}
    n_mrna = 0
    n_main = 0

    with open_func(gff_path, mode, encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, ftype, start_str, end_str, strand, attr = (
                parts[0],
                parts[2].lower(),
                parts[3],
                parts[4],
                parts[6],
                parts[8],
            )
            if ftype not in target_types:
                continue
            n_mrna += 1
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue

            best_id = parse_attributes_for_best_id(attr)
            if not best_id:
                continue

            id_core = extract_id_core(best_id)

            if seqid not in chr_map:
                continue
            new_chr, is_chr = chr_map[seqid]
            if not (is_chr and new_chr):
                continue

            n_main += 1
            gene_map[id_core] = {
                "species_id": species_id,
                "gene_id": best_id,
                "id_core": id_core,
                "seqid_raw": seqid,
                "chr": new_chr,
                "start": start,
                "end": end,
                "strand": strand,
            }

    logger.info(
        "物种 %s：mRNA/transcript 条目 = %d，其中主染色体上的条目 = %d，id_core 数量 = %d",
        species_id,
        n_mrna,
        n_main,
        len(gene_map),
    )
    return gene_map


def load_single_copy_ogs(sc_file: Path, logger: logging.Logger) -> Set[str]:
    sc_set: Set[str] = set()
    if not sc_file.exists():
        logger.warning("未找到单拷贝 OG 列表：%s", sc_file)
        return sc_set

    logger.info("读取单拷贝 OG 列表：%s", sc_file)
    if sc_file.suffix.lower() == ".txt":
        with sc_file.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                og_id = line.split()[0]
                sc_set.add(og_id)
    else:
        with sc_file.open("r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if not header:
                return sc_set
            for row in reader:
                if not row:
                    continue
                og_id = row[0].strip()
                if og_id:
                    sc_set.add(og_id)

    logger.info("单拷贝 OG 数量 = %d", len(sc_set))
    return sc_set


def write_tsv(path: Path, fieldnames: List[str], records: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in records:
            writer.writerow(r)

# =========================
# 主流程
# =========================

def main() -> None:
    logger = setup_logging(OUTPUT_ROOT / "logs")
    clean_output_root(OUTPUT_ROOT, logger)

    species_ids, ref_species_id = load_species_meta(SPECIES_META_FILE, logger)
    species_set = set(species_ids)

    species_chr_map: Dict[str, Dict[str, Tuple[str | None, bool]]] = {}
    species_gene_map: Dict[str, Dict[str, Dict[str, object]]] = {}

    for sid in species_ids:
        logger.info("------ 预处理物种：%s ------", sid)
        chr_map = load_chr_rename_for_species(sid, logger)
        species_chr_map[sid] = chr_map
        gene_map = build_gene_mapping_for_species(sid, chr_map, logger)
        species_gene_map[sid] = gene_map

    sc_ogs = load_single_copy_ogs(ORTHOGROUPS_SC_FILE, logger)

    if not ORTHOGROUPS_FILE.exists():
        logger.error("未找到 Orthogroups.tsv：%s", ORTHOGROUPS_FILE)
        sys.exit(1)

    qc_stats: Dict[str, Dict[str, object]] = {}
    for sid in species_ids:
        qc_stats[sid] = {
            "species_id": sid,
            "n_proteins_in_OG": 0,
            "n_id_core_extracted": 0,
            "n_id_core_matched_GFF": 0,
            "n_id_core_unmatched_GFF": 0,
            "unmatched_examples": set(),
        }

    anchors_records: List[Dict[str, object]] = []

    n_ogs_total = 0
    n_ogs_after_sc_filter = 0
    n_ogs_after_species_filter = 0
    n_ogs_written = 0

    with ORTHOGROUPS_FILE.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            logger.error("Orthogroups.tsv 为空")
            sys.exit(1)

        species_cols = header[1:]
        header_species_set = set(species_cols)

        missing = species_set - header_species_set
        extra = header_species_set - species_set

        if missing:
            logger.error("Orthogroups.tsv 缺少 meta 物种：%s", ",".join(sorted(missing)))
            sys.exit(1)
        if extra:
            logger.warning("Orthogroups.tsv 中存在 meta 未用物种列：%s", ",".join(sorted(extra)))

        species_to_col = {}
        for i, col in enumerate(header[1:], start=1):
            if col in species_set:
                species_to_col[col] = i

        logger.info("Orthogroups.tsv 中有效物种列数（与 meta 一致） = %d", len(species_to_col))

        for row in reader:
            if not row:
                continue
            og_id = row[0].strip()
            if not og_id:
                continue
            n_ogs_total += 1

            is_sc = og_id in sc_ogs
            if USE_SINGLE_COPY_ONLY and not is_sc:
                continue
            n_ogs_after_sc_filter += 1

            og_member_records: List[Dict[str, object]] = []
            og_species_mapped: Set[str] = set()

            for sid in species_ids:
                col_idx = species_to_col.get(sid)
                if col_idx is None or col_idx >= len(row):
                    continue
                cell = row[col_idx].strip()
                if not cell:
                    continue
                raw_ids = [x.strip() for x in cell.replace(",", " ").split() if x.strip()]
                if not raw_ids:
                    continue

                gene_map = species_gene_map.get(sid, {})

                for prot_id_raw in raw_ids:
                    qc_stats[sid]["n_proteins_in_OG"] += 1
                    id_core = extract_id_core(prot_id_raw)
                    if not id_core:
                        continue
                    qc_stats[sid]["n_id_core_extracted"] += 1

                    gene_rec = gene_map.get(id_core)
                    if gene_rec is None:
                        qc_stats[sid]["n_id_core_unmatched_GFF"] += 1
                        if len(qc_stats[sid]["unmatched_examples"]) < 5:
                            qc_stats[sid]["unmatched_examples"].add(prot_id_raw)
                        continue

                    qc_stats[sid]["n_id_core_matched_GFF"] += 1
                    og_species_mapped.add(sid)
                    og_member_records.append(
                        {
                            "og_id": og_id,
                            "prot_id_raw": prot_id_raw,
                            "id_core": id_core,
                            "species_id": sid,
                            "gene_id": gene_rec["gene_id"],
                            "chr": gene_rec["chr"],
                            "start": int(gene_rec["start"]),
                            "end": int(gene_rec["end"]),
                            "strand": gene_rec["strand"],
                        }
                    )

            if len(og_species_mapped) < MIN_SPECIES_PER_OG:
                continue
            n_ogs_after_species_filter += 1

            has_ref = any(rec["species_id"] == ref_species_id for rec in og_member_records)
            if not has_ref:
                continue

            ref_records = [rec for rec in og_member_records if rec["species_id"] == ref_species_id]
            ref_records.sort(key=lambda x: (x["chr"], x["start"]))
            ref_rec = ref_records[0]
            ref_chr = ref_rec["chr"]
            ref_start = ref_rec["start"]
            ref_end = ref_rec["end"]

            is_core = "yes" if og_species_mapped == species_set else "no"
            is_sc_flag = "yes" if is_sc else "no"

            for rec in og_member_records:
                anchors_records.append(
                    {
                        "og_id": og_id,
                        "is_single_copy": is_sc_flag,
                        "is_core": is_core,
                        "ref_species_id": ref_species_id,
                        "ref_chr": ref_chr,
                        "ref_start": ref_start,
                        "ref_end": ref_end,
                        "species_id": rec["species_id"],
                        "id_core": rec["id_core"],
                        "prot_id_raw": rec["prot_id_raw"],
                        "gene_id": rec["gene_id"],
                        "chr": rec["chr"],
                        "start": rec["start"],
                        "end": rec["end"],
                        "strand": rec["strand"],
                    }
                )
            n_ogs_written += 1

    logger.info("Orthogroups.tsv 总 OG 数量            = %d", n_ogs_total)
    logger.info("单拷贝筛选后 OG 数量                 = %d", n_ogs_after_sc_filter)
    logger.info("通过物种数 (>=MIN_SPECIES_PER_OG) OG 数 = %d", n_ogs_after_species_filter)
    logger.info("最终写入 anchors_global 的 OG 数量    = %d", n_ogs_written)
    logger.info("anchors_global 记录行数（基因级）     = %d", len(anchors_records))

    anchors_path = OUTPUT_ROOT / "anchors_global.tsv"
    anchors_fields = [
        "og_id",
        "is_single_copy",
        "is_core",
        "ref_species_id",
        "ref_chr",
        "ref_start",
        "ref_end",
        "species_id",
        "id_core",
        "prot_id_raw",
        "gene_id",
        "chr",
        "start",
        "end",
        "strand",
    ]
    write_tsv(anchors_path, anchors_fields, anchors_records)
    logger.info("写出 anchors_global.tsv：%s", anchors_path)

    qc_records: List[Dict[str, object]] = []
    for sid in species_ids:
        stat = qc_stats[sid]
        n_prot = int(stat["n_proteins_in_OG"])
        n_core = int(stat["n_id_core_extracted"])
        n_match = int(stat["n_id_core_matched_GFF"])
        n_unmatch = int(stat["n_id_core_unmatched_GFF"])
        examples = stat["unmatched_examples"]
        match_rate = n_match / n_core if n_core > 0 else 0.0

        qc_records.append(
            {
                "species_id": sid,
                "n_proteins_in_OG": n_prot,
                "n_id_core_extracted": n_core,
                "n_id_core_matched_GFF": n_match,
                "n_id_core_unmatched_GFF": n_unmatch,
                "match_rate": f"{match_rate:.4f}",
                "example_unmatched_ids": ",".join(sorted(examples)),
            }
        )

    qc_path = OUTPUT_ROOT / "id_mapping_qc.tsv"
    qc_fields = [
        "species_id",
        "n_proteins_in_OG",
        "n_id_core_extracted",
        "n_id_core_matched_GFF",
        "n_id_core_unmatched_GFF",
        "match_rate",
        "example_unmatched_ids",
    ]
    write_tsv(qc_path, qc_fields, qc_records)
    logger.info("写出 ID 匹配 QC 表：%s", qc_path)

    summary_records = [
        {"metric": "n_ogs_total", "value": n_ogs_total},
        {"metric": "n_ogs_after_single_copy_filter", "value": n_ogs_after_sc_filter},
        {"metric": "n_ogs_after_species_filter", "value": n_ogs_after_species_filter},
        {"metric": "n_ogs_written", "value": n_ogs_written},
        {"metric": "n_anchor_records", "value": len(anchors_records)},
        {"metric": "use_single_copy_only", "value": str(USE_SINGLE_COPY_ONLY)},
        {"metric": "min_species_per_og", "value": MIN_SPECIES_PER_OG},
        {"metric": "ref_species_id", "value": ref_species_id},
    ]
    summary_path = OUTPUT_ROOT / "anchors_global_summary.tsv"
    write_tsv(summary_path, ["metric", "value"], summary_records)
    logger.info("写出 anchors_global_summary.tsv：%s", summary_path)

    logger.info("synteny_01 完成。")


if __name__ == "__main__":
    main()

