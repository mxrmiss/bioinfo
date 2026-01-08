#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny01_mcscan_catalog.py
—— MCscan 线性共线性前置：geneorder + BLAST catalog（16 物种 · 多线程版）

当前职责：
  1) 读取 raw_data/synteny_species_meta.tsv，识别：
       - 全部 species_id
       - 唯一参考物种 ref_species（is_reference == "yes"）
  2) 对所有物种生成 MCscan 风格的 geneorder BED：
       - 只保留 synteny_00_chr_rename 中 is_chromosome == "yes" 的 Chr*
       - 从 GFF 中提取基因坐标，写出：
           chr  start  end  gene_id
  3) 基于参考物种蛋白质序列构建 BLASTP 数据库；
  4) 对每个“非参考物种”并行执行 BLASTP：
       - query = 该物种 protein/*.faa
       - db    = 参考物种 protein DB
       - 输出：blast/<ref>__vs__<qry>.blast
  5) 全程多线程控制：
       - 同时最多处理 PAIR_PARALLELISM 对物种
       - 每对物种的 BLAST 使用 BLAST_THREADS_PER_PAIR 线程

说明：
  - 本脚本当前版本只负责 geneorder + BLAST catalog，
    不在此处调用 jcvi 的 mcscan / screen。
"""

from __future__ import annotations

import os
import sys
import csv
import gzip
import shutil
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed


# =========================
# 参数区（皇上可在此修改）
# =========================

# 项目根目录（默认：当前脚本所在目录的上一级 magic/）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 原始数据目录
RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
GFF_DIR = RAW_DATA_DIR / "gff"
PROTEIN_DIR = RAW_DATA_DIR / "protein"
SPECIES_META_FILE = RAW_DATA_DIR / "synteny_species_meta.tsv"

# synteny_00 的输出（染色体重命名表）
CHR_RENAME_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"

# 本脚本输出目录
OUTPUT_ROOT = PROJECT_ROOT / "output" / "synteny_01_mcscan_catalog"

# 是否执行 BLAST 步骤（一般保持 True）
RUN_BLAST = True

# 物种对的并行度控制：
#   - 同时最多处理多少对 ref vs qry
PAIR_PARALLELISM = 3

# 每对物种 BLASTP 使用的线程数
BLAST_THREADS_PER_PAIR = 10

# BLAST 参数
BLAST_PROGRAM = "blastp"
BLAST_EVALUE = 1e-5
BLAST_MAX_TARGET_SEQS = 5

# 日志等级
LOG_LEVEL = "INFO"

# GFF 第 9 列中用于提取“蛋白 ID 对应的转录本 ID”的字段列表
# 按顺序优先级依次匹配；皇上可自行在此追加字段名
TRANSCRIPT_ATTR_KEYS: List[str] = [
    "transcript_id",
    "orig_transcript_id",
]


# =========================
# 工具函数
# =========================

def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "synteny01_mcscan_catalog.log"

    logger = logging.getLogger("synteny01")
    logger.setLevel(getattr(logging, LOG_LEVEL.upper(), logging.INFO))

    # 清空旧 handler，避免重复
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

    logger.info("========== synteny01 — MCscan catalog 构建 ==========")
    logger.info("PROJECT_ROOT = %s", PROJECT_ROOT)
    logger.info("RAW_DATA_DIR = %s", RAW_DATA_DIR)
    logger.info("OUTPUT_ROOT  = %s", OUTPUT_ROOT)
    logger.info("RUN_BLAST = %s", RUN_BLAST)
    logger.info("PAIR_PARALLELISM = %d, BLAST_THREADS_PER_PAIR = %d",
                PAIR_PARALLELISM, BLAST_THREADS_PER_PAIR)

    return logger


def clean_output_root(output_root: Path) -> None:
    """
    删除旧的输出目录并重建。
    不使用 logger，方便在初始化日志前调用。
    """
    if output_root.exists():
        shutil.rmtree(output_root)
    (output_root / "geneorder").mkdir(parents=True, exist_ok=True)
    (output_root / "blast").mkdir(parents=True, exist_ok=True)
    (output_root / "logs").mkdir(parents=True, exist_ok=True)


def load_species_meta(meta_file: Path, logger: logging.Logger) -> Tuple[List[str], str]:
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


def find_gff_for_species(species_id: str, logger: logging.Logger) -> Path:
    candidates = [
        GFF_DIR / f"{species_id}.gff3",
        GFF_DIR / f"{species_id}.gff",
        GFF_DIR / f"{species_id}.gff3.gz",
        GFF_DIR / f"{species_id}.gff.gz",
    ]
    for p in candidates:
        if p.exists():
            logger.info("物种 %s 使用 GFF 文件：%s", species_id, p)
            return p

    logger.error("未能在 %s 下找到物种 %s 的 GFF 文件。", GFF_DIR, species_id)
    sys.exit(1)


def load_chr_rename(species_id: str, logger: logging.Logger) -> Dict[str, str]:
    path = CHR_RENAME_DIR / f"chr_rename_{species_id}.tsv"
    if not path.exists():
        logger.error("找不到 chr_rename 表：%s", path)
        sys.exit(1)

    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("is_chromosome") or "") != "yes":
                continue
            raw = (row.get("seqid_raw") or "").strip()
            new = (row.get("new_chr_name") or "").strip()
            if raw and new:
                mapping[raw] = new

    logger.info("物种 %s：从 chr_rename 中读取到 %d 条主染色体映射。",
                species_id, len(mapping))
    return mapping


def normalize_id_core(raw_id: str) -> str:
    """
    对原始 ID 做标准化：
      1) 去掉首尾空白与引号；
      2) 若包含 '|'，取最后一个 '|' 后面的片段；
      3) 返回处理后的字符串。
    """
    val = raw_id.strip().strip('"').strip()
    if "|" in val:
        val = val.rsplit("|", 1)[1]
    return val


def extract_id_from_attr(attr: str) -> Optional[str]:
    """
    按皇上要求的优先级从 GFF 第 9 列中提取 gene_id：

    1) 依次尝试 TRANSCRIPT_ATTR_KEYS 中的字段（如 transcript_id、orig_transcript_id）：
         - 若找到对应值，取该值中“从后往前第一个 '|' 后面的片段”为 ID；
    2) 若上述字段均未命中：
         - 回退到 ID= 或 Parent= 的值；
         - 若值以 'rna-' 开头，则去掉 'rna-' 前缀；
         - 再同样取“最后一个 '|' 后面的片段”为 ID。
    """
    # 先拆成若干 key/value
    kv_list: List[Tuple[str, str]] = []
    for field in attr.split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        else:
            parts = field.split()
            if len(parts) >= 2:
                k, v = parts[0], parts[-1]
            elif len(parts) == 1:
                k, v = parts[0], ""
            else:
                continue
        k = k.strip()
        v = v.strip().strip('"').strip()
        if k:
            kv_list.append((k, v))

    # 1) 优先从 TRANSCRIPT_ATTR_KEYS 中取值
    for key in TRANSCRIPT_ATTR_KEYS:
        for k, v in kv_list:
            if k == key and v:
                return normalize_id_core(v)

    # 2) 回退到 ID= 或 Parent=
    id_val: Optional[str] = None
    for k, v in kv_list:
        if k == "ID" and v:
            id_val = v
            break
    if id_val is None:
        for k, v in kv_list:
            if k == "Parent" and v:
                id_val = v
                break

    if not id_val:
        return None

    # 去掉 mRNA 前缀 rna-
    if id_val.startswith("rna-"):
        id_val = id_val[len("rna-") :]

    return normalize_id_core(id_val)


def parse_gff_gene_intervals(
    gff_path: Path,
    chr_map: Dict[str, str],
    logger: logging.Logger,
) -> List[Tuple[str, int, int, str]]:
    """
    从 GFF 中抽取基因坐标，返回：
      [(chr, start, end, gene_id), ...]
    仅保留在 chr_map 中出现的 seqid。

    当前版本只使用 mRNA / transcript 层的记录，
    不再直接使用 feature_type == "gene"，避免同一坐标重复写入 gene-LOC 与转录本 ID。
    """
    intervals: List[Tuple[str, int, int, str]] = []

    is_gz = gff_path.suffix == ".gz"
    open_func = gzip.open if is_gz else open
    mode = "rt" if is_gz else "r"

    with open_func(gff_path, mode, encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, feature_type, start_str, end_str, attr = (
                parts[0],
                parts[2],
                parts[3],
                parts[4],
                parts[8],
            )

            if seqid not in chr_map:
                continue

            # 只保留 mRNA / transcript 层，丢弃 gene 层，避免重复
            if feature_type not in {"mRNA", "transcript"}:
                continue

            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue

            gene_id = extract_id_from_attr(attr)
            if not gene_id:
                continue

            chr_name = chr_map[seqid]
            intervals.append((chr_name, start - 1, end, gene_id))

    logger.info("从 GFF %s 中解析出 %d 个 gene 间隔。", gff_path.name, len(intervals))
    return intervals


def write_geneorder_bed(
    species_id: str,
    intervals: List[Tuple[str, int, int, str]],
    out_dir: Path,
    logger: logging.Logger,
) -> Path:
    out_path = out_dir / f"{species_id}.bed"
    with out_path.open("w", encoding="utf-8") as f:
        for chr_name, start, end, gene_id in intervals:
            f.write(f"{chr_name}\t{start}\t{end}\t{gene_id}\n")

    logger.info("写出 geneorder bed：%s", out_path)
    return out_path


def run_cmd(cmd: List[str], logger: logging.Logger, cwd: Optional[Path] = None) -> int:
    cmd_str = " ".join(str(x) for x in cmd)
    logger.info("运行命令：%s", cmd_str)
    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
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
# BLAST 相关
# =========================

def build_blast_db(ref_species: str, logger: logging.Logger) -> Path:
    faa = PROTEIN_DIR / f"{ref_species}.faa"
    if not faa.exists():
        logger.error("找不到参考物种蛋白文件：%s", faa)
        sys.exit(1)

    db_prefix = OUTPUT_ROOT / "blast" / f"{ref_species}_db"
    cmd = [
        "makeblastdb",
        "-in", str(faa),
        "-dbtype", "prot",
        "-out", str(db_prefix),
    ]
    rc = run_cmd(cmd, logger)
    if rc != 0:
        logger.error("makeblastdb 失败，无法继续。")
        sys.exit(1)

    return db_prefix


def run_blast_for_pair(
    ref_species: str,
    qry_species: str,
    db_prefix: Path,
    logger: logging.Logger,
) -> Path:
    qry_faa = PROTEIN_DIR / f"{qry_species}.faa"
    if not qry_faa.exists():
        logger.error("找不到物种 %s 的蛋白文件：%s", qry_species, qry_faa)
        raise FileNotFoundError(str(qry_faa))

    out_blast = OUTPUT_ROOT / "blast" / f"{ref_species}__vs__{qry_species}.blast"

    cmd = [
        BLAST_PROGRAM,
        "-query", str(qry_faa),
        "-db", str(db_prefix),
        "-evalue", str(BLAST_EVALUE),
        "-max_target_seqs", str(BLAST_MAX_TARGET_SEQS),
        "-num_threads", str(BLAST_THREADS_PER_PAIR),
        "-outfmt", "6",
        "-out", str(out_blast),
    ]
    rc = run_cmd(cmd, logger)
    if rc != 0:
        raise RuntimeError(f"BLAST 失败：{ref_species} vs {qry_species}")

    return out_blast


# =========================
# 物种对处理（目前只做 BLAST）
# =========================

def process_pair_worker(
    ref_species: str,
    qry_species: str,
    db_prefix: Path,
    logger: logging.Logger,
) -> None:
    logger.info("====== 处理物种对：%s (ref) vs %s (qry) ======",
                ref_species, qry_species)

    if RUN_BLAST:
        _ = run_blast_for_pair(ref_species, qry_species, db_prefix, logger)


# =========================
# 主流程
# =========================

def main() -> None:
    # 先清空输出目录，再初始化日志，避免删掉刚建好的 log 文件
    clean_output_root(OUTPUT_ROOT)
    logger = setup_logging(OUTPUT_ROOT / "logs")

    species_ids, ref_species = load_species_meta(SPECIES_META_FILE, logger)

    # 1) 为所有物种生成 geneorder
    for sid in species_ids:
        logger.info("====== 生成 geneorder：%s ======", sid)
        chr_map = load_chr_rename(sid, logger)
        gff_path = find_gff_for_species(sid, logger)
        intervals = parse_gff_gene_intervals(gff_path, chr_map, logger)
        write_geneorder_bed(sid, intervals, OUTPUT_ROOT / "geneorder", logger)

    # 2) 基于参考物种构建 BLAST DB
    db_prefix = build_blast_db(ref_species, logger)

    # 3) 并行处理 ref vs 其它物种的 BLAST
    qry_species_list = [s for s in species_ids if s != ref_species]
    logger.info("准备处理 %d 个物种对，最大并行数 = %d",
                len(qry_species_list), PAIR_PARALLELISM)

    with ThreadPoolExecutor(max_workers=PAIR_PARALLELISM) as executor:
        future_to_pair = {}
        for qry in qry_species_list:
            future = executor.submit(
                process_pair_worker,
                ref_species,
                qry,
                db_prefix,
                logger,
            )
            future_to_pair[future] = qry

        for future in as_completed(future_to_pair):
            qry = future_to_pair[future]
            try:
                future.result()
            except Exception as e:
                logger.error("处理物种对 %s vs %s 时发生异常：%s",
                             ref_species, qry, str(e))

    logger.info("synteny01 完成。geneorder + BLAST catalog 已生成。")


if __name__ == "__main__":
    main()

