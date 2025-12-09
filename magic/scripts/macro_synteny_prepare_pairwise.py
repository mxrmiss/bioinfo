#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
macro_synteny_prepare_pairwise.py

步骤 1：为参考物种 vs 其他所有物种运行 JCVI ortholog，
产出 anchors / anchors.simple，并写出物种信息与 anchors 统计表。

本版要点（ID 双保险修正版）：
1. 按“项目蓝图”约定目录结构和输出表头。
2. FASTA / GFF 后缀全面兼容（.fa/.fasta/.fna[.gz] 等）。
3. 从蛋白 FASTA 中同时记录：
   - full_id：原始头（Sco01g06570.1）
   - norm_id：去掉末尾 .数字 的 ID（Sco01g06570）
   并建立 full <-> norm 的双向映射。
4. GFF 属性里只要能匹配到 full_id 或 norm_id，就认作同一个 gene，
   坐标统一挂在 norm_id 上。
5. 写 BED 时，对每个 gene 输出“所有相关 ID 一起写”：
   - 行 4 列既有 Sco01g06570.1，又有 Sco01g06570（坐标相同），
     这样无论 anchors 用哪种 ID，都能在 BED 里找到。
6. 即便某物种 BED 行数为 0，也会写出一个空 .bed 文件，避免 JCVI FileNotFoundError。
"""

from __future__ import annotations

import os
import sys
import gzip
import logging
import subprocess
from pathlib import Path
from typing import Dict, Set, List, Tuple, Optional

# ====================== 用户参数区（皇上控制台） ======================

# 参考物种前缀（必须和 raw_data/ 下文件前缀一致）
REFERENCE_SPECIES = "Sinonovacula_constricta"

# 染色体长度阈值（Mb）
CHR_LENGTH_THRESHOLD_MBP = 10

# 宏观窗口大小（Mb），本脚本只记录到 species_info.tsv，真正用在 build_blocks 里
WINDOW_SIZE_MBP = 5

# 每条参考 chr × 每个物种最多保留的 links 数量（本脚本只记录，不实际裁剪）
MAX_LINKS_PER_CHR = 2000

# 挖掘组物种
GROUP_DIGGING = [
    "Sinonovacula_constricta",
    "Sinonovacula_rivularis",
    "Novaculina_chinensis",
    "Tegillarca_granosa",
    "Meretrix_meretrix",
    "Mya_arenaria",
    "Panopea_generosa",
    "Mercenaria_mercenaria",
]

# 非挖掘组物种
GROUP_NON_DIGGING = [
    "Crassostrea_gigas",
    "Ostrea_edulis",
    "Pecten_maximus",
    "Argopecten_irradians",
    "Mytilus_galloprovincialis",
    "Mytilus_edulis",
    "Pinctada_fucata",
    "Ctenoides_ales",
]

# synteny 的线程数说明（目前 jcvi 内部 lastal 默认 -P 32，这里仅作为记录/将来扩展用）
SYNTENY_THREADS = 32

# 日志等级
LOG_LEVEL = logging.INFO

# ====================== 路径约定 ======================

# 当前脚本：project/magic/scripts/macro_synteny_prepare_pairwise.py
PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_GENOME_DIR = PROJECT_ROOT / "raw_data" / "genome"
RAW_GFF_DIR = PROJECT_ROOT / "raw_data" / "gff"
RAW_PEP_DIR = PROJECT_ROOT / "raw_data" / "protein"

OUT_ROOT = PROJECT_ROOT / "output" / "macro_synteny"
LOG_DIR = OUT_ROOT / "logs"
SYNTENY_PAIRWISE_DIR = OUT_ROOT / "synteny_pairwise"
META_DIR = OUT_ROOT / "meta"

# JCVI 工作目录：每对物种一个子目录
JCVI_WORK_DIR = SYNTENY_PAIRWISE_DIR / "work"


# ====================== 通用工具函数 ======================

def ensure_dir(path: Path) -> None:
    """确保目录存在。"""
    path.mkdir(parents=True, exist_ok=True)


def open_maybe_gzip(path: Path, mode: str = "rt"):
    """透明打开 .gz 或普通文本。"""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def setup_logger() -> None:
    """配置日志：同时写文件和标准输出。"""
    ensure_dir(LOG_DIR)
    log_file = LOG_DIR / "macro_synteny_prepare_pairwise.log"

    root = logging.getLogger()
    root.setLevel(LOG_LEVEL)

    # 清理旧 handler，避免重复输出
    for h in list(root.handlers):
        root.removeHandler(h)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(LOG_LEVEL)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(LOG_LEVEL)
    ch.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))

    root.addHandler(fh)
    root.addHandler(ch)


def run_cmd(cmd: List[str], cwd: Optional[Path] = None, desc: Optional[str] = None) -> None:
    """
    在子进程中执行命令。
    - cmd: 命令及参数列表
    - cwd: 工作目录
    """
    msg = " ".join(cmd)
    if desc:
        logging.info("[CMD] %s: %s", desc, msg)
    else:
        logging.info("[CMD] %s", msg)

    try:
        subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        logging.error("命令执行失败：%s，退出码=%s", msg, e.returncode)
        raise


# ====================== 物种标签与分组 ======================

def infer_species_label(prefix: str) -> str:
    """
    prefix → 物种标签（S. constricta）。
    """
    parts = prefix.split("_")
    if len(parts) >= 2:
        g, s = parts[0], parts[1]
        if g:
            return f"{g[0]}. {s}"
    return prefix


def infer_group(prefix: str) -> str:
    """根据 GROUP_DIGGING / GROUP_NON_DIGGING 判定分组标签。"""
    if prefix in GROUP_DIGGING:
        return "digging"
    if prefix in GROUP_NON_DIGGING:
        return "nondigging"
    return "unknown"


# ====================== FASTA / GFF 扫描与统计 ======================

GENOME_EXTS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")
GFF_EXTS = (".gff3", ".gff", ".gff3.gz", ".gff.gz")
PEP_EXTS = (
    ".pep.fa", ".pep.fa.gz",
    ".faa", ".faa.gz",
    ".fa", ".fa.gz",
    ".fasta", ".fasta.gz",
)


def scan_prefixes_in_dir(root: Path, exts: Tuple[str, ...]) -> Set[str]:
    """扫描某个目录中具有指定后缀的前缀集合。"""
    prefixes: Set[str] = set()
    if not root.exists():
        return prefixes
    for p in root.iterdir():
        if not p.is_file():
            continue
        for ext in exts:
            if p.name.endswith(ext):
                prefixes.add(p.name[: -len(ext)])
                break
    return prefixes


def choose_genome_path(prefix: str) -> Path:
    """根据 prefix 尝试各种 genome 后缀。"""
    for ext in GENOME_EXTS:
        p = RAW_GENOME_DIR / f"{prefix}{ext}"
        if p.exists():
            return p
    raise FileNotFoundError(f"找不到 genome：{prefix} ({RAW_GENOME_DIR})")


def choose_gff_path(prefix: str) -> Path:
    """根据 prefix 尝试各种 GFF 后缀。"""
    for ext in GFF_EXTS:
        p = RAW_GFF_DIR / f"{prefix}{ext}"
        if p.exists():
            return p
    raise FileNotFoundError(f"找不到 GFF：{prefix} ({RAW_GFF_DIR})")


def choose_pep_path(prefix: str) -> Path:
    """根据 prefix 尝试各种蛋白 FASTA 后缀。"""
    for ext in PEP_EXTS:
        p = RAW_PEP_DIR / f"{prefix}{ext}"
        if p.exists():
            return p
    raise FileNotFoundError(f"找不到蛋白 FASTA：{prefix} ({RAW_PEP_DIR})")


def compute_genome_stats(genome_fa: Path) -> Tuple[int, int, int, int, List[str], List[int]]:
    """
    统计 genome 信息：
    - 总长度
    - seqid 总数
    - 染色体数（长度>=阈值）
    - scaffold 数
    - 染色体 ID 列表 & 对应长度
    """
    seq_lengths: Dict[str, int] = {}
    cur_id: Optional[str] = None

    with open_maybe_gzip(genome_fa, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                cur_id = line[1:].strip().split()[0]
                if cur_id not in seq_lengths:
                    seq_lengths[cur_id] = 0
            else:
                if cur_id is None:
                    continue
                seq_lengths[cur_id] += len(line.strip())

    total_len = sum(seq_lengths.values())
    n_seq_total = len(seq_lengths)

    chr_threshold = int(CHR_LENGTH_THRESHOLD_MBP * 1e6)
    chr_ids: List[str] = []
    chr_lengths: List[int] = []

    for sid, slen in seq_lengths.items():
        if slen >= chr_threshold:
            chr_ids.append(sid)
            chr_lengths.append(slen)

    n_chr = len(chr_ids)
    n_scaffold = n_seq_total - n_chr

    chr_pairs = sorted(zip(chr_ids, chr_lengths), key=lambda x: x[1], reverse=True)
    chr_ids_sorted = [cid for cid, _ in chr_pairs]
    chr_lengths_sorted = [clen for _, clen in chr_pairs]

    return total_len, n_seq_total, n_chr, n_scaffold, chr_ids_sorted, chr_lengths_sorted


def build_species_info(prefixes: List[str]) -> None:
    """统计所有物种 genome，写 species_info.tsv。"""
    ensure_dir(META_DIR)
    out_path = META_DIR / "species_info.tsv"

    logging.info("开始统计 genome 信息，共 %d 个物种。", len(prefixes))
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(
            "\t".join(
                [
                    "species_prefix",
                    "species_label",
                    "group",
                    "genome_size_bp",
                    "n_seq_total",
                    "n_chr",
                    "n_scaffold",
                    "chr_ids",
                    "chr_lengths_bp",
                    "notes",
                ]
            )
            + "\n"
        )

        for sp in prefixes:
            genome_fa = choose_genome_path(sp)
            logging.info("[GENOME] 统计物种 %s genome：%s", sp, genome_fa.name)
            total_len, n_seq, n_chr, n_scaffold, chr_ids, chr_lens = compute_genome_stats(genome_fa)

            species_label = infer_species_label(sp)
            group = infer_group(sp)
            chr_ids_str = ",".join(chr_ids)
            chr_lens_str = ",".join(str(x) for x in chr_lens)

            out.write(
                "\t".join(
                    [
                        sp,
                        species_label,
                        group,
                        str(total_len),
                        str(n_seq),
                        str(n_chr),
                        str(n_scaffold),
                        chr_ids_str,
                        chr_lens_str,
                        "",
                    ]
                )
                + "\n"
            )

    logging.info("写出物种信息表：%s", out_path)


# ====================== GFF / PEP → BED ======================

def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    """解析 GFF 第 9 列属性。"""
    attrs: Dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        if not part:
            continue
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key:
            attrs[key] = value
    return attrs


def normalize_version_suffix(raw_id: str) -> str:
    """
    更严格的版本号去除策略：
    - 仅当 ID 中恰好只有一个 '.' 且 '.' 之后全部是数字时，才去掉该尾部版本号。
      例如：
        XP_012345678.2 -> XP_012345678
        geneX.v1 -> geneX.v1  (不去)
        evm.model.HiC_scaffold_1.9 -> 保留（有多个点）
    - 这样避免把像 evm.model.HiC_scaffold_1.9 之类的结构化 ID 坍缩到 evm.model.HiC_scaffold_1
    """
    rid = raw_id.strip()
    if not rid:
        return rid
    # 只有恰好一个点时，可能是 NCBl-like 版本号
    if rid.count('.') == 1:
        left, right = rid.rsplit('.', 1)
        if right.isdigit():
            return left
    return rid


def cleanup_attr_id_for_match(raw: str) -> str:
    """
    对 GFF 属性值做最小清洗，让它更像 FASTA 里的 ID。

    主要处理 NCBI WGS 这种：
    - 'gnl|WGS:JARBDR|Teg00000001-RA:cds'
      -> 取最后一段 'Teg00000001-RA:cds'
      -> 再截到 ':' 之前，得到 'Teg00000001-RA'
    """
    val = raw.strip()
    if not val:
        return val

    if "gnl|WGS" in val and "|" in val:
        last = val.split("|")[-1]
        if ":" in last:
            last = last.split(":", 1)[0]
        val = last.strip()

    return val


def collect_ids_from_pep(pep_fasta: Path) -> Tuple[Set[str], Dict[str, str], Dict[str, Set[str]]]:
    """
    从蛋白 FASTA 中收集 ID 信息：

    返回：
    - pep_full_ids：所有 full_id（原始头第一个字段）
    - full_to_norm：full_id -> norm_id（去掉末尾 .数字）
    - norm_to_fulls：norm_id -> {所有对应的 full_id}
    """
    pep_full_ids: Set[str] = set()
    full_to_norm: Dict[str, str] = {}
    norm_to_fulls: Dict[str, Set[str]] = {}
    n_header = 0

    with open_maybe_gzip(pep_fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            n_header += 1
            header = line[1:].strip()
            if not header:
                continue
            full_id = header.split()[0]
            if not full_id:
                continue
            norm_id = normalize_version_suffix(full_id)

            pep_full_ids.add(full_id)
            full_to_norm[full_id] = norm_id
            norm_to_fulls.setdefault(norm_id, set()).add(full_id)

    logging.info(
        "[PEP] %s: 头行总数=%d, unique full_id=%d, unique norm_id=%d",
        pep_fasta.name,
        n_header,
        len(pep_full_ids),
        len(norm_to_fulls),
    )
    return pep_full_ids, full_to_norm, norm_to_fulls


def build_bed_from_gff_and_pep(
    prefix: str,
    gff_path: Path,
    pep_path: Path,
    bed_path: Path,
) -> Tuple[int, int]:
    """
    用蛋白 FASTA ID 集合作为约束，从 GFF 中构建 BED。

    逻辑：
    1. FASTA 中每个 full_id 归并到一个 norm_id（去版本号）。
    2. GFF 属性里若能匹配到 full_id 或 norm_id，就认为是这个 gene。
    3. coord 层面统一挂在 norm_id 上（合并所有 exon/CDS）。
    4. 写 BED 时，对该 norm_id 关联的所有 full_id 和 norm_id 本身都写一行 BED：
       -> 这样 anchors 用 full_id 或 norm_id，都能在 BED 里找到坐标。
    """
    logging.info(
        "[BED] 构建 BED：物种=%s, gff=%s, pep(用于BLAST)=%s, bed=%s",
        prefix,
        gff_path.name,
        pep_path.name,
        bed_path,
    )

    pep_full_ids, full_to_norm, norm_to_fulls = collect_ids_from_pep(pep_path)
    if not pep_full_ids:
        logging.warning("[BED] 物种=%s: 蛋白 FASTA 为空，写出空 BED：%s", prefix, bed_path)
        with open(bed_path, "w", encoding="utf-8") as out:
            pass
        return 0, 0

    # 收集每个 norm_id 的坐标：norm_id -> [chrom, start_min, end_max, strand]
    loci: Dict[str, List[object]] = {}

    with open_maybe_gzip(gff_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start_s, end_s, score, strand, phase, attr = parts

            attrs = parse_gff_attributes(attr)

            raw_candidates: List[str] = []
            if "protein_id" in attrs:
                raw_candidates.append(attrs["protein_id"])
            if "orig_protein_id" in attrs:
                raw_candidates.append(attrs["orig_protein_id"])
            for key in ("ID", "Parent", "gene_id", "Name"):
                val = attrs.get(key)
                if val:
                    raw_candidates.append(val)

            chosen_norm: Optional[str] = None

            for raw_val in raw_candidates:
                clean_val = cleanup_attr_id_for_match(raw_val)

                # 1) 直接看是否是 full_id
                if clean_val in pep_full_ids:
                    chosen_norm = full_to_norm[clean_val]
                    break

                # 2) 再看是否能归约成某个 norm_id
                nv = normalize_version_suffix(clean_val)
                if nv in norm_to_fulls:
                    chosen_norm = nv
                    break

            if chosen_norm is None:
                continue

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue

            if strand not in ("+", "-"):
                strand = "."

            if chosen_norm not in loci:
                loci[chosen_norm] = [chrom, start, end, strand]
            else:
                chrom0, s0, e0, strand0 = loci[chosen_norm]
                s0 = min(s0, start)
                e0 = max(e0, end)
                if strand0 == "." and strand in ("+", "-"):
                    strand0 = strand
                loci[chosen_norm] = [chrom0, s0, e0, strand0]

    # 写出 BED：对每个 norm_id，把所有 full_id + norm_id 自己都写出去
    written_lines = 0
    unique_genes = len(loci)

    with open(bed_path, "w", encoding="utf-8") as out:
        for norm_id, (chrom, start, end, strand) in loci.items():
            alias_ids: Set[str] = set()

            # 先挂上所有 full_id
            if norm_id in norm_to_fulls:
                alias_ids.update(norm_to_fulls[norm_id])

            # 再把 norm_id 自己也加进去（避免只有 full_id 没有 norm 的情况）
            alias_ids.add(norm_id)

            for gid in sorted(alias_ids):
                out.write(
                    "\t".join(
                        [
                            str(chrom),
                            str(start),
                            str(end),
                            gid,  # 可以是 full_id 也可以是 norm_id，两者全有
                            "0",
                            strand,
                        ]
                    )
                    + "\n"
                )
                written_lines += 1

    coverage = unique_genes / len(pep_full_ids) if pep_full_ids else 0.0
    logging.info(
        "[BED] 物种=%s: BED 基因数=%d, BED 行数=%d, 蛋白 full_id 总数=%d, 覆盖率(按gene) = %.3f",
        prefix,
        unique_genes,
        written_lines,
        len(pep_full_ids),
        coverage,
    )
    if coverage < 0.80:
        logging.warning(
            "[BED] 物种=%s: BED 覆盖率较低 (%.3f < 0.80)，后续 synteny 可能偏弱，请核查 GFF/FASTA 是否完全匹配。",
            prefix,
            coverage,
        )

    return written_lines, len(pep_full_ids)


def symlink_pep_in_workdir(prefix: str, pep_src: Path, workdir: Path) -> Path:
    """在 JCVI 工作目录下创建 pep 软链接，命名为 <prefix>.pep。"""
    ensure_dir(workdir)
    pep_link = workdir / f"{prefix}.pep"
    if pep_link.exists() or pep_link.is_symlink():
        return pep_link
    rel = os.path.relpath(pep_src, workdir)
    logging.info("[LINK] 物种=%s: 建立蛋白软链接 %s -> %s", prefix, pep_link, rel)
    os.symlink(rel, pep_link)
    return pep_link


def build_bed_for_species_in_workdir(prefix: str, workdir: Path) -> Path:
    """在指定 workdir 内为某个物种构建 BED；若已有则复用。"""
    bed_path = workdir / f"{prefix}.bed"
    if bed_path.exists():
        logging.info("[BED] 物种=%s: 检测到已有 BED，跳过重建：%s", prefix, bed_path)
        return bed_path

    gff_path = choose_gff_path(prefix)
    pep_path = choose_pep_path(prefix)
    build_bed_from_gff_and_pep(prefix, gff_path, pep_path, bed_path)
    return bed_path


# ====================== JCVI ortholog 执行与统计 ======================

def count_lines(path: Path) -> int:
    """统计文本文件行数，若文件不存在返回 0。"""
    if not path.is_file():
        return 0
    n = 0
    with open(path, "r", encoding="utf-8") as fh:
        for _ in fh:
            n += 1
    return n


def run_ortholog_for_pair(ref_prefix: str, sp_prefix: str) -> Tuple[str, int, int, int, str]:
    """
    为 ref vs sp 运行 jcvi.compara.catalog ortholog。
    返回：(sp_prefix, n_raw_pairs, n_filtered_pairs, n_pairs_ref_chr, comment)
    """
    pair_name = f"{ref_prefix}_vs_{sp_prefix}"
    workdir = JCVI_WORK_DIR / pair_name
    ensure_dir(workdir)

    logging.info("Preparing synteny: %s vs %s ...", ref_prefix, sp_prefix)
    logging.info("[PAIR] 开始处理物种对：%s vs %s", ref_prefix, sp_prefix)

    anchors_raw_path = SYNTENY_PAIRWISE_DIR / f"{ref_prefix}_vs_{sp_prefix}.anchors"
    anchors_simple_path = SYNTENY_PAIRWISE_DIR / f"{ref_prefix}_vs_{sp_prefix}.anchors.simple"

    # 建立 pep 软链接
    ref_pep_src = choose_pep_path(ref_prefix)
    sp_pep_src = choose_pep_path(sp_prefix)
    symlink_pep_in_workdir(ref_prefix, ref_pep_src, workdir)
    symlink_pep_in_workdir(sp_prefix, sp_pep_src, workdir)

    # 构建 BED
    build_bed_for_species_in_workdir(ref_prefix, workdir)
    build_bed_for_species_in_workdir(sp_prefix, workdir)

    # 已有 anchors.simple 则跳过
    if anchors_simple_path.exists():
        logging.info("[PAIR] %s: anchors.simple 已存在，跳过 ortholog。", pair_name)
        n_raw = count_lines(anchors_raw_path)
        n_simple = count_lines(anchors_simple_path)
        return sp_prefix, n_raw, n_simple, 0, "SKIP_RERUN"

    # 运行 ortholog
    cmd = [
        sys.executable,
        "-m",
        "jcvi.compara.catalog",
        "ortholog",
        ref_prefix,
        sp_prefix,
        "--dbtype",
        "prot",
    ]
    try:
        run_cmd(cmd, cwd=workdir, desc=f"ortholog:{pair_name}")
    except Exception as e:
        logging.error(
            "[ERROR] ortholog 失败：%s vs %s，命令执行失败（ortholog:%s），原因=%s",
            ref_prefix,
            sp_prefix,
            pair_name,
            e,
        )
        return sp_prefix, 0, 0, 0, "FAILED_ORTHOLOG"

    # 挪 anchors / anchors.simple
    anchors_name = f"{ref_prefix}.{sp_prefix}.anchors"
    anchors_simple_name = f"{anchors_name}.simple"
    src_anchors = workdir / anchors_name
    src_anchors_simple = workdir / anchors_simple_name

    ensure_dir(SYNTENY_PAIRWISE_DIR)

    if src_anchors.exists():
        if anchors_raw_path.exists():
            anchors_raw_path.unlink()
        src_anchors.rename(anchors_raw_path)

    if src_anchors_simple.exists():
        if anchors_simple_path.exists():
            anchors_simple_path.unlink()
        src_anchors_simple.rename(anchors_simple_path)

    n_raw_pairs = count_lines(anchors_raw_path)
    n_filtered_pairs = count_lines(anchors_simple_path)
    n_pairs_ref_chr = 0  # 细分留给后续步骤

    if n_filtered_pairs == 0:
        comment = "NO_SYNTENY"
    elif n_filtered_pairs < 50:
        comment = "LOW_SYNTENY"
    else:
        comment = "OK"

    logging.info(
        "[PAIR] %s: n_raw_pairs=%d, n_filtered_pairs=%d, comment=%s",
        pair_name,
        n_raw_pairs,
        n_filtered_pairs,
        comment,
    )

    return sp_prefix, n_raw_pairs, n_filtered_pairs, n_pairs_ref_chr, comment


def write_anchors_summary(ref_prefix: str, records: List[Tuple[str, int, int, int, str]]) -> None:
    """写出 anchors_summary.tsv。"""
    ensure_dir(META_DIR)
    out_path = META_DIR / "anchors_summary.tsv"
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(
            "\t".join(
                [
                    "ref_prefix",
                    "sp_prefix",
                    "n_raw_pairs",
                    "n_filtered_pairs",
                    "n_pairs_ref_chr",
                    "comment",
                ]
            )
            + "\n"
        )
        for sp_prefix, n_raw, n_filt, n_ref_chr, comment in records:
            out.write(
                "\t".join(
                    [
                        ref_prefix,
                        sp_prefix,
                        str(n_raw),
                        str(n_filt),
                        str(n_ref_chr),
                        comment,
                    ]
                )
                + "\n"
            )
    logging.info("写出 anchors 统计表：%s", out_path)


# ====================== 主流程 ======================

def main() -> None:
    setup_logger()

    logging.info("=== 宏观共线性步骤1：物种扫描 & JCVI pairwise ortholog 开始 ===")
    logging.info("项目根目录：%s", PROJECT_ROOT)
    logging.info("参考物种：%s", REFERENCE_SPECIES)
    logging.info(
        "参数：CHR_LENGTH_THRESHOLD_MBP=%s, WINDOW_SIZE_MBP=%s, MAX_LINKS_PER_CHR=%s",
        CHR_LENGTH_THRESHOLD_MBP,
        WINDOW_SIZE_MBP,
        MAX_LINKS_PER_CHR,
    )

    ensure_dir(OUT_ROOT)
    ensure_dir(SYNTENY_PAIRWISE_DIR)
    ensure_dir(JCVI_WORK_DIR)
    ensure_dir(META_DIR)

    genome_prefixes = scan_prefixes_in_dir(RAW_GENOME_DIR, GENOME_EXTS)
    gff_prefixes = scan_prefixes_in_dir(RAW_GFF_DIR, GFF_EXTS)
    pep_prefixes = scan_prefixes_in_dir(RAW_PEP_DIR, PEP_EXTS)

    logging.info("genome 目录物种数=%d：%s", len(genome_prefixes), ", ".join(sorted(genome_prefixes)))
    logging.info("gff    目录物种数=%d：%s", len(gff_prefixes), ", ".join(sorted(gff_prefixes)))
    logging.info("protein目录物种数=%d：%s", len(pep_prefixes), ", ".join(sorted(pep_prefixes)))

    common_prefixes = sorted(genome_prefixes & gff_prefixes & pep_prefixes)
    logging.info("三者交集物种数=%d（将参与 synteny）：%s", len(common_prefixes), ", ".join(common_prefixes))

    if REFERENCE_SPECIES not in common_prefixes:
        logging.error("参考物种 %s 不在三者交集物种列表中，流程终止。", REFERENCE_SPECIES)
        sys.exit(1)

    # 写 species_info.tsv（所有交集物种）
    build_species_info(common_prefixes)

    # 逐物种跑 synteny（参考物种 vs 其他）
    records: List[Tuple[str, int, int, int, str]] = []
    n_fail = 0
    total_sp = len(common_prefixes)

    for idx, sp in enumerate(common_prefixes, start=1):
        if sp == REFERENCE_SPECIES:
            continue
        logging.info(
            "Preparing synteny: %d/%d %s vs %s ...",
            idx,
            total_sp,
            REFERENCE_SPECIES,
            sp,
        )
        result = run_ortholog_for_pair(REFERENCE_SPECIES, sp)
        records.append(result)
        if result[-1].startswith("FAILED"):
            n_fail += 1

    write_anchors_summary(REFERENCE_SPECIES, records)

    if n_fail > 0:
        logging.error(
            "=== 宏观共线性步骤1 结束：共有 %d 个物种对 synteny 失败，请查看日志 & anchors_summary.tsv ===",
            n_fail,
        )
        sys.exit(1)
    else:
        logging.info("=== 宏观共线性步骤1 完成：所有物种对 synteny 成功 ===")


if __name__ == "__main__":
    main()

