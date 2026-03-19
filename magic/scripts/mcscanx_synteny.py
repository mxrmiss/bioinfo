#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAGIC / scripts/synteny_prep_mcscanx.py
======================================
功能：利用 JCVI（Python 版 MCscan）做两物种染色体共线性预处理，生成圈图所需表格。

目录假定（magic 工程）：
  project/
    magic/
      raw_data/      # 放原始 gff / protein / 可选 genome
      input/         # 本脚本输出的 TSV
      output/
        synteny/     # 后续 R 画圈图输出
        temp_mcscanx/# JCVI 中间文件（包含 temp_blast 等）

核心步骤：
  1. 扫描 raw_data/，识别物种及其 gff / protein（genome 可选）
  2. 用 GFF 解析基因坐标，生成 bed，按“Top N 染色体”过滤
  3. 清洗蛋白质 ID（移除物种前缀，保留和 GFF 一致的 ID）
  4. 调用 JCVI：
       python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names ...
       python -m jcvi.compara.synteny screen --simple ...
  5. 解析 anchors.simple → block 坐标表与汇总表（full/main 统一基于 anchors.simple）

注意：
  * 本脚本仅支持 2 个物种；若 raw_data 有多物种，请在 CONFIG 中显式指定 SPECIES_PAIR。
  * JCVI 环境需已安装：`pip/conda install jcvi`，以及 LAST 或 BLAST 等依赖。
"""

import os
import sys
import gzip
import subprocess
import multiprocessing
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

# =============================================================================
# 1. 用户配置区（所有参数都在这里改）
# =============================================================================

# 1.1 路径设置（按 magic 工程默认布局）
DIR_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_RAW = os.path.join(DIR_ROOT, "raw_data")
DIR_IN = os.path.join(DIR_ROOT, "input")

# 统一 output 根目录，所有输出与中间件集中到这里
DIR_OUTPUT_ROOT = os.path.join(DIR_ROOT, "output")
DIR_OUT = os.path.join(DIR_OUTPUT_ROOT, "synteny")        # 目前不直接写入图，预留
TEMP_DIR = os.path.join(DIR_OUTPUT_ROOT, "temp_mcscanx")  # JCVI 中间文件（含 temp_blast）

# 1.2 物种选择
# 如果 raw_data 中正好只有 2 个物种，可保持为 None 自动识别；
# 若有多个物种，请改成例如：["S_constricta", "S_rivularis"]
SPECIES_PAIR = None  # 例如 ["S_constricta", "S_rivularis"]

# 1.3 染色体过滤策略
#   MODE = "top_n" ：每个物种按“基因覆盖长度”排序，保留前 N 条（推荐，去掉碎 scaffolds）
#   MODE = "none"  ：不过滤，保留 GFF 中所有染色体 / scaffolds
CHR_FILTER_MODE = "top_n"  # "top_n" 或 "none"

# 全局默认 Top N（大部分物种用这个）
DEFAULT_TOP_N = 19

# 针对特定物种的例外（键是物种前缀 = 文件名前缀）
# 若物种名包含这些 key（子串匹配），则用指定值覆盖 DEFAULT_TOP_N
SPECIES_SPECIFIC_TOPN = {
    # "S_constricta": 19,
    # "S_rivularis": 19,
}

# 1.4 GFF 解析策略
# 自动探测 feature type：优先 mRNA，其次 transcript，最后 gene
GFF_DETECT_LINES = 5000  # 探测多少行用于判断类型

# 1.5 蛋白 ID 清洗
# 蛋白质 header 里是否带 "物种名|" 的前缀，例如 "S_constricta|Sco01g00010.1"
# 若为 True，则自动保留 '|' 后面的部分作为真正 ID
PROTEIN_STRIP_SPECIES_PREFIX = True

# 1.6 JCVI / 比对参数
# 只做蛋白模式（顶刊通用做法）
JCV_ORTHO_DBTYPE = "prot"  # JCVI 的 --dbtype 参数
JCV_ORTHO_CCPU = 1         # jcvi.compara.catalog ortholog 的 --cpu（很多 LAST 仅支持单线程，保险用 1）
JCV_ORTHO_CSCORE = 0.7    # --cscore，用于 RBH 过滤（常见写法 0.7~0.99）

# anchors.simple 的过滤
ANCHOR_MINSPAN = 5  # --minspan，决定 block 至少跨多少锚点

# 1.7 线程数（本脚本自身用）
CPU_THREADS = max(1, multiprocessing.cpu_count() - 1)

# =============================================================================
# 2. 工具函数
# =============================================================================

def smart_open(path, mode="rt"):
    """自动识别是否为 .gz 压缩文件并打开"""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def log(msg: str):
    """简单打印日志，统一前缀"""
    print(msg, flush=True)


def check_dependencies():
    """检查 Python 库与 JCVI 是否可用"""
    log(f">> [Init] 检查环境（可用 CPU 核心：{multiprocessing.cpu_count()}）")
    missing_py = []
    try:
        import Bio  # noqa: F401
    except ImportError:
        missing_py.append("biopython")
    try:
        import pandas  # noqa: F401
    except ImportError:
        missing_py.append("pandas")
    try:
        import jcvi  # noqa: F401
    except ImportError:
        missing_py.append("jcvi")

    if missing_py:
        log(f"[错误] 缺少 Python 库：{', '.join(missing_py)}")
        log("       请先在 magic 环境中安装上述依赖后再运行。")
        sys.exit(1)

    # JCVI 内部会再去调用 LAST / BLAST 等，这里不再一层层检查，交给 jcvi 自己报错
    log(">> [Init] 环境检查通过。")


def classify_file_type(filename: str) -> str:
    """
    根据文件名/后缀粗略判断类型：
      - 返回 'gff' / 'prot' / 'genome' / 'unknown'
    """
    lower = filename.lower()
    # gff / gff3 / gtf
    if lower.endswith((".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")):
        return "gff"
    # 明确的蛋白质后缀
    if lower.endswith((".pep", ".pep.gz", ".faa", ".faa.gz", ".prot.fa", ".prot.fasta", ".protein.fa", ".protein.fasta")):
        return "prot"
    # 文件名中包含关键词也判作蛋白
    if any(k in lower for k in ("protein", "pep", "aa")):
        # 避免 genome.fa 之类被误判
        if not any(lower.endswith(ext) for ext in (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")):
            return "prot"
    # 基因组 FASTA（且不带 pep/prot 等字样）
    if lower.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")):
        if not any(k in lower for k in ("pep", "prot", "protein")):
            return "genome"
    return "unknown"


def scan_species_in_raw():
    """
    从 raw_data/ 自动识别物种及其 gff/protein/genome 文件。

    返回：
      species_order: [sp1, sp2, ...]
      species_files: {species_name: {"gff": path, "prot": path, "genome": path or None}}
    """
    if not os.path.isdir(DIR_RAW):
        log(f"[错误] 原始数据目录不存在：{DIR_RAW}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # ✅ 修改处（仅修改扫描目录逻辑）：支持 raw_data/{gff,proteins,genome}/ 结构
    # -------------------------------------------------------------------------
    gff_dir = os.path.join(DIR_RAW, "gff")
    prot_dir = os.path.join(DIR_RAW, "proteins")
    genome_dir = os.path.join(DIR_RAW, "genome")

    # 至少需要 gff + proteins 目录存在
    if not os.path.isdir(gff_dir) or not os.path.isdir(prot_dir):
        # 兼容旧布局：若 raw_data/ 根目录下直接放文件，才允许继续旧扫描
        # 但你当前布局是分目录，因此这里优先提示目录缺失
        log(f"[错误] raw_data/ 需要包含子目录：gff/ 与 proteins/（genome/ 可选）")
        log(f"       当前检测到：{DIR_RAW}")
        log(f"       - gff_dir exists? {os.path.isdir(gff_dir)} : {gff_dir}")
        log(f"       - prot_dir exists? {os.path.isdir(prot_dir)} : {prot_dir}")
        sys.exit(1)

    species_files = defaultdict(dict)

    # 1) 扫描 gff/
    for fname in os.listdir(gff_dir):
        if fname.startswith("."):
            continue
        full_path = os.path.join(gff_dir, fname)
        if not os.path.isfile(full_path):
            continue
        ftype = classify_file_type(fname)
        if ftype != "gff":
            continue
        species = fname.split(".")[0]
        species_files[species]["gff"] = full_path

    # 2) 扫描 proteins/
    for fname in os.listdir(prot_dir):
        if fname.startswith("."):
            continue
        full_path = os.path.join(prot_dir, fname)
        if not os.path.isfile(full_path):
            continue
        ftype = classify_file_type(fname)
        if ftype != "prot":
            continue
        species = fname.split(".")[0]
        species_files[species]["prot"] = full_path

    # 3) 扫描 genome/（可选）
    if os.path.isdir(genome_dir):
        for fname in os.listdir(genome_dir):
            if fname.startswith("."):
                continue
            full_path = os.path.join(genome_dir, fname)
            if not os.path.isfile(full_path):
                continue
            ftype = classify_file_type(fname)
            if ftype != "genome":
                continue
            species = fname.split(".")[0]
            species_files[species]["genome"] = full_path

    # -------------------------------------------------------------------------
    # ✅ 修改结束：后续逻辑保持原样
    # -------------------------------------------------------------------------

    # 过滤掉缺少必需文件的
    valid = {
        sp: info
        for sp, info in species_files.items()
        if "gff" in info and "prot" in info
    }

    if len(valid) < 2:
        log(f"[错误] 至少需要 2 个物种（且每个物种至少有 gff + protein），当前有效物种数：{len(valid)}")
        log(f"       当前检测到的物种及文件：{species_files}")
        sys.exit(1)

    # 若用户显式指定物种对，则按指定顺序返回
    if SPECIES_PAIR is not None:
        pair = []
        for tag in SPECIES_PAIR:
            if tag in valid:
                pair.append(tag)
            else:
                # 允许“模糊匹配”：SPECIES_PAIR 中给短名，raw_data 中是长名，尝试包含关系
                matched = [sp for sp in valid.keys() if tag in sp]
                if len(matched) == 1:
                    pair.append(matched[0])
                else:
                    log(f"[错误] SPECIES_PAIR 中的 '{tag}' 在 raw_data/ 中找不到唯一匹配。")
                    log(f"       有效物种：{list(valid.keys())}")
                    sys.exit(1)
        if len(pair) != 2:
            log("[错误] SPECIES_PAIR 必须指定两个物种名称。")
            sys.exit(1)
        species_order = pair
    else:
        # 自动模式：若刚好 2 个物种，就用它们；否则报错让用户手动在 SPECIES_PAIR 指定
        if len(valid) != 2:
            log(f"[错误] 当前识别到 {len(valid)} 个物种：{list(valid.keys())}")
            log("       请在脚本顶部 SPECIES_PAIR 中显式指定要分析的两个物种。")
            sys.exit(1)
        species_order = sorted(valid.keys())

    log(f">> [Init] 识别到物种：{species_order}")
    # 只返回用到的那两个物种
    reduced = {sp: valid[sp] for sp in species_order}
    return species_order, reduced


def detect_gff_feature_type(gff_path: str) -> str:
    """扫描前若干行，自动判断使用 mRNA / transcript / gene"""
    feature_types = set()
    count = 0
    with smart_open(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            feature_types.add(parts[2])
            count += 1
            if count >= GFF_DETECT_LINES:
                break
    if "mRNA" in feature_types:
        ftype = "mRNA"
    elif "transcript" in feature_types:
        ftype = "transcript"
    else:
        ftype = "gene"
    return ftype


def decide_topn_for_species(species: str) -> int:
    """根据物种名称和 SPECIES_SPECIFIC_TOPN 选择 topN"""
    if CHR_FILTER_MODE != "top_n":
        return -1
    for key, val in SPECIES_SPECIFIC_TOPN.items():
        if key in species:
            return val
    return DEFAULT_TOP_N


def parse_gff_to_bed_and_coords(species: str, gff_path: str):
    """
    解析 GFF，提取目标 feature（mRNA/transcript/gene），生成：
      - bed 文件：TEMP_DIR/{species}.bed
      - coords 映射：gene_id -> (chr, start, end)
      - chr 长度统计：[{chr, size, species}, ...]，size 用该 chr 上基因 max_end 近似

    会根据 CHR_FILTER_MODE 和 TopN 决定保留哪些染色体（基于“基因 max_end”排序）。
    """
    log(f"\n>> [{species}] 解析注释：{os.path.basename(gff_path)}")
    ftype = detect_gff_feature_type(gff_path)
    log(f"   - [GFF] 锁定特征类型：{ftype}")

    # 先读取所有候选 feature
    tmp_features = []  # (chr, start, end, gene_id)
    chr_len_map = defaultdict(int)

    with smart_open(gff_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != ftype:
                continue
            chrom = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            attr = parts[8]
            gene_id = None
            # ✅ 修改处：优先使用 transcript_id=，兼容 NCBI Gnomon / RefSeq 风格
            if "transcript_id=" in attr:
                gene_id = attr.split("transcript_id=")[1].split(";")[0].strip().strip('"')
            elif "ID=" in attr:
                gene_id = attr.split("ID=")[1].split(";")[0]
            elif "gene_id" in attr:
                # 兼容部分 gtf 风格
                gene_id = (
                    attr.split("gene_id")[1]
                    .strip()
                    .lstrip("=")
                    .strip('"; ')
                    .split(";")[0]
                )
            if not gene_id:
                continue
            tmp_features.append((chrom, start, end, gene_id))
            if end > chr_len_map[chrom]:
                chr_len_map[chrom] = end

    if not tmp_features:
        log(f"[错误] 在 {gff_path} 中没有解析到任何 {ftype} 记录。")
        sys.exit(1)

    # 根据 chr_len_map 决定保留哪些染色体
    chr_records = []
    for chrom, size in chr_len_map.items():
        chr_records.append({"chr": chrom, "size": int(size), "species": species})

    # 过滤染色体
    if CHR_FILTER_MODE == "top_n":
        topn = decide_topn_for_species(species)
        chr_records_sorted = sorted(chr_records, key=lambda d: d["size"], reverse=True)
        keep_chrs = {d["chr"] for d in chr_records_sorted[:topn]}
        log(f"   - [Genome] 过滤策略：Top {topn}（保留 {len(keep_chrs)}/{len(chr_records)} 条 chr/scaffold）")
        chr_records = [d for d in chr_records_sorted if d["chr"] in keep_chrs]
    else:
        keep_chrs = set(chr_len_map.keys())
        log(f"   - [Genome] 过滤策略：不过滤（保留 {len(keep_chrs)} 条 chr/scaffold）")

    # 生成 bed 和 coords
    coords = {}
    bed_rows = []
    kept_gene = 0
    skipped_gene = 0

    for chrom, start, end, gene_id in tmp_features:
        if chrom not in keep_chrs:
            skipped_gene += 1
            continue
        coords[gene_id] = (chrom, start, end)
        # standard bed: chr, start, end, gene_id
        bed_rows.append((chrom, start, end, gene_id))
        kept_gene += 1

    log(f"   - [GFF] 总 feature 数：{len(tmp_features)}，保留：{kept_gene}，因 chr 过滤跳过：{skipped_gene}")

    # bed 排序：按 chr + start
    bed_rows.sort(key=lambda x: (x[0], x[1]))
    bed_path = os.path.join(TEMP_DIR, f"{species}.bed")
    with open(bed_path, "w") as bed_fh:
        for chrom, start, end, gene_id in bed_rows:
            bed_fh.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

    return bed_path, coords, chr_records


def clean_protein_fasta(species: str, prot_path: str, valid_gene_ids: set):
    """
    清洗蛋白质 FASTA：
      - 如 PROTEIN_STRIP_SPECIES_PREFIX=True，则将 header 里 'xxx|' 前缀去掉，只保留后半段
      - 仅保留 ID 能在 valid_gene_ids 中找到的条目
      - 输出至 TEMP_DIR/{species}.pep
    """
    log(f"   - [Prot] 处理蛋白质文件：{os.path.basename(prot_path)}")
    total = 0
    kept = 0
    out_path = os.path.join(TEMP_DIR, f"{species}.pep")

    cleaned_records = []

    with smart_open(prot_path, "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            total += 1
            raw_id = rec.id
            norm_id = raw_id
            if PROTEIN_STRIP_SPECIES_PREFIX and "|" in norm_id:
                norm_id = norm_id.split("|")[-1]
            # 不做版本截断，让 ID 与 GFF 完全一致（若有需要可再加规则）
            if norm_id not in valid_gene_ids:
                continue
            rec.id = norm_id
            rec.description = ""
            cleaned_records.append(rec)
            kept += 1

    if not cleaned_records:
        log(f"[警告] 物种 {species} 的蛋白文件中没有任何 ID 能匹配到 GFF 基因，请检查命名规则。")

    with open(out_path, "w") as out_fh:
        SeqIO.write(cleaned_records, out_fh, "fasta")

    log(f"   - [Prot] 总序列数：{total}，成功匹配到 GFF 的：{kept}")
    return out_path


def run_jcvi_ortholog(species1: str, species2: str):
    """
    调用 JCVI 进行同线性分析：
      python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names --cpu=N sp1 sp2
      python -m jcvi.compara.synteny screen --minspan=ANCHOR_MINSPAN --simple sp1.sp2.anchors sp1.sp2.anchors.simple
    """
    log("\n>> [JCVI] 调用 jcvi.compara.catalog ortholog 做蛋白同线性（MCscan）")
    cwd = os.getcwd()
    try:
        os.chdir(TEMP_DIR)
        anchors_prefix = f"{species1}.{species2}"

        # Step 1: ortholog
        cmd_orth = [
            sys.executable, "-m", "jcvi.compara.catalog", "ortholog",
            f"--dbtype={JCV_ORTHO_DBTYPE}",
            "--no_strip_names",
            f"--cpu={JCV_ORTHO_CCPU}",
            f"--cscore={JCV_ORTHO_CSCORE}",
            species1, species2,
        ]
        log("   - 运行： " + " ".join(cmd_orth))
        subprocess.run(cmd_orth, check=True)

        anchors_file = anchors_prefix + ".anchors"
        if not os.path.exists(anchors_file):
            log(f"[错误] JCVI 未生成 {anchors_file}，请检查 jcvi 报错信息。")
            sys.exit(1)

        # Step 2: anchors → anchors.simple
        anchors_simple = anchors_prefix + ".anchors.simple"
        cmd_screen = [
            sys.executable, "-m", "jcvi.compara.synteny", "screen",
            "--minspan", str(ANCHOR_MINSPAN),
            "--simple", anchors_file, anchors_simple,
        ]
        log("   - 运行： " + " ".join(cmd_screen))
        subprocess.run(cmd_screen, check=True)

        if not os.path.exists(anchors_simple):
            log(f"[错误] JCVI 未生成 {anchors_simple}，请检查 jcvi 报错信息。")
            sys.exit(1)

        log(f"   - [JCVI] anchors 文件：{os.path.abspath(anchors_file)}")
        log(f"   - [JCVI] anchors.simple 文件：{os.path.abspath(anchors_simple)}")
        return os.path.abspath(anchors_file), os.path.abspath(anchors_simple)

    finally:
        os.chdir(cwd)


def load_bed_coords(bed_path: str):
    """
    从 bed 文件读取坐标映射：
      gene_id -> (chr, start, end)
    """
    coords = {}
    with open(bed_path, "r") as fh:
        for line in fh:
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
            coords[gene_id] = (chrom, start, end)
    return coords


def parse_anchors_to_blocks(anchors_path: str,
                            coords1: dict,
                            coords2: dict,
                            label_prefix: str):
    """
    解析 JCVI 生成的 anchors.simple 文件，转换为 block 坐标：
      anchors.simple 每行格式（6 列）：
        col1: 物种1的起始基因 ID
        col2: 物种1的结束基因 ID
        col3: 物种2的起始基因 ID
        col4: 物种2的结束基因 ID
        col5: 评分（score）
        col6: 方向（strand）

    输出：
      blocks: [[block_id, chr1, start1, end1, chr2, start2, end2], ...]
    """
    blocks = []
    if not os.path.exists(anchors_path):
        log(f"[警告] anchors.simple 文件不存在：{anchors_path}")
        return blocks

    with open(anchors_path, "r") as fh:
        idx = 0
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            g1_start, g1_end, g2_start, g2_end = parts[:4]
            if (g1_start not in coords1 or g1_end not in coords1 or
                    g2_start not in coords2 or g2_end not in coords2):
                continue
            chr1, s1a, e1a = coords1[g1_start]
            chr1b, s1b, e1b = coords1[g1_end]
            chr2, s2a, e2a = coords2[g2_start]
            chr2b, s2b, e2b = coords2[g2_end]
            # 若同一物种内部出现跨染色体 block，这里仍然保留，方便分析复杂重排
            start1 = min(s1a, s1b, e1a, e1b)
            end1 = max(s1a, s1b, e1a, e1b)
            start2 = min(s2a, s2b, e2a, e2b)
            end2 = max(s2a, s2b, e2a, e2b)
            idx += 1
            block_id = f"{label_prefix}_{idx}"
            blocks.append([block_id, chr1, start1, end1, chr2, start2, end2])
    return blocks


def summarize_blocks(blocks, species1: str, species2: str):
    """
    对 block 进行 chr1-chr2 统计，生成 synteny_summary.tsv 的数据：
      species1, chr1, species2, chr2, n_blocks
    """
    combo_counts = defaultdict(int)
    for _bid, chr1, _s1, _e1, chr2, _s2, _e2 in blocks:
        combo_counts[(chr1, chr2)] += 1

    rows = []
    for (chr1, chr2), n in sorted(combo_counts.items(), key=lambda x: (-x[1], x[0][0], x[0][1])):
        rows.append({
            "species1": species1,
            "chr1": chr1,
            "species2": species2,
            "chr2": chr2,
            "n_blocks": n,
        })
    return rows


# =============================================================================
# 3. 主流程
# =============================================================================

def main():
    check_dependencies()

    # 准备目录
    os.makedirs(DIR_IN, exist_ok=True)
    os.makedirs(DIR_OUTPUT_ROOT, exist_ok=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    os.makedirs(TEMP_DIR, exist_ok=True)

    # 识别物种 & 文件
    species_order, species_files = scan_species_in_raw()
    sp1, sp2 = species_order[0], species_order[1]

    all_chr_records = []
    bed_paths = {}
    coords_by_species = {}

    # 逐物种处理：GFF → bed / coords，protein → pep
    for sp in species_order:
        info = species_files[sp]
        gff_path = info["gff"]
        prot_path = info["prot"]

        # GFF → bed/coords + chr 长度
        bed_path, coords, chr_records = parse_gff_to_bed_and_coords(sp, gff_path)
        bed_paths[sp] = bed_path
        coords_by_species[sp] = coords
        all_chr_records.extend(chr_records)

        # protein → pep（ID 清洗 + 与 coords 交集）
        _pep_path = clean_protein_fasta(sp, prot_path, set(coords.keys()))

    # 输出 genome_len.tsv（给 R 画圈图用）
    genome_len_path = os.path.join(DIR_IN, "genome_len.tsv")
    pd.DataFrame(all_chr_records)[["chr", "size", "species"]].to_csv(
        genome_len_path, sep="\t", index=False
    )
    log(f"\n>> [Output] 染色体长度表：{genome_len_path}")

    # 调用 JCVI 做 ortholog + anchors.simple
    anchors_path, anchors_simple_path = run_jcvi_ortholog(sp1, sp2)

    # 重新从 bed 读 coords（防止未来 bed 手动再加工时 coords 失配，这里以 bed 为准）
    coords1 = load_bed_coords(bed_paths[sp1])
    coords2 = load_bed_coords(bed_paths[sp2])

    # 方案 A：仅解析 anchors.simple，用于 main 与 full（复杂重排也基于筛选后 block）
    log(f"\n>> [Parse] 解析 anchors.simple（筛选后 block）：{anchors_simple_path}")
    blocks_main = parse_anchors_to_blocks(
        anchors_simple_path, coords1, coords2, label_prefix="main"
    )
    log(f"   - main blocks 数量：{len(blocks_main)}")

    # full 直接等于 main（即 anchors.simple 的全部 block）
    blocks_full = list(blocks_main)
    log(f"   - full blocks（复用 anchors.simple）数量：{len(blocks_full)}")

    # 写出 synteny_blocks.tsv（主图）与 synteny_blocks_full.tsv（全量=同 main）
    cols = ["block_id", "chr1", "start1", "end1", "chr2", "start2", "end2"]

    blocks_main_path = os.path.join(DIR_IN, "synteny_blocks.tsv")
    pd.DataFrame(blocks_main, columns=cols).to_csv(
        blocks_main_path, sep="\t", index=False
    )
    log(f">> [Output] 主图 block 表：{blocks_main_path}（共 {len(blocks_main)} 条）")

    blocks_full_path = os.path.join(DIR_IN, "synteny_blocks_full.tsv")
    pd.DataFrame(blocks_full, columns=cols).to_csv(
        blocks_full_path, sep="\t", index=False
    )
    log(f">> [Output] 全量 block 表：{blocks_full_path}（共 {len(blocks_full)} 条）")

    # 生成 synteny_summary.tsv（基于 full blocks = anchors.simple）
    if blocks_full:
        summary_rows = summarize_blocks(blocks_full, sp1, sp2)
        summary_path = os.path.join(DIR_IN, "synteny_summary.tsv")
        pd.DataFrame(summary_rows)[
            ["species1", "chr1", "species2", "chr2", "n_blocks"]
        ].to_csv(summary_path, sep="\t", index=False)
        log(f">> [Output] 染色体对汇总表：{summary_path}")
    else:
        log("[警告] full blocks 为空，无法生成 synteny_summary.tsv")

    log("\n>> [完成] JCVI 同线性预处理流程结束。")
    log(f"   - 染色体长度表：{genome_len_path}")
    log(f"   - 主图 block 表：{blocks_main_path}")
    log(f"   - 全量 block 表：{blocks_full_path}")
    log(f"   - 汇总表（若有）：input/synteny_summary.tsv")
    log("   - 下一步：用 R / circlize 读取 genome_len.tsv + synteny_blocks.tsv 画圈图；")
    log("             若要分析复杂重排，可使用 synteny_blocks_full.tsv + synteny_summary.tsv。")


if __name__ == "__main__":
    main()

