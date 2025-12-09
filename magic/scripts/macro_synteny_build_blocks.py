#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
宏观染色体共线性 - 步骤2：根据 anchors 构建宏观 blocks & links

功能：
1. 读取 meta/species_info.tsv，确定参考物种染色体及其长度；
2. 构建参考染色体的累计坐标 & 颜色表，写出：
   - meta/ref_chr_color.tsv
   - meta/species_order.tsv
3. 从各物种 GFF 中解析基因位点，建立 gene_id -> (seqid, mid_pos) 映射；
4. 从 synteny_pairwise 中读取 <ref>_vs_<sp>.anchors.simple，
   将 anchor 映射到 (ref_chr, ref_pos_bp) 和 (sp_chr, sp_pos_bp)；
5. 将参考染色体按窗口大小划分，统计每物种、每窗口的 anchor 数，
   写出：
   - blocks/blocks.tsv  （条带用）
   - links/links.tsv    （丝带用）

兼容：
- genome：.fa/.fasta/.fa.gz/.fasta.gz/.fna/.fna.gz
- gff   ：.gff3/.gff/.gff.gz/.gff3.gz
- GFF feature：gene + mRNA 均可
"""

import sys
import csv
import gzip
import logging
from math import ceil
from pathlib import Path
from typing import Dict, Tuple, List

# ====================== 用户参数区 ======================

REFERENCE_SPECIES = "Sinonovacula_constricta"

CHR_LENGTH_THRESHOLD_MBP = 10   # 仅用于 sanity check
WINDOW_SIZE_MBP = 5
MAX_LINKS_PER_CHR = 2000        # 每个物种、每条参考染色体最多保留多少条丝带

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

# ====================== 通用函数 ======================


def setup_paths() -> dict:
    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]

    raw_data_dir = project_root / "raw_data"
    genome_dir = raw_data_dir / "genome"
    gff_dir = raw_data_dir / "gff"

    output_root = project_root / "output" / "macro_synteny"
    log_dir = output_root / "logs"
    meta_dir = output_root / "meta"
    synteny_dir = output_root / "synteny_pairwise"
    blocks_dir = output_root / "blocks"
    links_dir = output_root / "links"

    for d in [output_root, log_dir, meta_dir, synteny_dir, blocks_dir, links_dir]:
        d.mkdir(parents=True, exist_ok=True)

    return {
        "project_root": project_root,
        "raw_data": raw_data_dir,
        "genome_dir": genome_dir,
        "gff_dir": gff_dir,
        "output_root": output_root,
        "log_dir": log_dir,
        "meta_dir": meta_dir,
        "synteny_dir": synteny_dir,
        "blocks_dir": blocks_dir,
        "links_dir": links_dir,
    }


def setup_logger(log_dir: Path, name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    log_file = log_dir / f"{name}.log"
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh.setFormatter(fmt)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def prefix_to_group(prefix: str) -> str:
    if prefix in GROUP_DIGGING:
        return "digging"
    if prefix in GROUP_NON_DIGGING:
        return "nondigging"
    return "unknown"


def load_species_info(path: Path, logger: logging.Logger):
    """
    读取 species_info.tsv，返回 dict: species_info[prefix] = {...}
    """
    if not path.exists():
        logger.error(f"缺少 species_info.tsv：{path}，请先运行步骤1。")
        sys.exit(1)

    species_info = {}
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            prefix = parts[idx["species_prefix"]]
            species_label = parts[idx["species_label"]]
            group = parts[idx["group"]]
            genome_size = int(parts[idx["genome_size_bp"]])
            n_seq_total = int(parts[idx["n_seq_total"]])
            n_chr = int(parts[idx["n_chr"]])
            n_scaffold = int(parts[idx["n_scaffold"]])
            chr_ids = parts[idx["chr_ids"]].split(",") if parts[idx["chr_ids"]] else []
            chr_lengths = (
                [int(x) for x in parts[idx["chr_lengths_bp"]].split(",")]
                if parts[idx["chr_lengths_bp"]]
                else []
            )
            species_info[prefix] = {
                "species_prefix": prefix,
                "species_label": species_label,
                "group": group,
                "genome_size_bp": genome_size,
                "n_seq_total": n_seq_total,
                "n_chr": n_chr,
                "n_scaffold": n_scaffold,
                "chr_ids": chr_ids,
                "chr_lengths_bp": chr_lengths,
            }

    logger.info(f"已加载 {len(species_info)} 个物种的 species_info。")
    return species_info


def read_gff_gene_positions(path: Path, logger: logging.Logger) -> Dict[str, Tuple[str, int]]:
    """
    从 gff/gff3(.gz) 文件中读取基因（或 mRNA）的位置信息：
    返回 gene_pos[gene_id] = (seqid, mid_pos_bp)

    兼容：
      - feature type 为 "gene" 或 "mRNA"
      - attributes 中包含 ID=xxx
    """
    gene_pos: Dict[str, Tuple[str, int]] = {}

    if not path.exists():
        logger.error(f"GFF 文件不存在：{path}")
        return gene_pos

    logger.info(f"解析 GFF：{path.name}")
    if str(path).endswith(".gz"):
        fh = gzip.open(path, "rt")
    else:
        fh = open(path, "r")

    with fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr = parts
            if ftype not in ("gene", "mRNA"):
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            mid = (start_i + end_i) // 2

            gene_id = None
            for field in attr.split(";"):
                if field.startswith("ID="):
                    gene_id = field.split("=", 1)[1]
                    break
            if gene_id is None:
                continue

            gene_pos[gene_id] = (seqid, mid)

    logger.info(f"GFF 中共识别到 {len(gene_pos)} 个 gene/mRNA feature。")
    return gene_pos


def read_anchors_simple(path: Path, logger: logging.Logger) -> List[Tuple[str, str]]:
    """
    读取 anchors.simple：假定前两列分别为 ref_gene_id 和 sp_gene_id。
    """
    pairs: List[Tuple[str, str]] = []
    if not path.exists():
        logger.warning(f"anchors.simple 文件不存在：{path}")
        return pairs

    logger.info(f"读取 anchors.simple：{path.name}")
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split()
            if len(cols) < 2:
                continue
            ref_gene, sp_gene = cols[0], cols[1]
            pairs.append((ref_gene, sp_gene))

    logger.info(f"共读取 {len(pairs)} 条 anchors。")
    return pairs


def generate_ref_chr_layout(ref_info: dict) -> List[dict]:
    """
    根据参考物种的 chr_ids 和 chr_lengths，生成累计 X 坐标和颜色编号。
    返回列表中每个元素：
    {
        "ref_chr_id", "ref_chr_rank", "length_bp",
        "x_start", "x_end", "color_id"
    }
    """
    layout = []
    offset = 0.0
    chr_ids = ref_info["chr_ids"]
    chr_lengths = ref_info["chr_lengths_bp"]

    for i, (cid, ln) in enumerate(zip(chr_ids, chr_lengths), start=1):
        x_start = offset
        x_end = offset + ln
        layout.append({
            "ref_chr_id": cid,
            "ref_chr_rank": i,
            "length_bp": ln,
            "x_start": x_start,
            "x_end": x_end,
            "color_id": i,
        })
        offset = x_end

    return layout


def hsv_color_wheel(n: int) -> List[str]:
    """
    简单生成 n 个分布均匀的 HSV 颜色，返回 hex 列表。
    """
    import colorsys
    colors = []
    for i in range(n):
        h = i / float(max(1, n))
        s = 0.6
        v = 0.9
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        colors.append("#{0:02X}{1:02X}{2:02X}".format(int(r * 255), int(g * 255), int(b * 255)))
    return colors


# ====================== 主流程 ======================


def main():
    # 复用步骤1中的 read_fasta_lengths 函数
    from macro_synteny_prepare_pairwise import read_fasta_lengths

    paths = setup_paths()
    logger = setup_logger(paths["log_dir"], "macro_synteny_build_blocks")

    logger.info("=== 宏观共线性步骤2：构建 blocks & links 开始 ===")
    window_bp = int(WINDOW_SIZE_MBP * 1_000_000)
    logger.info(f"参考物种：{REFERENCE_SPECIES}")
    logger.info(f"窗口大小：{WINDOW_SIZE_MBP} Mb ({window_bp} bp)")

    species_info_path = paths["meta_dir"] / "species_info.tsv"
    species_info = load_species_info(species_info_path, logger)

    if REFERENCE_SPECIES not in species_info:
        logger.error(f"species_info.tsv 中找不到参考物种 {REFERENCE_SPECIES}")
        sys.exit(1)

    ref_info = species_info[REFERENCE_SPECIES]
    if not ref_info["chr_ids"]:
        logger.error(f"参考物种 {REFERENCE_SPECIES} 没有识别到任何染色体（chr_ids 为空），请检查步骤1阈值。")
        sys.exit(1)

    # 1. 生成参考染色体布局 & 颜色表
    ref_layout = generate_ref_chr_layout(ref_info)
    n_chr = len(ref_layout)
    logger.info(f"参考物种共有 {n_chr} 条染色体参与宏观图。")

    color_hex_list = hsv_color_wheel(n_chr)
    for i, row in enumerate(ref_layout):
        row["color_hex"] = color_hex_list[i]

    ref_chr_color_path = paths["meta_dir"] / "ref_chr_color.tsv"
    logger.info(f"写入参考染色体颜色表 -> {ref_chr_color_path}")
    with open(ref_chr_color_path, "w", encoding="utf-8") as fo:
        fo.write("\t".join([
            "ref_chr_id",
            "ref_chr_rank",
            "length_bp",
            "x_start",
            "x_end",
            "color_id",
            "color_hex",
        ]) + "\n")
        for row in ref_layout:
            fo.write("\t".join([
                row["ref_chr_id"],
                str(row["ref_chr_rank"]),
                str(row["length_bp"]),
                str(row["x_start"]),
                str(row["x_end"]),
                str(row["color_id"]),
                row["color_hex"],
            ]) + "\n")

    # 2. 生成物种顺序表（row_index）
    species_order_path = paths["meta_dir"] / "species_order.tsv"
    prefixes = sorted(species_info.keys())

    def sort_key(p):
        g = prefix_to_group(p)
        order = {"digging": 0, "nondigging": 1}.get(g, 2)
        return (order, p)

    sorted_prefixes = sorted(prefixes, key=sort_key)

    logger.info(f"为宏观图生成物种行顺序 -> {species_order_path}")
    with open(species_order_path, "w", encoding="utf-8") as fo:
        fo.write("\t".join(["row_index", "species_prefix", "species_label", "group"]) + "\n")
        # row_index=0 保留给参考物种
        fo.write("\t".join([
            "0",
            REFERENCE_SPECIES,
            ref_info["species_label"],
            ref_info["group"],
        ]) + "\n")
        row_index = 1
        for p in sorted_prefixes:
            if p == REFERENCE_SPECIES:
                continue
            info = species_info[p]
            fo.write("\t".join([
                str(row_index),
                p,
                info["species_label"],
                info["group"],
            ]) + "\n")
            row_index += 1

    # 映射 prefix -> row_index
    prefix_to_row = {}
    with open(species_order_path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            r = int(parts[idx["row_index"]])
            p = parts[idx["species_prefix"]]
            prefix_to_row[p] = r

    # 3. 为参考物种加载 genome 长度 & GFF gene 位置
    logger.info("为参考物种加载 genome 序列长度信息 …")
    ref_genome_file = None
    for ext in [".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna", ".fna.gz"]:
        cand = paths["genome_dir"] / f"{REFERENCE_SPECIES}{ext}"
        if cand.exists():
            ref_genome_file = cand
            break
    if ref_genome_file is None:
        logger.error(f"找不到参考物种 {REFERENCE_SPECIES} 的 genome 文件。")
        sys.exit(1)

    seq_len_cache: Dict[str, Dict[str, int]] = {
        REFERENCE_SPECIES: read_fasta_lengths(ref_genome_file)
    }

    logger.info("为参考物种加载 GFF gene/mRNA 位置信息 …")
    ref_gff_file = None
    for ext in [".gff3", ".gff", ".gff.gz", ".gff3.gz"]:
        cand = paths["gff_dir"] / f"{REFERENCE_SPECIES}{ext}"
        if cand.exists():
            ref_gff_file = cand
            break
    if ref_gff_file is None:
        logger.error(f"找不到参考物种 {REFERENCE_SPECIES} 的 GFF 文件。")
        sys.exit(1)

    gene_pos_cache: Dict[str, Dict[str, Tuple[str, int]]] = {}
    gene_pos_cache[REFERENCE_SPECIES] = read_gff_gene_positions(ref_gff_file, logger)

    ref_chr_by_id = {row["ref_chr_id"]: row for row in ref_layout}

    # 4. 遍历其它物种的 anchors.simple，构建 blocks_stats 和 links_raw
    blocks_stats = {}   # (sp_prefix, ref_chr_rank, window_id) -> {"n_anchors", "sp_chr_counts"}
    links_raw = {}      # (sp_prefix, ref_chr_rank) -> [(ref_x, sp_seqid, sp_mid), ...]

    for sp_prefix in sorted_prefixes:
        if sp_prefix == REFERENCE_SPECIES:
            continue

        anchors_file = paths["synteny_dir"] / f"{REFERENCE_SPECIES}_vs_{sp_prefix}.anchors.simple"
        if not anchors_file.exists():
            logger.warning(f"缺少 anchors.simple：{anchors_file}，该物种将不会出现在 blocks/links 中。")
            continue

        logger.info(f"处理物种 {sp_prefix} 的 anchors.simple …")

        # genome 长度
        if sp_prefix not in seq_len_cache:
            genome_file = None
            for ext in [".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna", ".fna.gz"]:
                cand = paths["genome_dir"] / f"{sp_prefix}{ext}"
                if cand.exists():
                    genome_file = cand
                    break
            if genome_file is None:
                logger.error(f"物种 {sp_prefix} 找不到 genome 文件，跳过。")
                continue
            seq_len_cache[sp_prefix] = read_fasta_lengths(genome_file)

        # GFF gene 位置
        if sp_prefix not in gene_pos_cache:
            gff_file = None
            for ext in [".gff3", ".gff", ".gff.gz", ".gff3.gz"]:
                cand = paths["gff_dir"] / f"{sp_prefix}{ext}"
                if cand.exists():
                    gff_file = cand
                    break
            if gff_file is None:
                logger.error(f"物种 {sp_prefix} 找不到 GFF 文件，跳过。")
                continue
            gene_pos_cache[sp_prefix] = read_gff_gene_positions(gff_file, logger)

        ref_gene_pos = gene_pos_cache[REFERENCE_SPECIES]
        sp_gene_pos = gene_pos_cache[sp_prefix]
        seq_len_sp = seq_len_cache[sp_prefix]

        sp_info = species_info[sp_prefix]
        sp_chr_set = set(sp_info["chr_ids"])
        sp_chr_rank_map = {cid: i + 1 for i, cid in enumerate(sp_info["chr_ids"])}

        pairs = read_anchors_simple(anchors_file, logger)
        if not pairs:
            logger.warning(f"{sp_prefix} 与 {REFERENCE_SPECIES} 的 anchors.simple 为空或读取失败。")
            continue

        # 初始化 links_raw
        for ref_row in ref_layout:
            key = (sp_prefix, ref_row["ref_chr_rank"])
            links_raw.setdefault(key, [])

        skipped = 0
        for ref_gene, sp_gene in pairs:
            if ref_gene not in ref_gene_pos or sp_gene not in sp_gene_pos:
                skipped += 1
                continue
            ref_seqid, ref_mid = ref_gene_pos[ref_gene]
            sp_seqid, sp_mid = sp_gene_pos[sp_gene]

            if ref_seqid not in ref_chr_by_id:
                continue
            ref_chr_row = ref_chr_by_id[ref_seqid]
            ref_chr_rank = ref_chr_row["ref_chr_rank"]
            ref_x = ref_chr_row["x_start"] + ref_mid

            key = (sp_prefix, ref_chr_rank)
            links_raw[key].append((ref_x, sp_seqid, sp_mid))

            # 统计 blocks
            window_id = int(ref_mid // window_bp) + 1
            stats_key = (sp_prefix, ref_chr_rank, window_id)
            if stats_key not in blocks_stats:
                blocks_stats[stats_key] = {"n_anchors": 0, "sp_chr_counts": {}}
            blocks_stats[stats_key]["n_anchors"] += 1
            d = blocks_stats[stats_key]["sp_chr_counts"]
            d[sp_seqid] = d.get(sp_seqid, 0) + 1

        logger.info(f"{sp_prefix}：anchors 总数={len(pairs)}，缺位基因={skipped}")

    # 5. 写 blocks.tsv
    blocks_path = paths["blocks_dir"] / "blocks.tsv"
    logger.info(f"写出宏观条带表 -> {blocks_path}")

    # 读取 species_order，方便取 row_index / label / group
    order_records = {}
    with open(species_order_path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            row_index = int(parts[idx["row_index"]])
            sp_prefix = parts[idx["species_prefix"]]
            sp_label = parts[idx["species_label"]]
            group = parts[idx["group"]]
            order_records[sp_prefix] = {
                "row_index": row_index,
                "species_label": sp_label,
                "group": group,
            }

    with open(blocks_path, "w", encoding="utf-8") as fo:
        fo.write("\t".join([
            "species_prefix",
            "species_label",
            "group",
            "row_index",
            "ref_chr_id",
            "ref_chr_rank",
            "window_id",
            "window_start_bp",
            "window_end_bp",
            "x_start",
            "x_end",
            "sp_chr_id",
            "sp_chr_is_chr",
            "sp_chr_rank",
            "n_anchors",
            "color_id",
        ]) + "\n")

        for sp_prefix in sorted_prefixes:
            if sp_prefix == REFERENCE_SPECIES:
                continue
            if sp_prefix not in order_records:
                continue

            order = order_records[sp_prefix]
            row_index = order["row_index"]
            species_label = order["species_label"]
            group = order["group"]

            sp_info = species_info[sp_prefix]
            sp_chr_set = set(sp_info["chr_ids"])
            sp_chr_rank_map = {cid: i + 1 for i, cid in enumerate(sp_info["chr_ids"])}

            for ref_row in ref_layout:
                ref_chr_id = ref_row["ref_chr_id"]
                ref_chr_rank = ref_row["ref_chr_rank"]
                chr_len = ref_row["length_bp"]
                chr_x_start = ref_row["x_start"]
                color_id = ref_row["color_id"]

                n_windows = ceil(chr_len / window_bp)
                for window_id in range(1, n_windows + 1):
                    window_start_bp = (window_id - 1) * window_bp
                    window_end_bp = min(window_id * window_bp, chr_len)
                    x_start = chr_x_start + window_start_bp
                    x_end = chr_x_start + window_end_bp

                    stats_key = (sp_prefix, ref_chr_rank, window_id)
                    if stats_key in blocks_stats:
                        st = blocks_stats[stats_key]
                        n_anchors = st["n_anchors"]
                        if st["sp_chr_counts"]:
                            sp_chr_id = max(st["sp_chr_counts"].items(), key=lambda x: x[1])[0]
                        else:
                            sp_chr_id = ""
                    else:
                        n_anchors = 0
                        sp_chr_id = ""

                    if sp_chr_id and sp_chr_id in sp_chr_set:
                        sp_chr_is_chr = 1
                        sp_chr_rank = sp_chr_rank_map.get(sp_chr_id, 0)
                    elif sp_chr_id:
                        sp_chr_is_chr = 0
                        sp_chr_rank = 0
                    else:
                        sp_chr_is_chr = 0
                        sp_chr_rank = 0

                    fo.write("\t".join([
                        sp_prefix,
                        species_label,
                        group,
                        str(row_index),
                        ref_chr_id,
                        str(ref_chr_rank),
                        str(window_id),
                        str(window_start_bp),
                        str(window_end_bp),
                        str(x_start),
                        str(x_end),
                        sp_chr_id,
                        str(sp_chr_is_chr),
                        str(sp_chr_rank),
                        str(n_anchors),
                        str(color_id),
                    ]) + "\n")

    # 6. 写 links.tsv（抽样后的丝带坐标）
    links_path = paths["links_dir"] / "links.tsv"
    logger.info(f"写出丝带坐标表 -> {links_path}")

    with open(links_path, "w", encoding="utf-8") as fo:
        fo.write("\t".join([
            "species_prefix",
            "species_label",
            "group",
            "row_index",
            "ref_chr_id",
            "ref_chr_rank",
            "ref_x",
            "sp_chr_id",
            "sp_chr_rank",
            "sp_local_y",
            "color_id",
        ]) + "\n")

        for sp_prefix in sorted_prefixes:
            if sp_prefix == REFERENCE_SPECIES:
                continue
            if sp_prefix not in order_records:
                continue

            order = order_records[sp_prefix]
            row_index = order["row_index"]
            species_label = order["species_label"]
            group = order["group"]

            sp_info = species_info[sp_prefix]
            sp_chr_set = set(sp_info["chr_ids"])
            sp_chr_rank_map = {cid: i + 1 for i, cid in enumerate(sp_info["chr_ids"])}
            seq_len_sp = seq_len_cache.get(sp_prefix, {})

            for ref_row in ref_layout:
                ref_chr_id = ref_row["ref_chr_id"]
                ref_chr_rank = ref_row["ref_chr_rank"]
                color_id = ref_row["color_id"]
                key = (sp_prefix, ref_chr_rank)
                if key not in links_raw:
                    continue
                anchors_list = links_raw[key]
                if not anchors_list:
                    continue

                total = len(anchors_list)
                if total > MAX_LINKS_PER_CHR:
                    step = total / MAX_LINKS_PER_CHR
                    sampled = [anchors_list[int(i * step)] for i in range(MAX_LINKS_PER_CHR)]
                else:
                    sampled = anchors_list

                for ref_x, sp_chr_id, sp_pos_bp in sampled:
                    if sp_chr_id in sp_chr_set:
                        sp_chr_rank = sp_chr_rank_map.get(sp_chr_id, 0)
                    else:
                        sp_chr_rank = 0

                    chr_len_sp = seq_len_sp.get(sp_chr_id, None)
                    if chr_len_sp and chr_len_sp > 0:
                        frac = sp_pos_bp / chr_len_sp
                        frac = max(0.0, min(1.0, frac))
                    else:
                        frac = 0.5
                    # 映射到 0.1~0.9 的局部 y
                    sp_local_y = 0.1 + 0.8 * frac

                    fo.write("\t".join([
                        sp_prefix,
                        species_label,
                        group,
                        str(row_index),
                        ref_chr_id,
                        str(ref_chr_rank),
                        str(ref_x),
                        sp_chr_id,
                        str(sp_chr_rank),
                        f"{sp_local_y:.4f}",
                        str(color_id),
                    ]) + "\n")

    logger.info("=== 宏观共线性步骤2 完成 ===")


if __name__ == "__main__":
    main()

