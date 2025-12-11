#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny05_plot_karyotype.py —— JCVI karyotype 线性宏观共线性终版脚本

职责：
  1) 从 Step03 的 FASTA 风格 seqids 中解析物种顺序；
  2) 为每个物种构建 JCVI 需要的 bed_noheader（基因级别，每行一个基因，统一灰色）；
  3) 从 Step02 的 anchors_for_links 中构建 simple/ 目录，并且：
       - 仅修正 score 列，把 "62.0" 这种浮点数改成 "62" 这种整数字符串；
       - 其他字段保持不动，避免破坏 JCVI 原有解析逻辑；
  4) 按物种顺序构建 layout 文件，并将参考物种放在靠近底部；
  5) 仅调用一次 JCVI 生成 PDF：synteny_linear.pdf（不再尝试 SVG）。

注意：
  - 所有路径、参数集中在脚本顶部 PARAM 区域配置；
  - 不依赖命令行参数。
"""

import os
import sys
import shutil
import logging
import subprocess
from typing import Dict, List, Tuple

# =========================
# PARAM 区：皇上可按需修改
# =========================

# 项目根目录
PROJECT_ROOT = "/home/mxrmiss/project/magic"

# 各步骤输出目录
OUTPUT_SUBDIR = "output/synteny_05_plot"
STEP01_DIR = os.path.join(PROJECT_ROOT, "output/synteny_01_mcscan_catalog")
STEP02_DIR = os.path.join(PROJECT_ROOT, "output/synteny_02_blocks_macro")
STEP03_DIR = os.path.join(PROJECT_ROOT, "output/synteny_03_chr_order")
STEP04_DIR = os.path.join(PROJECT_ROOT, "output/synteny_04_layout")

# 参考物种 species_id（必须与 seqids 中的物种名一致）
REFERENCE_SPECIES = "Sinonovacula_constricta"

# 图像参数（JCVI figsize / dpi）
FIGSIZE = "18x10"
DPI = 300

# 统一的基因底色（灰色）
DEFAULT_GENE_COLOR = "#CCCCCC"

# 物种 → 短代号映射（作为 meta 缺失时的兜底）
# 必须与 anchors_for_links 中已经存在的文件名前缀一致
DEFAULT_SPECIES_SHORT: Dict[str, str] = {
    "Argopecten_irradians": "Air",
    "Pecten_maximus": "Pma",
    "Ctenoides_ales": "Cal",
    "Crassostrea_gigas": "Cgi",
    "Ostrea_edulis": "Oed",
    "Pinctada_fucata": "Pfu",
    "Mytilus_edulis": "Med",
    "Mytilus_galloprovincialis": "Mga",
    "Tegillarca_granosa": "Tgr",
    "Mercenaria_mercenaria": "Mmc",
    "Meretrix_meretrix": "Mme",
    "Mya_arenaria": "Mar",
    "Panopea_generosa": "Pge",
    "Novaculina_chinensis": "Nch",
    "Sinonovacula_constricta": "Sco",
    "Sinonovacula_rivularis": "Sri",
}


# =========================
# 日志与工具函数
# =========================

def setup_logger() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def run_cmd(cmd: List[str], cwd: str = None) -> Tuple[int, str, str]:
    """
    运行外部命令，返回 (returncode, stdout, stderr)
    """
    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    out, err = proc.communicate()
    return proc.returncode, out, err


# =========================
# 核心函数
# =========================

def get_paths() -> Dict[str, str]:
    """
    构建本脚本用到的关键路径。
    """
    out_root = os.path.join(PROJECT_ROOT, OUTPUT_SUBDIR)

    paths = {
        "OUTPUT_ROOT": out_root,
        "SEQIDS_SRC": os.path.join(STEP03_DIR, "seqids"),
        "SEQIDS_BAK": os.path.join(out_root, "seqids_from_step03"),
        "SEQIDS_JCVI": os.path.join(out_root, "seqids"),
        "LAYOUT": os.path.join(out_root, "layout"),
        "PLOT_META": os.path.join(out_root, "plot_meta_used_files.tsv"),
        "BED_DIR": os.path.join(out_root, "bed_noheader"),
        "SIMPLE_DIR": os.path.join(out_root, "simple"),
        "META_TSV": os.path.join(PROJECT_ROOT, "raw_data", "synteny_species_meta.tsv"),
        "ANCHORS_LINKS_DIR": os.path.join(STEP02_DIR, "anchors_for_links"),
        "GENEORDER_DIR": os.path.join(STEP01_DIR, "geneorder"),
        "PDF_OUT": os.path.join(out_root, "synteny_linear.pdf"),
    }
    return paths


def parse_species_order_from_seqids_fasta(seqids_fasta: str) -> List[str]:
    """
    从 Step03 的 FASTA 风格 seqids 中解析物种顺序：
      >Argopecten_irradians
      Chr5
      Chr1
      ...

    返回 species_id 列表（不带 ">"）。
    """
    species: List[str] = []
    with open(seqids_fasta) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sid = line[1:].strip()
                if sid:
                    species.append(sid)
    logging.info("从 Step03 seqids 中解析出 %d 个物种：%s", len(species), ", ".join(species))
    return species


def load_species_short_labels(meta_tsv: str, species_list: List[str]) -> Dict[str, str]:
    """
    从 meta 文件加载 species_id → short_label 映射。
    如果 meta 文件丢失必要列，则回退到 DEFAULT_SPECIES_SHORT。
    """
    mapping: Dict[str, str] = {}
    if os.path.exists(meta_tsv):
        with open(meta_tsv) as fh:
            header = fh.readline().strip().split("\t")
            if "species_id" in header and "short_label" in header:
                col_sp = header.index("species_id")
                col_sh = header.index("short_label")
                for line in fh:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= max(col_sp, col_sh):
                        continue
                    sp = parts[col_sp]
                    sh = parts[col_sh]
                    if sp and sh:
                        mapping[sp] = sh
            else:
                logging.warning("meta 文件缺少 species_id 或 short_label 列，将使用内置 DEFAULT_SPECIES_SHORT 映射")
    else:
        logging.warning("未找到 meta 文件：%s，将使用内置 DEFAULT_SPECIES_SHORT 映射", meta_tsv)

    # 若 meta 为空或不完整，补充 DEFAULT_SPECIES_SHORT
    for sp in species_list:
        if sp not in mapping:
            if sp in DEFAULT_SPECIES_SHORT:
                mapping[sp] = DEFAULT_SPECIES_SHORT[sp]
            else:
                # 最兜底：按名称自动生成（首字母 + 后两个字母），尽量不冲突
                parts = sp.split("_")
                short = (parts[0][0] + parts[1][0:2]).capitalize() if len(parts) >= 2 else sp[:3].capitalize()
                mapping[sp] = short
                logging.warning("物种 %s 未在 DEFAULT_SPECIES_SHORT 中，自动生成短代号 %s", sp, short)

    logging.info(
        "物种短代号映射：%s",
        ", ".join(f"{sp}→{sh}" for sp, sh in mapping.items())
    )
    return mapping


def write_jcvi_seqids(seqids_src: str, seqids_out: str, species_list: List[str]) -> None:
    """
    将 Step03 的 FASTA 风格 seqids 转换为 JCVI 需要的纯染色体列格式：

    原始：
      >Argopecten_irradians
      Chr5
      Chr1
      ...

    JCVI：
      Chr5
      Chr1
      ...
      ChrX
      ChrY
      ...

    各物种按 species_list 顺序依次写出，中间不加空行。
    """
    # 先解析每个物种的染色体顺序
    species_to_chrs: Dict[str, List[str]] = {sp: [] for sp in species_list}
    current_sp = None

    with open(seqids_src) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sid = line[1:].strip()
                current_sp = sid if sid in species_to_chrs else None
                continue
            if current_sp is not None:
                species_to_chrs[current_sp].append(line)

    # 写出 JCVI 格式 seqids
    with open(seqids_out, "w") as out:
        for sp in species_list:
            chrs = species_to_chrs.get(sp, [])
            for c in chrs:
                out.write(f"{c}\n")

    logging.info("已重写 JCVI 使用的 seqids：%s", seqids_out)


def build_bed_noheader_for_species(
    species: str,
    geneorder_dir: str,
    bed_dir: str,
    color: str = DEFAULT_GENE_COLOR,
) -> None:
    """
    为单个物种构建 bed_noheader：
      输入：Step01 geneorder/{species}.geneorder.bed
        约定四列：chr, start, end, gene_id
      输出：bed_noheader/{species}.bed
        五列：chr, start, end, gene_id, color
    """
    os.makedirs(bed_dir, exist_ok=True)
    geneorder_file = os.path.join(geneorder_dir, f"{species}.geneorder.bed")
    if not os.path.exists(geneorder_file):
        logging.error("geneorder 文件不存在：%s", geneorder_file)
        raise FileNotFoundError(geneorder_file)

    out_bed = os.path.join(bed_dir, f"{species}.bed")
    n = 0
    with open(geneorder_file) as fh_in, open(out_bed, "w") as fh_out:
        for line in fh_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chr_id, start, end, gid = parts[0], parts[1], parts[2], parts[3]
            fh_out.write("\t".join([chr_id, start, end, gid, color]) + "\n")
            n += 1

    logging.info("物种 %s 的 bed_noheader 已写出：%s（共 %d 行）", species, out_bed, n)


def build_all_bed_noheader(species_list: List[str], paths: Dict[str, str]) -> None:
    """
    为所有物种批量构建 bed_noheader。
    """
    bed_dir = paths["BED_DIR"]
    geneorder_dir = paths["GENEORDER_DIR"]
    os.makedirs(bed_dir, exist_ok=True)

    for sp in species_list:
        build_bed_noheader_for_species(sp, geneorder_dir, bed_dir, DEFAULT_GENE_COLOR)


def build_simple_from_anchors(
    ref_short: str,
    species_list: List[str],
    species_short: Dict[str, str],
    anchors_dir: str,
    simple_dir: str,
) -> None:
    """
    从 anchors_for_links 中构建 simple/ 目录，同时修正 score 列格式。

    输入示例（anchors_for_links/Sco__vs__Air.anchors.simple.filtered）可能为：
      Chr1  Sco01g00010.1  Chr3  XM_...  62.0

    JCVI 的 SimpleFile 在读取时会对 score 调用 int()，因此：
      - 将最后一列 "62.0" 转成 "62"
      - 保持前面所有列不变

    输出：simple/{ref_short}__vs__{target_short}.simple
    """
    os.makedirs(simple_dir, exist_ok=True)

    for sp in species_list:
        if sp == REFERENCE_SPECIES:
            continue
        target_short = species_short[sp]
        src = os.path.join(anchors_dir, f"{ref_short}__vs__{target_short}.anchors.simple.filtered")
        dst = os.path.join(simple_dir, f"{ref_short}__vs__{target_short}.simple")

        if not os.path.exists(src):
            logging.warning("缺少 anchors simple：%s，跳过该物种", src)
            continue

        n = 0
        with open(src) as fh_in, open(dst, "w") as fh_out:
            for line in fh_in:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 5:
                    # 结构不完整时直接原样写出
                    fh_out.write(line + "\n")
                    n += 1
                    continue

                # 只修正最后一列 score：把 "62.0" -> "62"
                raw_score = parts[-1]
                try:
                    score_int = int(float(raw_score))
                    parts[-1] = str(score_int)
                except ValueError:
                    # 如果解析失败，保持原样写回，避免破坏结构
                    pass

                fh_out.write("\t".join(parts) + "\n")
                n += 1

        logging.info("simple 文件已写出：%s（源自 %s，共 %d 行）", dst, src, n)


def build_layout(
    species_list: List[str],
    species_short: Dict[str, str],
    ref_species: str,
    paths: Dict[str, str],
) -> None:
    """
    构建 JCVI layout 文件：
      - 每个物种一行 track；
      - 参考物种使用 label_va = top，其余 bottom；
      - edges 区域：参考物种与其他物种之间的 simple 连接。
    """
    layout_file = paths["LAYOUT"]
    bed_dir = paths["BED_DIR"]

    if ref_species not in species_list:
        raise ValueError(f"REFERENCE_SPECIES={ref_species} 不在 seqids 物种列表中")

    ref_index = species_list.index(ref_species)
    ref_short = species_short[ref_species]

    # 轨道纵向位置：从上到下等间距
    n_sp = len(species_list)
    y_start = 0.9125  # 第一条轨道的 y
    y_step = 0.055    # 两条轨道间隔

    with open(layout_file, "w") as fh:
        # 1) tracks 部分
        for i, sp in enumerate(species_list):
            y = y_start - i * y_step
            label = sp.replace("_", " ")
            bed_path = os.path.join(bed_dir, f"{sp}.bed")
            va = "top" if sp == ref_species else "bottom"
            # JCVI layout format:
            # y, xstart, xend, rotation, color, label, va, bed
            fh.write(f"{y:.4f}, 0.05, 0.95, 0, , {label}, {va}, {bed_path}\n")

        # 2) edges 部分
        fh.write("# edges\n")
        for i, sp in enumerate(species_list):
            if sp == ref_species:
                continue
            tgt_short = species_short[sp]
            simple_rel = os.path.join("simple", f"{ref_short}__vs__{tgt_short}.simple")
            fh.write(f"e, {ref_index}, {i}, {simple_rel}\n")

    logging.info("layout 文件构建完成：%s", layout_file)


def write_plot_meta(paths: Dict[str, str], species_list: List[str]) -> None:
    """
    写出一个简单的元数据文件，记录本次出图使用到的关键输入。
    """
    with open(paths["PLOT_META"], "w") as fh:
        fh.write("key\tvalue\n")
        fh.write(f"PROJECT_ROOT\t{PROJECT_ROOT}\n")
        fh.write(f"STEP01_DIR\t{STEP01_DIR}\n")
        fh.write(f"STEP02_DIR\t{STEP02_DIR}\n")
        fh.write(f"STEP03_DIR\t{STEP03_DIR}\n")
        fh.write(f"STEP04_DIR\t{STEP04_DIR}\n")
        fh.write(f"species_order\t{','.join(species_list)}\n")
        fh.write(f"reference_species\t{REFERENCE_SPECIES}\n")
        fh.write(f"seqids_jcvi\t{paths['SEQIDS_JCVI']}\n")
        fh.write(f"layout\t{paths['LAYOUT']}\n")
        fh.write(f"bed_dir\t{paths['BED_DIR']}\n")
        fh.write(f"simple_dir\t{paths['SIMPLE_DIR']}\n")
        fh.write(f"pdf_out\t{paths['PDF_OUT']}\n")

    logging.info("plot_meta_used_files.tsv 已写出：%s", paths["PLOT_META"])


def run_jcvi_karyotype(paths: Dict[str, str]) -> None:
    """
    只调用一次 JCVI 生成 PDF。
    SVG 完全不再尝试。
    """
    out_pdf = paths["PDF_OUT"]
    cmd = [
        sys.executable,
        "-m",
        "jcvi.graphics.karyotype",
        "seqids",
        "layout",
        "--figsize",
        FIGSIZE,
        "--dpi",
        str(DPI),
        "-o",
        os.path.basename(out_pdf),
    ]
    logging.info("运行命令（PDF）：%s", " ".join(cmd))
    rc, out, err = run_cmd(cmd, cwd=paths["OUTPUT_ROOT"])
    if out:
        logging.info("[JCVI STDOUT]\n%s", out)
    if err:
        logging.error("[JCVI STDERR]\n%s", err)
    if rc != 0:
        logging.error("PDF 绘图失败（退出码=%d），请根据上方 JCVI 日志排查 simple / bed / layout / seqids。", rc)
    else:
        logging.info("PDF 绘图完成：%s", out_pdf)


# =========================
# 主流程
# =========================

def main() -> None:
    setup_logger()
    logging.info("========== synteny05 — JCVI karyotype 出图（终版） ==========")

    paths = get_paths()
    out_root = paths["OUTPUT_ROOT"]
    os.makedirs(out_root, exist_ok=True)

    # 1) 备份 Step03 seqids，并从中解析物种顺序
    seqids_src = paths["SEQIDS_SRC"]
    if not os.path.exists(seqids_src):
        logging.error("Step03 seqids 不存在：%s", seqids_src)
        sys.exit(1)

    shutil.copyfile(seqids_src, paths["SEQIDS_BAK"])
    logging.info("已备份 Step03 seqids 至：%s", paths["SEQIDS_BAK"])

    species_list = parse_species_order_from_seqids_fasta(seqids_src)
    logging.info("总计 %d 个物种，参考物种 = %s", len(species_list), REFERENCE_SPECIES)

    # 2) 物种短代号映射
    species_short = load_species_short_labels(paths["META_TSV"], species_list)
    ref_short = species_short[REFERENCE_SPECIES]

    # 3) 重写 JCVI 使用的 seqids
    write_jcvi_seqids(seqids_src, paths["SEQIDS_JCVI"], species_list)

    # 4) 构建 bed_noheader
    build_all_bed_noheader(species_list, paths)

    # 5) 构建 simple/，修正 score=62.0 → 62
    build_simple_from_anchors(
        ref_short=ref_short,
        species_list=species_list,
        species_short=species_short,
        anchors_dir=paths["ANCHORS_LINKS_DIR"],
        simple_dir=paths["SIMPLE_DIR"],
    )

    # 6) 构建 layout
    build_layout(species_list, species_short, REFERENCE_SPECIES, paths)

    # 7) 写 meta
    write_plot_meta(paths, species_list)

    # 8) 调用 JCVI 画 PDF（不再尝试 SVG）
    run_jcvi_karyotype(paths)

    logging.info("synteny05 运行结束，请检查输出目录：%s", paths["OUTPUT_ROOT"])


if __name__ == "__main__":
    main()

