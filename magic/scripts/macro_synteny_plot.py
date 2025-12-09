#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
宏观染色体共线性 - 步骤3：绘制线性多物种宏观共线性图

输入：
- meta/species_order.tsv
- meta/ref_chr_color.tsv
- blocks/blocks.tsv
- links/links.tsv

输出：
- figures/macro_synteny_16species.pdf/png/svg

特点：
- 顶部一行：参考物种染色体条带 + 编号；
- 下方多行：各物种的彩色窗口条带；
- 参考行 → 各物种行之间绘制彩色丝带（按参考染色体颜色）；
- 行背景按 digging / nondigging 上色。
"""

import sys
import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MplPath

# ====================== 用户参数区 ======================

REFERENCE_SPECIES = "Sinonovacula_constricta"

CHR_LENGTH_THRESHOLD_MBP = 10
WINDOW_SIZE_MBP = 5
MAX_LINKS_PER_CHR = 2000

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

BACKGROUND_COLOR_DIGGING = "#E7F5E6"     # 挖掘组背景（浅绿）
BACKGROUND_COLOR_NONDIGGING = "#FBE5EC"  # 非挖掘组背景（浅粉）
BACKGROUND_COLOR_REF = "#F2F2F2"         # 参考物种行背景


# ====================== 通用函数 ======================


def setup_paths() -> dict:
    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]

    output_root = project_root / "output" / "macro_synteny"
    log_dir = output_root / "logs"
    meta_dir = output_root / "meta"
    blocks_dir = output_root / "blocks"
    links_dir = output_root / "links"
    figures_dir = output_root / "figures"

    for d in [output_root, log_dir, meta_dir, blocks_dir, links_dir, figures_dir]:
        d.mkdir(parents=True, exist_ok=True)

    return {
        "project_root": project_root,
        "output_root": output_root,
        "log_dir": log_dir,
        "meta_dir": meta_dir,
        "blocks_dir": blocks_dir,
        "links_dir": links_dir,
        "figures_dir": figures_dir,
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


def load_species_order(path: Path, logger: logging.Logger):
    """
    读取 species_order.tsv，返回：
    - ordered：按 row_index 排序的列表
    - prefix2info：prefix -> {row_index, species_label, group}
    """
    ordered = []
    prefix2info = {}

    if not path.exists():
        logger.error(f"缺少 species_order.tsv：{path}")
        sys.exit(1)

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row_index = int(row["row_index"])
            prefix = row["species_prefix"]
            label = row["species_label"]
            group = row["group"]
            rec = {
                "row_index": row_index,
                "species_prefix": prefix,
                "species_label": label,
                "group": group,
            }
            ordered.append(rec)
            prefix2info[prefix] = rec

    ordered.sort(key=lambda x: x["row_index"])
    logger.info(f"已读取物种行顺序，共 {len(ordered)} 行。")
    return ordered, prefix2info


def load_ref_chr_color(path: Path, logger: logging.Logger):
    """
    读取 ref_chr_color.tsv，返回：
    - ref_rows：按 ref_chr_rank 排序列表
    - total_x_bp：参考基因组总长度（bp）
    """
    if not path.exists():
        logger.error(f"缺少 ref_chr_color.tsv：{path}")
        sys.exit(1)

    ref_rows = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ref_rows.append({
                "ref_chr_id": row["ref_chr_id"],
                "ref_chr_rank": int(row["ref_chr_rank"]),
                "length_bp": int(row["length_bp"]),
                "x_start": float(row["x_start"]),
                "x_end": float(row["x_end"]),
                "color_id": int(row["color_id"]),
                "color_hex": row["color_hex"],
            })

    ref_rows.sort(key=lambda x: x["ref_chr_rank"])
    total_x_bp = max(r["x_end"] for r in ref_rows) if ref_rows else 0
    logger.info(f"参考染色体信息读取完毕，chr数={len(ref_rows)}, total_length_bp={total_x_bp}")
    return ref_rows, total_x_bp


def load_blocks(path: Path, logger: logging.Logger):
    """
    读取 blocks.tsv，按 species_prefix 分组。
    """
    blocks_by_species = {}
    if not path.exists():
        logger.error(f"缺少 blocks.tsv：{path}")
        sys.exit(1)

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sp = row["species_prefix"]
            rec = {
                "species_prefix": sp,
                "species_label": row["species_label"],
                "group": row["group"],
                "row_index": int(row["row_index"]),
                "ref_chr_id": row["ref_chr_id"],
                "ref_chr_rank": int(row["ref_chr_rank"]),
                "x_start": float(row["x_start"]),
                "x_end": float(row["x_end"]),
                "color_id": int(row["color_id"]),
                "n_anchors": int(row["n_anchors"]),
            }
            blocks_by_species.setdefault(sp, []).append(rec)

    for sp in blocks_by_species:
        blocks_by_species[sp].sort(key=lambda r: (r["row_index"], r["ref_chr_rank"], r["x_start"]))

    logger.info(f"blocks.tsv 已加载，共 {len(blocks_by_species)} 个物种有 blocks 信息。")
    return blocks_by_species


def load_links(path: Path, logger: logging.Logger):
    """
    读取 links.tsv，按 species_prefix 分组。
    """
    links_by_species = {}
    if not path.exists():
        logger.error(f"缺少 links.tsv：{path}")
        sys.exit(1)

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sp = row["species_prefix"]
            rec = {
                "species_prefix": sp,
                "species_label": row["species_label"],
                "group": row["group"],
                "row_index": int(row["row_index"]),
                "ref_chr_id": row["ref_chr_id"],
                "ref_chr_rank": int(row["ref_chr_rank"]),
                "ref_x": float(row["ref_x"]),
                "sp_chr_id": row["sp_chr_id"],
                "sp_chr_rank": int(row["sp_chr_rank"]),
                "sp_local_y": float(row["sp_local_y"]),
                "color_id": int(row["color_id"]),
            }
            links_by_species.setdefault(sp, []).append(rec)

    logger.info(f"links.tsv 已加载，共 {len(links_by_species)} 个物种有 links 信息。")
    return links_by_species


# ====================== 主绘图流程 ======================


def main():
    paths = setup_paths()
    logger = setup_logger(paths["log_dir"], "macro_synteny_plot")

    logger.info("=== 宏观共线性步骤3：绘图开始 ===")
    logger.info(f"参考物种：{REFERENCE_SPECIES}")
    logger.info(f"窗口大小：{WINDOW_SIZE_MBP} Mb")

    species_order_path = paths["meta_dir"] / "species_order.tsv"
    ref_chr_color_path = paths["meta_dir"] / "ref_chr_color.tsv"
    blocks_path = paths["blocks_dir"] / "blocks.tsv"
    links_path = paths["links_dir"] / "links.tsv"

    ordered_species, prefix2info = load_species_order(species_order_path, logger)
    ref_rows, total_x_bp = load_ref_chr_color(ref_chr_color_path, logger)
    blocks_by_species = load_blocks(blocks_path, logger)
    links_by_species = load_links(links_path, logger)

    if total_x_bp <= 0:
        logger.error("参考染色体总长度为 0，无法绘图。")
        sys.exit(1)

    # X 轴单位换成 Mb，便于控制宽度
    scale = 1.0 / 1_000_000.0
    total_x_mb = total_x_bp * scale
    logger.info(f"X 轴总长度约为 {total_x_mb:.2f} Mb")

    # 行布局：一行参考物种 + 每个物种一行
    n_rows = len(ordered_species)
    row_height = 1.0
    top_margin = 1.5
    bottom_margin = 1.0
    total_height = top_margin + n_rows * row_height + bottom_margin

    fig_width = 8.0
    fig_height = max(4.0, total_height * 0.35)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_xlim(0, total_x_mb)
    ax.set_ylim(-bottom_margin, top_margin + n_rows * row_height)
    ax.set_facecolor("white")

    # 1. 画行背景（参考行 + 各物种行）
    for rec in ordered_species:
        row_index = rec["row_index"]
        prefix = rec["species_prefix"]
        group = rec["group"]

        if row_index == 0:
            # 参考物种行背景（画在最上方行）
            y0 = top_margin + n_rows * row_height
            y1 = y0 + 0.8 * row_height
            ax.add_patch(Rectangle(
                (0, y0), total_x_mb, y1 - y0,
                facecolor=BACKGROUND_COLOR_REF, edgecolor="none", zorder=0
            ))
        else:
            # 其它物种行：从上往下排
            y_center = top_margin + (n_rows - row_index) * row_height
            y0 = y_center - 0.4 * row_height
            y1 = y_center + 0.4 * row_height

            if group == "digging":
                bg = BACKGROUND_COLOR_DIGGING
            elif group == "nondigging":
                bg = BACKGROUND_COLOR_NONDIGGING
            else:
                bg = "#FFFFFF"

            ax.add_patch(Rectangle(
                (0, y0), total_x_mb, y1 - y0,
                facecolor=bg, edgecolor="none", zorder=0
            ))

    # 2. 画参考物种染色体条带 + 编号
    y_ref_center = top_margin + n_rows * row_height + 0.4 * row_height
    y_ref_height = 0.4 * row_height
    y_ref0 = y_ref_center - y_ref_height / 2
    y_ref1 = y_ref_center + y_ref_height / 2

    for row in ref_rows:
        x0_mb = row["x_start"] * scale
        x1_mb = row["x_end"] * scale
        color = row["color_hex"]
        ax.add_patch(Rectangle(
            (x0_mb, y_ref0), x1_mb - x0_mb, y_ref_height,
            facecolor=color, edgecolor="none", zorder=2
        ))
        x_mid = (x0_mb + x1_mb) / 2
        ax.text(
            x_mid, y_ref1 + 0.1 * row_height,
            str(row["ref_chr_rank"]),
            ha="center", va="bottom", fontsize=7
        )

    # 3. 各物种的 blocks 条带
    for sp_prefix, blist in blocks_by_species.items():
        if sp_prefix not in prefix2info:
            continue
        row_index = prefix2info[sp_prefix]["row_index"]
        if row_index == 0:
            continue

        y_center = top_margin + (n_rows - row_index) * row_height
        bar_height = 0.3 * row_height
        y0 = y_center - bar_height / 2

        for block in blist:
            x0_mb = block["x_start"] * scale
            x1_mb = block["x_end"] * scale
            ref_rank = block["ref_chr_rank"]

            # 根据 ref_chr_rank 找对应颜色
            color_hex = None
            for rr in ref_rows:
                if rr["ref_chr_rank"] == ref_rank:
                    color_hex = rr["color_hex"]
                    break
            if color_hex is None:
                color_hex = "#CCCCCC"

            ax.add_patch(Rectangle(
                (x0_mb, y0), x1_mb - x0_mb, bar_height,
                facecolor=color_hex, edgecolor="none", zorder=2
            ))

    # 4. 画 links（参考 -> 各物种）
    for sp_prefix, llist in links_by_species.items():
        if sp_prefix not in prefix2info:
            continue
        row_index = prefix2info[sp_prefix]["row_index"]
        if row_index == 0:
            continue

        y_target_center = top_margin + (n_rows - row_index) * row_height
        for rec in llist:
            ref_x_mb = rec["ref_x"] * scale
            sp_local_y = rec["sp_local_y"]   # 0~1 之间
            y_source = y_ref_center
            # 映射到该行内部区间
            y_target = y_target_center - 0.2 * row_height + sp_local_y * 0.4 * row_height
            ctrl_y = (y_source + y_target) / 2

            verts = [
                (ref_x_mb, y_source),
                (ref_x_mb, ctrl_y),
                (ref_x_mb, y_target),
            ]
            codes = [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]
            path = MplPath(verts, codes)

            color_hex = None
            for rr in ref_rows:
                if rr["ref_chr_rank"] == rec["ref_chr_rank"]:
                    color_hex = rr["color_hex"]
                    break
            if color_hex is None:
                color_hex = "#BBBBBB"

            patch = PathPatch(
                path,
                facecolor="none",
                edgecolor=color_hex,
                linewidth=0.3,
                alpha=0.15,
                zorder=1,
            )
            ax.add_patch(patch)

    # 5. 左侧物种标签
    for rec in ordered_species:
        row_index = rec["row_index"]
        label = rec["species_label"]
        if row_index == 0:
            y_label = y_ref_center
        else:
            y_label = top_margin + (n_rows - row_index) * row_height
        ax.text(
            -0.02 * total_x_mb,
            y_label,
            label,
            ha="right",
            va="center",
            fontsize=8,
        )

    # 6. 轴样式优化
    ax.set_xlabel("Genomic position (Mb) along reference chromosomes", fontsize=8)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()

    # 7. 保存图像
    fig_basename = "macro_synteny_16species"
    pdf_path = paths["figures_dir"] / f"{fig_basename}.pdf"
    png_path = paths["figures_dir"] / f"{fig_basename}.png"
    svg_path = paths["figures_dir"] / f"{fig_basename}.svg"

    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path, dpi=300)

    logger.info(f"已保存图像：{pdf_path}")
    logger.info(f"已保存图像：{png_path}")
    logger.info(f"已保存图像：{svg_path}")
    logger.info("=== 宏观共线性步骤3 完成 ===")


if __name__ == "__main__":
    main()

