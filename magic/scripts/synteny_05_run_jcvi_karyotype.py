#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_05_run_jcvi_karyotype.py —— JCVI 线性多基因组共线性堆叠图驱动脚本

职责：
  1) 读取 layout_species_tracks.tsv，确定物种顺序与轨道位置；
  2) 基于 02 步的 bed（带表头），为 JCVI 生成无表头 bed；
  3) 基于 02 步的 anchors.simple（基因对），生成 JCVI 需要的 6 列 simple 文件；
  4) 写出符合 jcvi.graphics.karyotype 语法的 seqids 和 layout 文件；
  5) 可选：在脚本内直接调用 JCVI 画图（默认关闭，由皇上手动执行）。
"""

import sys
import logging
from pathlib import Path
import shutil
import subprocess

# =========================== 顶部可调参数 ===========================

# 是否在脚本内直接调用 JCVI 作图
RUN_JCVI = False

# 图像输出前缀
OUTPUT_PREFIX = "synteny_linear_stack"

# DP I 设置（仅在 RUN_JCVI=True 时使用）
FIG_DPI = 300

# 日志级别
LOG_LEVEL = "INFO"

# =========================== 日志配置 ===============================


def setup_logger(log_file: Path, level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger("synteny_05")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    fh = logging.FileHandler(log_file)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# =========================== 核心函数 ===============================


def read_species_tracks(layout_species_tracks: Path, logger: logging.Logger):
    """
    读取 layout_species_tracks.tsv，返回物种顺序列表和每个物种的轨道参数字典。
    """
    species_order = []
    track_info = {}

    with layout_species_tracks.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            species_id = parts[col_idx["species_id"]]
            include = parts[col_idx["include_in_synteny"]].strip().lower()
            if include not in ("yes", "true", "1"):
                continue

            short_label = parts[col_idx["short_label"]]
            plot_label = parts[col_idx["plot_label"]]
            order_index = int(parts[col_idx["order_index"]])
            group = parts[col_idx["group"]]
            y_center = float(parts[col_idx["y_center"]])
            track_height = float(parts[col_idx["track_height"]])

            species_order.append((order_index, species_id))
            track_info[species_id] = {
                "short_label": short_label,
                "plot_label": plot_label,
                "group": group,
                "y_center": y_center,
                "track_height": track_height,
            }

    species_order.sort(key=lambda x: x[0])
    ordered_species_ids = [sid for _, sid in species_order]

    logger.info(
        "参与 synteny 的物种数量 = %d", len(ordered_species_ids)
    )
    logger.info(
        "物种顺序为：%s", " -> ".join(ordered_species_ids)
    )

    return ordered_species_ids, track_info


def write_seqids_file(
    seqids_dir: Path,
    chr_order_dir: Path,
    species_ids: list,
    logger: logging.Logger,
):
    """
    为每个物种写出 seqids/<species>.seqids，并汇总到 karyotype 的 seqids 文件。
    """
    seqids_dir.mkdir(parents=True, exist_ok=True)
    karyotype_seqids = seqids_dir.parent / "seqids"

    with karyotype_seqids.open("w") as out_all:
        for sid in species_ids:
            chr_order_file = chr_order_dir / f"chr_order_{sid}.tsv"
            if not chr_order_file.exists():
                logger.error("找不到 chr_order 文件：%s", chr_order_file)
                raise FileNotFoundError(chr_order_file)

            species_seqids_file = seqids_dir / f"{sid}.seqids"
            ranks = []
            with chr_order_file.open() as fh:
                header = fh.readline().rstrip("\n").split("\t")
                col_idx = {name: i for i, name in enumerate(header)}
                for line in fh:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    chr_name = parts[col_idx["chr"]]
                    orientation = parts[col_idx["orientation"]]
                    rank = int(parts[col_idx["rank"]])
                    ranks.append((rank, chr_name, orientation))
            ranks.sort(key=lambda x: x[0])

            with species_seqids_file.open("w") as fh2:
                for _, chr_name, orientation in ranks:
                    fh2.write(f"{chr_name}{orientation}\n")

            out_all.write(f"{sid}\t{species_seqids_file.name}\n")

    logger.info("seqids 写入完成，共 %d 行（物种）", len(species_ids))


def build_bed_noheader(
    bed_dir: Path,
    bed_noheader_dir: Path,
    species_ids: list,
    logger: logging.Logger,
):
    """
    从 02 步的 bed（带表头）生成 JCVI 需要的“无表头” bed 文件。
    JCVI 的 BedLine 只认：
      col1 = chr
      col2 = start
      col3 = end
      col4 = name
    其他列不用写。
    """
    bed_noheader_dir.mkdir(parents=True, exist_ok=True)

    for sid in species_ids:
        in_bed = bed_dir / f"{sid}.bed"
        if not in_bed.exists():
            logger.error("找不到 bed 文件：%s", in_bed)
            raise FileNotFoundError(in_bed)

        out_bed = bed_noheader_dir / f"{sid}.bed"
        n = 0
        with in_bed.open() as fi, out_bed.open("w") as fo:
            header = fi.readline()
            for line in fi:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                chr_name = parts[0]
                start = parts[1]
                end = parts[2]
                gene_id = parts[3]
                fo.write(f"{chr_name}\t{start}\t{end}\t{gene_id}\n")
                n += 1

        logger.info(
            "物种 %s：生成无表头 bed 文件 %s，记录数 = %d",
            sid,
            out_bed,
            n,
        )


def build_simple_files(
    anchors_dir: Path,
    simple_dir: Path,
    species_ids: list,
    logger: logging.Logger,
):
    """
    从 02 步的 anchors.simple（两列基因配对）生成 JCVI 需要的 6 列 simple 文件。
    暂时采用简化策略：每个基因配对视为一个“最小 synteny block”：
      geneA geneA geneB geneB 1 +
    后续若需要可以升级为真正的区段聚类。
    """
    simple_dir.mkdir(parents=True, exist_ok=True)

    pairs = []
    for i in range(len(species_ids) - 1):
        a = species_ids[i]
        b = species_ids[i + 1]
        pairs.append((a, b))

    for a, b in pairs:
        anchors_file = anchors_dir / f"{a}__vs__{b}.anchors.simple"
        if not anchors_file.exists():
            logger.error("找不到 anchors.simple 文件：%s", anchors_file)
            raise FileNotFoundError(anchors_file)

        simple_file = simple_dir / f"{a}__vs__{b}.simple"
        n = 0
        with anchors_file.open() as fi, simple_file.open("w") as fo:
            for line in fi:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split()
                if len(parts) < 2:
                    continue
                gene_a, gene_b = parts[0], parts[1]
                # 6 列：startA, endA, startB, endB, score, orientation
                fo.write(f"{gene_a}\t{gene_a}\t{gene_b}\t{gene_b}\t1\t+\n")
                n += 1

        logger.info(
            "物种对 %s__vs__%s：生成 simple 文件 %s，block 数 = %d",
            a,
            b,
            simple_file,
            n,
        )


def write_layout_file(
    layout_file: Path,
    bed_noheader_dir: Path,
    species_ids: list,
    track_info: dict,
    simple_dir: Path,
    logger: logging.Logger,
):
    """
    写出 jcvi.graphics.karyotype 使用的 layout 文件。
    轨道行：
      y, xstart, xend, rotation, color, label, va, bed
    edges 行：
      e, track_idx_A, track_idx_B, simple_file
    """
    with layout_file.open("w") as fh:
        fh.write("# y, xstart, xend, rotation, color, label, va, bed\n")

        y_min = min(track_info[sid]["y_center"] for sid in species_ids)
        y_max = max(track_info[sid]["y_center"] for sid in species_ids)
        y_range = y_max - y_min if y_max > y_min else 1.0

        for idx, sid in enumerate(species_ids):
            info = track_info[sid]
            y = info["y_center"]
            y_norm = (y - y_min) / y_range
            xstart = 0.1
            xend = 0.9
            rotation = 0
            color = ""
            label = info["plot_label"] or sid
            va = "top"
            bed_rel = f"{bed_noheader_dir.name}/{sid}.bed"
            fh.write(
                f"{y_norm:.4f}, {xstart:.3f}, {xend:.3f}, {rotation}, {color}, {label}, {va}, {bed_rel}\n"
            )

        fh.write("# edges\n")

        for idx in range(len(species_ids) - 1):
            a = species_ids[idx]
            b = species_ids[idx + 1]
            simple_rel = f"{simple_dir.name}/{a}__vs__{b}.simple"
            fh.write(f"e, {idx}, {idx + 1}, {simple_rel}\n")

    logger.info(
        "layout 写入完成：包含 %d 条轨道，%d 条 edges",
        len(species_ids),
        len(species_ids) - 1,
    )


def main():
    project_root = Path(__file__).resolve().parents[1]
    output_dir = project_root / "output" / "synteny_05_plot"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_file = output_dir / "synteny_05_run_jcvi_karyotype.log"
    logger = setup_logger(log_file, LOG_LEVEL)

    logger.info("========== Synteny 05 — JCVI karyotype 驱动脚本启动 ==========")
    logger.info("项目根目录：%s", project_root)

    layout_species_tracks = (
        project_root / "output" / "synteny_04_layout" / "layout_species_tracks.tsv"
    )
    chr_order_dir = project_root / "output" / "synteny_03_chr_order"
    bed_dir = project_root / "output" / "synteny_02_pairwise_bed" / "bed"
    anchors_dir = project_root / "output" / "synteny_02_pairwise_bed" / "anchors"

    bed_noheader_dir = output_dir / "bed_noheader"
    seqids_dir = output_dir / "seqids_dir"
    simple_dir = output_dir / "simple"
    layout_file = output_dir / "layout"

    if not layout_species_tracks.exists():
        raise FileNotFoundError(layout_species_tracks)

    species_ids, track_info = read_species_tracks(
        layout_species_tracks, logger
    )

    write_seqids_file(seqids_dir, chr_order_dir, species_ids, logger)
    build_bed_noheader(bed_dir, bed_noheader_dir, species_ids, logger)
    build_simple_files(anchors_dir, simple_dir, species_ids, logger)
    write_layout_file(
        layout_file, bed_noheader_dir, species_ids, track_info, simple_dir, logger
    )

    plot_meta = output_dir / "plot_meta_used_files.tsv"
    with plot_meta.open("w") as pm:
        pm.write("file_role\tpath\n")
        pm.write(f"layout_species_tracks\t{layout_species_tracks}\n")
        pm.write(f"chr_order_dir\t{chr_order_dir}\n")
        pm.write(f"bed_dir\t{bed_dir}\n")
        pm.write(f"anchors_dir\t{anchors_dir}\n")
        pm.write(f"bed_noheader_dir\t{bed_noheader_dir}\n")
        pm.write(f"simple_dir\t{simple_dir}\n")
        pm.write(f"seqids_file\t{(seqids_dir.parent / 'seqids')}\n")
        pm.write(f"layout_file\t{layout_file}\n")

    logger.info("写出绘图元信息文件：%s", plot_meta)

    if RUN_JCVI:
        logger.info("尝试在脚本内直接调用 JCVI 作图 ...")
        cmd = [
            sys.executable,
            "-m",
            "jcvi.graphics.karyotype",
            "seqids",
            "layout",
            "--format=pdf",
            "--shadestyle=line",
            "--notex",
        ]
        logger.info("调用命令：%s", " ".join(cmd))
        subprocess.run(cmd, cwd=output_dir, check=True)
        logger.info("JCVI 作图完成，输出位于：%s", output_dir)
    else:
        logger.info("RUN_JCVI=False，本次仅生成 seqids + 无表头 bed + simple + layout，不自动作图。")
        logger.info("皇上可在有 JCVI 的环境中手动执行：")
        logger.info("  cd %s", output_dir)
        logger.info("  python -m jcvi.graphics.karyotype seqids layout --format=pdf --shadestyle=line --notex")

    logger.info("========== synteny_05_run_jcvi_karyotype 结束 ==========")


if __name__ == "__main__":
    main()

