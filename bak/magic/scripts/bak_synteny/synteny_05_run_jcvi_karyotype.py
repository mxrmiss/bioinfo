#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny_05_run_jcvi_karyotype.py —— JCVI 线性多基因组共线性堆叠图驱动脚本（seqids 修正版）

职责（v4 版本）：
  1) 读取 layout_species_tracks.tsv，确定物种顺序与轨道位置；
  2) 基于 synteny_02 的带表头 bed，为 JCVI 生成无表头 bed（chr / start / end / gene_id）；
  3) 直接复用 synteny_02 中的 block-level simple（simple_blocks，6 列），
     拷贝到本步骤的 simple/ 目录，供 JCVI 使用；
  4) 正确写出 JCVI 需要的 seqids 文件（核心修正点）：
       - 不再有 “>species” 头行；
       - 不再出现 “Chr5?” 之类带 ? 的名字；
       - 每个物种一行，逗号分隔：Chr1,Chr2,Chr3,...
  5) 写出符合 jcvi.graphics.karyotype 语法的 layout 文件（官方 8 列格式）；
  6) 可选：在脚本内直接调用 JCVI 画图（默认关闭，由皇上手动执行）。

使用方式：
  1) 在 magic 项目根目录下运行本脚本；
  2) 运行完成后：
       cd output/synteny_05_plot
       python -m jcvi.graphics.karyotype seqids layout --format=pdf --shadestyle=line --notex
"""

import sys
import logging
from pathlib import Path
import shutil
import subprocess
import csv

# =========================== 顶部可调参数 ===========================

# 是否在脚本内直接调用 JCVI 作图
RUN_JCVI = False

# 图像输出前缀
OUTPUT_PREFIX = "synteny_linear_stack"

# DPI 设置（仅在 RUN_JCVI=True 时使用）
FIG_DPI = 300

# 日志级别
LOG_LEVEL = "INFO"


# =========================== 日志配置 ===============================


def setup_logger(log_file: Path, level: str = "INFO") -> logging.Logger:
    """
    配置日志：屏幕 + 文件双通道输出。
    """
    logger = logging.getLogger("synteny_05")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # 清掉旧 handler，避免重复输出
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


# =========================== 工具函数 ===============================


def read_species_tracks(layout_species_tracks: Path, logger: logging.Logger):
    """
    读取 layout_species_tracks.tsv，返回：
      - ordered_species: 参与 synteny 的物种 ID 列表（按 order_index 排序）
      - track_info: 每个物种的轨道信息（y_center / plot_label 等）
    要求 layout_species_tracks.tsv 至少包含列：
      - species_id
      - short_label
      - plot_label
      - order_index
      - group
      - y_center
      - track_height
      - include_in_synteny （yes/no）
    """
    species_order = []
    track_info = {}

    with layout_species_tracks.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        required_cols = [
            "species_id",
            "short_label",
            "plot_label",
            "order_index",
            "group",
            "y_center",
            "track_height",
        ]
        missing = [c for c in required_cols if c not in col_idx]
        if missing:
            logger.error(
                "layout_species_tracks.tsv 缺少必要列：%s", ", ".join(missing)
            )
            raise SystemExit(1)

        # include_in_synteny 是后加列，不一定都有：没有就默认都用
        has_include = "include_in_synteny" in col_idx

        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            species_id = parts[col_idx["species_id"]]

            if has_include:
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
    ordered_species = [sid for _, sid in species_order]

    logger.info("layout 中参与 synteny 的物种顺序：%s", ", ".join(ordered_species))
    return ordered_species, track_info


def write_seqids_files(
    seqids_file: Path,
    seqids_species_dir: Path,
    chr_order_dir: Path,
    species_ids: list,
    logger: logging.Logger,
):
    """
    ==== 核心修正点：正确写出 JCVI 的 seqids 文件 ====

    JCVI 期望的 seqids 文件格式（以 3 条轨道为例）：
        Chr1,Chr2,Chr3,Chr4
        Chr1,Chr2,Chr3
        Chr1,Chr2,Chr3,Chr4,Chr5,Chr6

    特点：
      - 每个物种对应一行（顺序要与 layout 中轨道顺序一致）；
      - 行内用逗号分隔；
      - 不要写物种名；
      - 不要写 “Chr1?”、“Chr1+” 之类后缀；
      - 名字要能在 bed 文件的第 1 列中找到（例如 Chr1, Chr2, ...）。

    我们这里的逻辑：
      1) 读取 synteny_03 的 chr_order_<species>.tsv；
      2) 按 rank 排序，收集 chr 列；
      3) 写 seqids_species/<species>.seqids（每行一个 Chr，便于皇上排查）；
      4) 写 seqids（每个物种一行，逗号分隔，不带正负号、不带 ?）。
    """
    seqids_species_dir.mkdir(parents=True, exist_ok=True)

    lines_for_global = []

    for sid in species_ids:
        chr_order_file = chr_order_dir / f"chr_order_{sid}.tsv"
        if not chr_order_file.exists():
            logger.error("找不到 chr_order 文件：%s", chr_order_file)
            raise FileNotFoundError(chr_order_file)

        ranks = []
        with chr_order_file.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col_idx = {name: i for i, name in enumerate(header)}
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                chr_name = parts[col_idx["chr"]]
                rank = int(parts[col_idx["rank"]])
                # orientation 可以用于 debug，但 seqids 不再写 + 或 ?，只保留 Chr 名
                ranks.append((rank, chr_name))

        ranks.sort(key=lambda x: x[0])
        chr_list = [chr_name for _, chr_name in ranks]

        # 写每个物种单独的 seqids，方便皇上检查
        species_seq_file = seqids_species_dir / f"{sid}.seqids"
        with species_seq_file.open("w") as fo:
            for chr_name in chr_list:
                fo.write(f"{chr_name}\n")

        logger.info(
            "物种 %s：seqids 染色体顺序 = %s",
            sid,
            ", ".join(chr_list),
        )

        # 为全局 seqids 收集一行（逗号分隔）
        lines_for_global.append(",".join(chr_list))

    # 写全局 seqids 文件（JCVI 实际读取的就是这个）
    with seqids_file.open("w") as out_all:
        for line in lines_for_global:
            out_all.write(line + "\n")

    logger.info("综合 seqids 文件写出完成：%s", seqids_file)


def build_bed_noheader(
    bed_dir: Path,
    bed_noheader_dir: Path,
    species_ids: list,
    logger: logging.Logger,
):
    """
    从 synteny_02 的带表头 bed 文件生成 JCVI 需要的无表头 bed：
      输入：bed/<species>.bed（含 header）
      输出：bed_noheader/<species>.bed（仅保留 chr, start, end, gene_id）
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
            reader = csv.DictReader(fi, delimiter="\t")
            for row in reader:
                chr_name = row["chr"]
                start = row["start"]
                end = row["end"]
                gene_id = row["gene_id"]
                fo.write(f"{chr_name}\t{start}\t{end}\t{gene_id}\n")
                n += 1

        logger.info(
            "物种 %s：生成无表头 bed 文件 %s，记录数 = %d",
            sid,
            out_bed,
            n,
        )


def use_block_simple_files(
    block_simple_src_dir: Path,
    simple_dir: Path,
    species_ids: list,
    logger: logging.Logger,
):
    """
    复用 synteny_02 中已经聚类好的区段级 simple_blocks：

      输入：synteny_02_pairwise_bed/simple_blocks/<A>__vs__<B>.anchors.simple
              （6 列 block simple）
      输出：当前输出目录 simple/<A>__vs__<B>.simple
              （内容完全拷贝，仅文件名后缀不同）

    注意：
      - 物种对仍然采用链式：
          (S1, S2), (S2, S3), ..., (S_{n-1}, S_n)
      - 这一步只做复制 + 简单统计，不改内容。
    """
    simple_dir.mkdir(parents=True, exist_ok=True)

    pairs = []
    for i in range(len(species_ids) - 1):
        a = species_ids[i]
        b = species_ids[i + 1]
        pairs.append((a, b))

    for a, b in pairs:
        src = block_simple_src_dir / f"{a}__vs__{b}.anchors.simple"
        if not src.exists():
            logger.error("找不到 block-level simple 文件：%s", src)
            raise FileNotFoundError(src)

        dst = simple_dir / f"{a}__vs__{b}.simple"
        shutil.copyfile(src, dst)

        # 简单统计一下 block 数量（行数）
        n_blocks = 0
        with dst.open() as fh:
            for line in fh:
                if line.strip():
                    n_blocks += 1

        logger.info(
            "物种对 %s vs %s：复用 block-level simple -> %s，blocks = %d",
            a,
            b,
            dst,
            n_blocks,
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
    写出 jcvi.graphics.karyotype 所需的 layout 配置文件。

    官方格式（参考 jcvi 文档）：
      # y, xstart, xend, rotation, color, label, va, bed
      .7, .1, .8, 0, , ZM, top, ZM.bed
      .5, .1, .8, 0, , SB, top, SB.bed
      # edges
      e, 0, 1, ZM.SB.anchors.simple

    我们这里：
      - y 使用 layout_species_tracks.tsv 中的 y_center，经归一化到 [0,1]；
      - xstart/xend 统一用 0.1 / 0.9；
      - label 使用 plot_label；
      - bed 使用相对路径：bed_noheader/<species>.bed；
      - edges 按轨道 index 链式连接，对应 simple/<A>__vs__<B>.simple。
    """
    with layout_file.open("w") as fh:
        fh.write("# karyotype layout for synteny linear stack\n")
        fh.write("# y, xstart, xend, rotation, color, label, va, bed\n")

        # 归一化 y，使不同物种的纵向间距大致合理
        y_vals = [track_info[sid]["y_center"] for sid in species_ids]
        y_min = min(y_vals)
        y_max = max(y_vals)
        y_range = y_max - y_min if y_max > y_min else 1.0

        for idx, sid in enumerate(species_ids):
            info = track_info[sid]
            y = info["y_center"]
            y_norm = (y - y_min) / y_range
            xstart = 0.10
            xend = 0.90
            rotation = 0
            color = ""  # 颜色后续由 karyotype 内部/主题决定
            label = info["plot_label"] or sid
            va = "top"
            bed_rel = f"{bed_noheader_dir.name}/{sid}.bed"

            fh.write(
                f"{y_norm:.4f}, {xstart:.3f}, {xend:.3f}, "
                f"{rotation}, {color}, {label}, {va}, {bed_rel}\n"
            )

        fh.write("# edges\n")

        # 按轨道 index 链式连接：0-1, 1-2, ..., (n-2)-(n-1)
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


# =========================== 主流程 ===============================


def main():
    # 假定脚本位于 magic/scripts/ 下
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "synteny_05_plot"

    # 皇上要求：每次运行 05 前清空自己的输出目录
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
    block_simple_src_dir = (
        project_root / "output" / "synteny_02_pairwise_bed" / "simple_blocks"
    )

    bed_noheader_dir = output_dir / "bed_noheader"
    seqids_species_dir = output_dir / "seqids_species"
    simple_dir = output_dir / "simple"
    layout_file = output_dir / "layout"
    seqids_file = output_dir / "seqids"

    if not layout_species_tracks.exists():
        logger.error("找不到 layout_species_tracks.tsv：%s", layout_species_tracks)
        raise FileNotFoundError(layout_species_tracks)

    # 1) 读取物种顺序 + 轨道信息
    species_ids, track_info = read_species_tracks(layout_species_tracks, logger)

    # 2) 写 seqids（修正部分）
    write_seqids_files(
        seqids_file,
        seqids_species_dir,
        chr_order_dir,
        species_ids,
        logger,
    )

    # 3) 生成无表头 bed
    build_bed_noheader(bed_dir, bed_noheader_dir, species_ids, logger)

    # 4) 复用 synteny_02 的 block-level simple
    use_block_simple_files(block_simple_src_dir, simple_dir, species_ids, logger)

    # 5) 写 layout（官方 8 列格式）
    write_layout_file(
        layout_file, bed_noheader_dir, species_ids, track_info, simple_dir, logger
    )

    # 6) 记录本次作图所用文件
    plot_meta = output_dir / "plot_meta_used_files.tsv"
    with plot_meta.open("w") as pm:
        pm.write("file_role\tpath\n")
        pm.write(f"layout_species_tracks\t{layout_species_tracks}\n")
        pm.write(f"chr_order_dir\t{chr_order_dir}\n")
        pm.write(f"bed_dir\t{bed_dir}\n")
        pm.write(f"anchors_dir\t{anchors_dir}\n")
        pm.write(f"simple_blocks_src_dir\t{block_simple_src_dir}\n")
        pm.write(f"bed_noheader_dir\t{bed_noheader_dir}\n")
        pm.write(f"simple_dir\t{simple_dir}\n")
        pm.write(f"seqids_species_dir\t{seqids_species_dir}\n")
        pm.write(f"seqids_file\t{seqids_file}\n")
        pm.write(f"layout_file\t{layout_file}\n")

    # 7) 可选：在脚本内直接调用 JCVI
    if RUN_JCVI:
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
        logger.info(
            "  python -m jcvi.graphics.karyotype seqids layout "
            "--format=pdf --shadestyle=line --notex"
        )

    logger.info("========== synteny_05_run_jcvi_karyotype 结束 ==========")


if __name__ == "__main__":
    main()

