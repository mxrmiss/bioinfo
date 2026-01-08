#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny04_layout_tracks.py
—— 朴素模式 Step04：生成 tracks 布局表（y_center 等），不再强做苛刻校验

输入（硬位置）：
- raw_data/synteny_species_meta.tsv

输出（硬位置）：
output/synteny_04_layout_tracks/
  layout_species_tracks.tsv
  summaries/step04.summary.tsv

说明：
- 物种名颜色不管：本表仍保留 label_color_hex 字段，但 Step05 不再把颜色写进 label 前缀
- chr_bar_color_hex：你可以继续用“参考红、其它灰”
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

REFERENCE_SPECIES_ID = "   "

TRACK_HEIGHT = 0.040
TRACK_SPACING = 0.018

REF_RED_HEX = "#E64B35"
OTHER_BAR_HEX = "#B3B3B3"
OTHER_LABEL_HEX = "#000000"

META_TSV_REL = "raw_data/synteny_species_meta.tsv"
OUTPUT_DIR_REL = "output/synteny_04_layout_tracks"
LOG_DIR_REL = "logs"

# ============================================================
# 实现区（皇上勿改）
# ============================================================

import csv
import gzip
import time
import shutil
import traceback
from pathlib import Path
from typing import List, Tuple, Iterable


def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


class Logger:
    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

    def _write(self, level: str, msg: str) -> None:
        line = f"{now_ts()} [{level}] {msg}"
        print(line, flush=True)
        with self.log_file.open("a", encoding="utf-8") as fw:
            fw.write(line + "\n")

    def info(self, msg: str) -> None:
        self._write("INFO", msg)

    def warn(self, msg: str) -> None:
        self._write("WARN", msg)

    def error(self, msg: str) -> None:
        self._write("ERROR", msg)


def die(logger: Logger, msg: str, code: int = 1) -> None:
    logger.error(msg)
    raise SystemExit(code)


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def clean_dir(path: Path, clean: bool, logger: Logger) -> None:
    if clean and path.exists():
        logger.info(f"Clean output dir: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fw:
        w = csv.writer(fw, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def read_meta(meta_tsv: Path) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with open_text_auto(meta_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != ["species_id", "group"]:
            raise ValueError(f"meta header must be: species_id<TAB>group, got: {header}")
        for ln, line in enumerate(fr, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"meta line {ln} must have 2 columns: {line}")
            sid = parts[0].strip()
            grp = parts[1].strip()
            if not sid:
                raise ValueError(f"meta line {ln} species_id empty")
            rows.append((sid, grp))
    if len(rows) < 2:
        raise ValueError("meta must contain at least 2 species")
    return rows


def parse_display_label(species_id: str) -> str:
    parts = species_id.split("_")
    if len(parts) < 2:
        return species_id
    genus = parts[0]
    sp = parts[1]
    if not genus or not sp:
        return species_id
    return f"{genus[0].upper()}. {sp.lower()}"


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        project = Path(__file__).resolve().parents[1]
    else:
        project = Path(str(pr)).expanduser().resolve()

    meta_tsv = project / META_TSV_REL
    out_dir = project / OUTPUT_DIR_REL
    log_dir = project / LOG_DIR_REL

    logger = Logger(log_dir / "synteny04_layout_tracks.log")
    t0 = time.time()

    logger.info("========== synteny04 — plain layout ==========")
    logger.info(f"PROJECT_ROOT={project}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"REFERENCE_SPECIES_ID={REFERENCE_SPECIES_ID}")

    if not meta_tsv.exists():
        die(logger, f"Missing meta: {meta_tsv}")

    meta_rows = read_meta(meta_tsv)
    species = [sid for sid, _ in meta_rows]
    logger.info("META species order=" + " | ".join(species))

    clean_dir(out_dir, CLEAN_OUTPUT, logger)
    d_sum = out_dir / "summaries"
    d_sum.mkdir(parents=True, exist_ok=True)

    N = len(species)
    top_margin = 0.95
    step = TRACK_HEIGHT + TRACK_SPACING

    layout_rows: List[List[str]] = []
    sum_rows: List[List[str]] = []

    for i, (sid, grp) in enumerate(meta_rows, start=1):
        y_center = top_margin - (i - 1) * step
        display_label = parse_display_label(sid)

        if sid == REFERENCE_SPECIES_ID:
            bar_hex = REF_RED_HEX
            label_hex = REF_RED_HEX
        else:
            bar_hex = OTHER_BAR_HEX
            label_hex = OTHER_LABEL_HEX

        layout_rows.append([
            str(i),
            sid,
            display_label,
            grp,
            f"{y_center:.6f}",
            f"{TRACK_HEIGHT:.6f}",
            bar_hex,
            label_hex,
        ])
        sum_rows.append([sid, str(i), display_label])

    write_tsv(
        out_dir / "layout_species_tracks.tsv",
        ["order_index", "species_id", "display_label", "group", "y_center", "track_height", "chr_bar_color_hex", "label_color_hex"],
        layout_rows
    )
    write_tsv(d_sum / "step04.summary.tsv", ["species_id", "order_index", "display_label"], sum_rows)

    logger.info(f"Done. runtime_sec={int(time.time()-t0)}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        pr = PROJECT_ROOT
        if pr is None or str(pr).strip() == "":
            project = Path(__file__).resolve().parents[1]
        else:
            project = Path(str(pr)).expanduser().resolve()
        lg = Logger(project / LOG_DIR_REL / "synteny04_layout_tracks.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

