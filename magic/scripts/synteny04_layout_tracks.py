#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny04_layout_tracks.py
—— Step04：轨道布局 + 横条/标签颜色表（严格按蓝图合同 + 加硬防线）

输入（硬位置）：
- raw_data/synteny_species_meta.tsv
- output/synteny_03_order_seqids/seqids_species/<species>.seqids   （用于校验非空 + 合法性）

输出（硬位置）：
output/synteny_04_layout_tracks/
  layout_species_tracks.tsv
  summaries/step04.summary.tsv

layout_species_tracks.tsv（硬表头）：
  order_index
  species_id
  display_label
  group
  y_center
  track_height
  chr_bar_color_hex
  label_color_hex

硬合同要点：
- 顺序严格按 meta 行顺序（order_index=1..N）
- display_label：Genus_species_xxx -> "G. species"
- 参考物种：chr_bar_color_hex 与 label_color_hex = 红
- 其它物种：chr_bar_color_hex=浅灰，label_color_hex=黑
- 防线：seqids token 必须是 ChrNN 或 ChrNN-；严禁 ChrNNr；发现即报错退出
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

REFERENCE_SPECIES_ID = "Sinonovacula_constricta"

# 轨道几何（固定默认写在参数区，允许皇上手动改）
TRACK_HEIGHT = 0.040
TRACK_SPACING = 0.018

# 默认颜色（按蓝图）
REF_RED_HEX = "#E64B35"
OTHER_BAR_HEX = "#D9D9D9"
OTHER_LABEL_HEX = "#000000"

META_TSV_REL = "raw_data/synteny_species_meta.tsv"
STEP03_DIR_REL = "output/synteny_03_order_seqids"
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
import re
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


_SEQ_TOKEN_RE = re.compile(r"^Chr\d+(?:-)?$")


def validate_seqids_line(sid: str, line: str, logger: Logger) -> None:
    """
    防线：seqids 单行逗号分隔；每个 token 必须是 ChrNN 或 ChrNN-
    严禁 ChrNNr；严禁空 token；严禁重复 base Chr。
    """
    raw = (line or "").strip()
    if not raw:
        raise ValueError(f"[{sid}] seqids is empty")

    tokens = [t.strip() for t in raw.split(",")]
    if any(t == "" for t in tokens):
        raise ValueError(f"[{sid}] seqids contains empty token (check extra commas): {raw}")

    base_seen = set()
    for t in tokens:
        if t.endswith("r") or t.endswith("R"):
            raise ValueError(f"[{sid}] seqids token uses 'r' (FORBIDDEN). Use '-' for flip: {t}")

        if not _SEQ_TOKEN_RE.match(t):
            raise ValueError(f"[{sid}] illegal seqids token: {t}  (must be ChrNN or ChrNN-)")

        base = t[:-1] if t.endswith("-") else t
        if base in base_seen:
            raise ValueError(f"[{sid}] duplicated chromosome in seqids: {base}")
        base_seen.add(base)

    logger.info(f"[SEQIDS-OK] {sid}: n_tokens={len(tokens)} flip_tokens={sum(1 for x in tokens if x.endswith('-'))}")


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        PROJECT = Path(__file__).resolve().parents[1]
    else:
        PROJECT = Path(str(pr)).expanduser().resolve()

    meta_tsv = PROJECT / META_TSV_REL
    step03_dir = PROJECT / STEP03_DIR_REL
    out_dir = PROJECT / OUTPUT_DIR_REL
    log_dir = PROJECT / LOG_DIR_REL

    logger = Logger(log_dir / "synteny04_layout_tracks.log")
    t0 = time.time()

    logger.info("========== synteny04 — layout tracks + bar/label colors ==========")
    logger.info(f"PROJECT_ROOT={PROJECT}")
    logger.info(f"OUTPUT_DIR={out_dir}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"REFERENCE_SPECIES_ID={REFERENCE_SPECIES_ID}")
    logger.info(f"TRACK_HEIGHT={TRACK_HEIGHT}")
    logger.info(f"TRACK_SPACING={TRACK_SPACING}")
    logger.info(f"REF_RED_HEX={REF_RED_HEX}")
    logger.info(f"OTHER_BAR_HEX={OTHER_BAR_HEX}")
    logger.info(f"OTHER_LABEL_HEX={OTHER_LABEL_HEX}")
    logger.info(f"META={meta_tsv}")
    logger.info(f"STEP03_DIR={step03_dir}")
    logger.info("Seqids defense: allow ChrNN / ChrNN- ; forbid ChrNNr")

    if not meta_tsv.exists():
        die(logger, f"Missing meta: {meta_tsv}")
    if not step03_dir.exists():
        die(logger, f"Missing step03 dir: {step03_dir}")

    meta_rows = read_meta(meta_tsv)
    species = [sid for sid, grp in meta_rows]
    logger.info("META species order=" + " | ".join(species))

    if REFERENCE_SPECIES_ID not in species:
        die(logger, f"REFERENCE_SPECIES_ID not in meta: {REFERENCE_SPECIES_ID}")

    # 校验 Step03 seqids 存在、非空、合法（硬防线）
    seqids_dir = step03_dir / "seqids_species"
    if not seqids_dir.exists():
        die(logger, f"Missing step03 seqids_species dir: {seqids_dir}")

    for sid in species:
        p = seqids_dir / f"{sid}.seqids"
        if not p.exists():
            die(logger, f"Missing seqids for {sid}: {p}")
        txt = p.read_text(encoding="utf-8", errors="replace")
        try:
            validate_seqids_line(sid, txt, logger)
        except Exception as e:
            die(logger, f"Seqids validation failed: {e}")

    clean_dir(out_dir, CLEAN_OUTPUT, logger)
    d_sum = out_dir / "summaries"
    d_sum.mkdir(parents=True, exist_ok=True)

    # y_center 从上到下递减，范围大致落在 (0,1)
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

    layout_path = out_dir / "layout_species_tracks.tsv"
    write_tsv(
        layout_path,
        ["order_index", "species_id", "display_label", "group", "y_center", "track_height", "chr_bar_color_hex", "label_color_hex"],
        layout_rows
    )

    sum_path = d_sum / "step04.summary.tsv"
    write_tsv(sum_path, ["species_id", "order_index", "display_label"], sum_rows)

    logger.info(f"Layout: {layout_path}")
    logger.info(f"Summary: {sum_path}")
    logger.info(f"Done. runtime_sec={int(time.time()-t0)}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        pr = PROJECT_ROOT
        if pr is None or str(pr).strip() == "":
            PROJECT = Path(__file__).resolve().parents[1]
        else:
            PROJECT = Path(str(pr)).expanduser().resolve()
        lg = Logger(PROJECT / LOG_DIR_REL / "synteny04_layout_tracks.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

