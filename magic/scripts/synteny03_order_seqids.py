#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny03_order_seqids.py
—— Step03：颜色链式继承 + chr 排序/翻转 + seqids 生成（严格按蓝图合同）

输入（硬位置）：
- raw_data/synteny_species_meta.tsv
- raw_data/palette.tsv                           (rank, color_hex)
- output/synteny_00_chr_rename/chr_rename_<sp>.tsv
- output/synteny_02_blocks_macro/pair_chr_weight/<A>__vs__<B>.chr_weight.tsv
- output/synteny_02_blocks_macro/blocks_merged/<A>__vs__<B>.blocks.merged.tsv

输出（硬位置）：
output/synteny_03_order_seqids/
  global_chr_style.tsv
  seqids_species/<species>.seqids                (jcvi 使用：单行逗号分隔；ChrXX- 表示翻转)
  chr_order/<species>.chr_order.tsv              (便于检查；global 的子集视图)
  chr_colors/<species>.chr_colors.tsv            (便于检查)
  summaries/step03.summary.tsv
  summaries/palette.tsv                          (落盘备查)

硬合同要点：
- 颜色链式继承：从 meta 的第一个物种 s1 开始，向后传播到 sN（不以 reference 为锚点）
- dominant mapping：对 pair (Si-1,Si)，对每个 Si_chr 取 total_span_bp 最大的 prev_chr
- orientation：对每个 Si_chr，用 blocks_merged 中涉及该 chr 的所有块做投票，权重 = span_bp
  （固定实现：对 i>1 的物种 Si，仅使用 pair (Si-1,Si) 的 blocks_merged 来给 Si 投票，符合“相邻链式”）
- chr 排序键（写死）：
  1) color_id（1..K 在前，0 在后）
  2) dominant_prev_chr 在上一物种的 order_index（小的在前）
  3) Step00 的 rank（小的在前）
- 翻转编码（蓝图硬规则）：ChrNN- 表示翻转；严禁 ChrNNr
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

META_TSV_REL = "raw_data/synteny_species_meta.tsv"
PALETTE_TSV_REL = "raw_data/palette.tsv"

STEP00_DIR_REL = "output/synteny_00_chr_rename"
STEP02_DIR_REL = "output/synteny_02_blocks_macro"

OUTPUT_DIR_REL = "output/synteny_03_order_seqids"
LOG_DIR_REL = "logs"

# ============================================================
# 实现区（皇上勿改）
# ============================================================

import re
import csv
import gzip
import time
import shutil
import traceback
from pathlib import Path
from typing import Dict, List, Tuple, Iterable


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


def read_meta(meta_tsv: Path) -> List[str]:
    species: List[str] = []
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
            if not sid:
                raise ValueError(f"meta line {ln} species_id empty")
            species.append(sid)
    if len(species) < 2:
        raise ValueError("meta must contain at least 2 species")
    return species


def safe_int(x: str, default: int = 0) -> int:
    try:
        return int(x)
    except Exception:
        return default


def load_palette(palette_tsv: Path) -> List[str]:
    """
    raw_data/palette.tsv
    header: rank, color_hex
    返回：按 rank 升序的 color_hex 列表（K = len）
    """
    if not palette_tsv.exists():
        raise FileNotFoundError(f"Missing palette: {palette_tsv}")
    rows: List[Tuple[int, str]] = []
    with open_text_auto(palette_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != ["rank", "color_hex"]:
            raise ValueError(f"palette header must be: rank<TAB>color_hex, got: {header}")
        for ln, line in enumerate(fr, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"palette line {ln} must have 2 columns: {line}")
            rk = safe_int(parts[0], 0)
            hx = parts[1].strip()
            if rk <= 0 or not hx:
                continue
            rows.append((rk, hx))
    rows.sort(key=lambda x: x[0])
    pal = [hx for _, hx in rows]
    if len(pal) < 6:
        raise ValueError(f"palette must contain >=6 colors, got {len(pal)}")
    return pal


def load_chr_list_and_rank(step00_dir: Path, sid: str) -> Tuple[List[str], Dict[str, int]]:
    """
    Step00 chr_rename_<sid>.tsv
    header:
      species_id seqid_raw seqid_renamed rank length_bp is_chromosome
    返回：
      chr_list：按 rank 升序的 ChrNN 列表（仅 is_chromosome==yes）
      rank_map：ChrNN -> rank
    """
    p = step00_dir / f"chr_rename_{sid}.tsv"
    if not p.exists():
        raise FileNotFoundError(f"Missing step00 file: {p}")
    rows: List[Tuple[int, str]] = []
    rank_map: Dict[str, int] = {}
    with open_text_auto(p) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        need = ["species_id", "seqid_raw", "seqid_renamed", "rank", "length_bp", "is_chromosome"]
        if header != need:
            raise ValueError(f"Bad header in {p.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                continue
            if parts[5] != "yes":
                continue
            rk = safe_int(parts[3], 0)
            chr_ = parts[2]
            if rk > 0 and chr_ != "NA":
                rows.append((rk, chr_))
                rank_map[chr_] = rk
    rows.sort(key=lambda x: x[0])
    chr_list = [c for _, c in rows]
    if not chr_list:
        raise ValueError(f"No chromosomes found in {p.name}")
    return chr_list, rank_map


def load_pair_chr_weight(step02_dir: Path, pair_id: str) -> List[Tuple[str, str, int]]:
    """
    pair_chr_weight/<pair>.chr_weight.tsv
    header: pair_id a_chr b_chr total_span_bp n_blocks
    返回：(a_chr, b_chr, total_span_bp)
    """
    p = step02_dir / "pair_chr_weight" / f"{pair_id}.chr_weight.tsv"
    if not p.exists():
        raise FileNotFoundError(f"Missing chr_weight: {p}")
    out: List[Tuple[str, str, int]] = []
    with open_text_auto(p) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        need = ["pair_id", "a_chr", "b_chr", "total_span_bp", "n_blocks"]
        if header != need:
            raise ValueError(f"Bad header in {p.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 5:
                continue
            a_chr, b_chr = parts[1], parts[2]
            tot = safe_int(parts[3], 0)
            out.append((a_chr, b_chr, tot))
    return out


def build_dominant_prev_map(weight_rows: List[Tuple[str, str, int]]) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    对每个 b_chr：选择 total_span_bp 最大的 a_chr 作为 dominant_prev_chr
    返回：
      dom_prev: b_chr -> a_chr
      dom_w:    b_chr -> total_span_bp
    """
    best: Dict[str, Tuple[int, str]] = {}
    for a_chr, b_chr, tot in weight_rows:
        if b_chr not in best or tot > best[b_chr][0]:
            best[b_chr] = (tot, a_chr)
    dom_prev = {b: a for b, (tot, a) in best.items()}
    dom_w = {b: tot for b, (tot, a) in best.items()}
    return dom_prev, dom_w


def load_blocks_for_pair(step02_dir: Path, pair_id: str) -> List[Tuple[str, str, str, int, int, str, int]]:
    """
    blocks_merged/<pair>.blocks.merged.tsv
    header:
      pair_id a_species b_species a_chr a_start a_end b_chr b_start b_end orientation block_id n_anchors span_bp
    返回（仅保留需要字段）：
      (a_chr, b_chr, orientation, a_start, a_end, b_chr(占位), span_bp)
    """
    p = step02_dir / "blocks_merged" / f"{pair_id}.blocks.merged.tsv"
    if not p.exists():
        raise FileNotFoundError(f"Missing blocks_merged: {p}")
    out: List[Tuple[str, str, str, int, int, str, int]] = []
    with open_text_auto(p) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        need = [
            "pair_id", "a_species", "b_species",
            "a_chr", "a_start", "a_end",
            "b_chr", "b_start", "b_end",
            "orientation", "block_id", "n_anchors", "span_bp",
        ]
        if header != need:
            raise ValueError(f"Bad header in {p.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 13:
                continue
            a_chr = parts[3]
            b_chr = parts[6]
            orient = parts[9].strip()
            a_start = safe_int(parts[4], 0)
            a_end = safe_int(parts[5], 0)
            span_bp = safe_int(parts[12], 0)
            out.append((a_chr, b_chr, orient, a_start, a_end, b_chr, span_bp))
    return out


def orientation_vote_for_bchr(blocks: List[Tuple[str, str, str, int, int, str, int]]) -> Dict[str, str]:
    """
    对 b_chr 做投票：权重 = span_bp
    若 plus >= minus 则 '+', 否则 '-'
    """
    acc: Dict[str, List[int]] = {}  # b_chr -> [plus_w, minus_w]
    for a_chr, b_chr, orient, a_s, a_e, _b_chr2, span_bp in blocks:
        if b_chr not in acc:
            acc[b_chr] = [0, 0]
        w = max(span_bp, 1)
        if orient == "-":
            acc[b_chr][1] += w
        else:
            acc[b_chr][0] += w
    out: Dict[str, str] = {}
    for b_chr, (pw, mw) in acc.items():
        out[b_chr] = "+" if pw >= mw else "-"
    return out


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        PROJECT = Path(__file__).resolve().parents[1]
    else:
        PROJECT = Path(str(pr)).expanduser().resolve()

    meta_tsv = PROJECT / META_TSV_REL
    palette_tsv = PROJECT / PALETTE_TSV_REL
    step00_dir = PROJECT / STEP00_DIR_REL
    step02_dir = PROJECT / STEP02_DIR_REL
    out_dir = PROJECT / OUTPUT_DIR_REL
    log_dir = PROJECT / LOG_DIR_REL

    logger = Logger(log_dir / "synteny03_order_seqids.log")
    t0 = time.time()

    logger.info("========== synteny03 — chain colors + order + orientation ==========")
    logger.info(f"PROJECT_ROOT={PROJECT}")
    logger.info(f"OUTPUT_DIR={out_dir}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"META={meta_tsv}")
    logger.info(f"PALETTE={palette_tsv}")
    logger.info(f"STEP00_DIR={step00_dir}")
    logger.info(f"STEP02_DIR={step02_dir}")
    logger.info("Flip encoding: ChrNN- (strict); forbid ChrNNr")

    if not meta_tsv.exists():
        die(logger, f"Missing meta: {meta_tsv}")
    if not palette_tsv.exists():
        die(logger, f"Missing palette: {palette_tsv}")
    if not step00_dir.exists():
        die(logger, f"Missing step00 dir: {step00_dir}")
    if not step02_dir.exists():
        die(logger, f"Missing step02 dir: {step02_dir}")

    species = read_meta(meta_tsv)
    logger.info("META species order=" + " | ".join(species))
    if len(species) < 2:
        die(logger, "meta must contain at least 2 species")

    palette_hex = load_palette(palette_tsv)
    K = len(palette_hex)
    logger.info(f"PALETTE_N={K} (K=palette rows)")

    clean_dir(out_dir, CLEAN_OUTPUT, logger)
    d_seqids = out_dir / "seqids_species"
    d_chr_order = out_dir / "chr_order"
    d_chr_colors = out_dir / "chr_colors"
    d_sum = out_dir / "summaries"
    for d in (d_seqids, d_chr_order, d_chr_colors, d_sum):
        d.mkdir(parents=True, exist_ok=True)

    # palette 落盘备查（严格两列：rank,color_hex）
    write_tsv(d_sum / "palette.tsv", ["rank", "color_hex"], [[str(i + 1), palette_hex[i]] for i in range(K)])

    # Step00：chr 列表与 rank_map
    chr_list: Dict[str, List[str]] = {}
    chr_rank_map: Dict[str, Dict[str, int]] = {}
    for sid in species:
        cl, rm = load_chr_list_and_rank(step00_dir, sid)
        chr_list[sid] = cl
        chr_rank_map[sid] = rm
        logger.info(f"[{sid}] n_chr={len(cl)}")

    # 全局容器
    ordered: Dict[str, List[str]] = {}
    orient: Dict[str, Dict[str, str]] = {}
    dom_prev: Dict[str, Dict[str, str]] = {}
    dom_w: Dict[str, Dict[str, int]] = {}
    color_id: Dict[str, Dict[str, int]] = {}

    # s1 初始化：按 Step00 rank 顺序
    s1 = species[0]
    ordered[s1] = list(chr_list[s1])
    orient[s1] = {c: "+" for c in ordered[s1]}
    dom_prev[s1] = {c: "NA" for c in ordered[s1]}
    dom_w[s1] = {c: 0 for c in ordered[s1]}

    cid1: Dict[str, int] = {}
    for i, c in enumerate(ordered[s1], start=1):
        cid1[c] = i if i <= K else 0
    color_id[s1] = cid1

    def order_index_map(sid: str) -> Dict[str, int]:
        return {c: i for i, c in enumerate(ordered[sid], start=1)}

    # 向后链式传播：s2..sN
    for i in range(1, len(species)):
        prev = species[i - 1]
        cur = species[i]
        pair_id = f"{prev}__vs__{cur}"

        wrows = load_pair_chr_weight(step02_dir, pair_id)
        cur2prev, cur2w = build_dominant_prev_map(wrows)

        blocks = load_blocks_for_pair(step02_dir, pair_id)
        cur_orient_vote = orientation_vote_for_bchr(blocks)

        prev_order = order_index_map(prev)

        cur_dom_prev: Dict[str, str] = {}
        cur_dom_w: Dict[str, int] = {}
        cur_cid: Dict[str, int] = {}
        cur_ori: Dict[str, str] = {}

        for c in chr_list[cur]:
            dp = cur2prev.get(c, "NA")
            cur_dom_prev[c] = dp
            cur_dom_w[c] = cur2w.get(c, 0)

            if dp != "NA":
                cur_cid[c] = color_id[prev].get(dp, 0)
            else:
                cur_cid[c] = 0

            cur_ori[c] = cur_orient_vote.get(c, "+")  # 没证据则 '+'

        def sort_key(chr_c: str):
            cid = cur_cid.get(chr_c, 0)
            dp = cur_dom_prev.get(chr_c, "NA")
            oi = prev_order.get(dp, 10**9)
            rk = chr_rank_map[cur].get(chr_c, 10**9)
            return (999 if cid == 0 else cid, oi, rk, chr_c)

        ordered[cur] = sorted(chr_list[cur], key=sort_key)
        orient[cur] = cur_ori
        dom_prev[cur] = cur_dom_prev
        dom_w[cur] = cur_dom_w
        color_id[cur] = cur_cid

        mapped_n = sum(1 for x in cur_dom_prev.values() if x != "NA")
        colored_n = sum(1 for x in cur_cid.values() if x > 0)
        flipped_n = sum(1 for x in cur_ori.values() if x == "-")
        logger.info(f"[prop] {prev} -> {cur}: mapped={mapped_n}/{len(chr_list[cur])} colored={colored_n}/{len(chr_list[cur])} flipped={flipped_n}/{len(chr_list[cur])}")

    # 输出
    global_header = ["species_id", "chr", "order_index", "orientation", "dominant_prev_chr", "dominant_weight_span_bp", "color_id", "color_hex"]
    global_rows: List[List[str]] = []

    sum_header = ["species_id", "n_chr", "n_assigned_color", "n_flipped", "prev_species_id"]
    sum_rows: List[List[str]] = []

    for i, sid in enumerate(species):
        prev_sid = "NA" if i == 0 else species[i - 1]
        n_chr = len(ordered[sid])

        n_col = 0
        n_flip = 0

        chr_order_rows: List[List[str]] = []
        chr_colors_rows: List[List[str]] = []
        seqids_tokens: List[str] = []

        for j, c in enumerate(ordered[sid], start=1):
            o = orient[sid].get(c, "+")
            if o == "-":
                n_flip += 1

            cid = color_id[sid].get(c, 0)
            if cid > 0:
                n_col += 1
                chex = palette_hex[cid - 1]
            else:
                chex = "NA"

            dp = dom_prev[sid].get(c, "NA")
            dw = dom_w[sid].get(c, 0)

            global_rows.append([sid, c, str(j), o, dp, str(dw), str(cid), chex])

            # ====== 翻转编码：ChrNN-（硬规则）======
            seqid_for_plot = c + ("-" if o == "-" else "")
            seqids_tokens.append(seqid_for_plot)

            chr_order_rows.append([sid, c, str(j), o, seqid_for_plot, dp, str(dw), str(cid), chex])
            chr_colors_rows.append([sid, c, str(cid), chex])

        write_tsv(
            d_chr_order / f"{sid}.chr_order.tsv",
            ["species_id", "chr", "order_index", "orientation", "seqid_for_plot", "dominant_prev_chr", "dominant_weight_span_bp", "color_id", "color_hex"],
            chr_order_rows,
        )

        write_tsv(
            d_chr_colors / f"{sid}.chr_colors.tsv",
            ["species_id", "chr", "color_id", "color_hex"],
            chr_colors_rows,
        )

        # jcvi karyotype：每个 track 一行逗号分隔
        (d_seqids / f"{sid}.seqids").write_text(",".join(seqids_tokens) + "\n", encoding="utf-8")

        sum_rows.append([sid, str(n_chr), str(n_col), str(n_flip), prev_sid])
        logger.info(f"[{sid}] n_chr={n_chr} n_assigned_color={n_col} n_flipped={n_flip} prev_species_id={prev_sid}")

    write_tsv(out_dir / "global_chr_style.tsv", global_header, global_rows)
    write_tsv(d_sum / "step03.summary.tsv", sum_header, sum_rows)

    logger.info(f"Done. runtime_sec={int(time.time()-t0)}")
    logger.info(f"Global style: {out_dir / 'global_chr_style.tsv'}")


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
        lg = Logger(PROJECT / LOG_DIR_REL / "synteny03_order_seqids.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

