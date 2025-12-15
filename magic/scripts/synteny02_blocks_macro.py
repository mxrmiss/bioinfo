#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny02_blocks_macro.py
—— Step02：相邻物种对 anchors 生成 + 过滤 + macro blocks + chr 权重表（严格按蓝图合同）

输入（硬位置）：
- raw_data/synteny_species_meta.tsv
- output/synteny_01_mcscan_catalog/<species>.geneorder.bed
- output/synteny_01_mcscan_catalog/geneorder_index_<species>.tsv
- raw_data/proteins/<species>.(faa/fa/fasta/pep[.gz]) 自动探测

输出（硬位置）：
output/synteny_02_blocks_macro/
  raw_anchors/        <A>__vs__<B>.anchors.simple   (+ 同步保存 .anchors/.anchors.new 便于排查)
  anchors_filtered/   <A>__vs__<B>.anchors.simple.filtered
  blocks_merged/      <A>__vs__<B>.blocks.merged.tsv
  pair_chr_weight/    <A>__vs__<B>.chr_weight.tsv
  summaries/          step02.pairs.summary.tsv
  work_pairs/<pair>/  jcvi 工作目录（便于排查）

硬合同要点：
- 只连相邻 pairs： (s1,s2),(s2,s3)...（严格按 meta 行顺序）
- 稀疏旋钮只有 MINSPAN（jcvi synteny screen --minspan）
- anchors_filtered 只做“格式/合法性过滤”：必须 6 列；绝不把 #HEX*gene 当注释丢弃
- blocks_merged 与 chr_weight 的 span 使用“bp 级跨度定义”（来自 gene index 坐标）
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

THREADS = 25
MINSPAN = 5

# JCVI 运行入口：一般写 "python"；若需固定环境可写绝对路径
JCVI_PYTHON = "python"

META_TSV_REL = "raw_data/synteny_species_meta.tsv"
STEP01_DIR_REL = "output/synteny_01_mcscan_catalog"
PROTEIN_DIR_REL = "raw_data/proteins"
OUTPUT_DIR_REL = "output/synteny_02_blocks_macro"
LOG_DIR_REL = "logs"

PROTEIN_SUFFIX_CANDIDATES = [
    ".faa", ".fa", ".fasta", ".pep",
    ".faa.gz", ".fa.gz", ".fasta.gz", ".pep.gz",
]

# 下面这个是“固定实现细节”，不是旋钮：用于 blocks 合并（不建议改）
_FIXED_MERGE_GAP_BP = 2_000_000

# ============================================================
# 实现区（皇上勿改）
# ============================================================

import os
import sys
import csv
import gzip
import time
import shutil
import traceback
import subprocess
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterable


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


_HEX_PREFIX_RE = re.compile(r"^#[0-9a-fA-F]{6}\*")


def strip_color_prefix(token: str) -> str:
    token = token.strip()
    if _HEX_PREFIX_RE.match(token):
        return token.split("*", 1)[1]
    return token


def normalize_gene_id(raw: str) -> str:
    x = raw.strip().strip('"').strip("'").strip()
    if "|" in x:
        x = x.split("|")[-1]
    if x.startswith("rna-"):
        x = x[4:]
    return x


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


def detect_protein_file(protein_dir: Path, species_id: str) -> Path:
    for suf in PROTEIN_SUFFIX_CANDIDATES:
        p = protein_dir / f"{species_id}{suf}"
        if p.exists():
            return p
    raise FileNotFoundError(f"Protein FASTA not found for {species_id} in {protein_dir}")


def load_gene_index(step01_dir: Path, sid: str) -> Dict[str, Tuple[str, int, int]]:
    """
    geneorder_index_<sid>.tsv：gene_id chr start end strand
    返回：gene_id -> (chr, start, end)
    """
    p = step01_dir / f"geneorder_index_{sid}.tsv"
    if not p.exists():
        raise FileNotFoundError(f"Missing gene index: {p}")
    mp: Dict[str, Tuple[str, int, int]] = {}
    with open_text_auto(p) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header[:4] != ["gene_id", "chr", "start", "end"]:
            raise ValueError(f"Bad header in {p.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            gid = normalize_gene_id(strip_color_prefix(parts[0]))
            chr_ = parts[1]
            try:
                s = int(parts[2])
                e = int(parts[3])
            except Exception:
                continue
            if gid:
                mp[gid] = (chr_, s, e)
    if not mp:
        raise ValueError(f"Empty gene index: {p}")
    return mp


def run_cmd(logger: Logger, cmd: List[str], cwd: Path) -> None:
    cmd_str = " ".join(cmd)
    logger.info(f"[CMD] ({cwd}) {cmd_str}")
    r = subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = (r.stdout or "").rstrip("\n")
    if out:
        for ln in out.splitlines():
            logger.info(f"[JCVI] {ln}")
    if r.returncode != 0:
        raise RuntimeError(f"Command failed (exit={r.returncode}): {cmd_str}")


def filter_anchors_simple_6cols(in_path: Path, out_path: Path) -> Tuple[int, int]:
    """
    anchors.simple.filtered：只保留严格 6 列的行。
    注意：绝不因为行首 '#' 就丢弃（#HEX*gene 也是合法 token）。
    返回：(kept, dropped)
    """
    kept = 0
    dropped = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text_auto(in_path) as fr, out_path.open("w", encoding="utf-8") as fw:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                dropped += 1
                continue
            fw.write(line + "\n")
            kept += 1
    return kept, dropped


def _count_nonempty_lines(path: Path) -> int:
    n = 0
    with open_text_auto(path) as fr:
        for line in fr:
            if line and line.strip():
                n += 1
    return n


def build_blocks_from_simple(
    pair_id: str,
    a: str,
    b: str,
    simple_filtered: Path,
    idxA: Dict[str, Tuple[str, int, int]],
    idxB: Dict[str, Tuple[str, int, int]],
) -> List[Dict[str, object]]:
    """
    将 anchors.simple.filtered（6列：a1 a2 b1 b2 n orient）转成 block list（未合并）。
    span_bp 定义（固定实现）：span_bp = max(a_span, b_span)
      a_span = |a_end - a_start|, b_span = |b_end - b_start|
    """
    blocks: List[Dict[str, object]] = []
    with open_text_auto(simple_filtered) as fr:
        for ln, line in enumerate(fr, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                continue

            ga1 = normalize_gene_id(strip_color_prefix(parts[0]))
            ga2 = normalize_gene_id(strip_color_prefix(parts[1]))
            gb1 = normalize_gene_id(strip_color_prefix(parts[2]))
            gb2 = normalize_gene_id(strip_color_prefix(parts[3]))
            try:
                n_anchors = int(parts[4])
            except Exception:
                n_anchors = 0
            orient = parts[5].strip()
            if orient not in ("+", "-"):
                orient = "+"

            if ga1 not in idxA or ga2 not in idxA or gb1 not in idxB or gb2 not in idxB:
                continue

            a_chr1, a_s1, a_e1 = idxA[ga1]
            a_chr2, a_s2, a_e2 = idxA[ga2]
            b_chr1, b_s1, b_e1 = idxB[gb1]
            b_chr2, b_s2, b_e2 = idxB[gb2]

            # anchors.simple 理论上端点在同一 chr；若不一致，跳过（输入不自洽）
            if a_chr1 != a_chr2 or b_chr1 != b_chr2:
                continue

            a_pos1 = min(a_s1, a_e1)
            a_pos2 = min(a_s2, a_e2)
            b_pos1 = min(b_s1, b_e1)
            b_pos2 = min(b_s2, b_e2)

            a_start = min(a_pos1, a_pos2)
            a_end = max(a_pos1, a_pos2)
            b_start = min(b_pos1, b_pos2)
            b_end = max(b_pos1, b_pos2)

            a_span = max(0, a_end - a_start)
            b_span = max(0, b_end - b_start)
            span_bp = max(a_span, b_span)

            blocks.append({
                "pair_id": pair_id,
                "a_species": a,
                "b_species": b,
                "a_chr": a_chr1,
                "a_start": a_start,
                "a_end": a_end,
                "b_chr": b_chr1,
                "b_start": b_start,
                "b_end": b_end,
                "orientation": orient,
                "n_anchors": max(n_anchors, 0),
                "span_bp": span_bp,
            })
    return blocks


def merge_blocks(blocks: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """
    合并规则（固定实现，不作为旋钮）：
    - 必须同一 (a_chr,b_chr,orientation)
    - 且在两侧基因组上的间隔都 <= _FIXED_MERGE_GAP_BP
    """
    if not blocks:
        return []

    # 按 chrpair + 坐标排序
    blocks_sorted = sorted(
        blocks,
        key=lambda x: (x["a_chr"], x["b_chr"], x["orientation"], int(x["a_start"]), int(x["b_start"]))
    )

    merged: List[Dict[str, object]] = []
    cur = dict(blocks_sorted[0])

    for blk in blocks_sorted[1:]:
        same = (
            blk["a_chr"] == cur["a_chr"] and
            blk["b_chr"] == cur["b_chr"] and
            blk["orientation"] == cur["orientation"]
        )
        if not same:
            merged.append(cur)
            cur = dict(blk)
            continue

        gapA = int(blk["a_start"]) - int(cur["a_end"])
        gapB = int(blk["b_start"]) - int(cur["b_end"])

        if gapA <= _FIXED_MERGE_GAP_BP and gapB <= _FIXED_MERGE_GAP_BP:
            # 合并：扩展区间，累加 anchors，span_bp 重新按 max(a_span,b_span) 计算
            cur["a_start"] = min(int(cur["a_start"]), int(blk["a_start"]))
            cur["a_end"] = max(int(cur["a_end"]), int(blk["a_end"]))
            cur["b_start"] = min(int(cur["b_start"]), int(blk["b_start"]))
            cur["b_end"] = max(int(cur["b_end"]), int(blk["b_end"]))
            cur["n_anchors"] = int(cur["n_anchors"]) + int(blk["n_anchors"])

            a_span = max(0, int(cur["a_end"]) - int(cur["a_start"]))
            b_span = max(0, int(cur["b_end"]) - int(cur["b_start"]))
            cur["span_bp"] = max(a_span, b_span)
        else:
            merged.append(cur)
            cur = dict(blk)

    merged.append(cur)
    return merged


def build_chr_weight(
    pair_id: str,
    blocks_merged: List[Dict[str, object]]
) -> List[List[str]]:
    """
    输出 pair_chr_weight/<pair>.chr_weight.tsv
    Header:
      pair_id a_chr b_chr total_span_bp n_blocks
    """
    acc: Dict[Tuple[str, str], List[int]] = {}
    for blk in blocks_merged:
        k = (str(blk["a_chr"]), str(blk["b_chr"]))
        if k not in acc:
            acc[k] = [0, 0]  # total_span_bp, n_blocks
        acc[k][0] += int(blk["span_bp"])
        acc[k][1] += 1

    rows: List[List[str]] = []
    for (a_chr, b_chr), (tot, nb) in acc.items():
        rows.append([pair_id, a_chr, b_chr, str(tot), str(nb)])

    rows.sort(key=lambda r: (int(r[3]), int(r[4]), r[1], r[2]), reverse=True)
    return rows


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        PROJECT = Path(__file__).resolve().parents[1]
    else:
        PROJECT = Path(str(pr)).expanduser().resolve()

    meta_tsv = PROJECT / META_TSV_REL
    step01_dir = PROJECT / STEP01_DIR_REL
    protein_dir = PROJECT / PROTEIN_DIR_REL
    out_dir = PROJECT / OUTPUT_DIR_REL
    log_dir = PROJECT / LOG_DIR_REL

    logger = Logger(log_dir / "synteny02_blocks_macro.log")
    t0 = time.time()

    logger.info("========== synteny02 — anchors + filtered + blocks_merged + chr_weight ==========")
    logger.info(f"PROJECT_ROOT={PROJECT}")
    logger.info(f"OUTPUT_DIR={out_dir}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"THREADS={THREADS}")
    logger.info(f"MINSPAN={MINSPAN}")
    logger.info(f"JCVI_PYTHON={JCVI_PYTHON}")
    logger.info(f"MERGE_GAP_BP(fixed)={_FIXED_MERGE_GAP_BP}")
    logger.info(f"META={meta_tsv}")
    logger.info(f"STEP01_DIR={step01_dir}")
    logger.info(f"PROTEIN_DIR={protein_dir}")

    if not meta_tsv.exists():
        die(logger, f"Missing meta: {meta_tsv}")
    if not step01_dir.exists():
        die(logger, f"Missing step01 dir: {step01_dir}")
    if not protein_dir.exists():
        die(logger, f"Missing protein dir: {protein_dir}")

    species = read_meta(meta_tsv)
    pairs = [(species[i], species[i + 1]) for i in range(len(species) - 1)]
    logger.info("META species order=" + " | ".join(species))
    logger.info(f"N_species={len(species)} N_pairs={len(pairs)} (adjacent only)")

    clean_dir(out_dir, CLEAN_OUTPUT, logger)

    d_raw = out_dir / "raw_anchors"
    d_filtered = out_dir / "anchors_filtered"
    d_blocks = out_dir / "blocks_merged"
    d_chr_weight = out_dir / "pair_chr_weight"
    d_work = out_dir / "work_pairs"
    d_inputs = out_dir / "_jcvi_inputs"
    d_sum = out_dir / "summaries"
    for d in (d_raw, d_filtered, d_blocks, d_chr_weight, d_work, d_inputs, d_sum):
        d.mkdir(parents=True, exist_ok=True)

    # 为 jcvi 准备 <sid>.bed / <sid>.pep（软链接）
    for sid in species:
        bed_src = step01_dir / f"{sid}.geneorder.bed"
        if not bed_src.exists():
            die(logger, f"Missing geneorder bed (step01): {bed_src}")
        pep_src = detect_protein_file(protein_dir, sid)

        bed_dst = d_inputs / f"{sid}.bed"
        pep_dst = d_inputs / f"{sid}.pep"
        for src, dst in ((bed_src, bed_dst), (pep_src, pep_dst)):
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            dst.symlink_to(src)

    gene_idx_cache: Dict[str, Dict[str, Tuple[str, int, int]]] = {}

    def get_idx(sid: str) -> Dict[str, Tuple[str, int, int]]:
        if sid not in gene_idx_cache:
            gene_idx_cache[sid] = load_gene_index(step01_dir, sid)
        return gene_idx_cache[sid]

    sum_rows: List[List[str]] = []
    sum_header = [
        "pair_id", "a_species", "b_species",
        "anchors_raw", "anchors_kept",
        "blocks_merged", "chrpairs_total",
        "runtime_sec"
    ]

    for a, b in pairs:
        pair_id = f"{a}__vs__{b}"
        wdir = d_work / pair_id
        wdir.mkdir(parents=True, exist_ok=True)

        # 清理 pair 工作目录，避免旧残留
        for x in wdir.glob("*"):
            try:
                if x.is_dir():
                    shutil.rmtree(x)
                else:
                    x.unlink()
            except Exception:
                pass

        # 放入软链接：A.bed A.pep B.bed B.pep
        for sid in (a, b):
            for ext in (".bed", ".pep"):
                src = d_inputs / f"{sid}{ext}"
                dst = wdir / f"{sid}{ext}"
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                dst.symlink_to(src)

        t_pair = time.time()

        # 1) ortholog
        cmd1 = [
            JCVI_PYTHON, "-m", "jcvi.compara.catalog", "ortholog",
            "--dbtype", "prot",
            "--no_strip_names",
            a, b,
            f"--cpus={THREADS}",
        ]
        run_cmd(logger, cmd1, cwd=wdir)

        anchors_file = wdir / f"{a}.{b}.anchors"
        if not anchors_file.exists():
            die(logger, f"[{pair_id}] Missing anchors file after ortholog: {anchors_file}")

        # 2) screen（注意：screen 不传 --cpus）
        cmd2 = [
            JCVI_PYTHON, "-m", "jcvi.compara.synteny", "screen",
            f"--minspan={MINSPAN}",
            "--simple",
            f"{a}.{b}.anchors",
            f"{a}.{b}.anchors.new",
        ]
        run_cmd(logger, cmd2, cwd=wdir)

        anchors_new = wdir / f"{a}.{b}.anchors.new"
        anchors_simple = wdir / f"{a}.{b}.anchors.simple"
        if not anchors_new.exists():
            die(logger, f"[{pair_id}] Missing anchors.new after screen: {anchors_new}")
        if not anchors_simple.exists():
            die(logger, f"[{pair_id}] Missing anchors.simple after screen: {anchors_simple}")

        # 保存 raw
        out_anchors = d_raw / f"{pair_id}.anchors"
        out_new = d_raw / f"{pair_id}.anchors.new"
        out_simple = d_raw / f"{pair_id}.anchors.simple"
        shutil.copy2(anchors_file, out_anchors)
        shutil.copy2(anchors_new, out_new)
        shutil.copy2(anchors_simple, out_simple)

        # anchors_filtered（6列）
        out_filtered = d_filtered / f"{pair_id}.anchors.simple.filtered"
        kept, dropped = filter_anchors_simple_6cols(out_simple, out_filtered)
        logger.info(f"[{pair_id}] anchors.filtered kept={kept} dropped={dropped}")

        # blocks_merged
        idxA = get_idx(a)
        idxB = get_idx(b)
        raw_blocks = build_blocks_from_simple(pair_id, a, b, out_filtered, idxA, idxB)
        merged_blocks = merge_blocks(raw_blocks)

        # 输出 blocks_merged tsv
        blocks_path = d_blocks / f"{pair_id}.blocks.merged.tsv"
        blocks_header = [
            "pair_id", "a_species", "b_species",
            "a_chr", "a_start", "a_end",
            "b_chr", "b_start", "b_end",
            "orientation", "block_id", "n_anchors", "span_bp",
        ]
        block_rows: List[List[str]] = []
        for k, blk in enumerate(merged_blocks, start=1):
            block_rows.append([
                pair_id, a, b,
                str(blk["a_chr"]), str(int(blk["a_start"])), str(int(blk["a_end"])),
                str(blk["b_chr"]), str(int(blk["b_start"])), str(int(blk["b_end"])),
                str(blk["orientation"]), f"B{k:05d}",
                str(int(blk["n_anchors"])),
                str(int(blk["span_bp"])),
            ])
        write_tsv(blocks_path, blocks_header, block_rows)

        # chr_weight
        chrw_path = d_chr_weight / f"{pair_id}.chr_weight.tsv"
        chrw_rows = build_chr_weight(pair_id, merged_blocks)
        write_tsv(chrw_path, ["pair_id", "a_chr", "b_chr", "total_span_bp", "n_blocks"], chrw_rows)

        runtime_sec = int(time.time() - t_pair)

        anchors_raw_n = _count_nonempty_lines(out_simple)  # raw anchors.simple 行数（近似）
        anchors_kept_n = _count_nonempty_lines(out_filtered)

        sum_rows.append([
            pair_id, a, b,
            str(anchors_raw_n),
            str(anchors_kept_n),
            str(len(merged_blocks)),
            str(len(chrw_rows)),
            str(runtime_sec),
        ])

        logger.info(
            f"[{pair_id}] anchors_raw(simple)={anchors_raw_n} anchors_kept={anchors_kept_n} "
            f"blocks_merged={len(merged_blocks)} chrpairs_total={len(chrw_rows)} runtime_sec={runtime_sec}"
        )

        if anchors_kept_n == 0:
            logger.warn(f"[{pair_id}] anchors.simple.filtered is empty -> downstream will fail")
        if len(merged_blocks) == 0:
            logger.warn(f"[{pair_id}] blocks_merged is empty -> downstream will fail")

    out_sum = d_sum / "step02.pairs.summary.tsv"
    write_tsv(out_sum, sum_header, sum_rows)

    logger.info(f"Done. runtime_sec={int(time.time()-t0)}")
    logger.info(f"Summary: {out_sum}")


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
        lg = Logger(PROJECT / LOG_DIR_REL / "synteny02_blocks_macro.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

