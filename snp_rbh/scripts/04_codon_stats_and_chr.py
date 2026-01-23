#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
04_codon_stats_and_chr.py
================================================================================
[定位]
这是“统计层”脚本：只负责生产标准、可复用的统计产物（不为画图脚本迁就路径/口径）。

[输入]
1) CODON_DIR 目录下的 *.codon.fa
   - 每个文件应只有两条序列（两个物种的同源对），为 pal2nal 回译后的 codon alignment
2) REF_GFF（参考物种的 gff/gff3，可 .gz）
   - 用 transcript 的 ID= 映射到 chr + start + end（用于 chr 汇总与滑窗）
3) 可选：REF_RENAME_MAP（chr 重命名映射表 old_seqid -> chr01..）
   - 若留空，自动推断 results/step0_chr/<code>.chr_rename_map.tsv

[输出]（全部放在 results/step5_stats，避免目录混乱）
A) results/step5_stats/pair_stats.tsv
   - 每个同源对一行（包含 chr/坐标/差异统计/密度）
B) results/step5_stats/chr_summary.tsv
   - 以 chr 聚合（包含可复现的累计量 + 派生密度）
C) results/step5_stats/window_1Mb.tsv
   - 以 1Mb 窗口聚合（包含累计量 + 派生密度）
   - 仅对“在参考GFF能定位到坐标”的 ref_id 才进入滑窗统计

[统计口径]
- substitutions：双方均为 A/C/G/T 且不同的位点数
- indel_bases：一方为 '-' 另一方非 '-' 的位点数（按碱基位点计数）
- diff_sites = substitutions + indel_bases
- aligned_length：双方均为 A/C/G/T 的位点数
- *_per_kb：除以 aligned_length 后 *1000
- density_per_Mb_aligned：sum_diff_sites / (sum_aligned_nuc / 1e6)

[进度输出]
- 每 PROGRESS_EVERY 个文件打印一次（流式输出到 stderr）
"""

# ----------------------------- 用户自定义区（皇上只改这里即可） -----------------------------
# codon alignment 文件目录（每个文件 2 条序列）
CODON_DIR = "results/step4_codon/codon_aln"

# 参考物种 GFF（用于 ref_id -> chr + start + end）
REF_GFF = "data/annotation/Sinonovacula_rivularis.gff.gz"

# 参考物种 genome fasta（优先用于自动推断物种缩写 sco/sri）
REF_GENOME = "data/genomic/Sinonovacula_rivularis.fasta.gz"

# chr 重命名映射表（可选）
# - 留空 ""：自动推断 results/step0_chr/<code>.chr_rename_map.tsv
# - 填具体路径：优先使用你指定的
REF_RENAME_MAP = ""

# 滑窗参数
WINDOW_BP = 1_000_000

# 输出目录与文件名（统计层：统一放 step5_stats）
OUTDIR = "results/step5_stats"
OUT_PAIR_TSV = "pair_stats.tsv"
OUT_CHR_TSV = "chr_summary.tsv"
OUT_WIN_TSV = "window_1Mb.tsv"

# 进度输出
PROGRESS_EVERY = 100
# -----------------------------------------------------------------------------------------

import os
import re
import sys
import math
import gzip
from collections import defaultdict
from typing import Dict, Optional, Tuple, List

# ----------------------------- 基础工具 -----------------------------
def project_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def ts_now() -> str:
    import datetime
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log(msg: str) -> None:
    sys.stderr.write(f"{ts_now()} [04] {msg}\n")

def open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")

def infer_code_from_path(path: str) -> Optional[str]:
    """
    从文件名推断物种缩写：genus首字母 + species前两字母（全小写）
    例：Sinonovacula_constricta -> sco；Sinonovacula_rivularis -> sri
    """
    try:
        bn = os.path.basename(path)
        for suf in [".gz", ".gff3", ".gff", ".fasta", ".fa", ".fna", ".faa"]:
            if bn.endswith(suf):
                bn = bn[: -len(suf)]
        parts = bn.split("_")
        if len(parts) < 2:
            return None
        genus, species = parts[0], parts[1]
        if not genus or not species:
            return None
        return (genus[0] + species[:2]).lower()
    except Exception:
        return None

def auto_rename_map_path(ref_genome: str, ref_gff: str) -> Optional[str]:
    code = None
    if ref_genome and os.path.exists(ref_genome):
        code = infer_code_from_path(ref_genome)
    if code is None:
        code = infer_code_from_path(ref_gff)
    if code is None:
        return None
    return os.path.join("results", "step0_chr", f"{code}.chr_rename_map.tsv")

def load_chr_rename_map(path: str) -> Dict[str, str]:
    """
    读取 old_seqid -> new_chr
    文件必须带表头：old_seqid  new_chr ...
    """
    mp: Dict[str, str] = {}
    if not path:
        return mp
    if (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return mp
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline()
        if not header:
            return mp
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            mp[cols[0]] = cols[1]
    return mp

# ----------------------------- 解析 GFF：ref_id -> (chr_raw, start, end) -----------------------------
def load_tx_info_from_gff(gff_path: str) -> Dict[str, Tuple[str, int, int]]:
    """
    从 GFF 中读取 transcript 信息：
    key = transcript ID（属性 ID=xxx）
    val = (seqid, start, end)
    只认 feature = mRNA / transcript
    """
    tx: Dict[str, Tuple[str, int, int]] = {}
    pat_id = re.compile(r"(?:^|;)ID=([^;]+)")
    with open_maybe_gz(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid = parts[0]
            ftype = parts[2]
            if ftype not in ("mRNA", "transcript"):
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except Exception:
                continue
            attr = parts[8]
            m = pat_id.search(attr)
            if not m:
                continue
            tid = m.group(1)
            if tid in tx:
                continue
            if start <= end:
                tx[tid] = (seqid, start, end)
            else:
                tx[tid] = (seqid, end, start)
    return tx

# ----------------------------- 读取 codon alignment：每文件两条序列 -----------------------------
def read_fasta_two(path: str) -> Optional[Tuple[str, str, str, str]]:
    ids: List[str] = []
    seqs: List[str] = []
    cur_id = None
    buf: List[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    ids.append(cur_id)
                    seqs.append("".join(buf))
                cur_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur_id is not None:
            ids.append(cur_id)
            seqs.append("".join(buf))
    if len(ids) != 2:
        return None
    return ids[0], seqs[0], ids[1], seqs[1]

def is_acgt(c: str) -> bool:
    return c in ("A", "C", "G", "T", "a", "c", "g", "t")

def codon_stats(seq1: str, seq2: str) -> Tuple[int, int, int]:
    """
    substitutions: A/C/G/T 之间不同算 1
    indel_bases: 单边为 '-' 的位置算 1
    aligned_length: 双边都是 A/C/G/T 的位置数
    """
    n = min(len(seq1), len(seq2))
    subs = 0
    indel = 0
    aligned = 0
    for i in range(n):
        a = seq1[i]
        b = seq2[i]
        if a == "-" and b == "-":
            continue
        if a == "-" and b != "-":
            indel += 1
            continue
        if b == "-" and a != "-":
            indel += 1
            continue
        if is_acgt(a) and is_acgt(b):
            aligned += 1
            if a.upper() != b.upper():
                subs += 1
    return subs, indel, aligned

# ----------------------------- 滑窗分箱 -----------------------------
def win_bounds_from_mid(mid_bp: float, win_bp: int) -> Tuple[int, int]:
    """
    用 midpoint 将记录归入固定窗口（1-based inclusive）
    win_start = floor((mid-1)/win_bp)*win_bp + 1
    win_end   = win_start + win_bp - 1
    """
    if (not math.isfinite(mid_bp)) or mid_bp <= 0:
        return (0, 0)
    k = int((mid_bp - 1) // win_bp)
    win_start = k * win_bp + 1
    win_end = win_start + win_bp - 1
    return win_start, win_end

# ----------------------------- 主流程 -----------------------------
def main():
    os.chdir(project_root())

    out_dir = OUTDIR
    os.makedirs(out_dir, exist_ok=True)

    out_pair = os.path.join(out_dir, OUT_PAIR_TSV)
    out_chr = os.path.join(out_dir, OUT_CHR_TSV)
    out_win = os.path.join(out_dir, OUT_WIN_TSV)

    log(f"CODON_DIR={CODON_DIR}")
    log(f"REF_GFF={REF_GFF}")
    log(f"REF_GENOME={REF_GENOME}")
    log(f"WINDOW_BP={WINDOW_BP}")
    log(f"OUTDIR={OUTDIR}")

    # 读 GFF：ref_id -> (chr_raw, start, end)
    log("Loading transcript info from GFF...")
    tx_info = load_tx_info_from_gff(REF_GFF)
    log(f"tx_info={len(tx_info)}")

    # chr 重命名映射表
    rename_map_path = REF_RENAME_MAP.strip()
    if not rename_map_path:
        auto_path = auto_rename_map_path(REF_GENOME, REF_GFF)
        if auto_path:
            rename_map_path = auto_path

    chr_map = load_chr_rename_map(rename_map_path)
    if chr_map:
        log(f"chr_rename_map=ON path={rename_map_path} n={len(chr_map)}")
    else:
        if rename_map_path:
            log(f"chr_rename_map=OFF (not found/empty) path={rename_map_path}")
        else:
            log("chr_rename_map=OFF (no path)")

    # codon 文件
    if not os.path.exists(CODON_DIR):
        log(f"ERROR: CODON_DIR not found: {CODON_DIR}")
        sys.exit(1)

    files = [f for f in os.listdir(CODON_DIR) if f.endswith(".codon.fa")]
    files.sort()
    n_total = len(files)
    if n_total == 0:
        log(f"ERROR: no .codon.fa found in {CODON_DIR}")
        sys.exit(1)

    # 统计容器
    pair_rows = []  # 每对一行
    chr_acc = defaultdict(lambda: {
        "n_pairs": 0,
        "sum_aligned": 0,
        "sum_sub": 0,
        "sum_indel": 0,
        "sum_diff": 0,
    })
    win_acc = defaultdict(lambda: {
        "n_pairs": 0,
        "sum_aligned": 0,
        "sum_sub": 0,
        "sum_indel": 0,
        "sum_diff": 0,
    })

    n_bad = 0
    n_no_ref = 0
    n_no_coord = 0

    log("Processing codon alignments...")
    for i, fn in enumerate(files, start=1):
        path = os.path.join(CODON_DIR, fn)
        rec = read_fasta_two(path)
        if rec is None:
            n_bad += 1
            continue
        id1, s1, id2, s2 = rec

        # 判定参考端：在 GFF tx_info 里能找到的那个就是 ref
        if id1 in tx_info:
            ref_id, other_id = id1, id2
            ref_seq, other_seq = s1, s2
        elif id2 in tx_info:
            ref_id, other_id = id2, id1
            ref_seq, other_seq = s2, s1
        else:
            ref_id, other_id = id1, id2
            ref_seq, other_seq = s1, s2
            n_no_ref += 1

        subs, indel, aligned = codon_stats(ref_seq, other_seq)
        diff_sites = subs + indel

        if aligned > 0:
            sub_per_kb = subs / aligned * 1000.0
            indel_per_kb = indel / aligned * 1000.0
            total_per_kb = diff_sites / aligned * 1000.0
        else:
            sub_per_kb = math.nan
            indel_per_kb = math.nan
            total_per_kb = math.nan

        # chr + 坐标
        chr_raw = "NA"
        chr = "NA"
        start = "NA"
        end = "NA"
        mid_bp = math.nan

        if ref_id in tx_info:
            chr_raw0, st0, en0 = tx_info[ref_id]
            chr_raw = chr_raw0
            chr = chr_map.get(chr_raw0, chr_raw0) if chr_map else chr_raw0
            start = st0
            end = en0
            mid_bp = (st0 + en0) / 2.0
        else:
            n_no_coord += 1

        # pair 行
        pair_rows.append((
            chr, chr_raw,
            ref_id, other_id,
            start, end,
            mid_bp if math.isfinite(mid_bp) else "NA",
            subs, indel, diff_sites, aligned,
            sub_per_kb, indel_per_kb, total_per_kb
        ))

        # chr 聚合（只要能定位到 chr 且 aligned>0）
        if chr != "NA" and aligned > 0:
            a = chr_acc[chr]
            a["n_pairs"] += 1
            a["sum_aligned"] += aligned
            a["sum_sub"] += subs
            a["sum_indel"] += indel
            a["sum_diff"] += diff_sites

        # window 聚合（只要能定位到坐标 且 aligned>0）
        if chr != "NA" and aligned > 0 and math.isfinite(mid_bp):
            wst, wen = win_bounds_from_mid(mid_bp, WINDOW_BP)
            if wst > 0:
                key = (chr, wst, wen)
                a = win_acc[key]
                a["n_pairs"] += 1
                a["sum_aligned"] += aligned
                a["sum_sub"] += subs
                a["sum_indel"] += indel
                a["sum_diff"] += diff_sites

        if (i % PROGRESS_EVERY) == 0 or i == n_total:
            log(f"{i}/{n_total} processed | bad={n_bad} no_ref={n_no_ref} no_coord={n_no_coord}")

    # ----------------------------- 写 pair_stats.tsv -----------------------------
    log(f"Writing: {out_pair}")
    with open(out_pair, "w", encoding="utf-8") as w:
        w.write("\t".join([
            "chr", "chr_raw",
            "ref_id", "other_id",
            "ref_start", "ref_end", "ref_mid_bp",
            "substitutions", "indel_bases", "diff_sites",
            "aligned_length",
            "sub_per_kb", "indel_per_kb", "total_per_kb"
        ]) + "\n")
        for r in pair_rows:
            w.write("\t".join(map(str, r)) + "\n")

    # ----------------------------- 写 chr_summary.tsv -----------------------------
    # 同时给出可复现累计量 + 派生密度（per kb / per Mb）
    log(f"Writing: {out_chr}")
    with open(out_chr, "w", encoding="utf-8") as w:
        w.write("\t".join([
            "chr",
            "n_pairs",
            "sum_aligned_nuc",
            "sum_substitutions",
            "sum_indel_bases",
            "sum_diff_sites",
            "sub_per_kb",
            "indel_per_kb",
            "diff_per_kb",
            "density_per_Mb_aligned"
        ]) + "\n")
        for chr_name in sorted(chr_acc.keys()):
            d = chr_acc[chr_name]
            aln = d["sum_aligned"]
            if aln > 0:
                sub_kb = d["sum_sub"] / aln * 1000.0
                indel_kb = d["sum_indel"] / aln * 1000.0
                diff_kb = d["sum_diff"] / aln * 1000.0
                dens_mb = d["sum_diff"] / (aln / 1e6)
            else:
                sub_kb = math.nan
                indel_kb = math.nan
                diff_kb = math.nan
                dens_mb = math.nan
            w.write("\t".join([
                chr_name,
                str(d["n_pairs"]),
                str(aln),
                str(d["sum_sub"]),
                str(d["sum_indel"]),
                str(d["sum_diff"]),
                str(sub_kb),
                str(indel_kb),
                str(diff_kb),
                str(dens_mb),
            ]) + "\n")

    # ----------------------------- 写 window_1Mb.tsv -----------------------------
    log(f"Writing: {out_win}")
    with open(out_win, "w", encoding="utf-8") as w:
        w.write("\t".join([
            "chr",
            "win_start",
            "win_end",
            "win_mid_bp",
            "n_pairs",
            "sum_aligned_nuc",
            "sum_substitutions",
            "sum_indel_bases",
            "sum_diff_sites",
            "sub_per_kb",
            "indel_per_kb",
            "diff_per_kb",
            "density_per_Mb_aligned"
        ]) + "\n")

        keys = sorted(win_acc.keys(), key=lambda x: (x[0], x[1]))
        for (chr_name, wst, wen) in keys:
            d = win_acc[(chr_name, wst, wen)]
            aln = d["sum_aligned"]
            win_mid = int((wst + wen) / 2)
            if aln > 0:
                sub_kb = d["sum_sub"] / aln * 1000.0
                indel_kb = d["sum_indel"] / aln * 1000.0
                diff_kb = d["sum_diff"] / aln * 1000.0
                dens_mb = d["sum_diff"] / (aln / 1e6)
            else:
                sub_kb = math.nan
                indel_kb = math.nan
                diff_kb = math.nan
                dens_mb = math.nan

            w.write("\t".join([
                chr_name,
                str(wst),
                str(wen),
                str(win_mid),
                str(d["n_pairs"]),
                str(aln),
                str(d["sum_sub"]),
                str(d["sum_indel"]),
                str(d["sum_diff"]),
                str(sub_kb),
                str(indel_kb),
                str(diff_kb),
                str(dens_mb),
            ]) + "\n")

    log(f"DONE | pair_rows={len(pair_rows)} chr_bins={len(chr_acc)} win_bins={len(win_acc)}")
    log(f"Wrote: {out_pair}")
    log(f"Wrote: {out_chr}")
    log(f"Wrote: {out_win}")

if __name__ == "__main__":
    main()
