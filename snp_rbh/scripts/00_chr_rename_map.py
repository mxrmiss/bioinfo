#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
00_chr_rename_map.py
========================================================
[目标]
为 genome FASTA 的 seqid 生成“标准 chr01 风格”的重命名映射表（TSV带表头）。

[核心规则（皇上指定）]
1) 对长度 >= CHROM_MIN_BP（默认10Mb）的序列，优先从 seqid 末尾提取数字作为染色体编号：
   - 例如：chr17 / Chr17 / HiC_scaffold_17 / scaffold17 / LG17  ->  idx=17  ->  chr17（两位补零：chr17，1-9会变chr01..chr09）
   - 输出统一前缀为 chr，数字两位补零（CHR_PAD_WIDTH=2）

2) 如果候选序列里末尾数字提取失败（比如没有末尾数字，或存在冲突/歧义导致无法可靠使用），
   则回退到“按长度降序自动编号”：
   - 最长 -> chr01，次长 -> chr02 ...

3) 默认只重命名 chromosome candidates（>=阈值）；小于阈值的序列保持原名（更安全可追溯）。
   若 RENAME_NONCHROM=True，可将其命名为 unplaced_0001 ...（不建议默认开）

[输出]
results/step0_chr/<code>.chr_rename_map.tsv
列：
old_seqid  new_chr  is_chrom  length_bp

[进度]
每 100 条序列打印一次
"""

# ----------------------------- 用户自定义区（皇上可改） -----------------------------
# 需要为每个物种各跑一次：把 GENOME_FASTA 改成对应文件即可
GENOME_FASTA = "data/genomic/Sinonovacula_rivularis.fasta.gz"

# 阈值：默认 10 Mb（可手动改）
CHROM_MIN_BP = 10_000_000

# chr 两位补零：chr01 风格
CHR_PAD_WIDTH = 2

# 如果末尾数字模式存在冲突（比如两个不同 seqid 末尾都是 17），是否强制回退到长度排序
FALLBACK_IF_TRAILING_CONFLICT = True

# 是否也重命名非染色体序列（默认 False）
RENAME_NONCHROM = False
NONCHROM_PREFIX = "unplaced_"
NONCHROM_PAD_WIDTH = 4

PROGRESS_EVERY = 100
# ------------------------------------------------------------------------------

import os
import re
import gzip
import sys
from typing import Dict, Tuple, List, Optional

def project_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def infer_code_from_path(path: str) -> str:
    bn = os.path.basename(path)
    for suf in [".gz", ".fasta", ".fa", ".fna"]:
        if bn.endswith(suf):
            bn = bn[: -len(suf)]
    parts = bn.split("_")
    if len(parts) < 2:
        raise ValueError(f"无法从文件名推断物种缩写：{os.path.basename(path)}（需要 Genus_species.*）")
    genus, species = parts[0], parts[1]
    return (genus[0] + species[:2]).lower()

def open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")

def fasta_lengths(path: str, progress_every: int) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    cur_id = None
    cur_len = 0
    n_rec = 0

    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line:
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    lengths[cur_id] = cur_len
                    n_rec += 1
                    if n_rec % progress_every == 0:
                        sys.stderr.write(f"[00_chr_rename_map] parsed {n_rec} seqs...\n")
                cur_id = line[1:].split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())

        if cur_id is not None:
            lengths[cur_id] = cur_len
            n_rec += 1

    sys.stderr.write(f"[00_chr_rename_map] parsed total {n_rec} seqs.\n")
    return lengths

def trailing_int(seqid: str) -> Optional[int]:
    """
    提取末尾连续数字作为编号；无则返回 None
    例：chr17 -> 17；HiC_scaffold_03 -> 3；LG0009 -> 9
    """
    m = re.search(r"(\d+)$", seqid)
    if not m:
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None

def fmt_chr(idx: int, pad: int) -> str:
    return f"chr{idx:0{pad}d}"

def main():
    os.chdir(project_root())
    os.makedirs("results/step0_chr", exist_ok=True)

    code = infer_code_from_path(GENOME_FASTA)
    out_map = f"results/step0_chr/{code}.chr_rename_map.tsv"

    sys.stderr.write(f"[00_chr_rename_map] genome={GENOME_FASTA}\n")
    sys.stderr.write(f"[00_chr_rename_map] CHROM_MIN_BP={CHROM_MIN_BP} CHR_PAD_WIDTH={CHR_PAD_WIDTH}\n")

    lengths = fasta_lengths(GENOME_FASTA, PROGRESS_EVERY)

    chrom = [(sid, l) for sid, l in lengths.items() if l >= CHROM_MIN_BP]
    nonchrom = [(sid, l) for sid, l in lengths.items() if l < CHROM_MIN_BP]

    sys.stderr.write(f"[00_chr_rename_map] chrom_candidates={len(chrom)} total={len(lengths)}\n")

    mapping: Dict[str, str] = {}
    chrom_sorted: List[Tuple[str, int]] = []

    # 1) 尝试末尾数字模式（只针对 chromosome candidates）
    idx_to_seqids: Dict[int, List[str]] = {}
    ok_trailing = True
    for sid, l in chrom:
        idx = trailing_int(sid)
        if idx is None:
            ok_trailing = False
            break
        idx_to_seqids.setdefault(idx, []).append(sid)

    conflict_idxs = [idx for idx, sids in idx_to_seqids.items() if len(sids) > 1]
    if conflict_idxs:
        sys.stderr.write(f"[00_chr_rename_map] WARNING: trailing-number conflicts at idx={conflict_idxs[:10]}\n")
        if FALLBACK_IF_TRAILING_CONFLICT:
            ok_trailing = False

    if len(chrom) == 0:
        ok_trailing = False

    if ok_trailing:
        # 按 idx 升序输出（保证 chr01..chrNN 与末尾数字对应）
        sys.stderr.write("[00_chr_rename_map] mode=trailing_number -> use trailing integer as chromosome index\n")
        chrom_sorted = sorted(chrom, key=lambda x: trailing_int(x[0]) if trailing_int(x[0]) is not None else 10**18)
        for sid, _ in chrom_sorted:
            idx = trailing_int(sid)
            mapping[sid] = fmt_chr(idx, CHR_PAD_WIDTH)
    else:
        # 2) 回退：按长度降序编号
        sys.stderr.write("[00_chr_rename_map] mode=length_rank_fallback -> rank by length desc\n")
        chrom_sorted = sorted(chrom, key=lambda x: (-x[1], x[0]))
        for rank, (sid, _) in enumerate(chrom_sorted, start=1):
            mapping[sid] = fmt_chr(rank, CHR_PAD_WIDTH)

    # 3) 可选：非染色体序列改名
    if RENAME_NONCHROM:
        nonchrom_sorted = sorted(nonchrom, key=lambda x: (-x[1], x[0]))
        for idx, (sid, _) in enumerate(nonchrom_sorted, start=1):
            mapping[sid] = f"{NONCHROM_PREFIX}{idx:0{NONCHROM_PAD_WIDTH}d}"

    # 写表（带表头）
    with open(out_map, "w", encoding="utf-8") as w:
        w.write("old_seqid\tnew_chr\tis_chrom\tlength_bp\n")

        # 先写 chromosome candidates（按当前模式排序），再写非染色体
        chrom_set = set(sid for sid, _ in chrom_sorted)
        for sid, l in chrom_sorted:
            w.write(f"{sid}\t{mapping.get(sid, sid)}\t1\t{l}\n")

        # 非染色体：默认保持原名
        for sid, l in sorted(nonchrom, key=lambda x: x[0]):
            new_name = mapping.get(sid, sid)
            w.write(f"{sid}\t{new_name}\t0\t{l}\n")

    # 额外报告：看看候选chr的 idx 分布是否“看起来像 1..N”
    if ok_trailing:
        idxs = sorted(idx_to_seqids.keys())
        sys.stderr.write(f"[00_chr_rename_map] trailing_idx_min={idxs[0]} max={idxs[-1]} n_unique={len(idxs)}\n")

    sys.stderr.write(f"[00_chr_rename_map] wrote: {out_map}\n")
    sys.stderr.write("[00_chr_rename_map] DONE\n")

if __name__ == "__main__":
    main()
