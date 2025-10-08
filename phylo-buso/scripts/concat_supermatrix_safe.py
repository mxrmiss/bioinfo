#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
纯 Python 安全拼接超矩阵（不依赖 AMAS）
- 输入：--trim-dir 下的一堆对齐后的蛋白序列文件（支持 .faa/.fa/.fasta 及其 .gz）
- 约定：同一 locus 文件内所有序列对齐长度相同；物种可能缺失某些 locus
- 输出：
    <out_prefix>.aa.fas 或 <out_prefix>.dna.fas     # 拼接后的超矩阵
    <out_dir>/partitions.txt                         # IQ-TREE 兼容的分区文件（"AA," 或 "DNA," 语法）
用法示例：
    python scripts/concat_supermatrix_safe.py \
        --trim-dir results/trim \
        --out-prefix results/concat/supermatrix \
        --datatype aa
"""

import os
import re
import sys
import glob
import gzip
from collections import OrderedDict

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode, encoding=None if "b" in mode else "utf-8", errors="ignore")
    return open(path, mode=mode, encoding=None if "b" in mode else "utf-8", errors="ignore")

def read_fasta(path):
    """极简 FASTA 读取：返回 [(header, seq), ...]；忽略空行。"""
    recs, hdr, seq = [], None, []
    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if hdr is not None:
            recs.append((hdr, "".join(seq)))
    return recs

def strip_exts(name):
    """去掉 .trim 以及常见扩展，得到 locus 名。"""
    base = re.sub(r"\.gz$", "", name, flags=re.IGNORECASE)
    base = re.sub(r"\.(fa|faa|fasta)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"\.trim$", "", base, flags=re.IGNORECASE)
    return base

def main():
    import argparse
    ap = argparse.ArgumentParser(description="纯 Python 拼接超矩阵（不依赖 AMAS）")
    ap.add_argument("--trim-dir", required=True, help="修剪后的对齐目录（*.trim.fa* 或 .gz）")
    ap.add_argument("--out-prefix", required=True, help="输出前缀（会生成 <prefix>.aa.fas/.dna.fas 与 partitions.txt）")
    ap.add_argument("--datatype", choices=["aa", "dna"], default="aa", help="序列类型：aa 或 dna（决定分区文件行头与输出文件后缀）")
    args = ap.parse_args()

    trim_dir   = os.path.abspath(args.trim_dir)
    out_prefix = os.path.abspath(args.out_prefix)
    out_dir    = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)
    out_fas  = out_prefix + (".aa.fas" if args.datatype == "aa" else ".dna.fas")
    out_part = os.path.join(out_dir, "partitions.txt")

    # 1) 收集所有输入文件
    exts = (".trim.faa", ".trim.fa", ".trim.fasta",
            ".trim.faa.gz", ".trim.fa.gz", ".trim.fasta.gz")
    files = []
    for ext in exts:
        files.extend(glob.glob(os.path.join(trim_dir, f"*{ext}")))
    files = sorted(set(files))
    if not files:
        # 与 AMAS 行为对齐：写空产物，方便流程继续
        open(out_part, "w", encoding="utf-8").close()
        open(out_fas, "w", encoding="utf-8").close()
        sys.stderr.write(f"[WARN] 在 {trim_dir} 未找到 *.trim.fa*，已写出空文件：\n  {out_part}\n  {out_fas}\n")
        return 0

    # 2) 第一次遍历：读取每个 locus，记录长度与物种->序列
    loci = []  # [(locus_name, L, {species: seq}), ...]
    all_species = set()
    for fp in files:
        recs = read_fasta(fp)
        if not recs:
            continue
        # 允许 header 如 'Species|gene', 取第一个 token，再取 '|' 前缀当物种名
        seqs = OrderedDict()
        lens = set()
        for hdr, seq in recs:
            sp = hdr.split()[0].split("|")[0]
            seq = seq.strip()
            seqs[sp] = seq
            lens.add(len(seq))
            all_species.add(sp)
        if len(lens) != 1:
            # 若同一 locus 内不同长度，说明前序步骤没对齐好；给出清晰报错
            msg = f"[ERROR] {os.path.basename(fp)} 内序列长度不一致：{sorted(lens)}。请检查上游比对/修剪。"
            sys.stderr.write(msg + "\n")
            return 2
        L = lens.pop()
        loci.append((strip_exts(os.path.basename(fp)), L, seqs))

    if not loci:
        open(out_part, "w", encoding="utf-8").close()
        open(out_fas, "w", encoding="utf-8").close()
        sys.stderr.write(f"[WARN] 可读的对齐文件为空，已写出空文件：\n  {out_part}\n  {out_fas}\n")
        return 0

    # 3) 固定物种与 locus 顺序（确定性）
    species = sorted(all_species)
    # 4) 为每个物种拼接所有片段（缺失用 '-' 填充）
    gap_blocks = {}
    for _, L, _ in loci:
        if L not in gap_blocks:
            gap_blocks[L] = "-" * L

    concat = {sp: [] for sp in species}
    starts, ends, pos = [], [], 1  # 记录分区起止
    for locus_name, L, seqs in loci:
        for sp in species:
            concat[sp].append(seqs.get(sp, gap_blocks[L]))
        starts.append(pos)
        ends.append(pos + L - 1)
        pos += L

    # 5) 写出 partitions.txt（IQ-TREE 兼容的简单格式）
    part_type = "AA" if args.datatype == "aa" else "DNA"
    with open(out_part, "w", encoding="utf-8") as fh:
        for (locus_name, _, _), s, e in zip(loci, starts, ends):
            fh.write(f"{part_type}, {locus_name} = {s}-{e}\n")

    # 6) 写出超矩阵 FASTA
    with open(out_fas, "w", encoding="utf-8") as fh:
        for sp in species:
            fh.write(f">{sp}\n")
            seq = "".join(concat[sp])
            # 为可读性每 60 列换行（IQ-TREE 不强制，但更友好）
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    sys.stderr.write(
        f"[OK] 拼接完成：物种 {len(species)}，locus {len(loci)}，总位点 {ends[-1]}。\n"
        f"     输出：\n       {out_fas}\n       {out_part}\n"
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
