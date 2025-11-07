#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
纯 Python 的超矩阵拼接脚本（安全版；可替代 AMAS 的常用功能）
- 输入：修剪后的多序列比对文件目录（*.trim.fa / *.trim.faa / *.trim.fasta / *.trim.fas，亦支持 .gz）
- 逻辑：逐个分区读取对齐 → 统一物种集合与顺序 → 对缺失物种补 gap → 依序拼接 → 写出 supermatrix 与 partitions.txt
- 输出：
    <out_prefix>.aa.fas 或 <out_prefix>.dna.fas   （根据 --datatype）
    partitions.txt                                 （分区坐标，AA/DNA 行头）

新增功能：
- --min-align-len：最小对齐长度过滤
- --min-taxa-per-locus：最小占有度过滤
- --top-n：仅保留最长的前 N 个分区
- --progress-every / --quiet：进度提示与安静模式
"""

import sys
import os
import gzip
import argparse
from pathlib import Path
from collections import OrderedDict

# -----------------------
# I/O & 工具函数
# -----------------------

def open_any(path):
    """兼容 .gz 与普通文本"""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")

def read_fasta(path):
    """读取 FASTA，yield (header, seq)。header 取第一段（空格前），seq 为大写"""
    name, seqs = None, []
    with open_any(path) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seqs).upper()
                name = line.strip()[1:].split()[0]
                seqs = []
            else:
                seqs.append(line.strip())
        if name is not None:
            yield name, "".join(seqs).upper()

def strip_exts(fname):
    """去掉常见扩展（.gz/.fa/.faa/.fasta/.fas/.trim.* 的组合），保留核心基名作为 locus 名"""
    f = fname
    if f.endswith(".gz"):
        f = f[:-3]
    for ext in (".trim.faa", ".trim.fa", ".trim.fasta", ".trim.fas",
                ".faa", ".fa", ".fasta", ".fas"):
        if f.endswith(ext):
            f = f[: -len(ext)]
            break
    return f

def natural_key(s):
    """用于按文件名做自然排序的 key（确保分区顺序稳定可重复）"""
    import re
    return [int(t) if t.isdigit() else t for t in re.findall(r"\d+|\D+", s)]

# -----------------------
# 主流程
# -----------------------

def main():
    ap = argparse.ArgumentParser(
        description="拼接超矩阵并生成 partitions.txt（支持最小对齐长度/最小占有度/Top-N 过滤）"
    )
    ap.add_argument("--trim-dir", required=True, help="修剪后的对齐目录（*.trim.fa* 或 .gz）")
    ap.add_argument("--out-prefix", required=True, help="输出前缀（会生成 <prefix>.aa.fas/.dna.fas 与 partitions.txt）")
    ap.add_argument("--datatype", choices=["aa", "dna"], default="aa",
                    help="序列类型：aa 或 dna（决定输出文件后缀与分区文件行头）")
    ap.add_argument("--min-align-len", type=int, default=0,
                    help="最小对齐长度阈值（按对齐列数计；0=关闭）")
    ap.add_argument("--min-taxa-per-locus", type=int, default=0,
                    help="最小占有度：该对齐至少有这么多物种含非缺失字符（0=关闭）")
    ap.add_argument("--top-n", type=int, default=0,
                    help="仅保留最长的前 N 个分区（0=关闭）")
    ap.add_argument("--progress-every", type=int, default=1000,
                    help="每处理多少个文件打印一次进度（stderr；0=不打印）")
    ap.add_argument("--quiet", action="store_true",
                    help="安静模式：不逐条打印过滤/警告，仅输出摘要")
    args = ap.parse_args()

    trim_dir = Path(args.trim_dir)
    if not trim_dir.is_dir():
        sys.exit(f"[ERROR] --trim-dir 不存在或不是目录：{trim_dir}")

    # 搜索候选对齐文件（包含常见扩展名），保持稳定顺序
    exts = [".trim.faa", ".trim.fa", ".trim.fasta", ".trim.fas", ".faa", ".fa", ".fasta", ".fas"]
    files = []
    for ext in exts:
        files.extend(sorted(trim_dir.glob(f"*{ext}"), key=lambda p: natural_key(p.name)))
        files.extend(sorted(trim_dir.glob(f"*{ext}.gz"), key=lambda p: natural_key(p.name)))

    # 去重
    seen = set()
    uniq_files = []
    for fp in files:
        if fp.name in seen:
            continue
        seen.add(fp.name)
        uniq_files.append(fp)

    # 空集：写空输出，保证流水线可继续
    out_dir = os.path.dirname(args.out_prefix) or "."
    out_fas = f"{args.out_prefix}.{args.datatype}.fas"
    out_part = os.path.join(out_dir, "partitions.txt")
    if not uniq_files:
        os.makedirs(out_dir, exist_ok=True)
        open(out_fas, "w").close()
        with open(out_part, "w") as g:
            g.write("")
        if not args.quiet:
            print("[WARN] 未发现对齐文件，已写出空的 supermatrix 与 partitions.txt", file=sys.stderr)
        return

    # 先扫描，决定 taxa 集合与各 locus 的属性
    loci = []     # [(locus_name, L, {taxon:seq}, orig_idx), ...]
    all_taxa = set()
    processed = 0
    bad_files = 0
    short_filtered = 0
    occ_filtered = 0

    def has_real_char(s: str) -> bool:
        # 视 '-' 和 '?' 为缺失；'X' 仍算观测
        return any(c not in "-?" for c in s)

    for idx, fp in enumerate(uniq_files):
        processed += 1
        if args.progress_every and processed % args.progress_every == 0:
            print(f"[PROGRESS] {processed}/{len(uniq_files)} files", file=sys.stderr)

        locus = strip_exts(fp.name)
        seqs = OrderedDict()
        lengths = set()

        try:
            for h, s in read_fasta(fp):
                seqs[h] = s
                lengths.add(len(s))
                all_taxa.add(h)
        except Exception as e:
            bad_files += 1
            if not args.quiet:
                print(f"[WARN] 解析失败，跳过：{fp.name} :: {e}", file=sys.stderr)
            continue

        if not seqs:
            if not args.quiet:
                print(f"[INFO] 跳过空对齐：{fp.name}", file=sys.stderr)
            continue

        if len(lengths) != 1:
            sys.exit(f"[ERROR] 对齐列长不一致：{fp.name}，观测到长度集合={lengths}")

        L = lengths.pop()

        # ① 最小长度过滤
        if args.min_align_len > 0 and L < args.min_align_len:
            short_filtered += 1
            if not args.quiet:
                print(f"[INFO] 过滤分区（过短）：{fp.name}  L={L} < {args.min_align_len}", file=sys.stderr)
            continue

        # ② 最小占有度过滤
        if args.min_taxa_per_locus > 0:
            occ = sum(1 for s in seqs.values() if has_real_char(s))
            if occ < args.min_taxa_per_locus:
                occ_filtered += 1
                if not args.quiet:
                    print(f"[INFO] 过滤分区（占有度不足）：{fp.name}  occ={occ} < {args.min_taxa_per_locus}", file=sys.stderr)
                continue

        loci.append((locus, L, seqs, idx))

    if not loci:
        os.makedirs(out_dir, exist_ok=True)
        open(out_fas, "w").close()
        with open(out_part, "w") as g:
            g.write("")
        print(f"[WARN] 所有分区被过滤（min-align-len={args.min_align_len}, "
              f"min-taxa-per-locus={args.min_taxa_per_locus})，已写出空文件。", file=sys.stderr)
        return

    # ③ Top-N：仅保留最长的前 N 个分区（保留原始顺序以保证坐标稳定）
    if args.top_n and args.top_n > 0 and len(loci) > args.top_n:
        loci = sorted(loci, key=lambda x: x[1], reverse=True)[:args.top_n]
        loci = sorted(loci, key=lambda x: x[3])  # 按原始索引恢复原顺序

    # 统一物种顺序（稳定可重复）：按名称字典序
    taxa = sorted(all_taxa)

    # 为每个物种构建拼接序列
    concat = {t: [] for t in taxa}
    starts, ends, names = [], [], []
    cur = 1  # 1-based 坐标

    for locus, L, seqs, _ in loci:
        gap = "-" * L
        for t in taxa:
            concat[t].append(seqs.get(t, gap))
        starts.append(cur)
        ends.append(cur + L - 1)
        names.append(locus)
        cur += L

    # 写超矩阵
    os.makedirs(out_dir, exist_ok=True)
    with open(out_fas, "w") as f:
        for t in taxa:
            f.write(f">{t}\n")
            f.write("".join(concat[t]) + "\n")

    # 写 partitions.txt（IQ-TREE -p 可读）
    kind = "AA" if args.datatype == "aa" else "DNA"
    with open(out_part, "w") as g:
        for nm, s, e in zip(names, starts, ends):
            g.write(f"{kind}, {nm} = {s}-{e}\n")

    totalL = sum(e - s + 1 for s, e in zip(starts, ends))
    kept = len(loci)
    msg = (f"[DONE] taxa={len(taxa)} loci_kept={kept} length={totalL} "
           f"(filtered: short={short_filtered}, occ={occ_filtered}, bad={bad_files}) "
           f"-> {out_fas}\n[INFO] partitions: {out_part}")
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
