#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_proteins.py
对蛋白质 fasta 文件进行质控过滤：
- 去除过短序列
- 去除内部终止密码子（*）
- 去除含 X 过多的序列
- 去除低复杂度（香农熵低）的序列
"""

import sys, math, collections, argparse

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="输入 fasta 文件")
    ap.add_argument("output", help="输出过滤后的 fasta 文件")
    ap.add_argument("--min-len", type=int, default=50, help="最短长度")
    ap.add_argument("--min-entropy", type=float, default=2.2, help="最低香农熵")
    ap.add_argument("--max-x-frac", type=float, default=0.05, help="最大X比例")
    ap.add_argument("--allow-stop", action="store_true", help="是否允许内部*")
    return ap.parse_args()

def fasta_iter(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).strip()
                name, seq = line.strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq).strip()

def entropy(seq):
    counts = collections.Counter([c for c in seq if c not in ("X","*")])
    n = sum(counts.values())
    if n == 0: return 0.0
    ent = 0.0
    for k in counts.values():
        p = k / n
        ent -= p * math.log(p, 2)
    return ent

def main():
    args = parse_args()
    kept = 0
    with open(args.output, "w") as fo:
        for name, seq in fasta_iter(args.input):
            if len(seq) < args.min_len: 
                continue
            if not args.allow_stop and ("*" in seq[:-1]): 
                continue  # 允许末尾终止，不允许内部
            xfrac = (seq.count("X")+seq.count("B")+seq.count("Z")) / max(1,len(seq))
            if xfrac > args.max_x_frac: 
                continue
            if entropy(seq) < args.min_entropy: 
                continue
            fo.write(f"{name}\n{seq}\n")
            kept += 1
    sys.stderr.write(f"[filter_proteins.py] kept={kept}\n")

if __name__ == "__main__":
    main()
