#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proteome_qc_simple_pct.py —— OrthoFinder 前的蛋白质(最长转录本) QC：精简版（百分数显示）

============================================================
输出文件（TSV）默认列名与含义（每列都可直接用于表格/论文）
============================================================

1) File
   - 蛋白质FASTA文件名（通常对应一个物种）

2) ProteinCount
   - FASTA中序列条目数（即 '>' 行数）

3) LenMin
   - 最短蛋白长度（aa）

4) LenP05
   - 蛋白长度第5百分位（aa）：反映“短蛋白尾部”是否异常

5) LenMedian
   - 蛋白长度中位数（aa）

6) LenP95
   - 蛋白长度第95百分位（aa）：反映“长蛋白尾部”是否异常

7) LenMax
   - 最长蛋白长度（aa）

8) Short<100_%
   - 长度 <100 aa 的蛋白比例（Short<100_N / ProteinCount * 100）

9) Short<200_%
   - 长度 <200 aa 的蛋白比例（Short<200_N / ProteinCount * 100）

10) Long>5000_%
    - 长度 >5000 aa 的蛋白比例（Long>5000_N / ProteinCount * 100）
      该比例若明显高，常提示拼接/错误ORF/异常延伸

11) DotSuffix_%
    - FASTA header 第一个字段以 ".数字" 结尾的比例
      例：Gene123.1 / ENSXXX.2
      用于辅助判断ID风格差异（不等同于一定是isoform）

12) IDRule
    - smart 规则判定：是否建议“去掉末尾 .数字”
      * strip_dotnum：去掉末尾 .数字 后，唯一ID数不会大幅塌方（更像isoform/版本号）
      * raw：去掉末尾 .数字 会导致唯一ID数严重塌方（更像ID的一部分，不能乱剥）

============================================================
注意
============================================================
- 本脚本只统计，不改动任何输入文件。
- 脚本不接收命令行参数；所有参数集中在下面“参数区”，皇上直接改那里即可。
"""

from __future__ import annotations

import os
import re
import glob
import gzip
from typing import List


# ============================================================
# 参数区（皇上只需要改这里）
# ============================================================

# 蛋白质文件夹：默认当前目录（推荐先 cd 到 proteomes 再运行）
PROTEOME_DIR = "."

# 识别的蛋白文件（按需增减）
FILE_GLOBS = [
    "*.faa", "*.fa", "*.fasta", "*.pep",
    "*.faa.gz", "*.fa.gz", "*.fasta.gz", "*.pep.gz",
]

# smart规则阈值：
# ratio = (去掉末尾 .数字 后的唯一ID数) / (原始唯一ID数)
# 若 ratio >= SMART_RATIO => 认为“.数字”更像isoform/版本号，可剥离 => IDRule=strip_dotnum
# 否则 => IDRule=raw
SMART_RATIO = 0.70

# 百分比统计阈值
SHORT_1 = 100
SHORT_2 = 200
LONG_1 = 5000

# 输出文件名（写到 PROTEOME_DIR 下）
OUT_TSV = "qc_summary.simple.pct.tsv"


# ============================================================
# 内部函数（一般不需要改）
# ============================================================

_strip_dotnum_re = re.compile(r"\.\d+$")


def _open_text_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def fasta_ids(path: str) -> List[str]:
    ids: List[str] = []
    with _open_text_maybe_gz(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def fasta_lengths(path: str) -> List[int]:
    lens: List[int] = []
    L = 0
    with _open_text_maybe_gz(path) as f:
        for line in f:
            if line.startswith(">"):
                if L:
                    lens.append(L)
                L = 0
            else:
                L += len(line.strip())
    if L:
        lens.append(L)
    return lens


def quantile(xs: List[int], q: float) -> float:
    if not xs:
        return 0.0
    xs = sorted(xs)
    pos = (len(xs) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(xs) - 1)
    frac = pos - lo
    return xs[lo] * (1 - frac) + xs[hi] * frac


def strip_dotnum(s: str) -> str:
    return _strip_dotnum_re.sub("", s)


def smart_rule_mode(ids: List[str], smart_ratio: float) -> str:
    """只返回 smart 判定的模式：raw 或 strip_dotnum"""
    if not ids:
        return "raw"
    uniq_raw = len(set(ids))
    if uniq_raw == 0:
        return "raw"
    ids_strip = [strip_dotnum(x) for x in ids]
    uniq_strip = len(set(ids_strip))
    ratio = uniq_strip / uniq_raw
    return "strip_dotnum" if ratio >= smart_ratio else "raw"


def list_proteome_files(proteome_dir: str, globs: List[str]) -> List[str]:
    paths = []
    for pat in globs:
        paths.extend(glob.glob(os.path.join(proteome_dir, pat)))
    return sorted(set(paths))


# ============================================================
# 主程序
# ============================================================

def main():
    proteome_dir = os.path.abspath(PROTEOME_DIR)
    files = list_proteome_files(proteome_dir, FILE_GLOBS)

    if not files:
        raise SystemExit(
            f"[ERROR] 在目录中找不到符合通配符的蛋白文件：{proteome_dir}\n"
            f"        FILE_GLOBS={FILE_GLOBS}"
        )

    header = [
        "File",
        "ProteinCount",
        "LenMin",
        "LenP05",
        "LenMedian",
        "LenP95",
        "LenMax",
        f"Short<{SHORT_1}_%",
        f"Short<{SHORT_2}_%",
        f"Long>{LONG_1}_%",
        "DotSuffix_%",
        "IDRule",
    ]

    out_lines = ["\t".join(header)]

    for path in files:
        fn = os.path.basename(path)

        ids = fasta_ids(path)
        n = len(ids)

        lens = fasta_lengths(path)
        if lens:
            len_min = min(lens)
            len_max = max(lens)
            len_p05 = quantile(lens, 0.05)
            len_med = quantile(lens, 0.50)
            len_p95 = quantile(lens, 0.95)
        else:
            len_min = len_max = 0
            len_p05 = len_med = len_p95 = 0.0

        short1 = sum(1 for L in lens if L < SHORT_1)
        short2 = sum(1 for L in lens if L < SHORT_2)
        long1 = sum(1 for L in lens if L > LONG_1)

        pct_short1 = (short1 / n * 100.0) if n else 0.0
        pct_short2 = (short2 / n * 100.0) if n else 0.0
        pct_long1 = (long1 / n * 100.0) if n else 0.0

        dot_suffix = sum(1 for x in ids if _strip_dotnum_re.search(x) is not None)
        pct_dot_suffix = (dot_suffix / n * 100.0) if n else 0.0

        mode = smart_rule_mode(ids, SMART_RATIO)

        row = [
            fn,
            str(n),
            str(len_min),
            f"{len_p05:.1f}",
            f"{len_med:.1f}",
            f"{len_p95:.1f}",
            str(len_max),
            f"{pct_short1:.2f}",
            f"{pct_short2:.2f}",
            f"{pct_long1:.3f}",
            f"{pct_dot_suffix:.2f}",
            mode,
        ]
        out_lines.append("\t".join(row))

    print("\n".join(out_lines))

    out_path = os.path.join(proteome_dir, OUT_TSV)
    with open(out_path, "w", encoding="utf-8") as w:
        w.write("\n".join(out_lines) + "\n")

    print(f"\n[OK] wrote: {out_path}")


if __name__ == "__main__":
    main()

