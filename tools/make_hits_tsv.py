#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 先使用grep -Ff 基因列表.txt Orthogroups.tsv > hits.tsv
# 然后再使用这个脚本，注意基因列表.txt是一列基因
"""
make_hits_tsv.py
把 grep 得到的 hits.tsv 转成矩阵 TSV：
- 第一列：OG
- 第一行：物种名（取 "species|gene" 里 | 前，按首次出现顺序去重）
- 单元格：基因ID（取 | 后），多个用 ", " 拼
- 缺失：NULL
"""

from __future__ import annotations

import os
from collections import OrderedDict, defaultdict


# =========================
# 参数区（皇上只改这里）
# =========================
INPUT_TSV = "hits.tsv"
OUTPUT_TSV = "hits.matrix.tsv"
MISSING_VALUE = "NULL"


def iter_species_gene_tokens(line: str):
    """
    从一行 hits.tsv 中提取所有 "species|gene" token。
    兼容：tab 分列；单元格内多个条目用逗号分隔。
    """
    parts = line.rstrip("\n").split("\t")
    if not parts:
        return
    for cell in parts[1:]:
        if not cell:
            continue
        cell = cell.strip()
        if not cell or "|" not in cell:
            continue
        for item in cell.split(","):
            item = item.strip()
            if not item or "|" not in item:
                continue
            sp, gene = item.split("|", 1)
            sp = sp.strip()
            gene = gene.strip()
            if sp and gene:
                yield sp, gene


def main():
    if not os.path.exists(INPUT_TSV):
        raise FileNotFoundError(f"找不到输入文件：{INPUT_TSV}")

    # 第一遍：收集 OG 顺序与物种顺序
    og_order = []
    seen_og = set()
    species_order = []
    seen_species = set()

    with open(INPUT_TSV, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            og = line.split(None, 1)[0].strip()
            if og and og not in seen_og:
                seen_og.add(og)
                og_order.append(og)

            for sp, _ in iter_species_gene_tokens(line):
                if sp not in seen_species:
                    seen_species.add(sp)
                    species_order.append(sp)

    # 第二遍：构建 OG×species 的 gene 列表
    og2sp2genes = OrderedDict((og, defaultdict(list)) for og in og_order)

    with open(INPUT_TSV, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            og = line.split(None, 1)[0].strip()
            if not og:
                continue
            if og not in og2sp2genes:
                og2sp2genes[og] = defaultdict(list)

            for sp, gene in iter_species_gene_tokens(line):
                og2sp2genes[og][sp].append(gene)

    # 写出矩阵 TSV
    with open(OUTPUT_TSV, "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(["OG"] + species_order) + "\n")
        for og in og_order:
            row = [og]
            sp2genes = og2sp2genes[og]
            for sp in species_order:
                genes = sp2genes.get(sp, [])
                row.append(", ".join(genes) if genes else MISSING_VALUE)
            out.write("\t".join(row) + "\n")

    print(f"[OK] 写出：{OUTPUT_TSV}（rows={len(og_order)}, species={len(species_order)}）")


if __name__ == "__main__":
    main()

