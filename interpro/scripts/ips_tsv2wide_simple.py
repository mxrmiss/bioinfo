#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ips_tsv2wide_simple.py
================================================================================
功能（极简测试版）：
  - 读取 InterProScan 5 的 TSV（无表头，一条命中一行）
  - 仅保留 7 个常规库：Pfam、SMART、ProSitePatterns、ProSiteProfiles、SUPERFAMILY、PANTHER、Gene3D
  - 输出“宽表”：一行一个 protein_id，每个库一列，单元格内容为 accession|description 的去重列表

输出格式（简洁模式）：
  Pfam 列举例：PF00069|Protein kinase domain; PF07714|PKinase_Tyr

注意：
  - 不输出 GO/Pathway（你说后续不用 -goterms/-pa）
  - 默认按 accession 字母序排序
================================================================================
"""

import os
import sys
import csv
from collections import defaultdict

# ----------------------------- 参数区（皇上可改） -----------------------------
IN_TSV = "Sinonovacula_constricta.faa.tsv"

# 输出文件名（放在同目录）
OUT_WIDE = "annotation.wide.tsv"

# 你选的 7 个库（列名固定就用这些）
TARGET_ANALYSES = [
    "Pfam",
    "SMART",
    "ProSitePatterns",
    "ProSiteProfiles",
    "SUPERFAMILY",
    "PANTHER",
    "Gene3D",
]

# 单元格内多条注释的分隔符
JOIN_SEP = "; "

# 描述里如果包含分隔符，替换掉，避免破坏表格可读性
REPLACE_SEP_IN_DESC = True
# ---------------------------------------------------------------------------


def die(msg: str, code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def norm_desc(s: str) -> str:
    s = (s or "").strip()
    if REPLACE_SEP_IN_DESC:
        s = s.replace("\t", " ").replace(";", ",")
    # 压缩多余空格
    s = " ".join(s.split())
    return s


def main() -> None:
    if not os.path.isfile(IN_TSV) or os.path.getsize(IN_TSV) == 0:
        die(f"输入文件不存在或为空：{IN_TSV}")

    target_set = set(TARGET_ANALYSES)

    # 存储结构：
    # hits[protein_id][analysis][accession] = best_description
    hits = defaultdict(lambda: defaultdict(dict))
    proteins_seen = set()

    # InterProScan 5 TSV 常见列：0 protein, 3 analysis, 4 accession, 5 description
    # 这里用“固定列位”做极简版（你是测试用，先跑起来）
    with open(IN_TSV, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            # 跳过可能的注释行
            if row[0].startswith("#"):
                continue
            if len(row) < 6:
                # 极简版：列不足直接跳过（也可以 die）
                continue

            protein_id = (row[0] or "").strip()
            analysis = (row[3] or "").strip()
            accession = (row[4] or "").strip()
            desc = norm_desc(row[5])

            if not protein_id:
                continue

            proteins_seen.add(protein_id)

            if analysis not in target_set:
                continue
            if not accession:
                continue

            # 去重：同一 protein + analysis + accession 只保留一个描述
            # 若描述不同，保留更长的（信息更完整的概率更高）
            old = hits[protein_id][analysis].get(accession, "")
            if len(desc) > len(old):
                hits[protein_id][analysis][accession] = desc

    if not proteins_seen:
        die("没有读到任何 protein_id，输入 TSV 可能不是标准 InterProScan TSV？")

    # 写宽表
    with open(OUT_WIDE, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        header = ["protein_id"] + TARGET_ANALYSES
        w.writerow(header)

        for pid in sorted(proteins_seen):
            row_out = [pid]
            for ana in TARGET_ANALYSES:
                acc2desc = hits.get(pid, {}).get(ana, {})
                if not acc2desc:
                    row_out.append("")
                    continue

                items = []
                for acc in sorted(acc2desc.keys()):
                    d = acc2desc[acc]
                    if d:
                        items.append(f"{acc}|{d}")
                    else:
                        items.append(acc)

                row_out.append(JOIN_SEP.join(items))
            w.writerow(row_out)

    print(f"[OK] 输入：{IN_TSV}")
    print(f"[OK] 输出：{OUT_WIDE}")
    print(f"[OK] proteins：{len(proteins_seen)}")


if __name__ == "__main__":
    main()

