#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
功能：
1) 读取锚点物种 S_constricta_node_cafe/og_up.tsv 的 family 列，得到锚点 OG 集合（去表头、去空行、去重）
2) 遍历当前 node 目录下所有“以大写字母开头”的文件夹，读取各自 og_up.tsv 的 family 列
3) 逐物种与锚点集合取交集
4) 输出到一个统一文件夹 node/_isct/：
   - S_constricta.og.list.tsv   （锚点 OG 单列，无表头）
   - isct.long.tsv              （长表：species, OG）
   - isct.matrix.tsv            （矩阵表：OG + 每个物种 0/1）
   - summary.tsv                （统计：每物种 OG 总数与交集数）

新增（皇上要求）：
5) 基于 isct.matrix.tsv / 内存统计结果，按“排除缢蛏后，其它物种中≥阈值比例为0”筛选 OG：
   - threshold_summary.tsv
   - og_threshold_table.tsv
   - og_zero_ge{T}pct_others.list   (T=100/80/60...)
"""

import os
import sys
import math


# =========================
# 皇上在这里填参数（不走命令行）
# =========================
NODE_DIR = "."  # 你现在所在的 node 目录；一般不用改
ANCHOR_FOLDER = "S_constricta_node_cafe"
OG_UP_FILENAME = "og_up.tsv"
FAMILY_COLNAME = "family"
FOLDER_SUFFIX = "_node_cafe"  # 皇上确认：一定是这个结尾
OUTPUT_SUBDIR = "_isct"       # 输出统一放到 node/_isct/ 下

# 矩阵里是否包含锚点物种列（通常包含，且你已确认一定包含）
INCLUDE_ANCHOR_AS_TARGET_COLUMN = True

# 新增：阈值（表示“其它物种中至少多少比例为0”）
# 例如 100 表示：排除缢蛏列后，所有其它物种都为0
ZERO_PCT_THRESHOLDS = [100, 80, 60]

# 新增：是否输出每个阈值的 OG 列表文件
WRITE_THRESHOLD_LIST_FILES = True


# =========================
# 工具函数
# =========================
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def read_family_og_set(tsv_path: str, family_colname: str) -> set:
    """
    读取 TSV 文件，按表头找到 family 列，返回 OG 的集合（去空、去重）
    """
    ogs = set()
    with open(tsv_path, "r", encoding="utf-8") as f:
        header = f.readline()
        if not header:
            raise RuntimeError(f"文件为空：{tsv_path}")

        header_cols = header.rstrip("\n").split("\t")
        try:
            idx = header_cols.index(family_colname)
        except ValueError:
            raise RuntimeError(
                f"找不到列名 '{family_colname}'：{tsv_path}\n"
                f"表头列为：{header_cols}"
            )

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if idx >= len(cols):
                # 这里不做“复杂防御”，遇到不完整行直接跳过
                continue
            og = cols[idx].strip()
            if og:
                ogs.add(og)
    return ogs


def og_sort_key(og: str):
    """
    让 OG0000012 这种按数字排序更自然；不影响结果，只影响输出顺序。
    """
    if og.startswith("OG"):
        tail = og[2:]
        if tail.isdigit():
            return (0, int(tail))
    return (1, og)


def list_target_folders(node_dir: str) -> list:
    """
    找到所有以大写字母开头的文件夹
    """
    items = []
    for name in os.listdir(node_dir):
        full = os.path.join(node_dir, name)
        if os.path.isdir(full) and name and name[0].isupper():
            items.append(name)
    items.sort()
    return items


def strip_suffix(folder_name: str, suffix: str) -> str:
    """
    皇上要求：不做检查，直接去掉末尾 suffix
    """
    return folder_name[:-len(suffix)]


def ceil_int(x: float) -> int:
    """
    向上取整（返回 int）
    """
    return int(math.ceil(x))


# =========================
# 主流程
# =========================
def main() -> int:
    node_dir = os.path.abspath(NODE_DIR)

    anchor_path = os.path.join(node_dir, ANCHOR_FOLDER, OG_UP_FILENAME)
    if not os.path.exists(anchor_path):
        print(f"[ERROR] 找不到锚点文件：{anchor_path}", file=sys.stderr)
        return 1

    out_dir = os.path.join(node_dir, OUTPUT_SUBDIR)
    ensure_dir(out_dir)

    print(f"[CONF] NODE_DIR = {node_dir}")
    print(f"[CONF] ANCHOR   = {anchor_path}")
    print(f"[CONF] OUT_DIR  = {out_dir}")

    # 1) 锚点 OG 集合
    anchor_ogs = read_family_og_set(anchor_path, FAMILY_COLNAME)
    anchor_ogs_sorted = sorted(anchor_ogs, key=og_sort_key)
    print(f"[INFO] Anchor OG unique count = {len(anchor_ogs)}")

    # 输出锚点 OG 单列（无表头）
    anchor_list_path = os.path.join(out_dir, "S_constricta.og.list.tsv")
    with open(anchor_list_path, "w", encoding="utf-8") as w:
        for og in anchor_ogs_sorted:
            w.write(f"{og}\n")
    print(f"[OK] Write anchor OG list -> {anchor_list_path}")

    # 2) 目标文件夹
    target_folders = list_target_folders(node_dir)

    # 皇上已确认矩阵一定包含缢蛏列，这里保持逻辑兼容：
    if not INCLUDE_ANCHOR_AS_TARGET_COLUMN:
        target_folders = [x for x in target_folders if x != ANCHOR_FOLDER]

    # 3) 逐物种取交集
    species_names = []
    species_ogsets = {}          # species -> set(all OG)
    species_intersections = {}   # species -> set(intersection with anchor)

    for folder in target_folders:
        species = strip_suffix(folder, FOLDER_SUFFIX)
        species_names.append(species)

        og_up_path = os.path.join(node_dir, folder, OG_UP_FILENAME)
        if not os.path.exists(og_up_path):
            print(f"[WARN] 缺少文件，跳过：{og_up_path}")
            species_ogsets[species] = set()
            species_intersections[species] = set()
            continue

        ogs = read_family_og_set(og_up_path, FAMILY_COLNAME)
        inter = anchor_ogs.intersection(ogs)

        species_ogsets[species] = ogs
        species_intersections[species] = inter

        print(f"[INFO] {species}: OG={len(ogs)}  intersect={len(inter)}")

    # 4) 输出 long 表
    long_path = os.path.join(out_dir, "isct.long.tsv")
    with open(long_path, "w", encoding="utf-8") as w:
        w.write("species\tOG\n")
        for species in species_names:
            inter_sorted = sorted(species_intersections[species], key=og_sort_key)
            for og in inter_sorted:
                w.write(f"{species}\t{og}\n")
    print(f"[OK] Write long table -> {long_path}")

    # 5) 输出 matrix 表（行=锚点OG，列=物种，值=0/1）
    matrix_path = os.path.join(out_dir, "isct.matrix.tsv")
    with open(matrix_path, "w", encoding="utf-8") as w:
        w.write("OG\t" + "\t".join(species_names) + "\n")
        for og in anchor_ogs_sorted:
            row = [og]
            for species in species_names:
                row.append("1" if og in species_intersections[species] else "0")
            w.write("\t".join(row) + "\n")
    print(f"[OK] Write matrix table -> {matrix_path}")

    # 6) 输出 summary
    summary_path = os.path.join(out_dir, "summary.tsv")
    with open(summary_path, "w", encoding="utf-8") as w:
        w.write("species\tn_og_total\tn_og_intersect_with_anchor\tn_anchor\n")
        for species in species_names:
            w.write(
                f"{species}\t{len(species_ogsets[species])}\t{len(species_intersections[species])}\t{len(anchor_ogs)}\n"
            )
    print(f"[OK] Write summary -> {summary_path}")

    # =========================
    # 7) 新增：阈值筛选（排除缢蛏列后）
    # =========================
    anchor_species = strip_suffix(ANCHOR_FOLDER, FOLDER_SUFFIX)

    # 其它物种列表（排除缢蛏）
    other_species = [s for s in species_names if s != anchor_species]
    n_other = len(other_species)

    # 对每个 OG 统计：其它物种中 ones / zeros / zero_pct
    # 并同时计算各阈值是否通过
    thresholds = list(ZERO_PCT_THRESHOLDS)
    thresholds_sorted = sorted(thresholds, reverse=True)  # 100,80,60...

    # 预计算每个阈值需要的最少 zeros 数
    thr2_min_zeros = {}
    for t in thresholds_sorted:
        min_zeros = ceil_int((t / 100.0) * n_other) if n_other > 0 else 0
        thr2_min_zeros[t] = min_zeros

    # 输出：og_threshold_table.tsv
    thr_table_path = os.path.join(out_dir, "og_threshold_table.tsv")
    with open(thr_table_path, "w", encoding="utf-8") as w:
        # 表头
        cols = ["OG", "ones_other", "zeros_other", "zero_pct_other"]
        cols += [f"pass_zero_ge{t}pct_others" for t in thresholds_sorted]
        w.write("\t".join(cols) + "\n")

        # 统计每个阈值通过数量
        thr_pass_counts = {t: 0 for t in thresholds_sorted}

        # 如果还要写 list 文件，就先准备容器
        thr2_oglist = {t: [] for t in thresholds_sorted}

        for og in anchor_ogs_sorted:
            ones = 0
            for sp in other_species:
                if og in species_intersections.get(sp, set()):
                    ones += 1
            zeros = n_other - ones
            zero_pct = (zeros / n_other * 100.0) if n_other > 0 else 0.0

            pass_flags = []
            for t in thresholds_sorted:
                ok = 1 if zeros >= thr2_min_zeros[t] else 0
                pass_flags.append(ok)
                if ok:
                    thr_pass_counts[t] += 1
                    if WRITE_THRESHOLD_LIST_FILES:
                        thr2_oglist[t].append(og)

            row = [og, str(ones), str(zeros), f"{zero_pct:.2f}"] + [str(x) for x in pass_flags]
            w.write("\t".join(row) + "\n")

    print(f"[OK] Write OG threshold table -> {thr_table_path}")

    # 输出：threshold_summary.tsv
    thr_summary_path = os.path.join(out_dir, "threshold_summary.tsv")
    with open(thr_summary_path, "w", encoding="utf-8") as w:
        w.write("threshold_zero_pct\tother_species_count\tmin_zeros_required\tn_og_pass\n")
        for t in thresholds_sorted:
            w.write(f"{t}\t{n_other}\t{thr2_min_zeros[t]}\t{thr_pass_counts[t]}\n")
    print(f"[OK] Write threshold summary -> {thr_summary_path}")

    # 输出：每个阈值一个 OG 列表（一列，无表头）
    if WRITE_THRESHOLD_LIST_FILES:
        for t in thresholds_sorted:
            out_list = os.path.join(out_dir, f"og_zero_ge{t}pct_others.list")
            with open(out_list, "w", encoding="utf-8") as w:
                for og in thr2_oglist[t]:
                    w.write(f"{og}\n")
            print(f"[OK] Write threshold OG list -> {out_list}")

    print("[DONE] All outputs are in:", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

