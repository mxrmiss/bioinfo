#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
输入文件：sample.tsv，差异表达表格
功能：
- 从“foot 特异基因表”中按指定指标挑选 Top N 基因（默认 50，可手动改）
- 读取 samples.tsv 将 SRR 映射到组织名称（gill/foot/...）
- 画 6 组织表达热图，输出 PNG + PDF

特点：
- 不使用命令行参数，所有参数集中在脚本顶部，皇上直接改就行
- 支持 log2(TPM+1) 变换
- 支持按“每个基因行 z-score”突出组织特异性（强烈推荐）
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############################################
# ========== 皇上手动参数区（只改这里） ==========
############################################

# 1) 输入文件
EXPR_TSV = "foot_enriched_candidates.tsv"     # 你的“foot特异”结果表（含 gene_id + 6个SRR列 + fold等）
SAMPLES_TSV = "../../data/samples.tsv"       # 你的 samples.tsv（sample/group/fastq1/fastq2）

# 2) 输出前缀（会生成 .png 和 .pdf）
OUT_PREFIX = "results/heatmap/foot_top_genes"

# 3) 选择 Top N 基因（皇上想看多少行就改多少）
TOP_N = 50

# 4) Top N 的排序依据（建议用 fold_target_vs_bgmax）
#    可选： "fold_target_vs_bgmax"  或  "target_tpm"  或你表里任何数值列名
RANK_BY = "fold_target_vs_bgmax"

# 5) 目标组织名（必须和 samples.tsv 里的 group 一致）
TARGET_TISSUE = "foot"

# 6) 列（组织）显示顺序：想固定顺序就写；想用 samples.tsv 顺序就设为 None
#    注意：这里填的是 group（组织名），不是 SRR
COLUMN_ORDER = ["foot", "gill", "adductor", "visceral", "mantle", "siphon"]

# 7) 表达量变换
DO_LOG2 = True             # True: log2(TPM+1)
DO_ROW_ZSCORE = True       # True: 每行做 z-score（突出特异性）；False: 显示绝对表达

# 8) 视觉参数
SHOW_GENE_LABELS = True    # Top_N 很大时可改 False
FIG_WIDTH = 8.5            # 图宽（英寸）
ROW_HEIGHT = 0.18          # 每行高度（英寸），Top_N大可适当调小
DPI = 300

############################################
# ========== 不要改下面内容 ==========
############################################

def _ensure_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d and (not os.path.exists(d)):
        os.makedirs(d, exist_ok=True)

def _safe_log2(x: np.ndarray) -> np.ndarray:
    return np.log2(x + 1.0)

def _row_zscore(mat: np.ndarray) -> np.ndarray:
    # 对每行： (x - mean) / std；std=0 时置 0，避免 NaN
    mean = mat.mean(axis=1, keepdims=True)
    std = mat.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (mat - mean) / std

def main():
    # 读 samples.tsv：获取 SRR -> tissue(group)
    try:
        smp = pd.read_csv(SAMPLES_TSV, sep="\t", dtype=str)
    except Exception as e:
        print(f"[ERROR] 无法读取 samples.tsv: {SAMPLES_TSV}\n{e}", file=sys.stderr)
        sys.exit(1)

    for col in ["sample", "group"]:
        if col not in smp.columns:
            print(f"[ERROR] samples.tsv 缺少列: {col}", file=sys.stderr)
            sys.exit(1)

    sample_to_group = dict(zip(smp["sample"], smp["group"]))

    # 读表达表
    try:
        df = pd.read_csv(EXPR_TSV, sep="\t")
    except Exception as e:
        print(f"[ERROR] 无法读取表达表: {EXPR_TSV}\n{e}", file=sys.stderr)
        sys.exit(1)

    if "gene_id" not in df.columns:
        print("[ERROR] 表中必须有 gene_id 列", file=sys.stderr)
        sys.exit(1)

    if RANK_BY not in df.columns:
        print(f"[ERROR] 排序列 RANK_BY='{RANK_BY}' 不在表头中", file=sys.stderr)
        print(f"[HINT] 你当前表头有：{', '.join(df.columns[:30])} ...", file=sys.stderr)
        sys.exit(1)

    # 确定表达矩阵的 SRR 列：用 samples.tsv 中的 sample 去匹配表达表的列
    srr_cols = [s for s in smp["sample"].tolist() if s in df.columns]
    if len(srr_cols) == 0:
        print("[ERROR] 在表达表中找不到任何 SRR 列，请确认 samples.tsv 的 sample 名与表达表列名一致", file=sys.stderr)
        sys.exit(1)

    # 确定列（组织）显示顺序
    if COLUMN_ORDER is None:
        # 按 samples.tsv 出现顺序显示
        ordered_srr = srr_cols
    else:
        # 按指定组织顺序拼接 SRR
        ordered_srr = []
        for tissue in COLUMN_ORDER:
            for srr in srr_cols:
                if sample_to_group.get(srr, "") == tissue:
                    ordered_srr.append(srr)

        # 把未覆盖到的 SRR（如果有）追加在后面，避免漏列
        for srr in srr_cols:
            if srr not in ordered_srr:
                ordered_srr.append(srr)

    # 过滤：只保留 target_tissue 存在的列
    if TARGET_TISSUE not in set(sample_to_group.values()):
        print(f"[ERROR] TARGET_TISSUE='{TARGET_TISSUE}' 不在 samples.tsv 的 group 列中", file=sys.stderr)
        sys.exit(1)

    # 选 Top N 基因
    df_ranked = df.sort_values(by=RANK_BY, ascending=False).head(int(TOP_N)).copy()

    # 提取表达矩阵
    expr = df_ranked[ordered_srr].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat = expr.to_numpy(dtype=float)

    # 变换
    if DO_LOG2:
        mat = _safe_log2(mat)

    # 行标准化（突出组织特异性）
    if DO_ROW_ZSCORE:
        mat_plot = _row_zscore(mat)
        cmap = "RdBu_r"
        cbar_label = "Row Z-score"
        value_label = "log2(TPM+1) then row z-score" if DO_LOG2 else "row z-score"
    else:
        mat_plot = mat
        cmap = "viridis"
        cbar_label = "Expression"
        value_label = "log2(TPM+1)" if DO_LOG2 else "TPM"

    # 组织标签（列标签）用 tissue 名
    col_labels = [sample_to_group.get(srr, srr) for srr in ordered_srr]
    gene_labels = df_ranked["gene_id"].tolist()

    # 画图尺寸：按行数自适应高度
    fig_h = max(3.0, ROW_HEIGHT * len(gene_labels))
    fig_w = float(FIG_WIDTH)

    _ensure_dir(OUT_PREFIX)

    plt.figure(figsize=(fig_w, fig_h))
    ax = plt.gca()

    im = ax.imshow(mat_plot, aspect="auto", interpolation="nearest", cmap=cmap)

    ax.set_title(f"Top {len(gene_labels)} {TARGET_TISSUE}-enriched genes ({RANK_BY})\n{value_label}", fontsize=12)

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")

    if SHOW_GENE_LABELS:
        ax.set_yticks(np.arange(len(gene_labels)))
        ax.set_yticklabels(gene_labels, fontsize=7)
    else:
        ax.set_yticks([])
        ax.set_yticklabels([])

    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(cbar_label)

    plt.tight_layout()

    out_png = f"{OUT_PREFIX}.png"
    out_pdf = f"{OUT_PREFIX}.pdf"
    plt.savefig(out_png, dpi=DPI)
    plt.savefig(out_pdf)
    plt.close()

    print("[OK] Heatmap saved:")
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()

