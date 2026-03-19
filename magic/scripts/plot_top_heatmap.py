#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
功能（皇上最新三输入结构）：
- 输入三个文件：
  1) input/gene_tpm.tsv            ：全基因 TPM 矩阵（gene_id + 18个SRR列）
  2) input/samples.tsv             ：sample/group/fastq1/fastq2（用于 SRR->组织映射）
  3) input/top_gene_heatmap.tsv     ：要画热图的基因清单（GeneID/GeneSymbol/Theme/PutativeFunction）
- 从 gene_tpm.tsv 提取 top_gene_heatmap.tsv 里的基因表达
- 支持两种模式：
  - tissue：按组织聚合成 6 列（median/mean 可选）——主图推荐
  - sample：保留 18 个样本列——补图推荐
- 支持 log2(TPM+1) + 行 z-score
- y轴标签可选 GeneID 或 GeneSymbol（新增开关）
- 支持 3 套配色方案开关 + 反转开关
- 输出 PNG + PDF，并可选导出用于作图的矩阵 TSV

特别修改（皇上最新要求）：
- 当 AGG_MODE="sample" 时：横轴标签显示 foot-1, foot-2, foot-3（其他组织同理），编号按列出现顺序递增。
- 热图去掉横纵轴黑框（隐藏 spines）
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

############################################
# ========== 皇上手动参数区（只改这里） ==========
############################################

# 1) 输入文件（相对 magic/）
GENE_TPM_TSV = "input/gene_tpm.tsv"
SAMPLES_TSV  = "input/samples.tsv"
TOP_TSV      = "input/top_gene_heatmap.tsv"

# 2) 输出
OUT_PREFIX = "output/heatmap/foot34_heatmap"
DPI = 600

# 3) 模式：tissue（6组织聚合） or sample（18样本）
AGG_MODE = "sample"          # "tissue" / "sample"
AGG_METHOD = "median"        # "median" / "mean"（仅 AGG_MODE="tissue" 时生效）

# 4) 组织列顺序（仅 tissue 模式生效；sample 模式也用于 grouped 拼列）
COLUMN_ORDER = ["foot", "gill", "intest", "lips", "mantle", "siphon"]

# 5) sample 模式列顺序
# - "samples": 按 samples.tsv 出现顺序
# - "grouped": 按 COLUMN_ORDER 组织分块拼接（每块内部仍按 samples.tsv 顺序）
SAMPLE_ORDER_MODE = "grouped"  # "samples" / "grouped"

# 6) 行（基因）顺序
# - "top_file": 完全按 top_gene_heatmap.tsv 出现顺序（最推荐，等于你手工排的故事顺序）
# - "theme": 按 Theme 分块（THEME_ORDER）+ 块内保持 top 文件顺序
ROW_ORDER_MODE = "theme"
THEME_ORDER = [
    "Development/TF",
    "Muscle/Ca2+",
    "Transport/Osmoregulation",
    "Glycosylation/Mucus",
    "Neural/Signaling",
    "ECM/Adhesion/Remodeling",
]

# 7) y轴基因标识使用哪一列（新增开关）
# - "GeneID"     ：用 top_gene_heatmap.tsv 第一列 GeneID（如 Sco02g03420.1）
# - "GeneSymbol" ：用第二列 GeneSymbol（默认；会走大小写规范化与重复消歧）
Y_LABEL_SOURCE = "GeneID"   # "GeneID" / "GeneSymbol"

# GeneSymbol 大小写统一开关（仅当 Y_LABEL_SOURCE="GeneSymbol" 时生效）
# - "keep": 原样
# - "upper": 尽量全大写（-like 后缀保持小写）
# - "title": 首字母大写（全大写符号保持不变；-like 保持小写）
# - "smart": 推荐：尽量不破坏 CHRNA10/SLC16A1 这类符号，同时修正明显的小写符号（如 pax3-b/coe3）
SYMBOL_CASE_MODE = "smart"

# 若规范化后出现重复 GeneSymbol，是否自动加 #1/#2 消歧（只在重复时触发）
DISAMBIGUATE_DUP_SYMBOLS = True

# 8) 表达量变换
DO_LOG2 = True               # True: log2(TPM+1)
DO_ROW_ZSCORE = True         # True: 每行 z-score（突出组织特异性）

# z-score 是否截断（防止极少数极端值把色阶拉爆）
CLIP_ZSCORE = True
Z_CLIP = 3.0                 # 常用 2~4；3.0 很稳

# 9) 配色方案开关
# - "scheme1": 梦幻清新发散（薄荷绿 ↔ 白 ↔ 粉紫系；适合 z-score）
# - "scheme2": 薄荷绿 → 蓝 递进（更“清新脱俗”）
# - "scheme3": Mint(#7FD3C8) → White → Peach(#F79D93)（皇上截图指定）
PALETTE_SCHEME = "scheme3"   # "scheme1" / "scheme2" / "scheme3"
PALETTE_REVERSE = False      # True: 反转颜色顺序（低↔高对调）

# 10) 作图尺寸
SHOW_GENE_LABELS = True
FIG_WIDTH = 7.5
ROW_HEIGHT = 0.18
FONT_SIZE_Y = 9
FONT_SIZE_X = 11

# 11) 可选：导出用于作图的矩阵（便于复现/附录）
EXPORT_MATRIX_TSV = True
EXPORT_MATRIX_SUFFIX = "matrix_used_for_plot.tsv"

############################################
# ========== 不要改下面内容 ==========
############################################

TOP_REQUIRED_COLS = ["GeneID", "GeneSymbol", "Theme", "PutativeFunction"]

def _ensure_dir_for_file(path: str) -> None:
    d = os.path.dirname(path)
    if d and (not os.path.exists(d)):
        os.makedirs(d, exist_ok=True)

def _safe_log2(x: np.ndarray) -> np.ndarray:
    return np.log2(x + 1.0)

def _row_zscore(mat: np.ndarray) -> np.ndarray:
    mean = mat.mean(axis=1, keepdims=True)
    std = mat.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (mat - mean) / std

def _clip(mat: np.ndarray, v: float) -> np.ndarray:
    return np.clip(mat, -v, v)

def _make_palette(scheme: str, reverse: bool) -> LinearSegmentedColormap:
    scheme = (scheme or "scheme1").strip().lower()

    # 方案1（梦幻清新发散：薄荷绿 ↔ 白 ↔ 粉紫系）
    # 说明：比旧版 scheme1 更柔和、更“奶”，但仍保留发散对比，适合 row z-score。
    scheme1 = [
        "#7FD3C8",  # mint
        "#9FD5CB",  # softer mint
        "#D7F2EC",  # very light mint
        "#FFFFFF",  # white
        "#FDE0EF",  # very light pink
        "#F7B6D2",  # soft pink
        "#F6C6E7",  # pink-lilac
        "#9E90E6",  # soft purple
    ]

    # 方案2（薄荷绿 → 蓝 递进）：来自旧图片 hex（保留）
    scheme2 = [
        "#B0F2BCFF",
        "#8EEEAEFF",
        "#70DEA7FF",
        "#57D0A3FF",
        "#43BEA3FF",
        "#34AAA2FF",
        "#2A949EFF",
        "#257D98FF",
    ]

    # 方案3（Mint(负) - White(0) - Peach(正)）：皇上截图指定
    # 注意：这是发散三锚点，适合 z-score（负/0/正）
    scheme3 = [
        "#7FD3C8",  # Mint (low / negative)
        "#FFFFFF",  # White (mid / zero)
        "#F79D93",  # Peach (high / positive)
    ]

    if scheme == "scheme1":
        colors = scheme1
    elif scheme == "scheme2":
        colors = scheme2
    elif scheme == "scheme3":
        colors = scheme3
    else:
        raise ValueError("PALETTE_SCHEME 只能是 'scheme1'/'scheme2'/'scheme3'")

    if reverse:
        colors = list(reversed(colors))

    return LinearSegmentedColormap.from_list(f"heatmap_{scheme}", colors, N=256)

def _looks_allcaps_symbol(s: str) -> bool:
    if not s:
        return False
    allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-./")
    return all((c in allowed) for c in s)

def _normalize_symbol(sym: str, mode: str) -> str:
    if sym is None:
        return ""
    sym = str(sym).strip()
    if sym == "":
        return ""

    mode = (mode or "keep").strip().lower()

    like_suffix = None
    if sym.lower().endswith("-like"):
        like_suffix = "-like"
        base = sym[:-5]
    else:
        base = sym

    if mode == "keep":
        out = base
    elif mode == "upper":
        out = base.upper()
    elif mode == "title":
        if _looks_allcaps_symbol(base):
            out = base
        else:
            out_chars = list(base)
            for i, ch in enumerate(out_chars):
                if ch.isalpha():
                    out_chars[i] = ch.upper()
                    break
            out = "".join(out_chars)
    elif mode == "smart":
        if _looks_allcaps_symbol(base):
            out = base
        else:
            has_upper = any(ch.isupper() for ch in base)
            has_lower = any(ch.islower() for ch in base)
            if (not has_upper) and has_lower:
                out_chars = list(base)
                for i, ch in enumerate(out_chars):
                    if ch.isalpha():
                        out_chars[i] = ch.upper()
                        break
                out = "".join(out_chars)
            else:
                out = base
    else:
        raise ValueError("SYMBOL_CASE_MODE 不合法：只能 keep/upper/title/smart")

    if like_suffix is not None:
        out = out + like_suffix
    return out

def _disambiguate_labels(labels):
    counts = {}
    for lab in labels:
        counts[lab] = counts.get(lab, 0) + 1
    dup_set = {k for k, v in counts.items() if v > 1}
    if not dup_set:
        return labels

    seen = {}
    out2 = []
    for lab in labels:
        if lab in dup_set:
            seen[lab] = seen.get(lab, 0) + 1
            out2.append(f"{lab} #{seen[lab]}")
        else:
            out2.append(lab)
    return out2

def _read_samples(samples_path: str) -> pd.DataFrame:
    smp = pd.read_csv(samples_path, sep=r"\s+", engine="python", dtype=str)
    for col in ["sample", "group"]:
        if col not in smp.columns:
            raise ValueError(f"samples.tsv 缺少列：{col}；当前列：{list(smp.columns)}")
    smp["sample"] = smp["sample"].astype(str)
    smp["group"] = smp["group"].astype(str)
    return smp

def _read_top(top_path: str) -> pd.DataFrame:
    top = pd.read_csv(top_path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in TOP_REQUIRED_COLS if c not in top.columns]
    if missing:
        raise ValueError(f"top_gene_heatmap.tsv 缺少必要列：{missing}；当前列：{list(top.columns)}")
    return top[TOP_REQUIRED_COLS].copy()

def _read_gene_tpm(gene_tpm_path: str) -> pd.DataFrame:
    g = pd.read_csv(gene_tpm_path, sep="\t", dtype=str, keep_default_na=False)
    if "gene_id" not in g.columns:
        raise ValueError(f"gene_tpm.tsv 缺少 gene_id 列；当前表头前几列：{list(g.columns[:10])}")
    return g

def main():
    try:
        smp = _read_samples(SAMPLES_TSV)
    except Exception as e:
        print(f"[ERROR] 无法读取 samples.tsv: {SAMPLES_TSV}\n{e}", file=sys.stderr)
        sys.exit(1)

    try:
        top = _read_top(TOP_TSV)
    except Exception as e:
        print(f"[ERROR] 无法读取 top_gene_heatmap.tsv: {TOP_TSV}\n{e}", file=sys.stderr)
        sys.exit(1)

    try:
        gtp = _read_gene_tpm(GENE_TPM_TSV)
    except Exception as e:
        print(f"[ERROR] 无法读取 gene_tpm.tsv: {GENE_TPM_TSV}\n{e}", file=sys.stderr)
        sys.exit(1)

    sample_to_group = dict(zip(smp["sample"], smp["group"]))

    # SRR 列：以 samples.tsv 的 sample 为准
    srr_cols = [s for s in smp["sample"].tolist() if s in gtp.columns]
    if len(srr_cols) == 0:
        print("[ERROR] gene_tpm.tsv 中找不到任何 SRR 列（与 samples.tsv 的 sample 不匹配）。", file=sys.stderr)
        print(f"[HINT] samples.tsv 前几个 sample：{smp['sample'].tolist()[:5]}", file=sys.stderr)
        print(f"[HINT] gene_tpm.tsv 表头前几个：{list(gtp.columns[:10])}", file=sys.stderr)
        sys.exit(1)

    # 用 top 基因列表过滤 gene_tpm
    top_gene_ids = top["GeneID"].astype(str).tolist()
    gtp["gene_id"] = gtp["gene_id"].astype(str)

    gtp_index = set(gtp["gene_id"].tolist())
    present = [gid for gid in top_gene_ids if gid in gtp_index]
    missing = [gid for gid in top_gene_ids if gid not in gtp_index]

    if len(present) == 0:
        print("[ERROR] top_gene_heatmap.tsv 的 GeneID 在 gene_tpm.tsv 里一个都找不到。", file=sys.stderr)
        sys.exit(1)

    if len(missing) > 0:
        print(f"[WARN] 有 {len(missing)} 个 GeneID 在 gene_tpm.tsv 中找不到，将被跳过：", file=sys.stderr)
        print("       " + ", ".join(missing[:30]) + (" ..." if len(missing) > 30 else ""), file=sys.stderr)

    sub = gtp.loc[gtp["gene_id"].isin(present), ["gene_id"] + srr_cols].copy()

    # 按 top 文件顺序稳定排序
    order_map = {gid: i for i, gid in enumerate(top_gene_ids)}
    sub["_ord"] = sub["gene_id"].map(order_map)
    sub = sub.sort_values("_ord", kind="stable").drop(columns=["_ord"])

    # merge 注释四列
    merged = pd.merge(
        top,
        sub,
        left_on="GeneID",
        right_on="gene_id",
        how="inner",
        sort=False,
    ).drop(columns=["gene_id"])

    # 行顺序：top_file / theme
    row_mode = ROW_ORDER_MODE.strip().lower()
    if row_mode == "theme":
        cat = pd.Categorical(merged["Theme"], categories=THEME_ORDER, ordered=True)
        merged = merged.assign(_theme_cat=cat)
        merged["_theme_cat"] = merged["_theme_cat"].cat.add_categories(["__OTHER__"]).fillna("__OTHER__")
        merged = merged.sort_values(by=["_theme_cat"], kind="stable").drop(columns=["_theme_cat"])

    # ===== y轴标签：GeneID or GeneSymbol（新增）=====
    y_source = (Y_LABEL_SOURCE or "GeneSymbol").strip().lower()
    if y_source == "geneid":
        y_labels = merged["GeneID"].astype(str).tolist()
    elif y_source == "genesymbol":
        symbols_raw = merged["GeneSymbol"].tolist()
        symbols_norm = [_normalize_symbol(s, SYMBOL_CASE_MODE) for s in symbols_raw]

        # 如果某些 symbol 为空，稳妥起见用 GeneID 顶上（避免出现空标签）
        gene_ids_fallback = merged["GeneID"].astype(str).tolist()
        symbols_norm = [sym if (sym is not None and str(sym).strip() != "") else gene_ids_fallback[i]
                       for i, sym in enumerate(symbols_norm)]

        if DISAMBIGUATE_DUP_SYMBOLS:
            symbols_norm = _disambiguate_labels(symbols_norm)
        y_labels = symbols_norm
    else:
        print("[ERROR] Y_LABEL_SOURCE 只能是 'GeneID' 或 'GeneSymbol'", file=sys.stderr)
        sys.exit(1)

    # 表达矩阵
    expr = merged[srr_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # 决定列（tissue 或 sample）
    agg_mode = AGG_MODE.strip().lower()
    if agg_mode not in {"tissue", "sample"}:
        print("[ERROR] AGG_MODE 必须是 'tissue' 或 'sample'", file=sys.stderr)
        sys.exit(1)

    if agg_mode == "tissue":
        groups_present = set(sample_to_group.values())
        for t in COLUMN_ORDER:
            if t not in groups_present:
                print(f"[ERROR] COLUMN_ORDER 中的组织 '{t}' 不在 samples.tsv 的 group 中。", file=sys.stderr)
                print(f"[HINT] samples.tsv 中出现的 group 有：{sorted(list(groups_present))}", file=sys.stderr)
                sys.exit(1)

        tissue_cols = []
        tissue_labels = []
        for tissue in COLUMN_ORDER:
            cols = [s for s in srr_cols if sample_to_group.get(s, "") == tissue]
            if len(cols) == 0:
                print(f"[ERROR] 组织 '{tissue}' 在 gene_tpm.tsv 里没有任何 SRR 列。", file=sys.stderr)
                sys.exit(1)

            sub_mat = expr[cols].to_numpy(dtype=float)
            if AGG_METHOD.strip().lower() == "median":
                vec = np.median(sub_mat, axis=1)
            elif AGG_METHOD.strip().lower() == "mean":
                vec = np.mean(sub_mat, axis=1)
            else:
                print("[ERROR] AGG_METHOD 必须是 median 或 mean", file=sys.stderr)
                sys.exit(1)

            tissue_cols.append(vec)
            tissue_labels.append(tissue)

        mat = np.stack(tissue_cols, axis=1)
        col_labels = tissue_labels

    else:
        # sample 模式：18列，横轴显示组织名 + 编号（foot-1/2/3）
        if SAMPLE_ORDER_MODE.strip().lower() == "grouped":
            ordered_srr = []
            for tissue in COLUMN_ORDER:
                for s in srr_cols:
                    if sample_to_group.get(s, "") == tissue:
                        ordered_srr.append(s)
            for s in srr_cols:
                if s not in ordered_srr:
                    ordered_srr.append(s)
        else:
            ordered_srr = srr_cols[:]

        mat = expr[ordered_srr].to_numpy(dtype=float)

        # 横轴标签：tissue-1/tissue-2...（编号按出现顺序）
        rep_count = {}
        col_labels = []
        for s in ordered_srr:
            tissue = sample_to_group.get(s, "NA")
            rep_count[tissue] = rep_count.get(tissue, 0) + 1
            col_labels.append(f"{tissue}-{rep_count[tissue]}")

    # 变换
    if DO_LOG2:
        mat = _safe_log2(mat)

    if DO_ROW_ZSCORE:
        mat_plot = _row_zscore(mat)
        if CLIP_ZSCORE:
            mat_plot = _clip(mat_plot, float(Z_CLIP))
        cbar_label = "Row Z-score"
        value_label = "log2(TPM+1) then row z-score" if DO_LOG2 else "row z-score"
    else:
        mat_plot = mat
        cbar_label = "Expression"
        value_label = "log2(TPM+1)" if DO_LOG2 else "TPM"

    # 配色
    cmap = _make_palette(PALETTE_SCHEME, PALETTE_REVERSE)

    # z-score 对称色阶
    vmin = None
    vmax = None
    if DO_ROW_ZSCORE:
        v = float(np.nanmax(np.abs(mat_plot)))
        if v == 0:
            v = 1.0
        vmin, vmax = -v, v

    # 可选导出矩阵
    if EXPORT_MATRIX_TSV:
        _ensure_dir_for_file(OUT_PREFIX + ".png")
        out_mat_path = os.path.join(
            os.path.dirname(OUT_PREFIX),
            os.path.basename(OUT_PREFIX) + "." + EXPORT_MATRIX_SUFFIX
        )
        out_df = pd.DataFrame(mat_plot, columns=col_labels)

        # 导出时仍保留两列标识，方便复现/追溯
        out_df.insert(0, "GeneSymbol", merged["GeneSymbol"].astype(str).tolist())
        out_df.insert(0, "GeneID", merged["GeneID"].astype(str).tolist())
        out_df.insert(2, "Theme", merged["Theme"].tolist())
        out_df.insert(3, "PutativeFunction", merged["PutativeFunction"].tolist())

        # 额外记录：这次画图实际用的 y 轴标签是什么
        out_df.insert(0, "YLabelUsed", y_labels)

        out_df.to_csv(out_mat_path, sep="\t", index=False)

    # 画图
    n_rows = mat_plot.shape[0]
    n_cols = mat_plot.shape[1]
    fig_h = max(3.0, ROW_HEIGHT * n_rows)
    fig_w = float(FIG_WIDTH)

    _ensure_dir_for_file(OUT_PREFIX + ".png")

    plt.figure(figsize=(fig_w, fig_h))
    ax = plt.gca()

    im = ax.imshow(
        mat_plot,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    mode_tag = f"{AGG_MODE}"
    if AGG_MODE.strip().lower() == "tissue":
        mode_tag += f"({AGG_METHOD})"

    title = f"Heatmap ({mode_tag}) | {value_label} | {PALETTE_SCHEME}{'_rev' if PALETTE_REVERSE else ''} | y={Y_LABEL_SOURCE}"
    ax.set_title(title, fontsize=12)

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=FONT_SIZE_X)

    if SHOW_GENE_LABELS:
        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(y_labels, fontsize=FONT_SIZE_Y)
    else:
        ax.set_yticks([])
        ax.set_yticklabels([])

    # ===== 去掉横纵轴黑框（隐藏四条 spines）=====
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(False)

    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(cbar_label)
    cbar.outline.set_visible(False)

    plt.tight_layout()

    out_png = f"{OUT_PREFIX}.{AGG_MODE}"
    if AGG_MODE.strip().lower() == "tissue":
        out_png += f".{AGG_METHOD}"
    out_png += f".{PALETTE_SCHEME}{'_rev' if PALETTE_REVERSE else ''}.png"

    out_pdf = out_png.replace(".png", ".pdf")

    plt.savefig(out_png, dpi=DPI)
    plt.savefig(out_pdf)
    plt.close()

    print("[OK] Heatmap saved:")
    print(out_png)
    print(out_pdf)
    if EXPORT_MATRIX_TSV:
        print("[OK] Matrix exported (for reproducibility):")
        print(out_mat_path)

    if len(missing) > 0:
        print("[WARN] Missing genes were skipped (not found in gene_tpm.tsv).", file=sys.stderr)

if __name__ == "__main__":
    main()

