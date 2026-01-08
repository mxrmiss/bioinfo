#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
tissue_atlas_specific.py —— 通用版“组织富集 / 组织特异”候选基因筛选脚本

使用场景（通用设计）：
  - 输入：gene × sample 的 TPM 矩阵（例如：results/05_matrix/tpms/gene_tpm.tsv）
  - 目标：在多个组织样本中，筛选某一个“目标样本 / 目标组织”中高表达、
          而在其它样本中相对低表达的候选基因。

核心设计：
  - 不写死具体组织，例如 foot/gill，而是通过参数：
      TARGET_SAMPLE_ID  —— 目标样本 ID（表头中的一列）
      TARGET_LABEL      —— 对应的组织名称标签（只用于日志 & 文件名）
      BACKGROUND_SAMPLE_IDS —— 其他样本 ID 列表（作为背景组织）

输出文件（全部为通用命名，带上目标标签）：
  1) <TARGET_LABEL>_enriched_candidates.tsv
  2) <TARGET_LABEL>_specific_strict.tsv
  3) <TARGET_LABEL>_enriched_candidates.list
  4) <TARGET_LABEL>_specific_strict.list

注意：
  - 本脚本不依赖命令行参数，所有参数在脚本顶部配置，方便在不同项目中复用。
  - 逻辑对任意目标样本 / 组织适用，只需修改配置区即可。
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import List

import pandas as pd


# =============================================================================
# 0. 参数区 —— 皇上只需要在这里修改配置，就能在不同项目 / 不同组织中复用
# =============================================================================

# 1) 输入 TPM 矩阵路径（相对路径以脚本所在目录为基准）
#    默认假定项目结构：
#      project/
#        scripts/
#        results/05_matrix/tpms/gene_tpm.tsv
INPUT_TPM_FILE: str = "../results/05_matrix/tpms/gene_tpm.tsv"

# 2) 目标样本及背景样本设置
#    - TARGET_SAMPLE_ID：TPM 矩阵表头中目标样本的列名
#    - TARGET_LABEL    ：目标样本对应的“组织标签”，用于输出文件命名与日志展示
#    - BACKGROUND_SAMPLE_IDS：其余样本 ID 列表，用作“背景组织”
TARGET_SAMPLE_ID: str = "SRR2162887"   # 例如：foot 的样本 ID
TARGET_LABEL: str = "foot"            # 例如：foot（可改为 gill/mantle 等）
BACKGROUND_SAMPLE_IDS: List[str] = [
    "SRR2162883",  # gill
    "SRR2162892",  # adductor muscle
    "SRR2162895",  # visceral mass
    "SRR2162898",  # mantle
    "SRR2162902",  # siphon
]

# 3) “目标样本富集（宽松版）”筛选阈值
#    - 条件 1：目标样本 TPM 不低于该阈值
#    - 条件 2：目标样本 / 背景最大值 >= 该 fold（富集倍数）
TPM_MIN_ENRICHED: float = 5.0     # target TPM ≥ 5
FOLD_MIN_ENRICHED: float = 4.0    # target / max(background) ≥ 4

# 4) “目标样本高度特异（严格版）”筛选阈值
#    在满足“富集版”基础上，增加更严格的条件：
TPM_MIN_STRICT: float = 10.0          # target TPM ≥ 10
FOLD_MIN_STRICT: float = 8.0          # target / max(background) ≥ 8
TPM_MAX_BACKGROUND_STRICT: float = 2.0  # 所有背景组织 TPM ≤ 2

# 5) 输出目录（相对脚本目录；脚本会自动创建）
OUTPUT_DIR: str = "../results/tissue_atlas"

# 6) 计算 fold 时防止除以 0 的小常数
EPS: float = 1e-4


# =============================================================================
# 1. 小工具函数
# =============================================================================

def log(message: str) -> None:
    """简单日志打印函数，统一前缀，方便日后重定向输出。"""
    sys.stdout.write(f"[tissue_atlas] {message}\n")
    sys.stdout.flush()


# =============================================================================
# 2. 读取 TPM 矩阵并做基础检查
# =============================================================================

def load_tpm_matrix(tpm_path: Path) -> pd.DataFrame:
    """
    读取 TPM 矩阵：
      - 要求第一列为基因 ID（列名推荐为 'gene_id'）
      - 后面的列为各个样本的 TPM（列名为样本 ID，例如 SRR***）
      - 分隔符为制表符

    返回：
      DataFrame，索引为 gene_id，列为样本 ID
    """
    if not tpm_path.exists():
        raise FileNotFoundError(f"找不到 TPM 矩阵文件：{tpm_path}")

    log(f"读取 TPM 矩阵：{tpm_path}")
    df = pd.read_csv(tpm_path, sep="\t")

    # 尝试识别基因 ID 列
    if "gene_id" in df.columns:
        gene_col = "gene_id"
    elif "GeneID" in df.columns:
        gene_col = "GeneID"
    else:
        # 如果没有标准列名，则默认第一列为基因 ID
        gene_col = df.columns[0]
        log(f"警告：未找到 'gene_id' 列，默认使用第一列 '{gene_col}' 作为基因 ID。")

    # 将 gene_id 设为索引
    df = df.set_index(gene_col)

    # 将表达列统一转为 float，非数值全部视为 0.0
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    log(f"TPM 矩阵维度：{df.shape[0]} 个基因 × {df.shape[1]} 个样本")
    return df


# =============================================================================
# 3. 计算“目标样本富集”相关指标
# =============================================================================

def compute_target_enrichment(df: pd.DataFrame) -> pd.DataFrame:
    """
    在 TPM 矩阵上计算以下列：
      - target_tpm           ：目标样本表达量
      - background_max       ：背景样本中的最大 TPM
      - background_mean      ：背景样本的平均 TPM（供参考）
      - fold_target_vs_bgmax ：target_tpm / (background_max + EPS)
    """
    all_needed = [TARGET_SAMPLE_ID] + BACKGROUND_SAMPLE_IDS
    missing_samples = [s for s in all_needed if s not in df.columns]
    if missing_samples:
        raise ValueError(f"在 TPM 矩阵中找不到以下样本列：{missing_samples}")

    target = df[TARGET_SAMPLE_ID]
    background = df[BACKGROUND_SAMPLE_IDS]

    df_out = df.copy()
    df_out["target_tpm"] = target
    df_out["background_max"] = background.max(axis=1)
    df_out["background_mean"] = background.mean(axis=1)
    df_out["fold_target_vs_bgmax"] = (df_out["target_tpm"] + EPS) / (df_out["background_max"] + EPS)

    return df_out


# =============================================================================
# 4. 根据不同阈值筛选候选基因
# =============================================================================

def filter_enriched(df: pd.DataFrame) -> pd.DataFrame:
    """
    “富集（宽松版）”筛选：
      条件：
        - target_tpm ≥ TPM_MIN_ENRICHED
        - fold_target_vs_bgmax ≥ FOLD_MIN_ENRICHED
    """
    cond_target = df["target_tpm"] >= TPM_MIN_ENRICHED
    cond_fold = df["fold_target_vs_bgmax"] >= FOLD_MIN_ENRICHED
    mask = cond_target & cond_fold

    filtered = df.loc[mask].copy()
    log(f"目标样本富集（宽松版）候选基因数：{filtered.shape[0]}")
    return filtered


def filter_strict(df: pd.DataFrame) -> pd.DataFrame:
    """
    “高度特异（严格版）”筛选：
      基于富集候选进一步添加限制：
        - target_tpm ≥ TPM_MIN_STRICT
        - fold_target_vs_bgmax ≥ FOLD_MIN_STRICT
        - 所有背景样本 TPM ≤ TPM_MAX_BACKGROUND_STRICT
    """
    enriched = filter_enriched(df)

    cond_target_strict = enriched["target_tpm"] >= TPM_MIN_STRICT
    cond_fold_strict = enriched["fold_target_vs_bgmax"] >= FOLD_MIN_STRICT

    background = enriched[BACKGROUND_SAMPLE_IDS]
    cond_bg_low = (background <= TPM_MAX_BACKGROUND_STRICT).all(axis=1)

    mask = cond_target_strict & cond_fold_strict & cond_bg_low
    strict = enriched.loc[mask].copy()

    log(f"目标样本高度特异（严格版）候选基因数：{strict.shape[0]}")
    return strict


# =============================================================================
# 5. 输出结果到多个文件（全通用命名）
# =============================================================================

def write_outputs(df_all: pd.DataFrame,
                  enriched: pd.DataFrame,
                  strict: pd.DataFrame,
                  out_dir: Path) -> None:
    """
    写出以下结果：
      - <TARGET_LABEL>_enriched_candidates.tsv
      - <TARGET_LABEL>_specific_strict.tsv
      - <TARGET_LABEL>_enriched_candidates.list
      - <TARGET_LABEL>_specific_strict.list
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    enriched_tsv = out_dir / f"{TARGET_LABEL}_enriched_candidates.tsv"
    strict_tsv = out_dir / f"{TARGET_LABEL}_specific_strict.tsv"
    enriched_list = out_dir / f"{TARGET_LABEL}_enriched_candidates.list"
    strict_list = out_dir / f"{TARGET_LABEL}_specific_strict.list"

    # 写 TSV（保留所有样本列 + 计算列）
    log(f"写出富集候选表：{enriched_tsv}")
    enriched.to_csv(enriched_tsv, sep="\t", index=True, index_label="gene_id")

    log(f"写出高度特异候选表：{strict_tsv}")
    strict.to_csv(strict_tsv, sep="\t", index=True, index_label="gene_id")

    # 写 gene_id 列表（可直接用于富集分析 / 交集分析等）
    log(f"写出富集候选基因列表：{enriched_list}")
    with enriched_list.open("w") as f:
        for gid in enriched.index:
            f.write(f"{gid}\n")

    log(f"写出高度特异基因列表：{strict_list}")
    with strict_list.open("w") as f:
        for gid in strict.index:
            f.write(f"{gid}\n")

    # -------------------------------------------------------------------------
    # 新增功能：写出“背景基因集”列表（可检测基因宇宙）
    # 定义：在全部样本（TPM 表中的表达列）中，至少在一个样本中 TPM > 0 的基因
    # 文件名：<TARGET_LABEL>_background_detectable.list
    # -------------------------------------------------------------------------
    # 表中包含了原始 TPM 列 + 计算出的 4 个指标列，这里排除后者，仅用表达列判断“可检测”
    metric_cols = {"target_tpm", "background_max", "background_mean", "fold_target_vs_bgmax"}
    expr_cols = [c for c in df_all.columns if c not in metric_cols]

    if expr_cols:
        expr_max = df_all[expr_cols].max(axis=1)
        bg_mask = expr_max > 0.0
    else:
        # 极端容错：如果没有表达列（理论上不会发生），则全部基因都视为背景
        bg_mask = pd.Series(True, index=df_all.index)

    bg_genes = df_all.index[bg_mask]
    bg_list_path = out_dir / f"{TARGET_LABEL}_background_detectable.list"

    log(f"写出背景基因集（可检测基因）：{bg_list_path}")
    with bg_list_path.open("w") as f:
        for gid in bg_genes:
            f.write(f"{gid}\n")
    log(f"背景基因总数（可检测）：{bg_genes.shape[0]}")

    # 简单汇总信息
    log("========== 汇总信息 ==========")
    log(f"总基因数：{df_all.shape[0]}")
    log(f"富集候选基因数：{enriched.shape[0]}")
    log(f"高度特异候选基因数：{strict.shape[0]}")
    log(f"输出目录：{out_dir.resolve()}")


# =============================================================================
# 6. 主函数
# =============================================================================

def main() -> None:
    # 统一以脚本所在目录作为相对路径基准，方便在不同位置调用
    script_dir = Path(__file__).resolve().parent
    tpm_path = (script_dir / INPUT_TPM_FILE).resolve()
    out_dir = (script_dir / OUTPUT_DIR).resolve()

    log("============================================")
    log("  通用组织富集 / 特异基因筛选 —— 开始运行")
    log("============================================")
    log(f"脚本目录：{script_dir}")
    log(f"TPM 矩阵路径：{tpm_path}")
    log(f"输出目录：{out_dir}")
    log(f"目标样本 ID：{TARGET_SAMPLE_ID}")
    log(f"目标标签   ：{TARGET_LABEL}")
    log(f"背景样本 ID：{', '.join(BACKGROUND_SAMPLE_IDS)}")

    # 1) 读取 TPM 矩阵
    df_tpm = load_tpm_matrix(tpm_path)

    # 2) 计算目标样本富集指标
    df_with_metrics = compute_target_enrichment(df_tpm)

    # 3) 筛选富集 & 严格特异基因
    enriched = filter_enriched(df_with_metrics)
    strict = filter_strict(df_with_metrics)

    # 4) 写出结果
    write_outputs(df_with_metrics, enriched, strict, out_dir)

    log("============================================")
    log("  通用组织富集 / 特异基因筛选 —— 运行结束")
    log("============================================")


if __name__ == "__main__":
    main()

