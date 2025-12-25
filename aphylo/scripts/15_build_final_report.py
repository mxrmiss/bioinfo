#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
15_build_final_report.py —— 最终报告打包（加强版）

目标：把 APhylo 流水线中若干关键产物拷贝到
  paths.joint_dir/report/ 下，形成较完整的“可交付成品结构”。

主要包含：
  1) PSG (07_codeml_aggregate) 聚合结果：
       - D_fdr_genes.tsv
       - D_beb_sites.tsv
       - （若存在）D_fg_summary.tsv, D_fg_gene_counts.tsv
  2) CAFE5 (13_cafe5_aggregate, 13b_cafe_build_enrich_inputs) 结果：
       - cafe_significant_families.tsv
       - cafe_significant_families_no_highfail.tsv
       - cafe_branch_summary.tsv
       - cafe_family_universe.tsv
       - cafe_clade_summary.tsv
       - cafe_family_branch_status.tsv
       - inputs_used.tsv
  3) PSG×CAFE 联合整合 (14_joint_integration)：
       - integration_counts.tsv
       - integration_intersect.tsv
       - integration_union.tsv
  4) CAFE5 富集桥（若启用，来自 cafe5.enrich_bridge.outputs_dir）：
       - 各 tag 子目录下的 *.list / *.tsv 文件
       - 复制到 report/gene_sets/<tag>/ 目录
  5) MCMCTree 质控与 finetune 建议（与 mcmctree_postcheck/ess_report/finetune_suggest 对口）：
       - 从 mcmctree.work_dir / qc_dirname 目录读取：
           node_ages.tsv
           ess.tsv
           ess_summary.md
           ess_recommendation.txt
           finetune_suggestion.md
           finetune_new_line.txt
           summary.md
       - 复制到 report/qc/

只复制“存在的文件”；不做任何计算。
"""

from __future__ import annotations

# ============================ 可配置区 ============================
CONFIG_PATH: str = "config.yaml"
LOG_LEVEL: str = "INFO"
LOG_FILE_BASENAME: str = "15_build_final_report.log"
# ================================================================

import sys
import shutil
import logging
import traceback
from pathlib import Path
from datetime import datetime

import yaml


# -------------------- 基础工具 --------------------
def mkdir_p(p: Path) -> None:
    """类似 mkdir -p，目录已存在时不报错。"""
    p.mkdir(parents=True, exist_ok=True)


def load_yaml(p: Path) -> dict:
    """读取 YAML 配置。"""
    with p.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def copy_if_exists(src: Path, dst: Path, copied: list[str]) -> None:
    """
    若源文件存在，则复制到目标位置，并记录到 copied 列表中。

    参数：
      src    : 源文件路径
      dst    : 目标文件路径（可在子目录下）
      copied : 用于记录“已复制文件名/相对路径”的列表
    """
    if src.exists() and src.is_file():
        mkdir_p(dst.parent)
        shutil.copy2(src, dst)
        # README 中展示相对 report 目录的路径不方便，这里沿用原行为：只记录文件名
        copied.append(dst.name)
        logging.info(f"[OK] 复制：{src} -> {dst}")


# -------------------- 主流程 --------------------
def main() -> None:
    # 1) 配置与基础路径
    cfg_path = Path(CONFIG_PATH).resolve()
    if not cfg_path.exists():
        print(f"[ERR] 配置不存在：{cfg_path}", file=sys.stderr)
        sys.exit(2)

    cfg = load_yaml(cfg_path)
    paths = cfg.get("paths", {}) or {}

    logs_dir = Path(paths.get("logs_dir", "logs")).resolve()
    joint_dir = Path(paths.get("joint_dir", "results/08_joint")).resolve()
    codeml_agg_dir = Path(paths.get("codeml_agg_dir", "results/05_cmlagg")).resolve()
    cafe_agg_dir = Path(paths.get("cafe_agg_dir", "results/07_cafeagg")).resolve()

    report_dir = joint_dir / "report"
    tables_dir = report_dir / "tables"
    figs_dir = report_dir / "figs"
    qc_dir = report_dir / "qc"
    gene_sets_dir = report_dir / "gene_sets"

    mkdir_p(logs_dir)
    mkdir_p(report_dir)
    mkdir_p(tables_dir)
    mkdir_p(figs_dir)

    # 2) 日志初始化
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(logs_dir / LOG_FILE_BASENAME, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("========== APhylo 15 — 最终报告打包（加强版） ==========")
    logging.info(f"[init] 使用配置文件：{cfg_path}")
    logging.info(f"[init] joint_dir       = {joint_dir}")
    logging.info(f"[init] codeml_agg_dir  = {codeml_agg_dir}")
    logging.info(f"[init] cafe_agg_dir    = {cafe_agg_dir}")
    logging.info(f"[init] report_dir      = {report_dir}")

    # 用于 README 记录所有“实际复制”的文件名
    copied: list[str] = []

    # 3) PSG (07) 结果 —— codeml 聚合
    # 原有两个：
    copy_if_exists(codeml_agg_dir / "D_fdr_genes.tsv", tables_dir / "D_fdr_genes.tsv", copied)
    copy_if_exists(codeml_agg_dir / "D_beb_sites.tsv", tables_dir / "D_beb_sites.tsv", copied)
    # 预留一些可能存在的额外表：
    copy_if_exists(codeml_agg_dir / "D_fg_summary.tsv", tables_dir / "D_fg_summary.tsv", copied)
    copy_if_exists(codeml_agg_dir / "D_fg_gene_counts.tsv", tables_dir / "D_fg_gene_counts.tsv", copied)

    # 4) CAFE5 (13/13b) 结果
    # 13_cafe5_aggregate.py 产物
    copy_if_exists(cafe_agg_dir / "cafe_significant_families.tsv",
                   tables_dir / "cafe_significant_families.tsv", copied)
    copy_if_exists(cafe_agg_dir / "cafe_significant_families_no_highfail.tsv",
                   tables_dir / "cafe_significant_families_no_highfail.tsv", copied)
    copy_if_exists(cafe_agg_dir / "cafe_branch_summary.tsv",
                   tables_dir / "cafe_branch_summary.tsv", copied)
    copy_if_exists(cafe_agg_dir / "cafe_family_universe.tsv",
                   tables_dir / "cafe_family_universe.tsv", copied)
    copy_if_exists(cafe_agg_dir / "cafe_clade_summary.tsv",
                   tables_dir / "cafe_clade_summary.tsv", copied)
    copy_if_exists(cafe_agg_dir / "inputs_used.tsv",
                   tables_dir / "inputs_used.tsv", copied)

    # 13b_cafe_build_enrich_inputs.py 产物
    copy_if_exists(cafe_agg_dir / "cafe_family_branch_status.tsv",
                   tables_dir / "cafe_family_branch_status.tsv", copied)

    # 5) Integration (14) 结果
    copy_if_exists(joint_dir / "integration_counts.tsv",
                   tables_dir / "integration_counts.tsv", copied)
    copy_if_exists(joint_dir / "integration_intersect.tsv",
                   tables_dir / "integration_intersect.tsv", copied)
    copy_if_exists(joint_dir / "integration_union.tsv",
                   tables_dir / "integration_union.tsv", copied)

    # 6) CAFE5 富集桥（per-species gene sets，来自 13b，可选）
    cafe_cfg = cfg.get("cafe5", {}) or {}
    enrich_cfg = cafe_cfg.get("enrich_bridge", {}) or {}
    enrich_out_root = enrich_cfg.get("outputs_dir", "results/07_cafeagg/enrich_inputs")
    enrich_out_root = Path(enrich_out_root).resolve()

    if enrich_cfg.get("enable", False) and enrich_out_root.exists():
        logging.info(f"[enrich] 发现 CAFE 富集桥输出目录：{enrich_out_root}")
        mkdir_p(gene_sets_dir)

        for tag_dir in enrich_out_root.iterdir():
            if not tag_dir.is_dir():
                continue
            tag_name = tag_dir.name
            dst_tag_dir = gene_sets_dir / tag_name
            for f in tag_dir.iterdir():
                if not f.is_file():
                    continue
                # 只拷 .list / .tsv 文件，其余（如中间日志）忽略
                if not (f.suffix in (".list", ".tsv")):
                    continue
                dst_f = dst_tag_dir / f.name
                copy_if_exists(f, dst_f, copied)
        logging.info("[enrich] CAFE per-species gene sets 已复制（若存在）")
    else:
        logging.info("[enrich] 未启用或未找到 CAFE 富集桥产物，跳过 gene_sets/")

    # 7) MCMCTree 质控产物（与 mcmctree_postcheck 对口）
    mcmctree_cfg = cfg.get("mcmctree", {}) or {}
    # work_dir 优先 mcmctree.work_dir；若未配置，则回退到 paths.mcmctree_work_dir 或默认路径
    mcmc_work_dir = mcmctree_cfg.get("work_dir") or paths.get("mcmctree_work_dir") or "results/06_cafe/mcmctree"
    mcmc_work_dir = Path(mcmc_work_dir).resolve()
    qc_dirname = mcmctree_cfg.get("qc_dirname", "qc_report")
    mcmc_qc_src = mcmc_work_dir / qc_dirname

    if mcmc_qc_src.exists():
        logging.info(f"[mcmctree] 发现 MCMCTree QC 目录：{mcmc_qc_src}")
        mkdir_p(qc_dir)
        qc_files = [
            "node_ages.tsv",
            "ess.tsv",
            "ess_summary.md",
            "ess_recommendation.txt",
            "finetune_suggestion.md",
            "finetune_new_line.txt",
            "summary.md",
        ]
        for name in qc_files:
            copy_if_exists(mcmc_qc_src / name, qc_dir / name, copied)
    else:
        logging.info(f"[mcmctree] 未找到 MCMCTree QC 目录：{mcmc_qc_src}，跳过 qc/")

    # 8) 写 README
    readme_lines = [
        "# APhylo Report",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        "- Included files (by name):",
    ] + [f"  - {name}" for name in sorted(set(copied))] + [
        "",
        "Directories:",
        "  - tables/: key result tables from PSG (07), CAFE5 (13/13b), and Integration (14).",
        "  - figs/: placeholder for future figures.",
        "  - gene_sets/: per-species gene sets from CAFE5 enrichment bridge (if enabled).",
        "  - qc/: MCMCTree QC reports (node ages, ESS, finetune), if available.",
        "",
        "Notes:",
        "  - Only existing files were copied.",
        "  - Paths are defined in config.yaml under `paths:` and `mcmctree:` / `cafe5:`.",
    ]
    (report_dir / "README.txt").write_text("\n".join(readme_lines) + "\n", encoding="utf-8")
    logging.info(f"[OK] 写出：{report_dir / 'README.txt'}")

    # 9) 哨兵文件
    (joint_dir / ".report.done").write_text("ok\n", encoding="utf-8")
    logging.info("[OK] .report.done 写入完成")
    logging.info("========== APhylo 15 — 完成 ==========")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)