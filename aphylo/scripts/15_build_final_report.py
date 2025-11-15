#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
15_build_final_report.py —— 最终报告打包（极简对口版）

目标：把 PSG(07)、CAFE(13)、Integration(14) 的关键表格复制到
paths.joint_dir/report/ 下，生成可交付的最小成品结构。

只复制“存在的文件”；不做计算。
"""

from __future__ import annotations

# ============================ 可配置区 ============================
CONFIG_PATH: str = "config.yaml"
LOG_LEVEL: str = "INFO"
LOG_FILE_BASENAME: str = "15_build_final_report.log"
# ================================================================

import sys, yaml, shutil, logging, traceback
from pathlib import Path
from datetime import datetime

# -------------------- 基础工具 --------------------
def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def load_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def copy_if_exists(src: Path, dst: Path, copied: list[str]) -> None:
    if src.exists():
        shutil.copy2(src, dst)
        copied.append(dst.name)
        logging.info(f"[OK] 复制：{src.name} -> {dst}")

# -------------------- 主流程 --------------------
def main() -> None:
    # 1) 配置与路径
    cfg_path = Path(CONFIG_PATH).resolve()
    if not cfg_path.exists():
        print(f"[ERR] 配置不存在：{cfg_path}", file=sys.stderr); sys.exit(2)
    cfg = load_yaml(cfg_path)
    paths = cfg.get("paths", {})
    logs_dir = Path(paths.get("logs_dir", "logs")).resolve()
    joint_dir = Path(paths.get("joint_dir", "results/08_joint")).resolve()
    codeml_agg_dir = Path(paths.get("codeml_agg_dir", "results/05_cmlagg")).resolve()
    cafe_agg_dir   = Path(paths.get("cafe_agg_dir",   "results/07_cafeagg")).resolve()
    report_dir = joint_dir / "report"
    tables_dir = report_dir / "tables"
    figs_dir   = report_dir / "figs"

    mkdir_p(logs_dir); mkdir_p(report_dir); mkdir_p(tables_dir); mkdir_p(figs_dir)

    # 2) 日志
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(logs_dir / LOG_FILE_BASENAME, encoding="utf-8"),
                  logging.StreamHandler(sys.stdout)]
    )
    logging.info("========== APhylo 15 — 最终报告打包 ==========")

    # 3) 复制表格（存在即拷贝）
    copied: list[str] = []

    # PSG (07)
    copy_if_exists(codeml_agg_dir / "D_fdr_genes.tsv", tables_dir / "D_fdr_genes.tsv", copied)
    copy_if_exists(codeml_agg_dir / "D_beb_sites.tsv", tables_dir / "D_beb_sites.tsv", copied)

    # CAFE (13)
    copy_if_exists(cafe_agg_dir / "cafe_significant_families.tsv",        tables_dir / "cafe_significant_families.tsv", copied)
    copy_if_exists(cafe_agg_dir / "cafe_significant_families_no_highfail.tsv", tables_dir / "cafe_significant_families_no_highfail.tsv", copied)
    copy_if_exists(cafe_agg_dir / "cafe_branch_summary.tsv",               tables_dir / "cafe_branch_summary.tsv", copied)
    copy_if_exists(cafe_agg_dir / "inputs_used.tsv",                       tables_dir / "inputs_used.tsv", copied)

    # Integration (14)
    copy_if_exists(joint_dir / "integration_counts.tsv",   tables_dir / "integration_counts.tsv", copied)
    copy_if_exists(joint_dir / "integration_intersect.tsv",tables_dir / "integration_intersect.tsv", copied)
    copy_if_exists(joint_dir / "integration_union.tsv",    tables_dir / "integration_union.tsv", copied)

    # 4) 写 README
    readme = [
        "# APhylo Report",
        f"- Generated at: {datetime.now().isoformat(timespec='seconds')}",
        "- Included tables:",
    ] + [f"  - {name}" for name in sorted(copied)] + [
        "",
        "Directories:",
        "  - tables/: key result tables from PSG (07), CAFE5 (13), and Integration (14).",
        "  - figs/: placeholder for future figures.",
        "",
        "Notes:",
        "  - Only existing files were copied.",
        "  - Paths are defined in config.yaml under `paths:`.",
    ]
    (report_dir / "README.txt").write_text("\n".join(readme) + "\n", encoding="utf-8")
    logging.info(f"[OK] 写出：{report_dir / 'README.txt'}")

    # 5) 哨兵
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
