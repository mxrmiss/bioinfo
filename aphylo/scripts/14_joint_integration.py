#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
14_joint_integration.py —— PSG×CAFE 联合整合（极简对口版）

仅接受现役产物与表头：
- PSG(07):  paths.codeml_agg_dir/D_fdr_genes.tsv    必须含列：OG, Q
- CAFE(13): paths.cafe_agg_dir/cafe_significant_families.tsv  必须含列：family, Q, significant

输出到 paths.joint_dir：
- integration_counts.tsv
- integration_intersect.tsv
- integration_union.tsv
- .integration.done

说明：
- 显著阈值 alpha = config.report.fdr_alpha（默认0.05）
- PSG 显著：Q <= alpha
- CAFE 显著：significant == "yes"（来自13，已按同一alpha判定）
"""

from __future__ import annotations

# ============================ 可配置区 ============================
CONFIG_PATH: str = "config.yaml"                # 配置文件路径
LOG_LEVEL: str = "INFO"                         # DEBUG/INFO/WARNING
LOG_FILE_BASENAME: str = "14_joint_integration.log"
# ================================================================

import sys, yaml, logging, traceback
from pathlib import Path
from datetime import datetime

# -------------------- 基础工具 --------------------
def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def load_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def read_tsv(p: Path) -> tuple[list[str], list[list[str]]]:
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise ValueError(f"空文件：{p}")
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    return header, rows

def write_tsv(p: Path, header: list[str], rows: list[list[object]]) -> None:
    with p.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

# -------------------- 主流程 --------------------
def main() -> None:
    # 1) 配置与路径
    cfg_path = Path(CONFIG_PATH).resolve()
    if not cfg_path.exists():
        print(f"[ERR] 配置不存在：{cfg_path}", file=sys.stderr); sys.exit(2)
    cfg = load_yaml(cfg_path)
    paths = cfg.get("paths", {})
    logs_dir = Path(paths.get("logs_dir", "logs")).resolve()
    codeml_agg_dir = Path(paths.get("codeml_agg_dir", "results/05_cmlagg")).resolve()
    cafe_agg_dir   = Path(paths.get("cafe_agg_dir",   "results/07_cafeagg")).resolve()
    joint_dir      = Path(paths.get("joint_dir",      "results/08_joint")).resolve()
    mkdir_p(logs_dir); mkdir_p(joint_dir)

    # 2) 日志
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(logs_dir / LOG_FILE_BASENAME, encoding="utf-8"),
                  logging.StreamHandler(sys.stdout)]
    )
    logging.info("========== APhylo 14 — PSG×CAFE 联合整合 ==========")

    # 3) 读取 alpha
    alpha = float(cfg.get("report", {}).get("fdr_alpha", 0.05))
    if not (0 < alpha <= 1):
        logging.warning(f"[WARN] report.fdr_alpha 非法：{alpha}，重置为0.05")
        alpha = 0.05
    logging.info(f"[Info] 使用 FDR α = {alpha}")

    # 4) 定位输入
    psg_tsv  = codeml_agg_dir / "D_fdr_genes.tsv"
    cafe_tsv = cafe_agg_dir   / "cafe_significant_families.tsv"
    if not psg_tsv.exists():
        logging.error(f"[ERR] 缺少 PSG 显著表：{psg_tsv}"); sys.exit(2)
    if not cafe_tsv.exists():
        logging.error(f"[ERR] 缺少 CAFE 显著表：{cafe_tsv}"); sys.exit(2)

    # 5) 读取 PSG（严格要求列：OG, Q）
    h_psg, r_psg = read_tsv(psg_tsv)
    try:
        i_og = h_psg.index("OG")
        i_q  = h_psg.index("Q")
    except ValueError as e:
        logging.error(f"[ERR] PSG表头缺列（需 OG 与 Q）：{h_psg}")
        sys.exit(2)
    psg_sig = set()
    for row in r_psg:
        if len(row) <= max(i_og, i_q): continue
        try:
            q = float(row[i_q])
        except Exception:
            continue
        if q <= alpha:
            psg_sig.add(row[i_og])

    # 6) 读取 CAFE（严格要求列：family, Q, significant）
    h_cafe, r_cafe = read_tsv(cafe_tsv)
    try:
        i_fam = h_cafe.index("family")
        i_qc  = h_cafe.index("Q")
        i_sig = h_cafe.index("significant")
    except ValueError:
        logging.error(f"[ERR] CAFE表头缺列（需 family, Q, significant）：{h_cafe}")
        sys.exit(2)
    cafe_sig = set()
    cafe_q_map: dict[str, float] = {}
    cafe_subset_map: dict[str, str] = {}
    try:
        i_subset = h_cafe.index("subset")
    except ValueError:
        i_subset = None
    for row in r_cafe:
        if len(row) <= max(i_fam, i_qc, i_sig): continue
        fam = row[i_fam]
        cafe_subset_map[fam] = (row[i_subset] if i_subset is not None and i_subset < len(row) else "")
        try:
            qc = float(row[i_qc]); cafe_q_map[fam] = qc
        except Exception:
            cafe_q_map[fam] = float("nan")
        if row[i_sig].lower() == "yes":
            cafe_sig.add(fam)

    # 7) 交/并集
    inter = sorted(psg_sig & cafe_sig)
    union = sorted(psg_sig | cafe_sig)

    # 8) 写 integration_counts.tsv
    counts_path = joint_dir / "integration_counts.tsv"
    write_tsv(counts_path,
              ["alpha","psg_sig_n","cafe_sig_n","intersect_n","union_n","timestamp"],
              [[f"{alpha:.6g}", len(psg_sig), len(cafe_sig), len(inter), len(union),
                datetime.now().isoformat(timespec="seconds")]])
    logging.info(f"[OK] 写出：{counts_path}")

    # 9) 写 integration_intersect.tsv（交集明细）
    inter_path = joint_dir / "integration_intersect.tsv"
    # PSG 侧把 Q 值也带上
    psg_q_map: dict[str, float] = {}
    for row in r_psg:
        if len(row) <= max(i_og, i_q): continue
        og = row[i_og]
        try:
            psg_q_map[og] = float(row[i_q])
        except Exception:
            pass
    inter_rows = []
    for og in inter:
        inter_rows.append([
            og,                              # OG / family
            f"{psg_q_map.get(og, float('nan')):.6g}",
            f"{cafe_q_map.get(og, float('nan')):.6g}",
            cafe_subset_map.get(og, "")
        ])
    write_tsv(inter_path, ["OG","psg_Q","cafe_Q","cafe_subset"], inter_rows)
    logging.info(f"[OK] 写出：{inter_path}（{len(inter_rows)} 行）")

    # 10) 写 integration_union.tsv（并集明细）
    union_path = joint_dir / "integration_union.tsv"
    union_rows = []
    for og in union:
        if og in psg_sig and og in cafe_sig:
            origin = "both"
        elif og in psg_sig:
            origin = "codeml_only"
        else:
            origin = "cafe_only"
        union_rows.append([og, origin, f"{psg_q_map.get(og, float('nan')):.6g}",
                           f"{cafe_q_map.get(og, float('nan')):.6g}",
                           cafe_subset_map.get(og, "")])
    write_tsv(union_path, ["OG","origin","psg_Q","cafe_Q","cafe_subset"], union_rows)
    logging.info(f"[OK] 写出：{union_path}（{len(union_rows)} 行）")

    # 11) 哨兵
    (joint_dir / ".integration.done").write_text("ok\n", encoding="utf-8")
    logging.info("[OK] .integration.done 写入完成")
    logging.info("========== APhylo 14 — 完成 ==========")

if __name__ == "__main__":
    try:
        main()
    except SystemExit:  # 让显式退出代码冒泡
        raise
    except Exception as e:
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

