#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_cafe5_aggregate.py —— CAFE5 聚合（对口版｜primary: Gamma 优先，缺则 Base 回退）

与皇上当前 12 的产物强绑定：
  - models/<model>/primary_global/：
      * 优先使用 Gamma_family_results.txt
      * 若缺失则回退 Base_family_results.txt（日志明确提示）
      * 表头严格为：#FamilyID   pvalue   ...
  - models/<model>/large/：
      * 仅使用 Base_family_results.txt（若存在）
  - models/<model>/flags/high_fail_ogs.list （可选）
  - *_results.txt 仅用于分支 λ 摘要：
      * primary：Gamma_results.txt 优先，若无则回退 Base_results.txt
      * large：  仅 Base_results.txt

输出（全部写到 config.paths.cafe_agg_dir）：
  1) cafe_significant_families.tsv
  2) cafe_significant_families_no_highfail.tsv
  3) cafe_branch_summary.tsv
  4) inputs_used.tsv

重要口径：
  - 列名与判定规则与先前一致：严格要求 family结果表首列 #FamilyID/FamilyID、第二列 pvalue。
  - subset ∈ {primary, large}；source ∈ {primary:Gamma, primary:Base, large:Base}。
  - FDR 为全局一次 BH（α = config.report.fdr_alpha，默认 0.05）。
  - 本脚本不做历史/宽松列名兼容，不猜路径，全部取自 config.yaml。
"""

from __future__ import annotations

# ============================ 可配置区（皇上只需改这里） ============================
CONFIG_PATH: str = "config.yaml"        # 配置文件路径
LOG_LEVEL: str = "INFO"                 # 日志等级：DEBUG / INFO / WARNING
LOG_FILE_BASENAME: str = "13_cafe5_aggregate.log"
# ====================================================================================

import re
import sys
import math
import yaml
import logging
import traceback
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime

# --------------------------- 基础工具 ---------------------------

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def read_lines(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\n\r") for ln in f]

def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

def file_stat(p: Path) -> tuple[int, str]:
    try:
        st = p.stat()
        return st.st_size, datetime.fromtimestamp(st.st_mtime).isoformat(sep=" ", timespec="seconds")
    except Exception:
        return 0, ""

def human_bytes(n: int) -> str:
    units = ["B","KB","MB","GB","TB","PB","EB"]
    x = float(n)
    for u in units:
        if x < 1024:
            return f"{x:.1f}{u}"
        x /= 1024
    return f"{x:.1f}ZB"

def list_models(models_dir: Path) -> List[str]:
    if not models_dir.exists():
        return []
    return sorted([p.name for p in models_dir.iterdir() if p.is_dir()])

# --------------------------- BH-FDR ---------------------------

def bh_fdr(pvals: List[float]) -> List[float]:
    """Benjamini–Hochberg：返回与输入顺序一一对应的 q 值列表。"""
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q_sorted = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(order, start=1):
        q = pvals[i] * m / rank
        if q > 1.0: q = 1.0
        if q < 0.0: q = 0.0
        prev = min(prev, q)
        q_sorted[i] = prev
    return q_sorted

# --------------------------- 读取 family 结果（严格口径） ---------------------------

def read_family_pairs_strict(path: Path) -> List[Tuple[str, float]]:
    """
    严格读取 family 结果：
      - 首行必须包含 '#FamilyID' 或 'FamilyID' 作为第1列，'pvalue' 作为第2列
      - 数据行第1列为家族ID，第2列为p值；分隔符可为任意空白（制表符/空格）
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"空文件：{path}")
    header = re.split(r"\s+", lines[0].strip())
    if len(header) < 2:
        raise ValueError(f"表头列数不足（期望>=2）：{path}\n实际：{lines[0]}")
    if header[0] not in ("#FamilyID", "FamilyID"):
        raise ValueError(f"family 列名不符合严格口径（期望 '#FamilyID' 或 'FamilyID'）：{path}\n实际：{header[0]}")
    if header[1] != "pvalue":
        raise ValueError(f"p 值列名不符合严格口径（期望 'pvalue'）：{path}\n实际：{header[1]}")

    out: List[Tuple[str, float]] = []
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = re.split(r"\s+", ln.strip())
        if len(parts) < 2:
            continue
        fam = parts[0].strip()
        try:
            p = float(parts[1])
        except Exception:
            continue
        if math.isnan(p):
            continue
        if p < 0: p = 0.0
        if p > 1: p = 1.0
        out.append((fam, p))
    if not out:
        raise ValueError(f"未解析到任何 (family, pvalue)：{path}")
    return out

# --------------------------- λ 解析（results 优先，含回退） ---------------------------

def parse_lambda_from_results(path: Path) -> List[Tuple[str, float]]:
    """
    从 *results.txt 提取 (branch, lambda)：
      - 优先匹配 'branch <name> lambda = <val>'
      - 其次匹配 '<branch>\\t<lambda>' 或 '<branch> <lambda>'（最后一列为数字）
    文件不存在或未解析到则返回空列表（不报错）。
    """
    if not path.exists():
        return []
    lines = read_lines(path)
    out: List[Tuple[str, float]] = []

    pat1 = re.compile(r"branch\s+([^\s:=]+)\s+lambda\s*=\s*([0-9eE\+\-\.]+)")
    for ln in lines:
        m = pat1.search(ln)
        if m:
            b = m.group(1)
            try:
                lam = float(m.group(2))
                out.append((b, lam))
            except Exception:
                pass

    if out:
        # 去重（保留首次）
        seen = set()
        res = []
        for b, lam in out:
            if b not in seen:
                seen.add(b)
                res.append((b, lam))
        return res

    # 退一步：按列解析（末列数字）
    for ln in lines:
        if not ln.strip():
            continue
        parts = re.split(r"\s+", ln.strip())
        if len(parts) >= 2:
            try:
                lam = float(parts[-1])
                branch = " ".join(parts[:-1]).strip()
                if branch and branch.lower() not in ("family", "clade", "change"):
                    out.append((branch, lam))
            except Exception:
                pass

    # 去重
    seen = set()
    res = []
    for b, lam in out:
        if b not in seen:
            seen.add(b)
            res.append((b, lam))
    return res

def detect_error_model(runlog: Path) -> bool:
    """在 run.log 头部 '# CMD:' 行检测是否包含 '-e'；仅做注记。"""
    if not runlog.exists():
        return False
    head = read_lines(runlog)[:50]
    return any(ln.startswith("# CMD:") and "-e" in ln for ln in head)

# --------------------------- 主流程 ---------------------------

def main() -> None:
    # 1) 读取配置
    cfg_path = Path(CONFIG_PATH).resolve()
    if not cfg_path.exists():
        print(f"[ERR] 配置文件不存在：{cfg_path}", file=sys.stderr)
        sys.exit(2)
    cfg = load_yaml(cfg_path)

    # 2) 关键路径
    paths = cfg.get("paths", {})
    cafe_run_dir = Path(paths.get("cafe_run_dir", "results/06_cafe")).resolve()
    models_dir   = cafe_run_dir / "models"
    agg_dir      = Path(paths.get("cafe_agg_dir", str(cafe_run_dir / "agg"))).resolve()
    logs_dir     = Path(paths.get("logs_dir", "logs")).resolve()

    mkdir_p(logs_dir); mkdir_p(agg_dir)

    # 3) 日志
    logfile = logs_dir / LOG_FILE_BASENAME
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(logfile, encoding="utf-8"),
                  logging.StreamHandler(sys.stdout)]
    )
    logging.info("========== APhylo 13 — CAFE5 聚合（Gamma 优先，缺则 Base） ==========")
    logging.info(f"[Info] 使用配置：{cfg_path}")
    logging.info(f"[Info] cafe_run_dir = {cafe_run_dir}")
    logging.info(f"[Info] models_dir   = {models_dir}")
    logging.info(f"[Info] agg_dir      = {agg_dir}")
    logging.info(f"[Info] logs_dir     = {logs_dir}")

    # 4) 前置哨兵
    sentinel = cafe_run_dir / ".cafe.done"
    if not sentinel.exists():
        logging.error(f"[ERR] 未发现 12 步完成哨兵：{sentinel}")
        sys.exit(2)

    # 5) FDR α
    fdr_alpha = float(cfg.get("report", {}).get("fdr_alpha", 0.05))
    if not (0.0 < fdr_alpha <= 1.0):
        logging.warning(f"[WARN] report.fdr_alpha 非法：{fdr_alpha}，重置为 0.05")
        fdr_alpha = 0.05
    logging.info(f"[Info] 使用 BH-FDR α = {fdr_alpha}")

    # 6) 模型列表
    models = list_models(models_dir)
    if not models:
        logging.error(f"[ERR] 未在 {models_dir} 发现任何模型目录")
        sys.exit(2)
    logging.info(f"[Info] 发现模型：{', '.join(models)}")

    # 容器
    fam_all: List[Tuple[str, str, str, float, str]] = []  # (model, subset, family, p, source)
    inputs_used: List[List[object]] = []
    branch_rows: List[List[object]] = []

    # 7) 遍历模型
    for model in models:
        mdir = models_dir / model
        primary_dir = mdir / "primary_global"
        large_dir   = mdir / "large"
        flags_dir   = mdir / "flags"

        if not primary_dir.exists():
            logging.error(f"[ERR] 缺少目录：{primary_dir}")
            sys.exit(2)

        # ------- primary family：Gamma 优先，缺则 Base 回退 -------
        pf_gamma = primary_dir / "Gamma_family_results.txt"
        pf_base  = primary_dir / "Base_family_results.txt"

        primary_family_file = None
        primary_source = None
        if pf_gamma.exists():
            primary_family_file = pf_gamma
            primary_source = "Gamma"
        elif pf_base.exists():
            primary_family_file = pf_base
            primary_source = "Base"
            logging.info(f"[Info] {model}/primary_global 未发现 Gamma_family_results.txt，已回退 Base_family_results.txt")
        else:
            logging.error(f"[ERR] {model}/primary_global 未发现 family 结果（Gamma/Base），无法继续")
            sys.exit(2)

        # ------- large family：仅 Base（可选） -------
        large_family_file = None
        if large_dir.exists():
            lf = large_dir / "Base_family_results.txt"
            if lf.exists():
                large_family_file = lf
            else:
                logging.info(f"[Info] {model}/large 目录存在，但未发现 Base_family_results.txt（忽略 large 线）")
        else:
            logging.info(f"[Info] {model} 未启用 two-stage（无 large 目录）")

        # 高失败率列表（可选）
        hf_path = flags_dir / "high_fail_ogs.list"
        if hf_path.exists():
            sz, mt = file_stat(hf_path)
            inputs_used.append([model, "-", str(hf_path), "flags_high_fail", "N/A", sz, human_bytes(sz), mt, ""])
            logging.info(f"[Info] {model} 读取 high_fail_ogs：{len([x for x in read_lines(hf_path) if x.strip()])} 条")
        else:
            logging.info(f"[Info] {model} 未发现 high_fail_ogs.list（按空集处理）")

        # 误差模型注记（仅溯源）
        used_em = detect_error_model(primary_dir / "run.log")
        if used_em:
            logging.info(f"[Info] {model}/primary_global 检测到误差模型（-e）")

        # 读取 primary family
        try:
            pairs = read_family_pairs_strict(primary_family_file)
        except Exception as e:
            logging.error(f"[ERR] 读取 {primary_family_file} 失败：{e}")
            sys.exit(2)
        for fam, p in pairs:
            fam_all.append((model, "primary", fam, p, f"primary:{primary_source}"))
        sz, mt = file_stat(primary_family_file)
        inputs_used.append([model, "primary", str(primary_family_file), "family_results", primary_source, sz, human_bytes(sz), mt, "error_model" if used_em else ""])

        # 读取 large family（若有）
        if large_family_file is not None:
            try:
                lpairs = read_family_pairs_strict(large_family_file)
                for fam, p in lpairs:
                    fam_all.append((model, "large", fam, p, "large:Base"))
                sz, mt = file_stat(large_family_file)
                inputs_used.append([model, "large", str(large_family_file), "family_results", "Base", sz, human_bytes(sz), mt, ""])
            except Exception as e:
                logging.warning(f"[WARN] 读取 {large_family_file} 失败：{e}（忽略 large 线）")

        # ------- 分支 λ 摘要：primary 先 Gamma，缺则 Base；large 仅 Base -------
        # primary
        g_res = primary_dir / "Gamma_results.txt"
        b_res = primary_dir / "Base_results.txt"
        gl = parse_lambda_from_results(g_res)
        src_res = None
        if gl:
            src_res = "Gamma"
        else:
            bl = parse_lambda_from_results(b_res)
            if bl:
                gl = bl
                src_res = "Base"
        if gl:
            res_file = g_res if src_res == "Gamma" else b_res
            sz, mt = file_stat(res_file)
            inputs_used.append([model, "primary", str(res_file), "results", src_res, sz, human_bytes(sz), mt, ""])
            for b, lam in gl:
                branch_rows.append([model, "primary", b, lam, "results"])

        # large
        lb_res = large_dir / "Base_results.txt"
        lbl = parse_lambda_from_results(lb_res)
        if lbl:
            sz, mt = file_stat(lb_res)
            inputs_used.append([model, "large", str(lb_res), "results", "Base", sz, human_bytes(sz), mt, ""])
            for b, lam in lbl:
                branch_rows.append([model, "large", b, lam, "results"])

        # run.log 溯源
        for sub, d in (("primary", primary_dir), ("large", large_dir)):
            rlog = d / "run.log"
            if rlog.exists():
                sz, mt = file_stat(rlog)
                note = "error_model" if (sub == "primary" and used_em) else ""
                inputs_used.append([model, sub, str(rlog), "run_log", "-", sz, human_bytes(sz), mt, note])

    # 8) 全局 BH-FDR 并写表
    if not fam_all:
        logging.error("[ERR] 未收集到任何 family 记录")
        sys.exit(2)

    pvals = [x[3] for x in fam_all]
    qvals = bh_fdr(pvals)

    all_rows: List[List[object]] = []
    for (model, subset, family, p, source), q in zip(fam_all, qvals):
        sig = "yes" if q <= fdr_alpha else "no"
        all_rows.append([model, subset, family, f"{p:.6g}", f"{q:.6g}", sig, source])

    sig_path = agg_dir / "cafe_significant_families.tsv"
    write_tsv(sig_path, ["model", "subset", "family", "P", "Q", "significant", "source"], all_rows)
    logging.info(f"[OK] 写出：{sig_path}（{len(all_rows)} 行）")

    # 9) 剔除高失败率版本
    model_hf: Dict[str, set] = {}
    for model in models:
        hf = models_dir / model / "flags" / "high_fail_ogs.list"
        if hf.exists():
            model_hf[model] = set(x.strip() for x in read_lines(hf) if x.strip())
        else:
            model_hf[model] = set()
    nofail_rows: List[List[object]] = []
    for row in all_rows:
        model, subset, fam = row[0], row[1], row[2]
        if fam not in model_hf.get(model, set()):
            nofail_rows.append(row)

    sig_nofail_path = agg_dir / "cafe_significant_families_no_highfail.tsv"
    write_tsv(sig_nofail_path, ["model", "subset", "family", "P", "Q", "significant", "source"], nofail_rows)
    logging.info(f"[OK] 写出：{sig_nofail_path}（{len(nofail_rows)} 行）")

    # 10) 分支 λ 摘要
    branch_path = agg_dir / "cafe_branch_summary.tsv"
    write_tsv(branch_path, ["model", "subset", "branch", "lambda", "extra"], branch_rows)
    logging.info(f"[OK] 写出：{branch_path}（{len(branch_rows)} 行）")

    # 11) inputs 溯源清单
    inputs_path = agg_dir / "inputs_used.tsv"
    write_tsv(inputs_path,
              ["model","subset","file","role","source","bytes","size_h","mtime","note"],
              inputs_used)
    logging.info(f"[OK] 写出：{inputs_path}（{len(inputs_used)} 行）")

    logging.info("========== APhylo 13 — 完成 ==========")

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

