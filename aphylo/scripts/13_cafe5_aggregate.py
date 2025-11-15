#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_cafe5_aggregate.py —— CAFE5 聚合（与 12 脚本严丝合缝版）

设计要点（与皇上当前目录完全对口）：
  1) 仅解析 models/<model>/primary_global 与 models/<model>/large 两套目录
  2) primary_global：优先使用 Gamma_family_results.txt；若不存在回退 Base_family_results.txt
  3) large：仅使用 Base_family_results.txt（与 12 的二阶段产物一致）
  4) 只产一条主结果线（不生成 *_with_error 旁路），但会在 inputs_used.tsv 标注是否启用 -e
  5) 高失败率家族剔除来自：models/<model>/flags/high_fail_ogs.list（若无则当空集）
  6) 全局一次性 BH-FDR（alpha 来自 config.report.fdr_alpha，默认 0.05）
  7) 分支 λ 摘要同时扫描 primary_global 与 large；优先从 *_results.txt 提取，run.log 兜底
  8) 所有聚合产物写入 config.paths.cafe_agg_dir（与后续 14/15 对口）
  9) 运行前强校验：<cafe_run_dir>/.cafe.done 必须存在；primary_global 下需有 family 结果文件

输出文件（全部在 paths.cafe_agg_dir）：
  - cafe_significant_families.tsv
  - cafe_significant_families_no_highfail.tsv
  - cafe_branch_summary.tsv
  - inputs_used.tsv

注意：
  - 所有路径与参数仅来自 config.yaml；本脚本不接受命令行参数。
  - 兼顾多 model（例如 global、XYZ），自动遍历 models/*，列中保留 model 字段。
"""

from __future__ import annotations

# ============================ 可配置区（皇上只需改这里） ============================

# 配置文件路径（默认同项目根目录）
CONFIG_PATH: str = "config.yaml"

# 日志详细程度：DEBUG / INFO / WARNING
LOG_LEVEL: str = "INFO"

# 写入日志文件名（落在 config.paths.logs_dir）
LOG_FILE_BASENAME: str = "13_cafe5_aggregate.log"

# ====================================================================================

import os
import re
import sys
import time
import math
import yaml
import errno
import logging
import traceback
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime

# 为了避免额外依赖，尽量用标准库；仅在 pandas 存在时用于稳健读表（可无）
try:
    import pandas as pd
except Exception:
    pd = None

# ----------------------------- 基础工具函数 -----------------------------

def mkdir_p(path: Path) -> None:
    """安全创建目录（等价于 shell: mkdir -p）。"""
    path.mkdir(parents=True, exist_ok=True)

def load_yaml(path: Path) -> dict:
    """读取 YAML 配置。"""
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def human_bytes(n: int) -> str:
    """将字节数格式化为便读字符串。"""
    units = ["B","KB","MB","GB","TB"]
    x = float(n)
    for u in units:
        if x < 1024.0:
            return f"{x:.1f}{u}"
        x /= 1024.0
    return f"{x:.1f}PB"

def list_models(models_dir: Path) -> List[str]:
    """列出 models 目录下的 model 名称（子目录名）。"""
    if not models_dir.exists():
        return []
    return sorted([p.name for p in models_dir.iterdir() if p.is_dir()])

def file_stat(p: Path) -> Tuple[int, str]:
    """返回文件大小与修改时间（字符串）。"""
    try:
        st = p.stat()
        return st.st_size, datetime.fromtimestamp(st.st_mtime).isoformat(sep=" ", timespec="seconds")
    except Exception:
        return 0, ""

def read_lines(path: Path) -> List[str]:
    """读取文本文件每行，去掉换行。"""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\n\r") for ln in f]

def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    """写出 TSV（不依赖 pandas，保证环境最小化）。"""
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

# --------------------------- BH-FDR（Benjamini-Hochberg） ---------------------------

def bh_fdr(pvals: List[float]) -> List[float]:
    """
    纯标准库 BH-FDR；输入 p 值列表，输出对应的 q 值列表（与输入顺序一致）。
    算法：按升序排序、计算 q_i = p_(i) * (m / i)，再做向后累积最小化。
    """
    m = len(pvals)
    if m == 0:
        return []
    # 记录原索引
    idx_sorted = sorted(range(m), key=lambda i: (float("inf") if pvals[i] is None else pvals[i]))
    sorted_p = [1.0 if pvals[i] is None else max(0.0, min(1.0, pvals[i])) for i in idx_sorted]

    q_sorted = [0.0] * m
    prev = 1.0
    for rank, p in enumerate(sorted_p, start=1):
        q = p * m / rank
        if q < 0: q = 0.0
        if q > 1: q = 1.0
        prev = min(prev, q)
        q_sorted[-rank] = prev  # 注意是从后往前填充

    # 还原原顺序
    qvals = [0.0] * m
    for pos, i in enumerate(idx_sorted):
        qvals[i] = q_sorted[pos]
    return qvals

# ------------------------------ family 结果读取 ------------------------------

FAMILY_COL_CANDIDATES = ["Family", "Orthogroup", "OG", "family", "orthogroup"]
PVALUE_COL_CANDIDATES = ["P", "p", "Pvalue", "pvalue", "P-value", "p-value", "Prob", "prob"]

def try_read_family_table(path: Path) -> Tuple[List[str], List[List[object]]]:
    """
    读取 CAFE5 family 结果（Gamma_family_results.txt / Base_family_results.txt）。
    返回：header, rows（列表格式），其中至少包含 family 与 p 值列。
    - 优先用 pandas（若可用），否则手写解析器（TSV / 多空格都兼容）。
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"空文件：{path}")

    # 简单判断分隔符：若首行包含制表符则按 TSV，否则按任意空白切分
    if "\t" in lines[0]:
        sep_tab = True
    else:
        sep_tab = False

    # 用 pandas 更稳健（若存在）
    if pd is not None:
        try:
            if sep_tab:
                df = pd.read_csv(path, sep="\t", dtype=str, engine="python")
            else:
                # 以任意空白分隔
                df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")
            header = list(df.columns)
            rows = df.fillna("").values.tolist()
            return header, rows
        except Exception:
            # 退回手动解析
            pass

    # 手动解析
    if sep_tab:
        parts = [ln.split("\t") for ln in lines]
    else:
        parts = [re.split(r"\s+", ln.strip()) if ln.strip() else [] for ln in lines]

    if not parts or not parts[0]:
        raise ValueError(f"无法解析表头：{path}")

    header = parts[0]
    rows = parts[1:]
    # 对齐列数（避免某些行末尾空值导致列数不齐）
    ncol = len(header)
    fixed_rows: List[List[str]] = []
    for r in rows:
        if len(r) < ncol:
            r = r + [""] * (ncol - len(r))
        elif len(r) > ncol:
            r = r[:ncol]
        fixed_rows.append(r)
    return header, fixed_rows

def pick_col_index(header: List[str], cand: List[str]) -> Optional[int]:
    """在表头中选出第一个匹配候选的列索引（大小写宽松匹配）。"""
    lower = [h.lower() for h in header]
    for name in cand:
        name_l = name.lower()
        if name_l in lower:
            return lower.index(name_l)
    return None

def extract_family_pvalues(path: Path) -> Tuple[List[Tuple[str, float]], List[str]]:
    """
    从 family 结果文件中抽取 (family, p-value) 列表。
    返回：(pairs, problems)；problems 记录跳过行或格式问题，方便日志告警。
    """
    header, rows = try_read_family_table(path)
    fi = pick_col_index(header, FAMILY_COL_CANDIDATES)
    pi = pick_col_index(header, PVALUE_COL_CANDIDATES)
    if fi is None or pi is None:
        raise ValueError(f"未找到必要列：family={FAMILY_COL_CANDIDATES}, p={PVALUE_COL_CANDIDATES} in {path}")

    pairs: List[Tuple[str, float]] = []
    problems: List[str] = []
    for row in rows:
        try:
            fam = str(row[fi]).strip()
            if fam == "" or fam.lower() in ("nan", "none"):
                continue
            p = float(str(row[pi]).strip())
            if math.isnan(p):
                continue
            p = min(max(p, 0.0), 1.0)
            pairs.append((fam, p))
        except Exception as e:
            problems.append(f"[WARN] 跳过异常行 @ {path.name}: {e}")
    return pairs, problems

# ------------------------------ λ 解析（results 优先，log 兜底） ------------------------------

def parse_lambda_from_results(path: Path) -> List[Tuple[str, float, str]]:
    """
    尝试从 *_results.txt 提取每条分支的 lambda 值。
    不同版本/模型输出略有差异，这里尽量宽松匹配：
      - 区块标题如 "Rates per branch" / "Branch rates" 等
      - 行形如：BranchName <tab/space> <lambda> 或 "branch <name> lambda = <val>"
    返回：[(branch, lambda, "results"), ...]
    """
    if not path.exists():
        return []
    out = []
    lines = read_lines(path)

    # 先尝试“branch <name> lambda = <val>”
    pat1 = re.compile(r"branch\s+([^\s:=]+)\s+lambda\s*=\s*([0-9eE\+\-\.]+)")
    for ln in lines:
        m = pat1.search(ln)
        if m:
            b = m.group(1)
            try:
                lam = float(m.group(2))
                out.append((b, lam, "results"))
            except Exception:
                pass

    # 再尝试“Rates per branch”块：形式为 <branch>\t<lambda> 或 <branch> <lambda>
    # 通过检测包含 lambda 数字的行做宽松解析
    pat_num = re.compile(r"([0-9eE\+\-\.]+)")
    for ln in lines:
        # 跳过明显是标题的行
        if not ln or "branch" in ln.lower() and "lambda" in ln.lower():
            continue
        parts_tab = None
        if "\t" in ln:
            parts_tab = ln.split("\t")
        # 以制表符为主的键值对
        if parts_tab and len(parts_tab) >= 2:
            left = parts_tab[0].strip()
            right = parts_tab[1].strip()
            # right 应该是数字
            try:
                lam = float(right)
                if left and not left.lower().startswith(("base_", "gamma_", "family", "clade", "change")):
                    out.append((left, lam, "results"))
                continue
            except Exception:
                pass

        # 以空格分隔的 fallback：最后一个 token 是数字
        toks = re.split(r"\s+", ln.strip())
        if len(toks) >= 2:
            try:
                lam = float(toks[-1])
                branch = " ".join(toks[:-1]).strip()
                # 简单过滤非分支行
                if branch and len(branch) < 128 and not branch.lower().startswith(("base_", "gamma_", "family", "clade", "change")):
                    out.append((branch, lam, "results"))
            except Exception:
                pass

    # 去重（以 branch 为键保留首次）
    seen = set()
    uniq = []
    for b, lam, src in out:
        if b not in seen:
            seen.add(b)
            uniq.append((b, lam, src))
    return uniq

def parse_lambda_from_runlog(path: Path) -> List[Tuple[str, float, str]]:
    """
    从 run.log 兜底提取 lambda；匹配：
      - "branch <name> lambda = <val>"
      - 或 "# CMD:" 行仅用于标注 -e，不含 λ 直接值
    """
    if not path.exists():
        return []
    out = []
    lines = read_lines(path)
    pat = re.compile(r"branch\s+([^\s:=]+)\s+lambda\s*=\s*([0-9eE\+\-\.]+)")
    for ln in lines:
        m = pat.search(ln)
        if m:
            b = m.group(1)
            try:
                lam = float(m.group(2))
                out.append((b, lam, "runlog"))
            except Exception:
                pass
    # 去重
    seen = set()
    uniq = []
    for b, lam, src in out:
        if b not in seen:
            seen.add(b)
            uniq.append((b, lam, src))
    return uniq

def detect_error_model_in_runlog(path: Path) -> bool:
    """
    检测是否启用误差模型（-e），通过 run.log 头部的 "# CMD:" 记录。
    """
    if not path.exists():
        return False
    lines = read_lines(path)
    for ln in lines[:50]:  # 头部有限扫描即可
        if ln.startswith("# CMD:") and "-e" in ln:
            return True
    return False

# ------------------------------ 主流程 ------------------------------

def main() -> None:
    # 1) 读取配置
    cfg_path = Path(CONFIG_PATH).resolve()
    if not cfg_path.exists():
        print(f"[ERR] 配置文件不存在：{cfg_path}", file=sys.stderr)
        sys.exit(2)
    cfg = load_yaml(cfg_path)

    # 2) 解析关键路径（与 12 对口）
    paths = cfg.get("paths", {})
    cafe_run_dir = Path(paths.get("cafe_run_dir", "results/06_cafe")).resolve()
    cafe_models_dir = cafe_run_dir / "models"
    cafe_agg_dir = Path(paths.get("cafe_agg_dir", str(cafe_run_dir / "agg"))).resolve()
    logs_dir = Path(paths.get("logs_dir", "logs")).resolve()

    mkdir_p(logs_dir)
    mkdir_p(cafe_agg_dir)

    # 3) 设置日志
    logfile = logs_dir / LOG_FILE_BASENAME
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(logfile, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("========== APhylo 13 — CAFE5 聚合（严丝合缝版） ==========")
    logging.info(f"[Info] 使用配置：{cfg_path}")
    logging.info(f"[Info] cafe_run_dir = {cafe_run_dir}")
    logging.info(f"[Info] cafe_models_dir = {cafe_models_dir}")
    logging.info(f"[Info] cafe_agg_dir = {cafe_agg_dir}")
    logging.info(f"[Info] logs_dir = {logs_dir}")

    # 4) 前置哨兵校验
    sentinel = cafe_run_dir / ".cafe.done"
    if not sentinel.exists():
        logging.error(f"[ERR] 未发现 12 步完成哨兵：{sentinel}")
        sys.exit(2)

    # 5) BH-FDR alpha
    report_cfg = cfg.get("report", {})
    fdr_alpha = float(report_cfg.get("fdr_alpha", 0.05))
    if not (0.0 < fdr_alpha <= 1.0):
        logging.warning(f"[WARN] report.fdr_alpha 数值异常：{fdr_alpha}，重置为 0.05")
        fdr_alpha = 0.05
    logging.info(f"[Info] 使用 BH-FDR α = {fdr_alpha}")

    # 6) 枚举 models
    models = list_models(cafe_models_dir)
    if not models:
        logging.error(f"[ERR] 未在 {cafe_models_dir} 下发现任何模型目录")
        sys.exit(2)
    logging.info(f"[Info] 发现模型：{', '.join(models)}")

    # 7) 聚合容器
    # families_all: 收集所有 (model, subset, family, P, source)
    families_all: List[Tuple[str, str, str, float, str]] = []
    # inputs_used 记录
    inputs_used_rows: List[List[object]] = []
    # 分支 λ
    branch_rows: List[List[object]] = []

    # 8) 遍历各 model
    for model in models:
        model_dir = cafe_models_dir / model
        # 必备目录
        primary_dir = model_dir / "primary_global"
        large_dir = model_dir / "large"
        flags_dir = model_dir / "flags"

        # 8.1 primary_global 校验与 family 文件确定
        if not primary_dir.exists():
            logging.error(f"[ERR] 缺少目录：{primary_dir}")
            sys.exit(2)

        primary_gamma = primary_dir / "Gamma_family_results.txt"
        primary_base = primary_dir / "Base_family_results.txt"
        if primary_gamma.exists():
            primary_family_file = primary_gamma
            primary_source = "Gamma"
        elif primary_base.exists():
            primary_family_file = primary_base
            primary_source = "Base"
            logging.warning(f"[WARN] {model}/primary_global 缺少 Gamma_family_results.txt，已回退 Base_family_results.txt")
        else:
            logging.error(f"[ERR] {model}/primary_global 未发现 family 结果（Gamma/Base），无法继续")
            sys.exit(2)

        # 8.2 large（可选） family 文件
        large_family_file = None
        if large_dir.exists():
            lf = large_dir / "Base_family_results.txt"
            if lf.exists():
                large_family_file = lf
            else:
                logging.info(f"[Info] {model}/large 目录存在，但未发现 Base_family_results.txt（可能未完成 second stage）")
        else:
            logging.info(f"[Info] {model} 未启用 two-stage（无 large 目录）")

        # 8.3 读 high_fail 列表（可无）
        high_fail_list = []
        high_fail_path = flags_dir / "high_fail_ogs.list"
        if high_fail_path.exists():
            high_fail_list = [x.strip() for x in read_lines(high_fail_path) if x.strip()]
            logging.info(f"[Info] {model} 读取 high_fail_ogs：{len(high_fail_list)} 条")
            # 记录 inputs_used
            sz, mt = file_stat(high_fail_path)
            inputs_used_rows.append([model, "-", str(high_fail_path), "flags_high_fail", "N/A", sz, human_bytes(sz), mt, ""])
        else:
            logging.info(f"[Info] {model} 未发现 high_fail_ogs.list，按空集处理")

        # 8.4 标注是否启用误差模型（-e），仅用于 inputs_used 记录
        primary_runlog = primary_dir / "run.log"
        used_error_model = detect_error_model_in_runlog(primary_runlog)
        if used_error_model:
            logging.info(f"[Info] {model}/primary_global 检测到误差模型（-e）")

        # 8.5 读取 primary family 结果
        try:
            pairs, probs = extract_family_pvalues(primary_family_file)
            for m in probs:
                logging.warning(m)
            for fam, p in pairs:
                families_all.append((model, "primary", fam, p, f"primary:{primary_source}"))
            sz, mt = file_stat(primary_family_file)
            inputs_used_rows.append([model, "primary", str(primary_family_file), "family_results", primary_source, sz, human_bytes(sz), mt, "error_model" if used_error_model else ""])
        except Exception as e:
            logging.error(f"[ERR] 读取 {primary_family_file} 失败：{e}")
            sys.exit(2)

        # 8.6 读取 large family 结果（若存在）
        if large_family_file is not None:
            try:
                pairs, probs = extract_family_pvalues(large_family_file)
                for m in probs:
                    logging.warning(m)
                for fam, p in pairs:
                    families_all.append((model, "large", fam, p, "large:Base"))
                sz, mt = file_stat(large_family_file)
                inputs_used_rows.append([model, "large", str(large_family_file), "family_results", "Base", sz, human_bytes(sz), mt, ""])
            except Exception as e:
                logging.warning(f"[WARN] 读取 {large_family_file} 失败：{e}（此模型将仅聚合 primary）")

        # 8.7 分支 λ 摘要：results 优先，log 兜底
        def harvest_lambda(subset_dir: Path, subset_name: str) -> None:
            # 先选结果文件：subset=primary 优先 Gamma_results.txt 再 Base_results.txt；subset=large 仅 Base_results.txt
            results_file = None
            if subset_name == "primary":
                g = subset_dir / "Gamma_results.txt"
                b = subset_dir / "Base_results.txt"
                if g.exists():
                    results_file = g
                elif b.exists():
                    results_file = b
            else:
                b = subset_dir / "Base_results.txt"
                if b.exists():
                    results_file = b

            runlog = subset_dir / "run.log"

            lambda_items = []
            if results_file is not None:
                lambda_items = parse_lambda_from_results(results_file)

                # 记录 inputs_used
                sz, mt = file_stat(results_file)
                src = "Gamma" if results_file.name.startswith("Gamma_") else "Base"
                inputs_used_rows.append([model, subset_name, str(results_file), "results", src, sz, human_bytes(sz), mt, ""])

            # 若结果文件未取到任何 λ，则尝试 run.log
            if not lambda_items:
                lambda_items = parse_lambda_from_runlog(runlog)

            # 记录 run.log
            if runlog.exists():
                sz, mt = file_stat(runlog)
                inputs_used_rows.append([model, subset_name, str(runlog), "run_log", "-", sz, human_bytes(sz), mt, "error_model" if (subset_name == "primary" and used_error_model) else ""])

            # 写入到 branch_rows
            for bname, lam, src in lambda_items:
                branch_rows.append([model, subset_name, bname, lam, src])

        harvest_lambda(primary_dir, "primary")
        if large_dir.exists():
            harvest_lambda(large_dir, "large")

    # 9) 生成全量显著表（先计算全局 BH-FDR）
    if not families_all:
        logging.error("[ERR] 未收集到任何 family 记录，无法聚合")
        sys.exit(2)

    # 9.1 计算全局 q 值
    pvals = [x[3] for x in families_all]
    qvals = bh_fdr(pvals)

    # 9.2 写全量：cafe_significant_families.tsv
    all_rows: List[List[object]] = []
    for (model, subset, family, p, source), q in zip(families_all, qvals):
        sig = "yes" if (q <= fdr_alpha) else "no"
        all_rows.append([model, subset, family, f"{p:.6g}", f"{q:.6g}", sig, source])

    sig_path = cafe_agg_dir / "cafe_significant_families.tsv"
    write_tsv(sig_path, ["model", "subset", "family", "P", "Q", "significant", "source"], all_rows)
    logging.info(f"[OK] 写出：{sig_path}（共 {len(all_rows)} 行）")

    # 9.3 剔除高失败率版本：cafe_significant_families_no_highfail.tsv
    #    注意：各 model 的 high_fail 集不同，这里逐 model 过滤
    #    先构建 model -> set(high_fail) 映射（上面 8.3 中未显式保存，现重建一遍）
    model_highfail: Dict[str, set] = {}
    for model in models:
        hf = (cafe_models_dir / model / "flags" / "high_fail_ogs.list")
        if hf.exists():
            model_highfail[model] = set([x.strip() for x in read_lines(hf) if x.strip()])
        else:
            model_highfail[model] = set()

    nofail_rows: List[List[object]] = []
    for row in all_rows:
        model, subset, fam = row[0], row[1], row[2]
        if fam not in model_highfail.get(model, set()):
            nofail_rows.append(row)

    sig_nofail_path = cafe_agg_dir / "cafe_significant_families_no_highfail.tsv"
    write_tsv(sig_nofail_path, ["model", "subset", "family", "P", "Q", "significant", "source"], nofail_rows)
    logging.info(f"[OK] 写出：{sig_nofail_path}（共 {len(nofail_rows)} 行）")

    # 10) 写分支 λ 摘要
    branch_path = cafe_agg_dir / "cafe_branch_summary.tsv"
    write_tsv(branch_path, ["model", "subset", "branch", "lambda", "extra"], branch_rows)
    logging.info(f"[OK] 写出：{branch_path}（共 {len(branch_rows)} 行）")

    # 11) 写 inputs_used.tsv
    inputs_path = cafe_agg_dir / "inputs_used.tsv"
    write_tsv(inputs_path,
              ["model", "subset", "file", "role", "source", "bytes", "size_h", "mtime", "note"],
              inputs_used_rows)
    logging.info(f"[OK] 写出：{inputs_path}（共 {len(inputs_used_rows)} 行）")

    logging.info("========== APhylo 13 — 完成 ==========")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as se:
        raise
    except Exception as e:
        # 兜底异常报告，避免静默失败
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
