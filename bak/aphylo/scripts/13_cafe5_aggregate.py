#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_cafe5_aggregate.py —— CAFE5 聚合与汇总（修正版）

职责概述：
  1) 汇总 models/<model>/primary_global/ 与 large/ 下的 family 结果：
       - primary：Gamma_family_results.txt 优先，缺则 Base_family_results.txt
       - large：  仅 Base_family_results.txt（若存在）
       - 严格要求：首列 #FamilyID/FamilyID，第2列 pvalue
  2) 对所有 family 的 p 值做一次全局 BH-FDR（修正算法，保证单调性）：
       - α = config.report.fdr_alpha（默认 0.05）
       - 输出：
           cafe_significant_families.tsv
             model subset family P Q significant source
           cafe_significant_families_no_highfail.tsv
             剔除 flags/high_fail_ogs.list 中列出的 OG
  3) 分支 λ 摘要：
       - primary：Gamma_results.txt 优先，缺则 Base_results.txt
       - large：  Base_results.txt（若存在）
       - 输出：cafe_branch_summary.tsv
  4) 家族“宇宙”汇总：
       - 汇总每个 family 是否出现在 primary / large，是否被标记 high_fail
       - 输出：cafe_family_universe.tsv
             model family in_primary in_large is_high_fail
  5) Clade 汇总（物种/内部节点的扩张收缩计数）：
       - 使用 Base_clade_results.txt 或 Gamma_clade_results.txt
       - 标准格式（示例）：
           #Taxon_ID   Increase Decrease
           Sinonovacula_constricta<10>  541 441
       - 输出：cafe_clade_summary.tsv
             model taxon increase decrease source

重要口径：
  - 只从 config.yaml 中读取路径与模型名，不扫描、不猜测；
  - 不改变 14_joint_integration 等下游所依赖的表头与语义。
"""

from __future__ import annotations

# ============================ 顶部可调区 ============================
CONFIG_PATH: str = "config.yaml"                 # 配置文件路径
LOG_LEVEL: str = "INFO"                          # 日志等级：DEBUG / INFO / WARNING
LOG_FILE_BASENAME: str = "13_cafe5_aggregate.log"
# ====================================================================

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
    """类似 mkdir -p。"""
    p.mkdir(parents=True, exist_ok=True)

def load_yaml(path: Path) -> dict:
    """读取 YAML 配置。"""
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def read_lines(path: Path) -> List[str]:
    """读取文本文件所有行（去掉行尾换行符）。"""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\r\n") for ln in f]

def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    """简单 TSV 写出工具。"""
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

def file_stat(p: Path) -> Tuple[int, str]:
    """返回 (文件大小字节数, mtime ISO字符串)。"""
    try:
        st = p.stat()
        return st.st_size, datetime.fromtimestamp(st.st_mtime).isoformat(sep=" ", timespec="seconds")
    except Exception:
        return 0, ""

def human_bytes(n: int) -> str:
    """人类友好文件大小。"""
    units = ["B", "KB", "MB", "GB", "TB", "PB", "EB"]
    x = float(n)
    for u in units:
        if x < 1024:
            return f"{x:.1f}{u}"
        x /= 1024
    return f"{x:.1f}ZB"

def list_models(models_dir: Path) -> List[str]:
    """列出 models 目录下所有模型子目录名。"""
    if not models_dir.exists():
        return []
    return sorted([p.name for p in models_dir.iterdir() if p.is_dir()])

# --------------------------- BH-FDR（修正版） ---------------------------

def bh_fdr(pvals: List[float]) -> List[float]:
    """
    Benjamini–Hochberg FDR 校正（返回与输入顺序一一对应的 q 值列表）。

    正确算法要点：
      - 先按 p 从小到大排序，计算 q_raw = p_i * m / rank_i
      - 再在排序空间中从后往前扫描：
          q(i) = min(q_raw(i), q(i+1))
      - 最后把 q(i) 映射回原始顺序
    """
    m = len(pvals)
    if m == 0:
        return []
    # 按 p 从小到大排序，order 存的是原始索引
    order = sorted(range(m), key=lambda i: (float("inf") if math.isnan(pvals[i]) else pvals[i]))
    q = [1.0] * m
    prev = 1.0
    # 反向扫描，保证单调性
    for rank, i in reversed(list(enumerate(order, start=1))):
        p = pvals[i]
        if math.isnan(p) or p < 0:
            qi = 1.0
        else:
            qi = p * m / rank
        if qi > 1.0:
            qi = 1.0
        if qi < 0.0:
            qi = 0.0
        prev = min(prev, qi)
        q[i] = prev
    return q

# --------------------------- 读取 family 结果（严格口径） ---------------------------

def read_family_pairs_strict(path: Path) -> List[Tuple[str, float]]:
    """
    严格读取 CAFE5 的 *family_results.txt：
      - 首行必须包含 '#FamilyID' 或 'FamilyID'、'pvalue' 作为前两列
      - 每行：family_id  pvalue  ...
    返回：[(family, pvalue), ...]
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
        raise ValueError(f"第2列表头必须为 'pvalue'：{path}\n实际：{header[1]}")

    out: List[Tuple[str, float]] = []
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = re.split(r"\s+", ln.strip())
        if len(parts) < 2:
            continue
        fam = parts[0]
        try:
            p = float(parts[1])
        except Exception:
            continue
        out.append((fam, p))
    return out

# --------------------------- 读取 *results.txt 中的 branch λ ---------------------------

def parse_lambda_from_results(path: Path) -> List[Tuple[str, float]]:
    """
    从 *results.txt 提取 (branch, lambda)：
      - 优先匹配 "branch <name> lambda = <val>"
      - 若无则尝试按 "<branch>\\t<lambda>" 或 "<branch> <lambda>"（末列为数字）
    文件不存在或解析不到时返回空列表。
    """
    if not path.exists():
        return []
    lines = read_lines(path)
    out: List[Tuple[str, float]] = []

    # 模式1：branch xxx lambda = 0.00123
    pat1 = re.compile(r"branch\s+([^\s:=]+)\s+lambda\s*=\s*([0-9eE+\-\.]+)")
    for ln in lines:
        m = pat1.search(ln)
        if m:
            b = m.group(1)
            try:
                lam = float(m.group(2))
            except Exception:
                continue
            out.append((b, lam))

    if out:
        # 去重（保持首次出现）
        seen = set()
        res = []
        for b, lam in out:
            if b not in seen:
                seen.add(b)
                res.append((b, lam))
        return res

    # 模式2：列形式，末列为数值
    tmp: List[Tuple[str, float]] = []
    for ln in lines:
        if not ln.strip():
            continue
        parts = re.split(r"\s+", ln.strip())
        if len(parts) < 2:
            continue
        try:
            lam = float(parts[-1])
        except Exception:
            continue
        branch = " ".join(parts[:-1]).strip()
        if not branch:
            continue
        # 排除明显不是分支名的行
        if branch.lower() in ("family", "clade", "change"):
            continue
        tmp.append((branch, lam))

    seen = set()
    res = []
    for b, lam in tmp:
        if b not in seen:
            seen.add(b)
            res.append((b, lam))
    return res

# --------------------------- 解析 Base_clade_results.txt ---------------------------

def parse_clade_results(path: Path) -> List[Tuple[str, int, int]]:
    """
    解析 CAFE5 的 Base_clade_results.txt / Gamma_clade_results.txt

    兼容的表头示例：
      #Taxon_ID       Increase        Decrease
      Taxon           Increase        Decrease
      taxon           Increase        Decrease

    识别策略：
      - Increase / Decrease：严格匹配列名 "Increase" / "Decrease"
      - Taxon 列：
          * 对每个列名做规范化：
              去首位 '#'，去结尾的 '_ID' / '_id'，再 lower()
              例如：'#Taxon_ID' -> 'taxon'
          * 若规范化后包含子串 'taxon'，视为 Taxon 列
          * 若仍找不到，则兜底使用第 0 列

    返回：[(taxon, increase, decrease), ...]
    """
    if not path.exists():
        return []
    lines = read_lines(path)
    if not lines:
        return []
    header = re.split(r"\s+", lines[0].strip())
    if not header:
        return []

    # 1) 找 Increase / Decrease 列（严格匹配）
    try:
        i_inc = header.index("Increase")
        i_dec = header.index("Decrease")
    except ValueError:
        # 表头不含 Increase/Decrease，直接视为不支持
        return []

    # 2) 找 Taxon 列（宽松匹配）
    i_tax = None
    for idx, name in enumerate(header):
        # 去掉开头的 '#'，去掉结尾的 '_ID' 或 '_id'，再小写
        norm = name.lstrip("#")
        norm = re.sub(r"_id$", "", norm, flags=re.IGNORECASE)
        norm_lower = norm.lower()
        if "taxon" in norm_lower:
            i_tax = idx
            break

    # 若遍历后仍未找到，则兜底使用第 0 列
    if i_tax is None:
        if len(header) > 0:
            i_tax = 0
        else:
            return []

    out: List[Tuple[str, int, int]] = []
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = re.split(r"\s+", ln.strip())
        if len(parts) <= max(i_tax, i_inc, i_dec):
            continue
        tax = parts[i_tax]
        try:
            inc = int(parts[i_inc])
            dec = int(parts[i_dec])
        except Exception:
            continue
        out.append((tax, inc, dec))
    return out

# --------------------------- 日志初始化 ---------------------------

def setup_logging(logs_dir: Path) -> None:
    """初始化 logging：同时打屏幕和文件。"""
    mkdir_p(logs_dir)
    log_file = logs_dir / LOG_FILE_BASENAME

    fmt = "[%(asctime)s] [%(levelname)s] %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format=fmt,
        datefmt=datefmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, encoding="utf-8"),
        ],
    )
    logging.info(f"[init] 使用配置文件：{CONFIG_PATH}")
    logging.info(f"[init] 日志写入：{log_file}")

# --------------------------- 主流程 ---------------------------

def main() -> None:
    cfg = load_yaml(Path(CONFIG_PATH))

    # 1) 路径解析
    paths_cfg = cfg.get("paths", {})
    cafe_run_dir = Path(paths_cfg.get("cafe_run_dir", "results/06_cafe")).resolve()
    models_dir = cafe_run_dir / "models"
    agg_dir = Path(paths_cfg.get("cafe_agg_dir", "results/07_cafeagg")).resolve()
    logs_dir = Path(paths_cfg.get("logs_dir", "logs")).resolve()

    setup_logging(logs_dir)
    mkdir_p(agg_dir)

    logging.info(f"[init] cafe_run_dir = {cafe_run_dir}")
    logging.info(f"[init] models_dir   = {models_dir}")
    logging.info(f"[init] cafe_agg_dir = {agg_dir}")

    # 2) 模型列表（优先用 cafe5.models，其次 report.selected_cafe_model，最后默认为 ['global']）
    cafe5_cfg = cfg.get("cafe5", {})
    report_cfg = cfg.get("report", {})
    models = cafe5_cfg.get("models") or []
    if not models:
        sel = report_cfg.get("selected_cafe_model")
        if sel:
            models = [sel]
    if not models:
        models = ["global"]

    models = list(dict.fromkeys(models))   # 去重保持顺序
    logging.info(f"[init] 将处理模型：{models}")

    # 3) FDR 阈值
    alpha = float(report_cfg.get("fdr_alpha", 0.05))
    logging.info(f"[init] FDR α = {alpha}")

    if not models_dir.exists():
        logging.error(f"[ERR] 模型目录不存在：{models_dir}")
        sys.exit(2)

    # 汇总容器
    fam_records: List[Tuple[str, str, str, float, str]] = []  # (model, subset, family, p, source)
    branch_rows: List[List[object]] = []                      # [model, subset, branch, lambda, source]
    clade_rows: List[List[object]] = []                       # [model, taxon, increase, decrease, source]
    inputs_used: List[List[object]] = []                      # inputs 溯源
    model_high_fail: Dict[str, set] = {}                      # model -> {family}

    # 4) 遍历每个 model
    for model in models:
        mdir = models_dir / model
        primary_dir = mdir / "primary_global"
        large_dir = mdir / "large"
        flags_dir = mdir / "flags"

        if not primary_dir.exists():
            logging.error(f"[ERR] 缺少 primary 目录：{primary_dir}")
            sys.exit(2)

        logging.info(f"[model={model}] primary_dir = {primary_dir}")
        logging.info(f"[model={model}] large_dir   = {large_dir}")

        # -------- primary family：Gamma 优先，缺则 Base 回退 --------
        pf_gamma = primary_dir / "Gamma_family_results.txt"
        pf_base = primary_dir / "Base_family_results.txt"

        primary_family_file: Path | None = None
        primary_source = ""
        if pf_gamma.exists():
            primary_family_file = pf_gamma
            primary_source = "primary:Gamma"
        elif pf_base.exists():
            primary_family_file = pf_base
            primary_source = "primary:Base"
            logging.info(f"[model={model}] 未找到 Gamma_family_results.txt，已回退 Base_family_results.txt")
        else:
            logging.error(f"[ERR] {model}/primary_global 下找不到 Gamma/Base family 结果")
            sys.exit(2)

        # 读取 primary family 结果
        try:
            pairs = read_family_pairs_strict(primary_family_file)
        except Exception as e:
            logging.error(f"[ERR] 读取 primary family 结果失败：{primary_family_file} -> {e}")
            sys.exit(2)

        for fam, p in pairs:
            fam_records.append((model, "primary", fam, p, primary_source))

        sz, mt = file_stat(primary_family_file)
        inputs_used.append([
            model, "primary", str(primary_family_file), "family_results",
            primary_source, sz, human_bytes(sz), mt, "",
        ])

        # -------- large family：仅 Base_family_results.txt（若目录存在且文件存在） --------
        large_family_file: Path | None = None
        if large_dir.exists():
            lf_base = large_dir / "Base_family_results.txt"
            if lf_base.exists():
                large_family_file = lf_base
                logging.info(f"[model={model}] 发现 large Base_family_results.txt：{lf_base}")
            else:
                logging.info(f"[model={model}] large 目录存在，但缺少 Base_family_results.txt：{lf_base}")

        if large_family_file is not None:
            try:
                lpairs = read_family_pairs_strict(large_family_file)
                for fam, p in lpairs:
                    fam_records.append((model, "large", fam, p, "large:Base"))
                sz, mt = file_stat(large_family_file)
                inputs_used.append([
                    model, "large", str(large_family_file), "family_results",
                    "large:Base", sz, human_bytes(sz), mt, "",
                ])
            except Exception as e:
                logging.warning(f"[WARN] 读取 large family 结果失败（已忽略）：{large_family_file} -> {e}")

        # -------- 分支 λ 摘要：primary Gamma→Base，large Base --------
        # primary
        g_res = primary_dir / "Gamma_results.txt"
        b_res = primary_dir / "Base_results.txt"
        src = None
        lam_list: List[Tuple[str, float]] = []
        if g_res.exists():
            lam_list = parse_lambda_from_results(g_res)
            if lam_list:
                src = "primary:Gamma"
        if (not lam_list) and b_res.exists():
            lam_list = parse_lambda_from_results(b_res)
            if lam_list:
                src = "primary:Base"

        if lam_list and src is not None:
            res_file = g_res if "Gamma" in src else b_res
            sz, mt = file_stat(res_file)
            inputs_used.append([
                model, "primary", str(res_file), "results",
                src, sz, human_bytes(sz), mt, "",
            ])
            for br, lam in lam_list:
                branch_rows.append([model, "primary", br, lam, src])

        # large
        if large_dir.exists():
            l_res = large_dir / "Base_results.txt"
            if l_res.exists():
                lam_list_l = parse_lambda_from_results(l_res)
                if lam_list_l:
                    sz, mt = file_stat(l_res)
                    inputs_used.append([
                        model, "large", str(l_res), "results",
                        "large:Base", sz, human_bytes(sz), mt, "",
                    ])
                    for br, lam in lam_list_l:
                        branch_rows.append([model, "large", br, lam, "large:Base"])

        # -------- clade 结果：Base_clade_results.txt 或 Gamma_clade_results.txt --------
        base_clade = primary_dir / "Base_clade_results.txt"
        gamma_clade = primary_dir / "Gamma_clade_results.txt"
        clade_file: Path | None = None
        clade_source = ""
        if base_clade.exists():
            clade_file = base_clade
            clade_source = "Base"
        elif gamma_clade.exists():
            clade_file = gamma_clade
            clade_source = "Gamma"

        if clade_file is not None:
            cr = parse_clade_results(clade_file)
            if cr:
                sz, mt = file_stat(clade_file)
                inputs_used.append([
                    model, "primary", str(clade_file), "clade_results",
                    clade_source, sz, human_bytes(sz), mt, "",
                ])
                for tax, inc, dec in cr:
                    clade_rows.append([model, tax, inc, dec, clade_source])
            else:
                logging.info(f"[model={model}] clade 结果文件存在但未解析到有效数据：{clade_file}")

        # -------- high_fail_ogs.list --------
        hf_set: set = set()
        hf_file = flags_dir / "high_fail_ogs.list"
        if hf_file.exists():
            for ln in read_lines(hf_file):
                s = ln.strip()
                if s:
                    hf_set.add(s)
            sz, mt = file_stat(hf_file)
            inputs_used.append([
                model, "flags", str(hf_file), "high_fail_ogs",
                "", sz, human_bytes(sz), mt, "",
            ])
            logging.info(f"[model={model}] high_fail_ogs.list 中共有 {len(hf_set)} 个 OG")
        else:
            logging.info(f"[model={model}] 未发现 high_fail_ogs.list")

        model_high_fail[model] = hf_set

    # 5) 计算 BH-FDR 并写出显著表
    logging.info("[step] 开始 BH-FDR 计算与表格输出")

    if not fam_records:
        logging.error("[ERR] 没有任何 family 记录，检查 12 步是否成功产生 family_results.txt")
        sys.exit(2)

    pvals = [p for (_, _, _, p, _) in fam_records]
    qvals = bh_fdr(pvals)

    all_rows: List[List[object]] = []
    for (model, subset, fam, p, src), q in zip(fam_records, qvals):
        sig = "yes" if q <= alpha else "no"
        all_rows.append([model, subset, fam, f"{p:.6g}", f"{q:.6g}", sig, src])

    sig_path = agg_dir / "cafe_significant_families.tsv"
    write_tsv(sig_path,
              ["model", "subset", "family", "P", "Q", "significant", "source"],
              all_rows)
    logging.info(f"[OK] 写出：{sig_path}（{len(all_rows)} 行）")

    # 6) 剔除 high_fail OG 得到 no_highfail 版本
    nofail_rows: List[List[object]] = []
    for row in all_rows:
        model, subset, fam = row[0], row[1], row[2]
        if fam not in model_high_fail.get(model, set()):
            nofail_rows.append(row)

    sig_nofail_path = agg_dir / "cafe_significant_families_no_highfail.tsv"
    write_tsv(sig_nofail_path,
              ["model", "subset", "family", "P", "Q", "significant", "source"],
              nofail_rows)
    logging.info(f"[OK] 写出：{sig_nofail_path}（{len(nofail_rows)} 行）")

    # 7) 分支 λ 摘要
    branch_path = agg_dir / "cafe_branch_summary.tsv"
    write_tsv(branch_path,
              ["model", "subset", "branch", "lambda", "source"],
              branch_rows)
    logging.info(f"[OK] 写出：{branch_path}（{len(branch_rows)} 行）")

    # 8) 家族宇宙表：cafe_family_universe.tsv
    #    汇总每个 (model, family) 是否出现在 primary / large，是否 high_fail
    fam_universe: Dict[Tuple[str, str], Dict[str, bool]] = {}
    for model, subset, fam, p, src in fam_records:
        key = (model, fam)
        rec = fam_universe.setdefault(key, {"in_primary": False, "in_large": False})
        if subset == "primary":
            rec["in_primary"] = True
        elif subset == "large":
            rec["in_large"] = True

    family_univ_rows: List[List[object]] = []
    for (model, fam), info in sorted(fam_universe.items(), key=lambda x: (x[0][0], x[0][1])):
        hf = fam in model_high_fail.get(model, set())
        family_univ_rows.append([
            model,
            fam,
            "yes" if info["in_primary"] else "no",
            "yes" if info["in_large"] else "no",
            "yes" if hf else "no",
        ])

    fam_univ_path = agg_dir / "cafe_family_universe.tsv"
    write_tsv(fam_univ_path,
              ["model", "family", "in_primary", "in_large", "is_high_fail"],
              family_univ_rows)
    logging.info(f"[OK] 写出：{fam_univ_path}（{len(family_univ_rows)} 行）")

    # 9) Clade 汇总表：cafe_clade_summary.tsv
    clade_path = agg_dir / "cafe_clade_summary.tsv"
    write_tsv(clade_path,
              ["model", "taxon", "increase", "decrease", "source"],
              clade_rows)
    logging.info(f"[OK] 写出：{clade_path}（{len(clade_rows)} 行）")

    # 10) inputs 溯源表：inputs_used.tsv
    inputs_path = agg_dir / "inputs_used.tsv"
    write_tsv(inputs_path,
              ["model", "subset", "file", "role", "source", "bytes", "size_h", "mtime", "note"],
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

