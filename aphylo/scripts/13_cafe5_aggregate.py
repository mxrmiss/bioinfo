#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_cafe5_aggregate.py —— 汇总 CAFE5（“加固聚合版”，仅改 13，不动 12）

设计要点（与 12 自救增强 v2 适配）：
  1) 仅聚合“最后一轮”的干净产物：
     - 只扫描 models/<model>/ 根目录；忽略 _round*_tmp/ 等子目录；
     - 同名/同类结果若有多份，仅取“最新时间窗口”的那批（以 mtime 为准，窗口默认 ±600s）。
  2) 读取 models/<model>/run.log 的“最后一个运行块”（12 会追加头部：
     "===== MODEL=<name> ROUND=<N> ====="）：
     - 若该块含 "Failed to initialize any reasonable values" 或
       "You may want to try removing the top few families..."，直接报错拒绝聚合；
     - 解析该块的 “branch ... lambda ...” 行，生成分支 λ 摘要；
     - 若未找到任何标准 TSV，回退从该块提取 family 与 p 值。
  3) 产出（与旧版口径一致）：
     - cafe_significant_families.tsv：列 [model, family, pvalue, qvalue, signif, source]
     - cafe_branch_summary.tsv     ：列 [model, branch, lambda_estimate, source]
     - 附加生成 inputs_used.tsv：记录本次聚合采用/忽略的文件（可追溯）

依赖：
  - 读取 config.yaml -> paths.cafe_run_dir, paths.cafe_agg_dir, paths.logs_dir
  - 需要 12 成功后产生的 models/<model>/run.log；若无则降级为“按最新 mtime 选文件 + 无法做失败阻断”。

注意：
  - 若 paths.cafe_run_dir/.cafe.done 不存在，视为 12 未成功，直接拒绝聚合。
"""

from __future__ import annotations
import sys, io, logging, re, os
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import yaml
import datetime

DEFAULT_CONFIG = "config.yaml"

# ========== 基础工具 ==========
def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(config_path)
    if not p.exists(): raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub: cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name); lg.setLevel(level); lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt); fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout);     sh.setFormatter(fmt); sh.setLevel(level)
    lg.addHandler(fh); lg.addHandler(sh)
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout = _Flush(sys.stdout); sys.stderr = _Flush(sys.stderr)
    return lg

def banner(log: logging.Logger, text: str):
    bar = "=" * max(10, len(text)+2)
    log.info(bar); log.info(f" {text} "); log.info(bar)

# ========== 文本解析 ==========

FAILED_INIT_PAT = re.compile(r"Failed to initialize any reasonable values", re.I)
ADVICE_PAT = re.compile(r"You may want to try removing the top few families", re.I)
BRANCH_LAMBDA_PAT = re.compile(r"branch[:\s]+(\S+).+lambda\s*[:=]\s*([\d\.Ee+-]+)", re.I)

def read_last_run_block(runlog: Path) -> Optional[str]:
    """
    读取 models/<model>/run.log 的“最后一个运行块”文本。
    运行块用 12 写入的头部 "===== MODEL=... ROUND=... =====" 分割。
    """
    if not runlog.is_file():
        return None
    text = runlog.read_text(encoding="utf-8", errors="ignore")
    # 以头部分割，取最后一段
    parts = re.split(r"^=+\s*MODEL=.*?ROUND=.*?=+\s*$", text, flags=re.M)
    if not parts:
        # 整个文件都没有头，视为单块
        return text
    return parts[-1].strip()

def block_has_failure(block: str) -> bool:
    """判断最后块是否包含初始化失败/官方删家族建议。"""
    if not block: return False
    return bool(FAILED_INIT_PAT.search(block) or ADVICE_PAT.search(block))

def parse_branch_lambda_from_block(block: str) -> List[Tuple[str, str]]:
    """从最后块提取分支 λ 摘要 (branch, lambda_str)。"""
    out = []
    if not block: return out
    for line in block.splitlines():
        m = BRANCH_LAMBDA_PAT.search(line)
        if m:
            out.append((m.group(1), m.group(2)))
    return out

# ========== 文件收集与“最新窗口” ==========

def list_candidate_result_files(mdir: Path) -> List[Path]:
    """
    仅列出模型根目录下的候选结果文件（忽略 _round*_tmp/ 等子目录）：
      - *.tsv / *.csv
      - 文件名包含 'family' 且（包含 'result' 或 'signif' 或 'report'）
    """
    out: List[Path] = []
    for p in mdir.iterdir():
        if not p.is_file(): 
            continue
        name = p.name.lower()
        if (name.endswith(".tsv") or name.endswith(".csv")) and ("family" in name):
            if ("result" in name) or ("signif" in name) or ("report" in name):
                out.append(p)
    return out

def select_files_latest_window(files: List[Path], window_sec: int = 600) -> Tuple[List[Path], List[Path]]:
    """
    在候选文件中选取“最新时间窗口”的那批：
      - 以所有候选的 max(mtime) 为锚；选择 mtime >= (max_mtime - window_sec) 的文件；
      - 返回 (selected, ignored)。
    """
    if not files:
        return [], []
    mtimes = [p.stat().st_mtime for p in files]
    max_m = max(mtimes)
    cut = max_m - window_sec
    selected, ignored = [], []
    for p in files:
        if p.stat().st_mtime >= cut:
            selected.append(p)
        else:
            ignored.append(p)
    return selected, ignored

# ========== 读取 family 结果表 & FDR ==========

def parse_family_table(path: Path) -> Tuple[List[str], List[List[str]]]:
    """读取 TSV/CSV（按制表符优先，逗号次之），返回 (header_cols, rows)。"""
    text = path.read_text(encoding="utf-8")
    # 简单探测分隔符
    delim = "\t" if ("\t" in text.splitlines()[0]) else ","
    rows = [line for line in text.splitlines() if line.strip()]
    if not rows: return [], []
    header = rows[0].split(delim)
    body = [r.split(delim) for r in rows[1:]]
    return header, body

def bh_fdr(values: List[Optional[float]]) -> List[Optional[float]]:
    """Benjamini–Hochberg（step-up）+ 单调性修正。None 位置传回 None。"""
    idx = [i for i,v in enumerate(values) if (v is not None and v>=0.0)]
    m = len(idx)
    if m == 0: return [None]*len(values)
    idx_sorted = sorted(idx, key=lambda i: values[i])
    qtemp = [None]*len(values)
    min_q = 1.0
    # 逆序单调修正：q(i) = min(prev, p*m/rank)
    for rank1, i in reversed(list(enumerate(idx_sorted, start=1))):
        p = float(values[i])
        q = p * m / rank1
        if q < min_q: min_q = q
        qtemp[i] = min_q if min_q < 1.0 else 1.0
    return qtemp

# ========== 主流程 ==========

def main():
    cfg = load_config()
    paths = cfg["paths"]; report_cfg = cfg.get("report", {}) or {}
    fdr_alpha = float(report_cfg.get("fdr_alpha", 0.05))

    logs_dir = Path(paths["logs_dir"])
    log = get_logger("aphylo.13", logs_dir / "13_cafe5_aggregate.log")
    banner(log, "APhylo 13 — CAFE5 汇总（加固聚合版）")

    cafe_dir = Path(paths["cafe_run_dir"])
    agg_dir  = ensure_dir(Path(paths["cafe_agg_dir"]))

    # 保险：没有 .cafe.done 直接拒绝聚合（避免 12 失败后的误聚合）
    done_flag = cafe_dir / ".cafe.done"
    if not done_flag.is_file():
        log.error(f"[ERR] 未发现 {done_flag} —— 12 步未成功完成，拒绝聚合。")
        sys.exit(2)

    models_root = cafe_dir / "models"
    if not models_root.is_dir():
        log.error(f"[ERR] 缺少 CAFE 模型目录：{models_root}")
        sys.exit(2)

    # 汇总容器
    fam_rows: List[Tuple[str, str, Optional[float], str, str]] = []  # (model, family, p_val, p_str, source)
    branch_rows: List[Tuple[str, str, str, str]] = []               # (model, branch, lambda_estimate, source)
    inputs_used_lines = ["model\tfile\tmtime\tbytes\tselected\treason"]

    # 遍历每个模型目录
    for mdir in sorted(models_root.glob("*")):
        if not mdir.is_dir(): continue
        model = mdir.name
        runlog = mdir / "run.log"

        # 解析最后一个运行块
        last_block = read_last_run_block(runlog) if runlog.is_file() else None
        if runlog.is_file():
            # 失败阻断：最后块若含失败/建议语，拒绝聚合
            if block_has_failure(last_block or ""):
                log.error(f"[ERR] 模型 {model} 的最后运行块仍包含初始化失败/删家族建议，拒绝聚合。请先解决 12 步错误。")
                sys.exit(2)

        # 分支 λ：仅从最后块提取（避免历史干扰）
        for br, lam in parse_branch_lambda_from_block(last_block or ""):
            branch_rows.append((model, br, lam, "run.log:last_block"))

        # 收集候选结果文件（仅根目录）
        cand_files = list_candidate_result_files(mdir)

        # 过滤出“最新窗口”的那批文件
        sel_files, ign_files = select_files_latest_window(cand_files, window_sec=600)

        # 记录 inputs_used（无论选/不选都记，便于追溯）
        def fmt_ts(ts): 
            return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")
        for p in sorted(cand_files):
            stat = p.stat()
            tag = "yes" if p in sel_files else "no"
            reason = "latest_window" if tag=="yes" else "stale_or_old"
            inputs_used_lines.append("\t".join([
                model, str(p.name), fmt_ts(stat.st_mtime), str(stat.st_size), tag, reason
            ]))

        # 若“最新窗口”内没有任何标准文件，回退用 run.log 的最后块解析 family p 值
        if not sel_files:
            if last_block:
                # 解析 family ... pvalue ... 行（宽松匹配）
                for line in last_block.splitlines():
                    m = re.search(r"family\s+(\S+).+p[- ]?value\s*[:=]\s*([\d\.Ee+-]+)", line, re.I)
                    if m:
                        fam = m.group(1); p_str = m.group(2)
                        try:
                            p_val = float(p_str)
                        except ValueError:
                            p_val = None
                        fam_rows.append((model, fam, p_val, p_str, "run.log:last_block"))
                log.warning(f"[WARN] 模型 {model} 未找到标准结果文件，已从 run.log 最后块回退提取 family/p。")
            else:
                log.error(f"[ERR] 模型 {model} 既无标准结果文件，也无 run.log 可用于回退。")
                sys.exit(2)
        else:
            # 从选中文件解析 family/p
            for ff in sel_files:
                try:
                    hdr, body = parse_family_table(ff)
                except Exception as e:
                    log.warning(f"[WARN] 无法解析 {ff.name}：{e}；将跳过。")
                    continue
                # 列定位（尽量宽容）
                name_lc = [c.strip().lower() for c in hdr]
                # family 列一般在首列
                idx_family = 0
                # p 值列名可能是 pvalue 或 p
                idx_p = None
                for key in ("pvalue", "p"):
                    if key in name_lc:
                        idx_p = name_lc.index(key); break
                for row in body:
                    fam = row[idx_family] if idx_family < len(row) else ""
                    p_str = row[idx_p] if (idx_p is not None and idx_p < len(row)) else ""
                    p_val = None
                    try:
                        if p_str != "": p_val = float(p_str)
                    except ValueError:
                        p_val = None
                    fam_rows.append((model, fam, p_val, p_str, ff.name))

    # ===== 计算 BH-FDR 并写出 =====
    # family 表
    p_list = [r[2] for r in fam_rows]
    q_list = bh_fdr(p_list)
    out_lines = ["model\tfamily\tpvalue\tqvalue\tsignif\tsource"]
    for (model, fam, p_val, p_str, src), qv in zip(fam_rows, q_list):
        q_str = f"{qv:.3g}" if isinstance(qv, float) else ""
        signif = "yes" if (isinstance(qv, float) and qv <= fdr_alpha) else ("no" if isinstance(qv, float) else "")
        out_lines.append("\t".join([model, fam, p_str, q_str, signif, src]))
    (agg_dir/"cafe_significant_families.tsv").write_text("\n".join(out_lines)+"\n", encoding="utf-8")

    # branch 表
    br_lines = ["model\tbranch\tlambda_estimate\tsource"]
    for model, br, lam, src in branch_rows:
        br_lines.append("\t".join([model, br, lam, src]))
    (agg_dir/"cafe_branch_summary.tsv").write_text("\n".join(br_lines)+"\n", encoding="utf-8")

    # inputs_used 清单
    (agg_dir/"inputs_used.tsv").write_text("\n".join(inputs_used_lines)+"\n", encoding="utf-8")

    print(f"[OK] 写出：{agg_dir/'cafe_significant_families.tsv'}")
    print(f"[OK] 写出：{agg_dir/'cafe_branch_summary.tsv'}")
    print(f"[OK] 写出：{agg_dir/'inputs_used.tsv'}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)