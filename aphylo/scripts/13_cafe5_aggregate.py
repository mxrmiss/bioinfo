#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
10_cafe5_aggregate.py — 汇总 CAFE5（输出到 cafe_agg_dir）
增强：在不改变文件名/目录/下游契约的前提下，为 cafe_significant_families.tsv
      新增两列：qvalue（BH-FDR）与 signif（q<=report.fdr_alpha ? yes : no）
"""

# ===== APhylo utils (config, logging, sentinels) =====
from __future__ import annotations
import sys, io, logging, subprocess, re
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable, Optional
import yaml

DEFAULT_CONFIG = "config.yaml"

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

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True); return p

def need_dir(p: Path, what: str):
    p = Path(p)
    if not p.is_dir(): raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p

def need_file(p: Path, what: str):
    p = Path(p)
    if not p.is_file(): raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

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

def banner(logger: logging.Logger, text: str):
    bar = "=" * max(10, len(text)+2); logger.info(bar); logger.info(f" {text} "); logger.info(bar)

def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    recs=[]; name=None; seq=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None: recs.append((name, "".join(seq)))
            name=line[1:].strip(); seq=[]
        else:
            seq.append(line.strip())
    if name is not None: recs.append((name, "".join(seq)))
    return recs

def write_fasta(path: Path, recs: Iterable[Tuple[str,str]]):
    with path.open("w", encoding="utf-8") as w:
        for h, s in recs:
            w.write(f">{h}\n")
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + "\n")
# ===== utils end =====

from pathlib import Path
import re

# --------- 解析 CAFE family 结果（尽量温和，不做强假设） ---------
def parse_family_table(path: Path):
    rows=[]; hdr=None
    for i,line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        if i==0: hdr=line; continue
        if not line.strip(): continue
        rows.append(line.split("\t"))
    return hdr, rows

# --------- 纯 Python 的 BH-FDR 实现（保持单文件可运行） ---------
def bh_fdr(values: List[Optional[float]]) -> List[Optional[float]]:
    """
    输入：p 值列表（float 或 None）；输出：同长度的 q 值（None 表示该位无 p）
    算法：Benjamini–Hochberg（step-up），再做单调性修正
    """
    # 收集有效 p 的索引
    idx_valid = [i for i,p in enumerate(values) if (p is not None and p>=0.0)]
    m = len(idx_valid)
    if m == 0:
        return [None for _ in values]
    # 按 p 升序
    idx_sorted = sorted(idx_valid, key=lambda i: values[i])
    qtemp = [None for _ in values]
    min_q = 1.0
    # 逆序做单调性修正：q(i) = min( min_q, p(i)*m/i_rank )
    for rank_from_1, i in reversed(list(enumerate(idx_sorted, start=1))):
        p = float(values[i])
        q = p * m / rank_from_1
        if q < min_q:
            min_q = q
        qtemp[i] = min_q if min_q < 1.0 else 1.0
    return qtemp

def main():
    cfg = load_config()
    paths = cfg["paths"]
    report_cfg = cfg.get("report", {}) or {}
    fdr_alpha = float(report_cfg.get("fdr_alpha", 0.05))  # config.yaml 已提供该键
    logs_dir = Path(paths["logs_dir"]); LOG_FILE = logs_dir / "10_cafe5_aggregate.log"
    log = get_logger("aphylo.10", LOG_FILE)
    banner(log, "APhylo 10 — CAFE5 汇总")

    models = need_dir(Path(paths["cafe_run_dir"])/"models", "CAFE 模型目录")
    out_dir = ensure_dir(Path(paths["cafe_agg_dir"]))

    # 收集所有（model, family, p, source）
    rows: List[Tuple[str,str,Optional[float],str,str]] = []  # (model, family, p_float_or_None, p_str, source)
    lines_branch=["model\tbranch\tlambda_estimate\tsource"]

    for mdir in sorted(models.glob("*")):
        model = mdir.name

        # 1) 先尝试从标准 TSV 抓取
        fam_files = list(mdir.glob("*family*result*.tsv")) + list(mdir.glob("*famil*signif*.tsv"))
        parsed=False
        for ff in fam_files:
            try:
                hdr, rows_tsv = parse_family_table(ff)
                # 列定位：family 一般在第 1 列；pvalue 列名通常为 "pvalue"
                idx_family = 0
                idx_p = None
                if hdr:
                    cols = hdr.split("\t")
                    if "pvalue" in cols:
                        idx_p = cols.index("pvalue")
                    elif "p" in cols:
                        idx_p = cols.index("p")
                for r in rows_tsv:
                    fam = r[idx_family] if idx_family < len(r) else ""
                    p_str = r[idx_p] if (idx_p is not None and idx_p < len(r)) else ""
                    p_val = None
                    try:
                        p_val = float(p_str) if p_str != "" else None
                    except ValueError:
                        p_val = None
                    rows.append((model, fam, p_val, p_str, ff.name))
                parsed=True
            except Exception:
                pass

        # 2) 若没有 TSV，就从 run.log 的文本里兜底提取
        if not parsed:
            runlog = mdir/"run.log"
            if runlog.is_file():
                for line in runlog.read_text(encoding="utf-8").splitlines():
                    m = re.search(r"family\s+(\S+).+p[- ]?value\s*[:=]\s*([\d\.Ee-]+)", line, re.I)
                    if m:
                        fam = m.group(1)
                        p_str = m.group(2)
                        try:
                            p_val = float(p_str)
                        except ValueError:
                            p_val = None
                        rows.append((model, fam, p_val, p_str, "run.log"))
                    b = re.search(r"branch[:\s]+(\S+).+lambda\s*[:=]\s*([\d\.Ee-]+)", line, re.I)
                    if b:
                        lines_branch.append("\t".join([model, b.group(1), b.group(2), "run.log"]))

    # 计算 BH-FDR（对所有可解析 p 统一校正；不改变原有行数/顺序）
    p_list = [r[2] for r in rows]  # Optional[float]
    q_list = bh_fdr(p_list)

    # 写 families 表：在原三列基础上新增 qvalue 与 signif（q<=alpha）
    out_lines = ["model\tfamily\tpvalue\tqvalue\tsignif\tsource"]
    for (model, fam, p_val, p_str, src), qv in zip(rows, q_list):
        q_str = f"{qv:.3g}" if isinstance(qv, float) else ""
        signif = ""
        if isinstance(qv, float):
            signif = "yes" if qv <= fdr_alpha else "no"
        out_lines.append("\t".join([model, fam, p_str, q_str, signif, src]))

    (out_dir/"cafe_significant_families.tsv").write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    (out_dir/"cafe_branch_summary.tsv").write_text("\n".join(lines_branch) + "\n", encoding="utf-8")

    write_done(Path(paths["cafe_agg_dir"])/".agg.done")
    log.info(f"[DONE] CAFE5 汇总完成：families={len(rows)}，FDR α={fdr_alpha}")

if __name__ == "__main__":
    try: main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)

