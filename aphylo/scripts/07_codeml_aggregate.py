#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_codeml_aggregate.py — 汇总 codeml（显式 alt/null 配对；输出到 codeml_agg_dir）
修复要点：
1) lnL 解析仅取冒号/等号后的 lnL 数值，并取“最后一次出现”的 lnL（最终优化值）。
2) LRT 负值截断为 0；p = 0.5 * erfc(sqrt(LRT/2))（χ²₁一半混合），数值稳定。
3) BEB 表逐行解析（不跨行），自动剥离星标 *；兼容含/不含“Positively selected sites/BEB”表头。
4) BH-FDR 纯 Python 实现；输入/输出目录与文件名保持不变。
"""

from __future__ import annotations
import sys, io, logging, re, math
from pathlib import Path
from typing import Dict, Any, List, Tuple
import yaml

DEFAULT_CONFIG = "config.yaml"

# ---------- 通用工具 ----------
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

# ---------- 解析函数 ----------
NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"

def last_lnl(txt: str) -> float | None:
    """
    只抓 lnL 行里“冒号/等号后的第一个数值”，并取最后一次出现（最终 lnL）
      例：lnL(ntime: 27  np: 32):  -8539.712388  +0.000000
          lnL = -8539.712388
    """
    # 带括号的主流格式
    m1 = re.findall(rf"^lnL\s*\([^\n]*?\)\s*:\s*({NUM})", txt, flags=re.M)
    # 简写格式
    m2 = re.findall(rf"^lnL\s*=\s*({NUM})", txt, flags=re.M)
    m = (m1 if m1 else []) + (m2 if m2 else [])
    if not m: 
        return None
    try:
        return float(m[-1])
    except Exception:
        return None

def chi2_mix_half_df1_sf(x: float) -> float:
    """ 0.5*χ²₁ 的右尾概率：p = 0.5 * erfc(sqrt(x/2)) """
    if x <= 0: 
        return 0.5
    return 0.5 * math.erfc(math.sqrt(x / 2.0))

FLOAT_RE = re.compile(NUM)

def parse_beb_lines(mlc_text: str, cutoff: float = 0.95) -> List[Tuple[int,str,float]]:
    """
    逐行解析：匹配形如
        <site:int> <AA:单字母> ... <post:float[*]>
    的行；取该行“最后一个数值 token”（去掉尾部星标）作为 posterior。
    不强依赖固定表头；若检测到包含 “Positively selected sites” 或 “BEB” 的行，会在随后的非空行里优先解析。
    """
    rows: List[Tuple[int,str,float]] = []
    in_beb_block = False
    for raw in mlc_text.splitlines():
        line = raw.rstrip()
        # 进入/退出 BEB 区（尽量宽松）
        if re.search(r"Positively selected sites|Bayes Empirical Bayes|BEB", line, flags=re.I):
            in_beb_block = True
            continue
        if in_beb_block and not line.strip():
            in_beb_block = False
            continue

        parts = line.split()
        if len(parts) < 3: 
            continue
        if not parts[0].isdigit(): 
            continue
        if not re.fullmatch(r"[A-Z\*]", parts[1]): 
            continue

        # 从行尾向前找数值 token（去掉星标）
        post_tok = None
        for tok in reversed(parts[2:]):
            t = tok.rstrip("*")
            if FLOAT_RE.fullmatch(t):
                post_tok = t
                break
        if post_tok is None: 
            continue

        try:
            post = float(post_tok)
            if post >= cutoff:
                rows.append((int(parts[0]), parts[1], post))
        except Exception:
            continue
    return rows

def bh_fdr(pvals: List[float]) -> List[float]:
    """ Benjamini–Hochberg（单调回填） """
    m = len(pvals)
    if m == 0: return []
    order = sorted(range(m), key=lambda i: pvals[i])  # 升序
    q = [1.0]*m
    min_q = 1.0
    for rank_from_end, i in enumerate(reversed(order), start=1):
        j = m - rank_from_end + 1
        val = pvals[i] * m / j
        if val < min_q: min_q = val
        q[i] = min_q if min_q < 1.0 else 1.0
    return q

# ---------- 主流程 ----------
def main():
    cfg = load_config()
    paths = cfg["paths"]
    logs_dir = Path(paths["logs_dir"]); LOG_FILE = logs_dir / "07_codeml_aggregate.log"
    log = get_logger("aphylo.07", LOG_FILE)
    banner(log, "APhylo 07 — codeml 汇总")

    raw_root = need_dir(Path(paths["codeml_dir"]) / "raw", "codeml 原始结果目录")
    out_dir  = ensure_dir(Path(paths["codeml_agg_dir"]))

    rows_genes: List[List[str]] = []
    rows_sites: List[List[str]] = []

    for og_dir in sorted(raw_root.glob("OG*")):
        og = og_dir.name
        for fg_dir in sorted(og_dir.glob("*")):
            fg = fg_dir.name
            alt = need_file(fg_dir/"alt"/"mlc.txt",  "ALT mlc")
            nul = need_file(fg_dir/"null"/"mlc.txt", "NULL mlc")

            alt_txt = alt.read_text(encoding="utf-8", errors="ignore")
            nul_txt =  nul.read_text(encoding="utf-8", errors="ignore")

            la = last_lnl(alt_txt)
            ln = last_lnl(nul_txt)
            if la is None or ln is None:
                raise ValueError(f"[ERR] 无法解析 lnL：{alt} 或 {nul}")

            lrt = 2.0 * (la - ln)
            if lrt < 0: lrt = 0.0
            pval = chi2_mix_half_df1_sf(lrt)
            rows_genes.append([og, fg, f"{lrt:.6f}", f"{pval:.6g}"])

            for site, aa, post in parse_beb_lines(alt_txt, cutoff=0.95):
                rows_sites.append([og, fg, str(site), aa, f"{post:.3f}"])

    # 写基因表 + BH-FDR
    head_g = "OG\tforeground\tLRT\tP\tQ\n"
    if rows_genes:
        pvals = [float(r[3]) for r in rows_genes]
        qvals = bh_fdr(pvals)
        lines = []
        for r, q in zip(rows_genes, qvals):
            lines.append("\t".join([r[0], r[1], r[2], r[3], f"{q:.6g}"]))
        (out_dir/"D_fdr_genes.tsv").write_text(head_g + "\n".join(lines) + "\n", encoding="utf-8")
    else:
        (out_dir/"D_fdr_genes.tsv").write_text(head_g, encoding="utf-8")

    # 写 BEB 位点表
    head_s = "OG\tforeground\tsite\taa\tpost\n"
    if rows_sites:
        (out_dir/"D_beb_sites.tsv").write_text(head_s + "\n".join("\t".join(x) for x in rows_sites) + "\n", encoding="utf-8")
    else:
        (out_dir/"D_beb_sites.tsv").write_text(head_s, encoding="utf-8")

    write_done(out_dir/".aggregate.done")
    log.info(f"[DONE] codeml 汇总完成：基因 {len(rows_genes)}；BEB 位点 {len(rows_sites)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)

