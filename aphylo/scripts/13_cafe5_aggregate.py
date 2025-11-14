#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_cafe5_aggregate.py —— CAFE5 聚合（适配 12 教程增强版）
功能：
  1) 自动定位每个模型目录中的“primary_global”（作为主结果），并合并 “large”（若存在）；
     - 若存在 “*_e” 误差版本，则并行输出“with_error”结果集；
  2) 从 flags/high_fail_ogs.list 读取高失败率家族，输出：
     - all_families.tsv                 …… 全量  （primary + large）
     - no_highfail.tsv                  …… 剔除 high-fail 家族
     - all_families_with_error.tsv      …… 全量（误差模型）
     - no_highfail_with_error.tsv       …… 剔除 high-fail（误差模型）
  3) 统一做 BH-FDR（q 值），显著阈值 α 读 config.report.fdr_alpha（默认 0.05）。
  4) 产出分支 λ 摘要（从 primary_global 的最后块与/或文件读取），保存：
     - cafe_branch_summary.tsv
     - 若误差模型存在，则另存 cafe_branch_summary_with_error.tsv（若能获取）

输入要求：
  - 12 跑完后存在 .cafe.done；
  - 目录结构见 12 脚本说明（primary_global/、large/ 等）。
"""

from __future__ import annotations
import sys, io, logging, re, os
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import yaml
import datetime

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
    bar = "=" * max(10, len(text)+2); log.info(bar); log.info(f" {text} "); log.info(bar)

# ---------- 解析与工具 ----------

def read_table_auto(path: Path) -> Tuple[List[str], List[List[str]]]:
    txt = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not txt: return [], []
    delim = "\t" if ("\t" in txt[0]) else " "
    head = [c.strip() for c in re.split(r"\s+" if delim==" " else "\t", txt[0].strip())]
    body = []
    for ln in txt[1:]:
        if not ln.strip(): continue
        body.append([c.strip() for c in (re.split(r"\s+", ln.strip()) if delim==" " else ln.split("\t"))])
    return head, body

def find_family_results_file(d: Path) -> Optional[Path]:
    # 优先精确名，其次包含 family & result 的 tsv/txt
    for name in ("Gamma_family_results.txt", "Base_family_results.txt"):
        p = d / name
        if p.is_file(): return p
    cand = [p for p in d.iterdir() if p.is_file() and p.suffix in (".tsv",".txt") and ("family" in p.name.lower() and "result" in p.name.lower())]
    return sorted(cand)[0] if cand else None

def take_last_block(runlog: Path) -> str:
    if not runlog.is_file(): return ""
    text = runlog.read_text(encoding="utf-8", errors="ignore")
    parts = re.split(r"^=+\s*MODEL=.*?ROUND=.*?=+\s*$", text, flags=re.M)
    return (parts[-1].strip() if parts else text)

BRANCH_LAMBDA_PAT = re.compile(r"branch[:\s]+(\S+).+lambda\s*[:=]\s*([\d\.Ee+-]+)", re.I)

def parse_branch_lambda_from_block(block: str) -> List[Tuple[str, str]]:
    out = []
    if not block: return out
    for line in block.splitlines():
        m = BRANCH_LAMBDA_PAT.search(line)
        if m: out.append((m.group(1), m.group(2)))
    return out

def bh_fdr(values: List[Optional[float]]) -> List[Optional[float]]:
    idx = [i for i,v in enumerate(values) if (v is not None and v>=0.0)]
    m = len(idx)
    if m == 0: return [None]*len(values)
    idx_sorted = sorted(idx, key=lambda i: values[i])
    qtemp = [None]*len(values); min_q=1.0
    for rank1, i in reversed(list(enumerate(idx_sorted, start=1))):
        p = float(values[i]); q=p*m/rank1
        if q < min_q: min_q = q
        qtemp[i] = min_q if min_q<1.0 else 1.0
    return qtemp

# ---------- 主流程 ----------

def main():
    cfg = load_config()
    paths = cfg["paths"]; report_cfg = cfg.get("report", {}) or {}
    fdr_alpha = float(report_cfg.get("fdr_alpha", 0.05))

    logs_dir = Path(paths["logs_dir"])
    log = get_logger("aphylo.13", logs_dir / "13_cafe5_aggregate.log")
    banner(log, "APhylo 13 — CAFE5 聚合（large并回 / 高失败敏感性 / 误差并行）")

    cafe_dir = Path(paths["cafe_run_dir"])
    agg_dir  = ensure_dir(Path(paths["cafe_agg_dir"]))
    doneflag = cafe_dir / ".cafe.done"
    if not doneflag.is_file():
        log.error(f"[ERR] 未发现 {doneflag} —— 12 未成功完成"); sys.exit(2)

    models_root = cafe_dir / "models"
    if not models_root.is_dir():
        log.error(f"[ERR] 缺少模型目录：{models_root}"); sys.exit(2)

    # 聚合容器
    fam_rows_all:   List[Tuple[str,str,Optional[float],str,str,str]] = []  # (model, subset, family, p, pstr, source)
    fam_rows_all_e: List[Tuple[str,str,Optional[float],str,str,str]] = []  # with error
    branch_rows:    List[Tuple[str,str,str,str]] = []                      # (model, branch, lambda, source)
    branch_rows_e:  List[Tuple[str,str,str,str]] = []                      # with error (若能抓到)

    inputs_used = ["model\tdir\tfile\tmtime\tbytes\trole"]

    for mdir in sorted(models_root.glob("*")):
        if not mdir.is_dir(): continue
        model = mdir.name

        # high-fail 清单
        highfail = []
        hf_file = mdir / "flags" / "high_fail_ogs.list"
        if hf_file.is_file():
            highfail = [x.strip() for x in hf_file.read_text(encoding="utf-8").splitlines() if x.strip()]
            log.info(f"[{model}] high-fail 家族：{len(highfail)}")

        # 主集合：primary_global
        pg = mdir / "primary_global"
        if not pg.is_dir():
            log.error(f"[ERR] {model} 缺少 primary_global/"); sys.exit(2)
        fr_pg = find_family_results_file(pg)
        if not fr_pg:
            log.error(f"[ERR] {model} primary_global/ 未找到 family_results 文件"); sys.exit(2)
        hdr, body = read_table_auto(fr_pg)
        # 列定位
        name_lc = [c.lower() for c in hdr]
        idx_f = None; idx_p = None
        for i,c in enumerate(name_lc):
            if c in ("family","orthogroup","og"): idx_f = i
            if c in ("p","pvalue"): idx_p = i
        if idx_f is None or idx_p is None:
            log.error(f"[ERR] {model} 无法定位 family/p 列：{fr_pg}"); sys.exit(2)
        # 收集 primary_global
        for r in body:
            fam = r[idx_f] if idx_f < len(r) else ""
            pstr= r[idx_p] if idx_p < len(r) else ""
            pval= None
            try:
                if pstr!="": pval=float(pstr)
            except Exception: pval=None
            fam_rows_all.append((model, "primary", fam, pval, pstr, fr_pg.name))
        st = fr_pg.stat()
        inputs_used.append("\t".join([model, "primary_global", fr_pg.name,
                                      datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                                      str(st.st_size), "family_results"]))

        # 分支 λ（从 run.log 最后块抓）
        blk = take_last_block(pg / "run.log")
        for br, lam in parse_branch_lambda_from_block(blk):
            branch_rows.append((model, br, lam, "primary_global.runlog"))

        # large 并回（若存在）
        lg = mdir / "large"
        if lg.is_dir():
            fr_lg = find_family_results_file(lg)
            if fr_lg:
                h2,b2 = read_table_auto(fr_lg)
                nl2 = [c.lower() for c in h2]
                jf = idx_f if "family" in nl2 or "orthogroup" in nl2 or "og" in nl2 else None
                jp = idx_p if "p" in nl2 or "pvalue" in nl2 else None
                if jf is None or jp is None:
                    # 再次定位
                    jf = next((i for i,c in enumerate(nl2) if c in ("family","orthogroup","og")), None)
                    jp = next((i for i,c in enumerate(nl2) if c in ("p","pvalue")), None)
                if jf is None or jp is None:
                    log.warning(f"[WARN] {model} large/ 结果无法定位 family/p 列：{fr_lg}")
                else:
                    for r in b2:
                        fam = r[jf] if jf < len(r) else ""
                        pstr= r[jp] if jp < len(r) else ""
                        pval= None
                        try:
                            if pstr!="": pval=float(pstr)
                        except Exception: pval=None
                        fam_rows_all.append((model, "large", fam, pval, pstr, fr_lg.name))
                    st = fr_lg.stat()
                    inputs_used.append("\t".join([model, "large", fr_lg.name,
                                                  datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                                                  str(st.st_size), "family_results"]))

        # 误差模型并行集（若存在）
        pge = mdir / "primary_global_e"
        if pge.is_dir():
            fr_pge = find_family_results_file(pge)
            if fr_pge:
                h3,b3 = read_table_auto(fr_pge)
                n3 = [c.lower() for c in h3]
                jf = next((i for i,c in enumerate(n3) if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(n3) if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in b3:
                        fam = r[jf] if jf < len(r) else ""
                        pstr= r[jp] if jp < len(r) else ""
                        pval= None
                        try:
                            if pstr!="": pval=float(pstr)
                        except Exception: pval=None
                        fam_rows_all_e.append((model, "primary", fam, pval, pstr, fr_pge.name))
                    st = fr_pge.stat()
                    inputs_used.append("\t".join([model, "primary_global_e", fr_pge.name,
                                                  datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                                                  str(st.st_size), "family_results(e)"]))
            # 分支 λ（如果能抓到，以 _e 标记）
            blk_e = take_last_block(pge / "run.log")
            for br, lam in parse_branch_lambda_from_block(blk_e):
                branch_rows_e.append((model, br, lam, "primary_global_e.runlog"))

        lge = mdir / "large_e"
        if lge.is_dir():
            fr_lge = find_family_results_file(lge)
            if fr_lge:
                h4,b4 = read_table_auto(fr_lge)
                n4 = [c.lower() for c in h4]
                jf = next((i for i,c in enumerate(n4) if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(n4) if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in b4:
                        fam = r[jf] if jf < len(r) else ""
                        pstr= r[jp] if jp < len(r) else ""
                        pval= None
                        try:
                            if pstr!="": pval=float(pstr)
                        except Exception: pval=None
                        fam_rows_all_e.append((model, "large", fam, pval, pstr, fr_lge.name))
                    st = fr_lge.stat()
                    inputs_used.append("\t".join([model, "large_e", fr_lge.name,
                                                  datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                                                  str(st.st_size), "family_results(e)"]))

    # ---------- 输出 family 聚合（全量 / 剔除 high-fail） ----------
    def write_family_outputs(rows: List[Tuple[str,str,Optional[float],str,str,str]], suffix: str, highfail_all: Dict[str, set]):
        if not rows:
            return
        # rows: (model, subset, family, pval, pstr, source)
        p_list = [r[3] for r in rows]
        q_list = bh_fdr(p_list)
        out_all = ["model\tsubset\tfamily\tpvalue\tqvalue\tsignif\tsource"]
        out_noh= ["model\tsubset\tfamily\tpvalue\tqvalue\tsignif\tsource"]
        for (model, subset, fam, pval, pstr, src), qv in zip(rows, q_list):
            qstr = f"{qv:.3g}" if isinstance(qv, float) else ""
            signif = "yes" if (isinstance(qv, float) and qv <= fdr_alpha) else ("no" if isinstance(qv, float) else "")
            line = "\t".join([model, subset, fam, (pstr or ""), qstr, signif, src])
            out_all.append(line)
            if fam not in highfail_all.get(model, set()):
                out_noh.append(line)
        (agg_dir/f"cafe_significant_families{suffix}.tsv").write_text("\n".join(out_all)+"\n", encoding="utf-8")
        (agg_dir/f"cafe_significant_families_no_highfail{suffix}.tsv").write_text("\n".join(out_noh)+"\n", encoding="utf-8")

    # 汇总每模型的 high-fail 集合
    highfail_map: Dict[str,set] = {}
    for mdir in sorted(models_root.glob("*")):
        if not mdir.is_dir(): continue
        model = mdir.name
        hf_file = mdir / "flags" / "high_fail_ogs.list"
        if hf_file.is_file():
            s=set(x.strip() for x in hf_file.read_text(encoding="utf-8").splitlines() if x.strip())
            highfail_map[model]=s

    write_family_outputs(fam_rows_all, "", highfail_map)
    write_family_outputs(fam_rows_all_e, "_with_error", highfail_map)

    # ---------- 输出分支 λ 摘要 ----------
    def write_branch(rows: List[Tuple[str,str,str,str]], suffix: str):
        if not rows: return
        lines = ["model\tbranch\tlambda_estimate\tsource"]
        for model, br, lam, src in rows:
            lines.append("\t".join([model, br, lam, src]))
        (agg_dir/f"cafe_branch_summary{suffix}.tsv").write_text("\n".join(lines)+"\n", encoding="utf-8")

    write_branch(branch_rows, "")
    write_branch(branch_rows_e, "_with_error")

    # ---------- 记录 inputs_used ----------
    (agg_dir/"inputs_used.tsv").write_text("\n".join(inputs_used)+"\n", encoding="utf-8")

    print(f"[OK] 写出：{agg_dir/'cafe_significant_families.tsv'}")
    print(f"[OK] 写出：{agg_dir/'cafe_significant_families_no_highfail.tsv'}")
    if fam_rows_all_e:
        print(f"[OK] 写出：{agg_dir/'cafe_significant_families_with_error.tsv'}")
        print(f"[OK] 写出：{agg_dir/'cafe_significant_families_no_highfail_with_error.tsv'}")
    print(f"[OK] 写出：{agg_dir/'cafe_branch_summary.tsv'}")
    if branch_rows_e:
        print(f"[OK] 写出：{agg_dir/'cafe_branch_summary_with_error.tsv'}")
    print(f"[OK] 写出：{agg_dir/'inputs_used.tsv'}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)