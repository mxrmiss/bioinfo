#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, re, math
from pathlib import Path
import yaml

try:
    from statistics import median
except Exception:
    median = None

def bh_fdr(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0]*m
    prev = 1.0
    for rank,i in enumerate(reversed(order), start=1):
        pi = pvals[i]
        val = min(prev, pi * m / (m-rank+1))
        q[i] = val
        prev = val
    return q

def chi2_sf(x, k=1):
    z = math.sqrt(max(0.0, x))
    return math.erfc(z / math.sqrt(2.0))

CFG_PATH=str(Path(__file__).resolve().parent.parent / "config.yaml")  # 改
if not os.path.exists(CFG_PATH): sys.exit(f"[ERR] 缺少配置: {CFG_PATH}")
CFG=yaml.safe_load(open(CFG_PATH,encoding="utf-8")) or {}
POST=CFG.get("post",{}) or {}

OUTDIR = POST.get("outdir")
RAW_DIR = os.path.join(OUTDIR,"codeml","raw")
R_FDR  = os.path.join(OUTDIR,"codeml","D_fdr_genes.tsv")
R_BEB  = os.path.join(OUTDIR,"codeml","D_beb_sites.tsv")
if not os.path.isdir(RAW_DIR): sys.exit(f"[ERR] 缺少目录: {RAW_DIR}")

def parse_lnl(path):
    if not os.path.exists(path): return None
    txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"lnL\s*=\s*([-\d\.E+]+)", txt)
    return float(m.group(1)) if m else None

def parse_beb(path):
    if not os.path.exists(path): return []
    txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    if "Bayes Empirical Bayes" not in txt: return []
    lines = txt.splitlines()
    sites=[]; capture=False
    for ln in lines:
        if "Bayes Empirical Bayes" in ln:
            capture=True; continue
        if capture:
            ln=ln.strip()
            if not ln: continue
            parts = ln.split()
            if len(parts)>=3:
                try:
                    site = int(parts[0])
                    post = float(parts[-1]) if parts[-1].replace(".","",1).isdigit() else float(parts[-2])
                    if post>=0.95:
                        aa = parts[1]
                        sites.append((site, aa, post))
                except Exception:
                    pass
    return sites

rows=[]; pvals=[]; ogs=[]
for og_dir in Path(RAW_DIR).iterdir():
    if not og_dir.is_dir(): continue
    alt = og_dir/"alt"/"mlc.txt"
    null= og_dir/"null"/"mlc.txt"
    lnl_alt = parse_lnl(alt)
    lnl_null= parse_lnl(null)
    if lnl_alt is None or lnl_null is None:
        continue
    x = 2.0*(lnl_alt - lnl_null)
    p = 0.5 * chi2_sf(max(0.0,x), 1)
    rows.append((og_dir.name, x, p))
    pvals.append(p); ogs.append(og_dir.name)

if not rows:
    open(R_FDR,"w").write("OG\t2DeltaL\tp\tq\n")
    open(R_BEB,"w").write("OG\tsite\taa\tpostP\n")
    sys.exit("[WARN] 未找到可聚合的 codeml 结果")

qvals = bh_fdr(pvals)
with open(R_FDR,"w",encoding="utf-8") as w:
    w.write("OG\t2DeltaL\tp\tq\n")
    for (og, x, p), q in zip(rows, qvals):
        w.write(f"{og}\t{x:.6f}\t{p:.6e}\t{q:.6e}\n")

with open(R_BEB,"w",encoding="utf-8") as w:
    w.write("OG\tsite\taa\tpostP\n")
    for og in ogs:
        mlc = Path(RAW_DIR)/og/"alt"/"mlc.txt"
        for site, aa, post in parse_beb(mlc):
            w.write(f"{og}\t{site}\t{aa}\t{post:.3f}\n")

sys.stderr.write(f"[DONE] 聚合完成 → {R_FDR}, {R_BEB}\n")

