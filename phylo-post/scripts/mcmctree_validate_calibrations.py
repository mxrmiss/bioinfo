#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, csv, re
from pathlib import Path
import yaml

CFG_PATH = str(Path(__file__).resolve().parent.parent / "config.yaml")  # 改
if not os.path.exists(CFG_PATH):
    sys.exit(f"[ERR] 缺少配置文件: {CFG_PATH}")

CFG = yaml.safe_load(open(CFG_PATH, encoding="utf-8")) or {}
POST = CFG.get("post", {}) or {}

def need(d, key):
    cur = d
    for k in key.split("."):
        if not isinstance(cur, dict) or k not in cur:
            sys.exit(f"[ERR] config.yaml 缺少 post.{key}")
        cur = cur[k]
    return cur

OUTDIR = need(POST, "outdir")
CALIB  = os.path.normpath(need(POST, "mcmctree.calibrations"))  # 改：不再强拼 ../phylo
TREE   = os.path.join(OUTDIR, "tree", "supermatrix.rooted.contree")
DEST   = os.path.join(OUTDIR, "mcmctree", "calibrations.check.tsv")

for p in [CALIB, TREE]:
    if not os.path.exists(p):
        sys.exit(f"[ERR] 缺少输入路径: {p}")

txt = Path(TREE).read_text(encoding="utf-8").strip()
if not txt or "(" not in txt or ")" not in txt:
    sys.exit(f"[ERR] 树文件内容异常：{TREE}")

tips = set(re.findall(r"([A-Za-z0-9_.:-]+):[0-9eE\.\-\+]+", txt))
if not tips:
    tips = set(re.findall(r"[,\(\)]([A-Za-z0-9_.:-]+)[,\)\:]", txt))
if not tips:
    sys.exit("[ERR] 未从树中解析到物种名，请检查 newick 格式/叶名字符")

REQUIRED = {"node_id","clade_A","clade_B","median_Ma","lo_Ma","hi_Ma","source","notes"}

with open(CALIB, "r", encoding="utf-8") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    header = set(rdr.fieldnames or [])
    missing = REQUIRED - header
    if missing:
        sys.exit(f"[ERR] calibrations.tsv 表头不匹配，缺少列：{', '.join(sorted(missing))}")

    ok, bad = 0, 0
    rows_out = []
    for i, rec in enumerate(rdr, start=2):
        node_id  = (rec.get("node_id")  or "").strip()
        a        = (rec.get("clade_A")  or "").strip()
        b        = (rec.get("clade_B")  or "").strip()
        med      = (rec.get("median_Ma") or "").strip()
        lo_s     = (rec.get("lo_Ma")    or "").strip()
        hi_s     = (rec.get("hi_Ma")    or "").strip()
        source   = (rec.get("source")   or "").strip()
        notes    = (rec.get("notes")    or "").strip()

        reasons = []
        if not a:
            reasons.append("clade_A_empty")
        elif a not in tips:
            reasons.append("clade_A_not_in_tree")
        if b:
            if b not in tips:
                reasons.append("clade_B_not_in_tree")
        try:
            lo = float(lo_s); hi = float(hi_s)
            if lo > hi:
                reasons.append("lo_gt_hi")
        except Exception:
            reasons.append("age_non_numeric")

        status = "OK" if not reasons else ("FAIL:" + ",".join(reasons))
        if reasons: bad += 1
        else: ok += 1
        rows_out.append([node_id, a, b, med, lo_s, hi_s, source, notes, status])

Path(os.path.dirname(DEST)).mkdir(parents=True, exist_ok=True)
with open(DEST, "w", encoding="utf-8") as w:
    w.write("\t".join(["node_id","clade_A","clade_B","median_Ma","lo_Ma","hi_Ma","source","notes","status"]) + "\n")
    for r in rows_out:
        w.write("\t".join(r) + "\n")

sys.stderr.write(f"[DONE] 校准校验完成：OK={ok}, FAIL={bad} → {DEST}\n")

