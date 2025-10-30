#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys
from pathlib import Path
import yaml

CFG_PATH=str(Path(__file__).resolve().parent.parent / "config.yaml")  # 改
if not os.path.exists(CFG_PATH): sys.exit(f"[ERR] 缺少配置: {CFG_PATH}")
CFG=yaml.safe_load(open(CFG_PATH,encoding="utf-8")) or {}
POST=CFG.get("post",{}) or {}

def need(d,k):
    cur=d
    for p in k.split("."):
        if not isinstance(cur,dict) or p not in cur:
            sys.exit(f"[ERR] config 缺少 post.{k}")
        cur=cur[p]
    return cur

OUTDIR = need(POST,"outdir")
CODON_DIR = os.path.join(OUTDIR,"codon_sco")
TARGETS   = os.path.join(OUTDIR,"targets","og_targets.list")
TREE_IN   = os.path.join(OUTDIR,"tree","supermatrix.rooted.contree")
TREES_OUT = os.path.join(OUTDIR,"codeml","trees")
MODE = need(POST,"codeml.mode")
FGS  = need(POST,"codeml.foreground")

for p in [CODON_DIR, TARGETS, TREE_IN]:
    if not os.path.exists(p): sys.exit(f"[ERR] 缺少输入: {p}")
Path(TREES_OUT).mkdir(parents=True, exist_ok=True)

try:
    from newick import loads, Node
except Exception:
    sys.exit("[ERR] 需要 pip 安装 newick：pip install newick")

txt = Path(TREE_IN).read_text(encoding="utf-8").strip()
big = loads(txt)[0]

def tip_names(node):
    return [n.name for n in node.get_terminals()]

def all_tips_present(node, fg_list):
    names=set(tip_names(node))
    return all(x in names for x in fg_list)

def mark_stem(node, fg_list):
    s = node.newick
    best=None; best_depth=10**9
    for n in node.walk():
        if n.is_leaf: continue
        names=set(tip_names(n))
        if all(x in names for x in fg_list):
            d=0; cur=n
            while cur.parent is not None:
                d+=1; cur=cur.parent
            if d<best_depth:
                best=n; best_depth=d
    if best is None:
        return None
    sub = best.newick
    marked = s.replace(sub, sub+"#1", 1)
    return marked

def mark_terminal(node, fg_list):
    s = node.newick
    for fg in fg_list:
        s = s.replace(f"{fg}:", f"{fg}#1:", 1)
    return s

ok=0; skip=0; tot=0
for og in [ln.strip() for ln in open(TARGETS,encoding="utf-8") if ln.strip()]:
    tot+=1
    codon = Path(CODON_DIR)/f"{og}.codon.fasta"
    if not codon.exists(): 
        skip+=1; continue
    sub = big.clone()
    if MODE=="stem":
        if not all_tips_present(sub, FGS):
            skip+=1; sys.stderr.write(f"[INFO] OG {og} 缺少前景物种，跳过标记\n"); continue
        newick = mark_stem(sub, FGS)
    else:
        newick = mark_terminal(sub, FGS)
    if not newick:
        skip+=1; continue
    Path(TREES_OUT,f"{og}.nwk").write_text(newick+"\n",encoding="utf-8")
    ok+=1

sys.stderr.write(f"[DONE] 剪枝/标记完成：成功 {ok}，跳过 {skip}，总计 {tot} → {TREES_OUT}\n")

