#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys
from pathlib import Path
import yaml

def need(cfg,key):
    cur=cfg
    for k in key.split("."):
        if not isinstance(cur,dict) or k not in cur:
            sys.exit(f"[ERR] config 缺少 post.{key}")
        cur=cur[k]
    return cur

CFG_PATH=str(Path(__file__).resolve().parent.parent / "config.yaml")  # 改
if not os.path.exists(CFG_PATH): sys.exit(f"[ERR] 缺少配置: {CFG_PATH}")
CFG=yaml.safe_load(open(CFG_PATH,encoding="utf-8")) or {}
POST=CFG.get("post",{}) or {}

OUTDIR = need(POST,"outdir")
AA_DIR = os.path.normpath(need(POST,"inputs.aa_align_dir"))  # 改
TARGETS = os.path.join(OUTDIR,"targets","og_targets.list")
MAP_TSV = os.path.join(OUTDIR,"codon_sco","og_species_to_gene.tsv")
DEST_DIR= os.path.join(OUTDIR,"prot_msa_geneid")

for p in [AA_DIR, TARGETS, MAP_TSV]:
    if not os.path.exists(p): sys.exit(f"[ERR] 缺少输入: {p}")

from collections import defaultdict
og_sp2gene=defaultdict(dict)
for ln in open(MAP_TSV,encoding="utf-8"):
    if ln.startswith("OG\t"): continue
    og, sp, gid = [x.strip() for x in ln.rstrip("\n").split("\t")]
    og_sp2gene[og][sp]=gid

def relabel_one(og):
    src = Path(AA_DIR)/f"{og}.trim.faa"
    if not src.exists(): return False
    dest = Path(DEST_DIR)/src.name
    dest.parent.mkdir(parents=True, exist_ok=True)
    sp2gene = og_sp2gene.get(og, {})
    if not sp2gene:
        sys.stderr.write(f"[WARN] OG 无映射，跳过: {og}\n"); return False
    with src.open("r",encoding="utf-8",errors="replace") as r, \
         dest.open("w",encoding="utf-8") as w:
        cur_sp=""
        for line in r:
            if line.startswith(">"):
                sp = line[1:].strip().split()[0]
                gid = sp2gene.get(sp, "")
                if not gid:
                    cur_sp=""; continue
                w.write(f">{gid}\n"); cur_sp=sp
            else:
                if cur_sp: w.write(line)
    return True

ok=0; tot=0
for og in [ln.strip() for ln in open(TARGETS,encoding="utf-8") if ln.strip()]:
    tot+=1
    if relabel_one(og): ok+=1
sys.stderr.write(f"[DONE] 重标完成：{ok}/{tot} 写入 {DEST_DIR}\n")

