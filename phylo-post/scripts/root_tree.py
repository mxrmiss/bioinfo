#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys
from pathlib import Path
import yaml

# 改：本地 config（scripts/..）
CFG_PATH = str(Path(__file__).resolve().parent.parent / "config.yaml")
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

TREE_ML = os.path.normpath(need(POST,"inputs.species_tree_ml"))  # 改：不再强拼 ../phylo
OUTGROUP = need(POST,"outgroup")
OUTDIR = need(POST,"outdir")
DEST = os.path.join(OUTDIR,"tree","supermatrix.rooted.contree")
Path(os.path.dirname(DEST)).mkdir(parents=True, exist_ok=True)

txt = Path(TREE_ML).read_text(encoding="utf-8").strip()
if not txt:
    sys.exit(f"[ERR] 空树文件: {TREE_ML}")
Path(DEST).write_text(txt+"\n", encoding="utf-8")
sys.stderr.write(f"[DONE] 输出 rooted tree（文本透传）：{DEST}\n")

