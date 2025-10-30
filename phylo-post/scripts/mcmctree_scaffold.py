#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, shutil
from pathlib import Path
import yaml
from textwrap import dedent

def need(d,k):
    cur=d
    for p in k.split("."):
        if not isinstance(cur,dict) or p not in cur:
            sys.exit(f"[ERR] config 缺少 post.{k}")
        cur=cur[p]
    return cur

CFG_PATH=str(Path(__file__).resolve().parent.parent / "config.yaml")  # 改
if not os.path.exists(CFG_PATH): sys.exit(f"[ERR] 缺少配置: {CFG_PATH}")
CFG=yaml.safe_load(open(CFG_PATH,encoding="utf-8")) or {}
POST=CFG.get("post",{}) or {}
MCMC=need(POST,"mcmctree")

OUTDIR = need(POST,"outdir")
TREE_IN= os.path.join(OUTDIR,"tree","supermatrix.rooted.contree")
CALIB  = os.path.normpath(need(MCMC,"calibrations"))  # 改
CHECK  = os.path.join(OUTDIR,"mcmctree","calibrations.check.tsv")
DEST_DIR = os.path.join(OUTDIR,"mcmctree","scaffold")

for p in [TREE_IN, CALIB, CHECK]:
    if not os.path.exists(p): sys.exit(f"[ERR] 缺少输入: {p}")

Path(DEST_DIR).mkdir(parents=True, exist_ok=True)
shutil.copyfile(TREE_IN, os.path.join(DEST_DIR,"tree.nwk"))
shutil.copyfile(CALIB,   os.path.join(DEST_DIR,"calibrations.tsv"))

ctl = dedent(f"""
    seed = -1
    seqfile = NULL
    treefile = tree.nwk
    outfile = out.txt

    model = 7
    clock = {MCMC.get('clock')}
    RootAge = {MCMC.get('root_age_Ma')}

    BDparas = {MCMC.get('bd_parameters')[0]} {MCMC.get('bd_parameters')[1]} {MCMC.get('bd_parameters')[2]}
    rgene_gamma = {MCMC.get('rgene_gamma')[0]} {MCMC.get('rgene_gamma')[1]}
    sigma2_gamma = {MCMC.get('sigma2_gamma')[0]} {MCMC.get('sigma2_gamma')[1]}
    ncatG = {MCMC.get('ncatG')}

    cleandata = 1
    print = 1

    burnin = {MCMC.get('burnin')}
    sampfreq = {MCMC.get('sampfreq')}
    nsample = {MCMC.get('nsample')}
""").strip() + "\n"

Path(os.path.join(DEST_DIR,"mcmctree.ctl")).write_text(ctl, encoding="utf-8")

readme = dedent("""
    MCMCTree scaffold (auto-generated)
    ----------------------------------
    * 本目录包含：tree.nwk（根化树）、calibrations.tsv（校准表）、mcmctree.ctl（控制模板）
    * 当前流程只生成脚手架，不实际运行 MCMC。
    * 如需运行，请在合适位置准备 approximate likelihood 文件（in.BV）并修改 ctl，再手动运行：
        mcmctree mcmctree.ctl
    * calibrations.check.tsv 位于 ../ 目录，请先确认其中 status=OK。
""").strip()+"\n"
Path(os.path.join(DEST_DIR,"README.txt")).write_text(readme, encoding="utf-8")

sys.stderr.write(f"[DONE] MCMCTree 脚手架已生成：{DEST_DIR}\n")

