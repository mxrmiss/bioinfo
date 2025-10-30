#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, subprocess
from pathlib import Path
import yaml
from textwrap import dedent

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
CODEML_BIN = need(POST,"codeml.bin")
CODON_DIR = os.path.join(OUTDIR,"codon_sco")
TREES_DIR = os.path.join(OUTDIR,"codeml","trees")
RAW_DIR   = os.path.join(OUTDIR,"codeml","raw")
TARGETS   = os.path.join(OUTDIR,"targets","og_targets.list")
LOG_DIR   = os.path.join(OUTDIR,"codeml","logs")
Path(LOG_DIR).mkdir(parents=True, exist_ok=True)

for p in [CODON_DIR, TREES_DIR, TARGETS]:
    if not os.path.exists(p): sys.exit(f"[ERR] 缺少输入: {p}")

ALT_CTL = dedent("""
      seqfile = seq.phy
      treefile = tree.nwk
      outfile = mlc.txt
      noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
      CodonFreq = 2
      clock = 0
      aaDist = 0
      model = 2
      NSsites = 2
      icode = 0
      fix_kappa = 0
      kappa = 2
      fix_omega = 0
      omega = 1.5
      cleandata = 1
""").strip()

NULL_CTL = ALT_CTL.replace("fix_omega = 0","fix_omega = 1").replace("omega = 1.5","omega = 1")

def fasta_to_phylip(in_fa, out_phy):
    seqs=[]
    name=None; buf=[]
    for ln in open(in_fa,encoding="utf-8"):
        if ln.startswith(">"):
            if name is not None: seqs.append((name,"".join(buf)))
            name=ln[1:].strip().split()[0]; buf=[]
        else:
            buf.append(ln.strip())
    if name is not None: seqs.append((name,"".join(buf)))
    if not seqs: return False
    L=len(seqs[0][1])
    for _,s in seqs:
        if len(s)!=L: sys.exit(f"[ERR] 非同长对齐：{in_fa}")
    with open(out_phy,"w",encoding="utf-8") as w:
        w.write(f"{len(seqs)} {L}\n")
        for n,s in seqs:
            w.write(f"{n:<20} {s}\n")
    return True

ok=0; skip=0; tot=0
for og in [ln.strip() for ln in open(TARGETS,encoding="utf-8") if ln.strip()]:
    tot+=1
    codon = Path(CODON_DIR)/f"{og}.codon.fasta"
    tree  = Path(TREES_DIR)/f"{og}.nwk"
    if not codon.exists() or not tree.exists():
        skip+=1; continue
    for which, ctl_txt in (("alt",ALT_CTL),("null",NULL_CTL)):
        work = Path(RAW_DIR)/og/which
        work.mkdir(parents=True, exist_ok=True)
        if not fasta_to_phylip(codon, work/"seq.phy"):
            sys.stderr.write(f"[WARN] 空序列，跳过 {og}\n"); continue
        (work/"tree.nwk").write_text(tree.read_text(encoding="utf-8"), encoding="utf-8")
        (work/"codeml.ctl").write_text(ctl_txt+"\n",encoding="utf-8")
        try:
            subprocess.run([CODEML_BIN, str(work/"codeml.ctl")], cwd=work, check=True,
                           stdout=open(Path(LOG_DIR)/f"{og}.{which}.log","w"),
                           stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            sys.stderr.write(f"[WARN] codeml 失败: {og}/{which}\n")
            continue
    ok+=1

sys.stderr.write(f"[DONE] codeml 成对运行完成：{ok} 个OG，跳过 {skip}，目标 {tot}\n")

