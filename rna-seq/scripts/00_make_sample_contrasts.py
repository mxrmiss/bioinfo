#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_make_sample_contrasts.py — 可选：扫描 FASTQ 生成 samples.tsv / contrasts.tsv（可人工再修）
运行：python scripts/00_make_sample_contrasts.py
"""

from __future__ import annotations
import re, sys, json, time, subprocess
from pathlib import Path
from typing import Dict, Any, Optional
import yaml
import pandas as pd

# ============== 配置 ==============
LOCAL_CONFIG = {
  "config_yaml": "config.yaml",  # 默认读取根配置
  "paths": { "fastq_dir": "data/fastq", "logs_dir": "logs" },
  "tables": { "samples": "data/samples.tsv", "contrasts": "data/contrasts.tsv" }
}
SCRIPT_TAG="00_make_sample_contrasts"

# ============== 工具函数（YAML 优先合并） ==============
def cwd() -> Path: return Path.cwd().resolve()
def to_abs(p: str | Path) -> Path:
    if p is None or str(p).strip()=="": return None
    q=Path(p); return (q if q.is_absolute() else cwd()/q).resolve()
def load_yaml(path: str) -> Dict[str, Any]:
    p=to_abs(path); 
    if not p or not p.exists(): return {}
    with p.open('r', encoding='utf-8') as f: return yaml.safe_load(f) or {}
def merge_params(local_cfg: dict, yaml_cfg: dict) -> dict:
    if not isinstance(local_cfg, dict): local_cfg={}
    if not isinstance(yaml_cfg, dict):  yaml_cfg={}
    def is_set(v):
        return v is not None and not (isinstance(v,(str,list,dict)) and len(v)==0) and not (isinstance(v,str) and v.strip()=="")
    def merge(a,b):
        if isinstance(a,dict) and isinstance(b,dict):
            out={}
            for k in set(a.keys())|set(b.keys()):
                av=a.get(k); bv=b.get(k)
                if isinstance(av,dict) or isinstance(bv,dict):
                    out[k]=merge(av if isinstance(av,dict) else {}, bv if isinstance(bv,dict) else {})
                else:
                    out[k]=bv if is_set(bv) else av
            return out
        return b if is_set(b) else a
    return merge(local_cfg, yaml_cfg)
def err(msg:str): print(f"[ERR] {msg}", file=sys.stderr); sys.exit(1)
def ensure_dir(p): p=to_abs(p); p.mkdir(parents=True, exist_ok=True); return p
def require_exists(p, kind="file"):
    p=to_abs(p); 
    if p is None: err("必需路径为空")
    if kind=="file" and not p.is_file(): err(f"缺少文件：{p}")
    if kind=="dir" and not p.is_dir(): err(f"缺少目录：{p}")
    return p
def log_open(d, tag): d=ensure_dir(d); return open(d/f"{tag}.log","a",encoding="utf-8")
def log(fp, msg):
    ts=time.strftime("%Y-%m-%d %H:%M:%S"); fp.write(f"[{ts}] {msg}\n"); fp.flush(); print(msg, flush=True)
def write_params_snapshot(params, logs_dir, tag):
    d=ensure_dir(logs_dir); out=d/f"{tag}.params.tsv"
    rows=[]
    def flat(prefix,obj):
        if isinstance(obj,dict):
            for k,v in obj.items(): flat(f"{prefix}.{k}" if prefix else k, v)
        else: rows.append((prefix, json.dumps(obj, ensure_ascii=False)))
    flat("",params); pd.DataFrame(rows,columns=["key","value"]).to_csv(out,sep="\t",index=False)

# ============== 主逻辑 ==============
def main():
    cfg = merge_params(LOCAL_CONFIG, load_yaml(LOCAL_CONFIG["config_yaml"]))
    fastq_dir = require_exists(cfg["paths"]["fastq_dir"], "dir")
    samples_tsv = Path(cfg["tables"]["samples"])
    contrasts_tsv = Path(cfg["tables"]["contrasts"])
    logs_dir = cfg["paths"].get("logs_dir","logs")
    fp = log_open(logs_dir, SCRIPT_TAG); write_params_snapshot(cfg, logs_dir, SCRIPT_TAG)

    fq = sorted(fastq_dir.glob("*.f*q.gz"))
    if not fq: err("fastq_dir 下未发现 *.fq.gz/*.fastq.gz")

    def sample_name(p: Path):
        n=p.name
        n=re.sub(r"(_R[12]|_[12])(?=\.)","",n)
        n=re.sub(r"\.f(ast)?q\.gz$","",n)
        return n

    samples={}
    for p in fq:
        s=sample_name(p)
        mate = "R1" if re.search(r"_R?1\.", p.name) else ("R2" if re.search(r"_R?2\.", p.name) else "SE")
        d=samples.setdefault(s,{"sample":s,"group":"GroupA","fastq1":"","fastq2":""})
        if mate=="R1" or mate=="SE": d["fastq1"]=str(p)
        else: d["fastq2"]=str(p)

    df=pd.DataFrame(samples.values()).sort_values("sample")
    ensure_dir(samples_tsv.parent); df.to_csv(samples_tsv, sep="\t", index=False); log(fp,f"写出 {samples_tsv}")
    ct=pd.DataFrame([{"contrast":"GroupA_vs_GroupB","case":"GroupA","control":"GroupB"}])
    ensure_dir(contrasts_tsv.parent); ct.to_csv(contrasts_tsv, sep="\t", index=False); log(fp,f"写出 {contrasts_tsv}")

if __name__=="__main__":
    main()

