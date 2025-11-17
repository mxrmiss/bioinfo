#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_deg_module.py — 差异分析外壳（调用 06_f_deg.R；仅表，无绘图）
"""

from __future__ import annotations
import sys, json, time, subprocess, tempfile, os
from pathlib import Path
from typing import Dict, Any, Optional
import yaml
import pandas as pd

LOCAL_CONFIG = {
  "config_yaml": "config.yaml",
  "paths": { "matrix_dir":"results/matrix", "deg_dir":"results/deg", "logs_dir":"logs" },
  "tables": { "contrasts":"data/contrasts.tsv", "samples":"data/samples.tsv" },
  "binaries": { "rscript":"Rscript" },
  "deg": { "alpha":0.05, "lfc":1.0, "transform":"vst" }
}
SCRIPT_TAG="06_deg"

def cwd(): return Path.cwd().resolve()
def to_abs(p): 
    if p is None or str(p).strip()=="": return None
    q=Path(p); return (q if q.is_absolute() else cwd()/q).resolve()
def load_yaml(path): 
    p=to_abs(path)
    if not p or not p.exists(): return {}
    return yaml.safe_load(open(p,'r',encoding='utf-8')) or {}
def merge_params(a,b):
    if not isinstance(a,dict): a={}
    if not isinstance(b,dict): b={}
    def is_set(v): return v is not None and not (isinstance(v,(str,list,dict)) and len(v)==0) and not (isinstance(v,str) and v.strip()=="")
    def merge(x,y):
        if isinstance(x,dict) and isinstance(y,dict):
            out={}
            for k in set(x)|set(y):
                xv=x.get(k); yv=y.get(k)
                if isinstance(xv,dict) or isinstance(yv,dict):
                    out[k]=merge(xv if isinstance(xv,dict) else {}, yv if isinstance(yv,dict) else {})
                else:
                    out[k]= yv if is_set(yv) else xv
            return out
        return y if is_set(y) else x
    return merge(a,b)
def err(msg): print(f"[ERR] {msg}", file=sys.stderr); sys.exit(1)
def ensure_dir(p): p=to_abs(p); p.mkdir(parents=True, exist_ok=True); return p
def require_exists(p, kind="file"):
    p=to_abs(p); 
    if p is None: err("必需路径为空")
    if kind=="file" and not p.is_file(): err(f"缺少文件：{p}")
    if kind=="dir" and not p.is_dir(): err(f"缺少目录：{p}")
    return p
def log_open(d, tag): d=ensure_dir(d); return open(d/f"{tag}.log","a",encoding="utf-8")
def log(fp, msg): ts=time.strftime("%Y-%m-%d %H:%M:%S"); fp.write(f"[{ts}] {msg}\n"); fp.flush(); print(msg, flush=True)
def write_params_snapshot(params, logs_dir, tag):
    d=ensure_dir(logs_dir); out=d/f"{tag}.params.tsv"; rows=[]
    def flat(prefix,obj):
        if isinstance(obj,dict):
            for k,v in obj.items(): flat(f"{prefix}.{k}" if prefix else k, v)
        else: rows.append((prefix, json.dumps(obj, ensure_ascii=False)))
    flat("",params); pd.DataFrame(rows, columns=["key","value"]).to_csv(out, sep="\t", index=False)
def run_cmd(cmd, fp, cwd_path=None, env=None):
    log(fp,"CMD: "+" ".join(cmd)); r=subprocess.run(cmd, cwd=to_abs(cwd_path) if cwd_path else None, env=env)
    if r.returncode!=0: err(f"外部命令失败，退出码={r.returncode}")

def main():
    cfg = merge_params(LOCAL_CONFIG, load_yaml(LOCAL_CONFIG["config_yaml"]))
    for k in ("matrix_dir","deg_dir"): ensure_dir(cfg["paths"][k])
    for k in ("contrasts","samples"): require_exists(cfg["tables"][k], "file")
    logs_dir = cfg["paths"].get("logs_dir","logs")
    fp = log_open(logs_dir, SCRIPT_TAG); write_params_snapshot(cfg, logs_dir, SCRIPT_TAG)

    tmp = Path(tempfile.gettempdir())/f"deg_cfg.json"
    tmp.write_text(json.dumps(cfg, ensure_ascii=False), encoding="utf-8")
    env = dict(**os.environ, RNASEQ_CONFIG_JSON=str(tmp))
    r_script = str(Path(__file__).parent/"06_f_deg.R")
    run_cmd([cfg["binaries"]["rscript"], r_script], fp, env=env)
    log(fp, "DE 完成（仅表）")

if __name__=="__main__":
    main()

