#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_qc_fastp_cutadapt.py — 质控与去接头（fastp；无绘图，产原生报告+汇总表）
输出：
  - results/qc/clean/*.clean.fq.gz
  - results/qc/report/fastp_json/<sample>.json
  - results/qc/report/fastp_html/<sample>.html
  - results/qc/report/fastq_stats.tsv / clean_stats.tsv / fastq_pairing.tsv
"""

from __future__ import annotations
import re, sys, json, time, subprocess
from pathlib import Path
from typing import Dict, Any, Optional
import yaml
import pandas as pd

# ============== 配置（YAML 优先） ==============
LOCAL_CONFIG = {
  "config_yaml": "config.yaml",
  "paths": { "fastq_dir":"data/fastq", "clean_dir":"results/qc/clean", "logs_dir":"logs" },
  "binaries": { "fastp":"fastp" },
  "resources": { "threads": { "qc": 8 } },
  "qc": { "min_read_len": 30, "trim_quality": 20 }
}
SCRIPT_TAG="01_qc_fastp"

# ============== 工具函数（统一实现） ==============
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
def run_cmd(cmd: list[str], fp, cwd_path: Optional[str]=None):
    log(fp, "CMD: " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=to_abs(cwd_path) if cwd_path else None)
    if r.returncode != 0: err(f"外部命令失败，退出码={r.returncode}")

# ============== 主逻辑 ==============
def main():
    cfg = merge_params(LOCAL_CONFIG, load_yaml(LOCAL_CONFIG["config_yaml"]))
    fastq_dir = require_exists(cfg["paths"]["fastq_dir"], "dir")
    clean_dir = ensure_dir(cfg["paths"]["clean_dir"])
    logs_dir  = cfg["paths"].get("logs_dir","logs")
    fp = log_open(logs_dir, SCRIPT_TAG); write_params_snapshot(cfg, logs_dir, SCRIPT_TAG)

    rep_dir  = ensure_dir(Path(clean_dir).parent/"report")
    rep_json = ensure_dir(rep_dir/"fastp_json")
    rep_html = ensure_dir(rep_dir/"fastp_html")

    fq = sorted([p for p in fastq_dir.glob("*.f*q.gz")])
    if not fq: err("fastq_dir 下未发现 FASTQ.GZ")

    # 配对识别
    pairs={}
    for p in fq:
        root = re.sub(r"(_R[12]|_[12])\.f.*q\.gz$", "", p.name)
        mate = "1" if re.search(r"_R?1\.", p.name) else ("2" if re.search(r"_R?2\.", p.name) else "SE")
        d = pairs.setdefault(root, {"r1":None,"r2":None,"se":None})
        if mate=="1": d["r1"]=p
        elif mate=="2": d["r2"]=p
        else: d["se"]=p

    pairing_rows=[]
    for s,d in sorted(pairs.items()):
        mode = "SE" if d["se"] is not None else "PE"
        r1p = d["se"] if d["se"] is not None else d["r1"]
        r2p = None if d["se"] is not None else d["r2"]
        out1 = Path(clean_dir)/ (f"{s}.clean.fq.gz" if mode=="SE" else f"{s}_R1.clean.fq.gz")
        out2 = "" if mode=="SE" else str(Path(clean_dir)/f"{s}_R2.clean.fq.gz")
        pairing_rows.append({
            "sample_root": s, "mode": mode,
            "pair_complete": "yes" if mode=="SE" or (d["r1"] and d["r2"]) else "no",
            "r1_path": str(r1p) if r1p else "",
            "r2_path": str(r2p) if r2p else "",
            "clean_r1": str(out1), "clean_r2": out2
        })

    # fastp 执行（不剪接头，做质控与报告；如需剪接头可改为启用 fastp 的 adapter 模块或外接 cutadapt）
    for s,d in sorted(pairs.items()):
        mode = "SE" if d["se"] is not None else "PE"
        if mode=="PE" and (d["r1"] is None or d["r2"] is None): err(f"成对样本不完整：{s}")
        json_out = rep_json/f"{s}.json"; html_out = rep_html/f"{s}.html"
        if mode=="SE":
            cmd = [
                cfg["binaries"]["fastp"],
                "-i", str(d["se"]),
                "-o", str(Path(clean_dir)/f"{s}.clean.fq.gz"),
                "--disable_adapter_trimming",
                "-w", str(cfg["resources"]["threads"]["qc"]),
                "-l", str(cfg["qc"]["min_read_len"]),
                "-q", str(cfg["qc"]["trim_quality"]),
                "-j", str(json_out),
                "-h", str(html_out)
            ]
        else:
            cmd = [
                cfg["binaries"]["fastp"],
                "-i", str(d["r1"]),
                "-I", str(d["r2"]),
                "-o", str(Path(clean_dir)/f"{s}_R1.clean.fq.gz"),
                "-O", str(Path(clean_dir)/f"{s}_R2.clean.fq.gz"),
                "--disable_adapter_trimming",
                "-w", str(cfg["resources"]["threads"]["qc"]),
                "-l", str(cfg["qc"]["min_read_len"]),
                "-q", str(cfg["qc"]["trim_quality"]),
                "-j", str(json_out),
                "-h", str(html_out)
            ]
        run_cmd(cmd, fp)

    # JSON 汇总
    rows_before=[]; rows_after=[]
    for j in sorted(rep_json.glob("*.json")):
        data = json.loads(Path(j).read_text(encoding="utf-8"))
        summ = data.get("summary", {})
        bf = data.get("before_filtering", {})
        af = data.get("after_filtering", {})
        def pick(d,k,default=None): return d.get(k, default) if isinstance(d, dict) else default
        name = data.get("sample_name") or j.stem
        rows_before.append({
            "sample": name,
            "total_reads": pick(bf,"total_reads",0),
            "total_bases": pick(bf,"total_bases",0),
            "q20_rate": pick(bf,"q20_rate",0.0),
            "q30_rate": pick(bf,"q30_rate",0.0),
            "gc_content": pick(bf,"gc_content",0.0),
            "read_len_mean": pick(summ,"read1_mean_length", pick(bf,"mean_read_length", None)),
            "duplication_rate": data.get("duplication_rate", None)
        })
        rows_after.append({
            "sample": name,
            "total_reads": pick(af,"total_reads",0),
            "total_bases": pick(af,"total_bases",0),
            "q20_rate": pick(af,"q20_rate",0.0),
            "q30_rate": pick(af,"q30_rate",0.0),
            "gc_content": pick(af,"gc_content",0.0),
            "read_len_mean": pick(summ,"read1_mean_length", pick(af,"mean_read_length", None)),
            "duplication_rate": data.get("duplication_rate", None)
        })

    rep_tab = ensure_dir(Path(clean_dir).parent/"report")
    pd.DataFrame(rows_before).to_csv(rep_tab/"fastq_stats.tsv", sep="\t", index=False)
    pd.DataFrame(rows_after).to_csv(rep_tab/"clean_stats.tsv", sep="\t", index=False)
    pd.DataFrame(pairing_rows).to_csv(rep_tab/"fastq_pairing.tsv", sep="\t", index=False)
    log(fp, "QC 完成（原生 fastp 报告 + 等效汇总表已生成）")

if __name__=="__main__":
    main()

