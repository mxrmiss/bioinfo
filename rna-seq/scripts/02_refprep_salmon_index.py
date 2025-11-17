#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_refprep_salmon_index.py — 构建 Salmon 索引（强制 gentrome + decoys，官方推荐）
用法：在项目根目录运行
  python scripts/02_refprep_salmon_index.py

依赖（从 config.yaml 读取，YAML 优先于脚本默认）：
paths:
  ref_fasta:    ref/Sinonovacula_constricta_genome.fa      # 基因组 FASTA
  ref_gtf:      ref/Sinonovacula_constricta_genome.gff3    # GTF/GFF3，提取转录本用
  salmon_index: results/refprep/salmon_index               # 索引输出目录
  logs_dir:     logs
binaries:
  salmon:  salmon
  gffread: gffread
resources:
  threads:
    salmon_index: 20
"""

from __future__ import annotations
import sys, json, time, subprocess, shutil
from pathlib import Path
from typing import Dict, Any, Optional
import yaml
import pandas as pd

# ================== 配置（脚本默认；YAML 覆盖） ==================
LOCAL_CONFIG = {
  "config_yaml": "config.yaml",
  "paths": {
    "ref_fasta":    "ref/genome.fa",
    "ref_gtf":      "ref/genes.gff3",
    "salmon_index": "results/refprep/salmon_index",
    "logs_dir":     "logs"
  },
  "binaries": { "salmon":"salmon", "gffread":"gffread" },
  "resources": { "threads": { "salmon_index": 8 } }
}
SCRIPT_TAG = "02_salmon_index"

# ================== 工具函数 ==================
def cwd() -> Path: return Path.cwd().resolve()

def to_abs(p: str | Path | None) -> Optional[Path]:
    if p is None: return None
    s = str(p).strip()
    if s == "": return None
    q = Path(s)
    return (q if q.is_absolute() else cwd()/q).resolve()

def load_yaml(path: str) -> Dict[str, Any]:
    p = to_abs(path)
    if not p or not p.exists(): return {}
    with p.open('r', encoding='utf-8') as f:
        return yaml.safe_load(f) or {}

def merge_params(local_cfg: dict, yaml_cfg: dict) -> dict:
    if not isinstance(local_cfg, dict): local_cfg = {}
    if not isinstance(yaml_cfg, dict):  yaml_cfg  = {}
    def is_set(v):
        return v is not None and not (isinstance(v,(str,list,dict)) and len(v)==0) \
               and not (isinstance(v,str) and v.strip()=="")
    def merge(a,b):
        if isinstance(a,dict) and isinstance(b,dict):
            out={}
            for k in set(a.keys())|set(b.keys()):
                av=a.get(k); bv=b.get(k)
                if isinstance(av,dict) or isinstance(bv,dict):
                    out[k]=merge(av if isinstance(av,dict) else {}, bv if isinstance(bv,dict) else {})
                else:
                    out[k]= bv if is_set(bv) else av
            return out
        return b if is_set(b) else a
    return merge(local_cfg, yaml_cfg)

def err(msg: str):
    print(f"[ERR] {msg}", file=sys.stderr); sys.exit(1)

def ensure_dir(p) -> Path:
    p = to_abs(p); p.mkdir(parents=True, exist_ok=True); return p

def require_exists(p, kind="file") -> Path:
    p = to_abs(p)
    if p is None: err("必需路径为空")
    if kind=="file" and not p.is_file(): err(f"缺少文件：{p}")
    if kind=="dir" and not p.is_dir(): err(f"缺少目录：{p}")
    return p

def exe_exists(name: str) -> bool:
    return shutil.which(name) is not None

def log_open(d, tag):
    d = ensure_dir(d); return open(d/f"{tag}.log","a",encoding="utf-8")

def log(fp, msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    fp.write(f"[{ts}] {msg}\n"); fp.flush(); print(msg, flush=True)

def write_params_snapshot(params, logs_dir, tag):
    d = ensure_dir(logs_dir); out = d/f"{tag}.params.tsv"
    rows=[]
    def flat(prefix,obj):
        if isinstance(obj,dict):
            for k,v in obj.items(): flat(f"{prefix}.{k}" if prefix else k, v)
        else:
            rows.append((prefix, json.dumps(obj, ensure_ascii=False)))
    flat("", params)
    pd.DataFrame(rows, columns=["key","value"]).to_csv(out, sep="\t", index=False)

def run_cmd(cmd: list[str], fp, cwd_path: Optional[str]=None, env=None):
    log(fp, "CMD: " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=to_abs(cwd_path) if cwd_path else None, env=env)
    if r.returncode != 0: err(f"外部命令失败，退出码={r.returncode}")

# ================== 业务函数 ==================
def make_transcripts_with_gffread(ref_gtf: Path, ref_fasta: Path, out_transcripts: Path, fp):
    if not exe_exists(LOCAL_CONFIG["binaries"]["gffread"]):
        err("未找到 gffread 可执行（binaries.gffread）。请在 PATH 中提供或在 config.yaml 中指定。")
    run_cmd([
        LOCAL_CONFIG["binaries"]["gffread"],
        str(ref_gtf),
        "-g", str(ref_fasta),
        "-w", str(out_transcripts)
    ], fp)

def write_decoys_from_genome_fa(ref_fasta: Path, out_decoys: Path, fp):
    n = 0
    with open(ref_fasta, 'r', encoding='utf-8', errors='ignore') as fin, \
         open(out_decoys, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('>'):
                name = line[1:].strip().split()[0]
                if name:
                    fout.write(name + "\n"); n += 1
    log(fp, f"decoys.txt 写入 {n} 条序列名：{out_decoys}")

def concat_fa(fa1: Path, fa2: Path, out_fa: Path, fp):
    with open(out_fa, 'w', encoding='utf-8') as fout:
        for src in (fa1, fa2):
            with open(src, 'r', encoding='utf-8', errors='ignore') as fin:
                for line in fin: fout.write(line)
    log(fp, f"gentrome.fa 生成：{out_fa}")

# ================== 主逻辑（强制 gentrome） ==================
def main():
    cfg = merge_params(LOCAL_CONFIG, load_yaml(LOCAL_CONFIG["config_yaml"]))
    logs_dir = cfg["paths"].get("logs_dir","logs")
    fp = log_open(logs_dir, SCRIPT_TAG)
    write_params_snapshot(cfg, logs_dir, SCRIPT_TAG)

    ref_fa  = require_exists(cfg["paths"]["ref_fasta"], "file")
    ref_gtf = require_exists(cfg["paths"]["ref_gtf"],   "file")
    out_idx = ensure_dir(cfg["paths"]["salmon_index"])

    if not exe_exists(cfg["binaries"]["salmon"]):
        err("未找到 salmon 可执行（binaries.salmon）。请在 PATH 中提供或在 config.yaml 中指定。")

    ref_dir = to_abs(cfg["paths"]["ref_fasta"]).parent
    transcripts_fa = ref_dir / "transcripts.fa"
    decoys_txt     = ref_dir / "decoys.txt"
    gentrome_fa    = ref_dir / "gentrome.fa"

    log(fp, "[步骤] 构建 transcripts.fa（gffread）")
    make_transcripts_with_gffread(ref_gtf, ref_fa, transcripts_fa, fp)

    log(fp, "[步骤] 生成 decoys.txt（来自基因组 FASTA 头）")
    write_decoys_from_genome_fa(ref_fa, decoys_txt, fp)

    log(fp, "[步骤] 拼接 gentrome.fa = transcripts.fa + genome.fa")
    concat_fa(transcripts_fa, ref_fa, gentrome_fa, fp)

    log(fp, "[步骤] salmon index（gentrome + decoys）")
    run_cmd([
        cfg["binaries"]["salmon"], "index",
        "-t", str(gentrome_fa),
        "-d", str(decoys_txt),
        "-i", str(out_idx),
        "-p", str(cfg["resources"]["threads"]["salmon_index"])
    ], fp)

    log(fp, "Salmon 索引完成 ✅")
    log(fp, f"索引目录：{out_idx}")
    log(fp, f"中间文件：{transcripts_fa.name}, {decoys_txt.name}, {gentrome_fa.name}")

if __name__ == "__main__":
    main()

