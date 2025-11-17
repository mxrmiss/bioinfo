#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_salmon_quant.py — Salmon 准比对定量（PE/SE，自适应解析 FASTQ 路径，带回退）
"""

from __future__ import annotations
import sys, json, time, subprocess, shutil
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List
import yaml
import pandas as pd

LOCAL_CONFIG = {
  "config_yaml": "config.yaml",
  "paths": {
    "salmon_index": "results/refprep/salmon_index",
    "quant_dir":    "results/quant",
    "logs_dir":     "logs",
    "tx2gene":      "ref/tx2gene.geneMap.tsv",
    "fastq_dir":    "data/fastq"
  },
  "tables": { "samples": "data/samples.tsv" },
  "binaries": { "salmon": "salmon" },
  "resources": { "threads": { "quant": 8 } },
  "salmon": { "libType": "A", "use_geneMap": True }
}
SCRIPT_TAG = "04_salmon_quant"

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
    if kind=="dir"  and not p.is_dir() : err(f"缺少目录：{p}")
    return p

def exe_exists(name: str) -> bool: return shutil.which(name) is not None

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

# -------- 样本表与路径解析 --------
def parse_samples_table(samples_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(samples_tsv, sep="\t", dtype=str).fillna("")
    cols = {c.lower(): c for c in df.columns}
    def pick(*cands):
        for k in cands:
            if k in cols: return cols[k]
        return None
    c_sample = pick("sample","samples","id","sid")
    c_r1     = pick("r1","fastq1","fq1","read1","file1")
    c_r2     = pick("r2","fastq2","fq2","read2","file2")
    if c_sample is None or c_r1 is None:
        err("samples.tsv 至少需要列：sample 与 R1（大小写/同义名均可）")
    out = df[[c_sample]].rename(columns={c_sample:"sample"})
    out["R1"] = df[c_r1].astype(str)
    out["R2"] = df[c_r2].astype(str) if c_r2 else ""
    out = out[out["sample"].astype(str).str.len()>0].copy()
    if out.empty: err("samples.tsv 解析为空，请检查")
    return out.reset_index(drop=True)

def resolve_fastq_with_fallback(raw: str, fastq_dir: Path, tried: List[Path]) -> Optional[Path]:
    """
    解析顺序（并记录尝试过的候选路径）：
      1) 绝对路径（若存在）→ 使用
      2) 相对路径（含 '/'）→ 以项目根解析；若不存在，则回退到 fastq_dir / basename(raw)
      3) 仅文件名 → fastq_dir / 文件名
    """
    if raw is None or str(raw).strip()=="":
        return None
    s = str(raw).strip()
    # 1) 绝对路径
    p = Path(s)
    if p.is_absolute():
        tried.append(p)
        return p if p.exists() else p  # 让上层统一校验
    # 2) 相对路径（含 '/'）
    if "/" in s:
        p1 = to_abs(p)
        tried.append(p1)
        if p1.exists():
            return p1
        # 回退到 fastq_dir / basename
        p2 = to_abs(fastq_dir / p.name)
        tried.append(p2)
        return p2
    # 3) 仅文件名
    p3 = to_abs(fastq_dir / s)
    tried.append(p3)
    return p3

def build_salmon_cmd(sample: str, r1: Path, r2: Optional[Path], outdir: Path, cfg: dict) -> list[str]:
    cmd = [
        cfg["binaries"]["salmon"], "quant",
        "-i", str(require_exists(cfg["paths"]["salmon_index"], "dir")),
        "-l", str(cfg["salmon"]["libType"]),
        "-p", str(cfg["resources"]["threads"]["quant"]),
        "-o", str(outdir),
        "--validateMappings",
        "--gcBias"
    ]
    if r2 and str(r2).strip()!="":
        cmd += ["-1", str(require_exists(r1, "file")), "-2", str(require_exists(r2, "file"))]
    else:
        cmd += ["-r", str(require_exists(r1, "file"))]
    if bool(cfg["salmon"].get("use_geneMap", True)):
        gm = cfg["paths"].get("tx2gene", "")
        if gm:
            gm_p = require_exists(gm, "file")
            cmd += ["--geneMap", str(gm_p)]
    return cmd

# -------- 主逻辑 --------
def main():
    cfg = merge_params(LOCAL_CONFIG, load_yaml(LOCAL_CONFIG["config_yaml"]))
    logs_dir = cfg["paths"].get("logs_dir","logs")
    fp = log_open(logs_dir, SCRIPT_TAG)
    write_params_snapshot(cfg, logs_dir, SCRIPT_TAG)

    if not exe_exists(cfg["binaries"]["salmon"]):
        err("未找到 salmon 可执行（binaries.salmon）。请在 PATH 中提供或在 config.yaml 中指定。")
    require_exists(cfg["paths"]["salmon_index"], "dir")
    samples_tsv = require_exists(cfg["tables"]["samples"], "file")
    fastq_dir   = ensure_dir(cfg["paths"].get("fastq_dir","data/fastq"))
    quant_dir   = ensure_dir(cfg["paths"]["quant_dir"])

    df = parse_samples_table(samples_tsv)
    log(fp, f"共解析到 {len(df)} 个样本。FASTQ 根目录：{fastq_dir}")

    for _, row in df.iterrows():
        sample = row["sample"].strip()
        r1_raw = row["R1"].strip()
        r2_raw = row["R2"].strip()

        tried1: List[Path] = []
        tried2: List[Path] = []

        r1 = resolve_fastq_with_fallback(r1_raw, fastq_dir, tried1)
        r2 = resolve_fastq_with_fallback(r2_raw, fastq_dir, tried2) if r2_raw else None

        log(fp, f"[{sample}] R1 = {r1}")
        if r2: log(fp, f"[{sample}] R2 = {r2}")

        # 明确的存在性校验 + 失败时打印全部候选
        if not (r1 and Path(r1).is_file()):
            hints = "\n  - ".join(str(p) for p in tried1)
            err(f"[{sample}] 找不到 R1：{r1}\n已尝试：\n  - {hints}")
        if r2 and not Path(r2).is_file():
            hints = "\n  - ".join(str(p) for p in tried2)
            err(f"[{sample}] 找不到 R2：{r2}\n已尝试：\n  - {hints}")

        outdir = ensure_dir(Path(quant_dir) / sample)
        cmd = build_salmon_cmd(sample, r1, r2, outdir, cfg)
        run_cmd(cmd, fp)

        qsf = outdir/"quant.sf"
        if not qsf.exists():
            err(f"[{sample}] 运行结束但未发现 {qsf}，请检查 Salmon 日志。")
        log(fp, f"[{sample}] 完成，输出目录：{outdir}")

    log(fp, "全部样本定量完成 ✅")
    log(fp, f"结果根目录：{quant_dir}")

if __name__ == "__main__":
    main()

