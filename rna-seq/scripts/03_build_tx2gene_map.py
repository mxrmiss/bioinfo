#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_build_tx2gene_map.py — 从 GFF3/GTF 生成 tx2gene（原始 + 清洗版）
产出：
  - ref/tx2gene.raw.tsv        （列：TX,GENE）
  - ref/tx2gene.geneMap.tsv    （清洗后，供 --geneMap 与 tximport 使用）
"""

from __future__ import annotations
import sys, json, time, re
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import yaml
import pandas as pd

# ================== 脚本默认（YAML 可覆盖） ==================
LOCAL_CONFIG = {
  "config_yaml": "config.yaml",
  "paths": {
    "ref_gtf":   "ref/genes.gff3",
    "anno_dir":  "ref/annotations",
    "logs_dir":  "logs"
  }
}
SCRIPT_TAG = "03_tx2gene"

# ================== 通用工具 ==================
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

# ================== GTF/GFF3 解析 ==================
def parse_attr_kv(attr_raw: str) -> Dict[str,str]:
    """兼容 GTF（k "v";）与 GFF3（k=v;）属性解析"""
    s = attr_raw.strip()
    kv = {}
    # 判断格式
    if "=" in s and ";" in s and not ('"' in s):
        # GFF3
        for item in s.split(";"):
            item=item.strip()
            if not item: continue
            if "=" in item:
                k,v = item.split("=",1)
                kv[k.strip()] = v.strip()
    else:
        # GTF 风格：key "value";
        for item in s.split(";"):
            item=item.strip()
            if not item: continue
            parts = item.split()
            if len(parts)>=2:
                k = parts[0].strip()
                v = " ".join(parts[1:]).strip().strip('"')
                kv[k]=v
    return kv

def iter_tx_gene_pairs(gtf_path: Path):
    """
    仅抓取转录本层面的记录：
      - GFF3: 第3列为 mRNA / transcript，属性 ID=<tx>; Parent=<gene>
      - GTF : 第3列任意（exon/CDS/transcript），属性 transcript_id / gene_id 存在即可
    """
    with open(gtf_path, "r", encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            if not line or line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9: continue
            feature = cols[2].lower()
            attrs = parse_attr_kv(cols[8])
            tx = None; gene = None

            # 优先 GTF 常用键
            if "transcript_id" in attrs and "gene_id" in attrs:
                tx   = attrs.get("transcript_id","").strip()
                gene = attrs.get("gene_id","").strip()
            else:
                # GFF3：mRNA/transcript 层
                if feature in ("mrna","transcript"):
                    tx   = attrs.get("ID","").strip()
                    gene = attrs.get("Parent","").strip()

            if tx and gene:
                yield tx, gene

# ================== 清洗规则（与老口径一致） ==================
def clean_tx_id(tx: str) -> str:
    """
    归一化转录本 ID，保证与 quant.sf 的 Name 列对齐：
      1) 去前缀：cds\d+.
      2) 去后缀：.exon\d+
      3) 去版本号：.<数字>（只去最后一段小数版本）
      4) 去竖线后缀：| 之后全部删除
    """
    x = tx
    x = re.sub(r"^cds\d+\.", "", x)
    x = re.sub(r"\.exon\d+$", "", x)
    x = re.sub(r"\|.*$", "", x)          # 忽略 | 之后
    x = re.sub(r"\.(\d+)$", "", x)       # 去结尾 .1/.2
    return x

# ================== 主逻辑 ==================
def main():
    cfg = merge_params(LOCAL_CONFIG, load_yaml(LOCAL_CONFIG["config_yaml"]))
    logs_dir = cfg["paths"].get("logs_dir","logs")
    fp = log_open(logs_dir, SCRIPT_TAG)
    write_params_snapshot(cfg, logs_dir, SCRIPT_TAG)

    gtf = require_exists(cfg["paths"]["ref_gtf"], "file")
    anno_dir = ensure_dir(cfg["paths"].get("anno_dir","ref/annotations"))
    raw_out  = to_abs("ref/tx2gene.raw.tsv")
    map_out  = to_abs("ref/tx2gene.geneMap.tsv")

    log(fp, f"读取注释：{gtf}")
    pairs = list(iter_tx_gene_pairs(gtf))
    if not pairs:
        err("未在注释中解析到任何 transcript→gene 对，应检查 GFF3/GTF 与属性名。")

    df = pd.DataFrame(pairs, columns=["TX","GENE"]).drop_duplicates()
    # 过滤掉明显异常的空白
    df = df[(df["TX"].astype(str)!="") & (df["GENE"].astype(str)!="")].copy()
    n_raw = len(df)
    df.to_csv(raw_out, sep="\t", index=False)
    log(fp, f"原始 tx2gene 写出：{raw_out}（{n_raw} 行）")

    # 清洗 ID
    dfc = df.copy()
    dfc["TX"] = dfc["TX"].astype(str).map(clean_tx_id)
    dfc = dfc.drop_duplicates()
    n_map = len(dfc)
    dfc.to_csv(map_out, sep="\t", index=False)
    log(fp, f"清洗版 geneMap 写出：{map_out}（{n_map} 行）")

    # 生成用于富集/注释的目录（若后续需要）
    ensure_dir(anno_dir)
    log(fp, "完成 ✅")

if __name__ == "__main__":
    main()

