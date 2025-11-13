#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_pal2nal_like_and_colmask_sync.py —— 与 phylo 发布“无缝对接”的最终版（路径按 PDF）
- 仅路径改动：BT_DIR 改为 results/03_codon
"""
import os, re, sys
from pathlib import Path
from collections import OrderedDict
from typing import Dict, List, Tuple

CONFIG_FILE   = "config.yaml"

# —— phylo 发布的输入（AA 与 AA 位串 colmask）——
MSA_DIR        = "../phylo/results/publish/aphylo_ready/strict_sco_msa"
MSA_SUFFIX     = ".trim.faa"
COLMASK_DIR    = "../phylo/results/publish/aphylo_ready/colmask"
COLMASK_SUFFIX = ".colmask"   # 0/1 位串

# —— 02 步产物 —— 
ORDER_DIR      = "results/02_bt/order"
CDS_DIR        = "results/02_bt/tmp"

# —— 03 输出目录：按 PDF 命名 ——（仅此处改动）
BT_DIR         = "results/03_codon"
CODON_MSA_DIR  = f"{BT_DIR}/codon_msa"
NTMASK_DIR     = f"{BT_DIR}/colmask"
SENTINEL_PATH  = f"{BT_DIR}/.pal2nal.done"

# —— 日志目录（保持不变）——
LOG_DIR        = "logs/03_pal2nal"

LINE_WRAP      = 60
STRICT_CHECK   = True
MKDIR_MODE     = 0o755

def time_str():
    import datetime as _dt
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
def log(msg: str):  print(f"[{time_str()}] INFO  - {msg}", flush=True)
def warn(msg: str): print(f"[{time_str()}] WARN  - {msg}", flush=True)
def err(msg: str):  print(f"[{time_str()}] ERROR - {msg}", flush=True)

def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    try: os.chmod(p, MKDIR_MODE)
    except Exception: pass

def read_yaml_if_present(path: Path) -> dict:
    if not path.exists(): return {}
    try:
        import yaml
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}

def coalesce_from_config(cfg: dict):
    # 允许存在但不强制：保持原逻辑（只改默认路径）
    global MSA_DIR, MSA_SUFFIX, COLMASK_DIR, COLMASK_SUFFIX
    global ORDER_DIR, CDS_DIR, BT_DIR, CODON_MSA_DIR, NTMASK_DIR, SENTINEL_PATH, LOG_DIR
    MSA_DIR        = cfg.get("msa_dir", MSA_DIR)
    MSA_SUFFIX     = cfg.get("msa_suffix", MSA_SUFFIX)
    COLMASK_DIR    = cfg.get("colmask_dir", COLMASK_DIR)
    COLMASK_SUFFIX = cfg.get("colmask_suffix", COLMASK_SUFFIX)
    ORDER_DIR      = cfg.get("order_dir", ORDER_DIR)
    CDS_DIR        = cfg.get("cds_dir_02tmp", cfg.get("cds_dir", CDS_DIR))
    # PDF 命名：若未覆盖，仍写到 results/03_codon/*
    BT_DIR         = cfg.get("bt_dir", BT_DIR)
    CODON_MSA_DIR  = cfg.get("codon_msa_dir", f"{BT_DIR}/codon_msa")
    NTMASK_DIR     = cfg.get("ntmask_dir",    f"{BT_DIR}/colmask")
    SENTINEL_PATH  = cfg.get("sentinel_03",   f"{BT_DIR}/.pal2nal.done")
    LOG_DIR        = cfg.get("log_dir_03", LOG_DIR)

    step = cfg.get("step03", {})
    sp = step.get("paths", {}) if isinstance(step, dict) else {}
    if sp.get("msa_dir"):        MSA_DIR = sp["msa_dir"]
    if sp.get("colmask_dir"):    COLMASK_DIR = sp["colmask_dir"]
    if sp.get("order_dir"):      ORDER_DIR = sp["order_dir"]
    if sp.get("cds_dir"):        CDS_DIR = sp["cds_dir"]
    if sp.get("bt_dir"):         BT_DIR = sp["bt_dir"]
    if sp.get("codon_msa_dir"):  CODON_MSA_DIR = sp["codon_msa_dir"]
    else:                        CODON_MSA_DIR = f"{BT_DIR}/codon_msa"
    if sp.get("ntmask_dir"):     NTMASK_DIR = sp["ntmask_dir"]
    else:                        NTMASK_DIR = f"{BT_DIR}/colmask"
    if sp.get("sentinel"):       SENTINEL_PATH = sp["sentinel"]
    else:                        SENTINEL_PATH = f"{BT_DIR}/.pal2nal.done"
    if sp.get("log_dir"):        LOG_DIR = sp["log_dir"]

def read_fasta_dict(path: Path) -> "OrderedDict[str,str]":
    d: "OrderedDict[str,str]" = OrderedDict()
    with open(path, "r", encoding="utf-8") as f:
        name=None; buf=[]
        for line in f:
            if line.startswith(">"):
                if name is not None: d[name] = re.sub(r"\s+","", "".join(buf)).upper()
                name=line[1:].strip(); buf=[]
            else:
                buf.append(line.strip())
        if name is not None: d[name] = re.sub(r"\s+","", "".join(buf)).upper()
    return d

def write_fasta_dict(path: Path, d: "OrderedDict[str,str]", linewrap: int = 60):
    with open(path, "w", encoding="utf-8") as f:
        for h, s in d.items():
            f.write(f">{h}\n")
            if linewrap and linewrap>0:
                for i in range(0, len(s), linewrap): f.write(s[i:i+linewrap] + "\n")
            else:
                f.write(s + "\n")

def is_og_file(p: Path) -> bool: return p.is_file() and p.name.startswith("OG") and p.name.endswith(MSA_SUFFIX)
def og_from_filename(p: Path) -> str: return p.name.split(MSA_SUFFIX, 1)[0]

def load_colmask_bits(cm_path: Path, n_cols: int) -> str:
    if not cm_path.exists(): raise FileNotFoundError(f"[ERR] 缺少 colmask：{cm_path}")
    bits = re.sub(r"\s+","", cm_path.read_text(encoding="utf-8"))
    if not bits: raise ValueError(f"[ERR] {cm_path.name} 为空")
    if re.search(r"[^01]", bits): raise ValueError(f"[ERR] {cm_path.name} 含非 0/1 字符")
    if len(bits)!=n_cols: raise ValueError(f"[ERR] {cm_path.name} 长度≠AA 列数：mask={len(bits)} vs AA={n_cols}")
    return bits

def _looks_like_header(parts: List[str]) -> bool:
    cols = [re.sub(r"\s+","", x).lower() for x in parts]
    return len(cols)>=4 and cols[:4]==["og","species","protein_id","cds_id"]

def load_order_map(order_path: Path, og: str) -> Dict[str, str]:
    if not order_path.exists(): raise FileNotFoundError(f"[ERR] 缺少 order：{order_path}")
    sp2cds: Dict[str, str] = {}
    seen=False
    with open(order_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f,1):
            s=line.lstrip("\ufeff").strip()
            if not s or s.startswith("#"): continue
            parts=re.split(r"[ \t]+", s)
            if not seen and _looks_like_header(parts): continue
            if len(parts)<4: raise ValueError(f"[ERR] {order_path.name} 第{ln}列数<{4}：{s}")
            og_in, sp, _pid, cds = parts[0], parts[1], parts[2], parts[3]
            if og_in!=og: raise ValueError(f"[ERR] {order_path.name} 第{ln}行 OG 不匹配：'{og_in}'≠'{og}'")
            if sp in sp2cds: raise ValueError(f"[ERR] {order_path.name} 重复物种键：{sp}（第{ln}行）")
            sp2cds[sp]=cds; seen=True
    if not sp2cds: raise ValueError(f"[ERR] {order_path.name} 为空")
    return sp2cds

def aa_to_codon_msa(aa_msa: "OrderedDict[str,str]", cds_map: Dict[str,str], strict: bool=True):
    out: "OrderedDict[str,str]" = OrderedDict()
    stats={}
    for sp, aaseq in aa_msa.items():
        cds = cds_map.get(sp,""); total=len(cds)
        if strict and (total%3!=0): raise ValueError(f"[pal2nal] {sp} 的 CDS 长度非 3 的倍数：total_nt={total}")
        pos=0; nts=[]
        for aa in aaseq:
            if aa=="-": nts.append("---")
            else:
                trip = cds[pos:pos+3]
                if len(trip)!=3: raise ValueError(f"[pal2nal] {sp} 的 CDS 不足：pos={pos}, got={len(trip)}")
                nts.append(trip); pos+=3
        out[sp]="".join(nts); used=pos
        if strict and used>total: raise ValueError(f"[pal2nal] {sp} used_nt({used}) > total_nt({total})")
        stats[sp]=(used,total,total-used)
    return out, stats

def main():
    cfg = read_yaml_if_present(Path(CONFIG_FILE))
    if isinstance(cfg, dict):
        step03 = cfg.get("step03", {})
        if isinstance(step03, dict):
            for k,v in step03.items():
                if k in ("paths",): continue
                cfg[k]=v
        coalesce_from_config(cfg)

    log("APhylo 03 —— pal2nal-like（输出到 results/03_codon/*）")
    log(f"msa_dir      = {MSA_DIR}")
    log(f"order_dir    = {ORDER_DIR}")
    log(f"cds_dir      = {CDS_DIR}")
    log(f"out_codon    = {CODON_MSA_DIR}")
    log(f"out_ntmask   = {NTMASK_DIR}")
    log(f"log_dir      = {LOG_DIR}")

    for p in (Path(CODON_MSA_DIR), Path(NTMASK_DIR), Path(LOG_DIR), Path(SENTINEL_PATH).parent):
        safe_mkdir(p)

    msa_dir = Path(MSA_DIR)
    og_files = sorted([p for p in msa_dir.iterdir() if is_og_file(p)])
    log(f"发现 OG 数量：{len(og_files)}")
    if not og_files:
        warn("未发现 OG*.trim.faa；请检查发布目录。")
        Path(SENTINEL_PATH).write_text("EMPTY\n", encoding="utf-8")
        return

    processed=0
    for msa_path in og_files:
        og = og_from_filename(msa_path)
        oglog_path = Path(LOG_DIR) / f"pal2nal_{og}.log"
        try:
            aa = read_fasta_dict(msa_path)
            headers=list(aa.keys())
            bad=[h for h in headers if "|" in h]
            if bad: raise ValueError(f"[契约违规] AA 表头必须仅为物种名；发现含 '|' 的 header：{bad[:3]}")

            order_path = Path(ORDER_DIR) / f"{og}.order.tsv"
            cds_path   = Path(CDS_DIR)   / f"{og}.ordered_cds.fna"
            sp2cds = load_order_map(order_path, og)
            cds_fa = read_fasta_dict(cds_path)
            cds_map={}
            for sp in headers:
                cds_id = sp2cds.get(sp)
                if cds_id is None: raise KeyError(f"[ERR] order.tsv 缺少物种 {sp} 的 cds_id")
                if cds_id not in cds_fa: raise KeyError(f"[ERR] ordered_cds.fna 缺少条目：{cds_id}")
                cds_map[sp]=cds_fa[cds_id]

            n_cols = len(next(iter(aa.values()))) if aa else 0
            codon_out, stats = aa_to_codon_msa(aa, cds_map, strict=STRICT_CHECK)

            out_fa = Path(CODON_MSA_DIR) / f"{og}.codon.fna"
            write_fasta_dict(out_fa, codon_out, linewrap=LINE_WRAP)

            cm_path = Path(COLMASK_DIR) / f"{og}{COLMASK_SUFFIX}"
            if cm_path.exists():
                bits = load_colmask_bits(cm_path, n_cols)
                ntmask = "".join(ch*3 for ch in bits)
                (Path(NTMASK_DIR)/f"{og}.ntmask").write_text(ntmask+"\n", encoding="utf-8")
                have_bits=True
            else:
                warn(f"OG={og} 缺少发布的 colmask；跳过 ntmask 生成。")
                have_bits=False

            with open(oglog_path,"w",encoding="utf-8") as lf:
                lf.write(f"# OG={og}\n")
                lf.write(f"aa_cols={n_cols}\n")
                lf.write(f"n_taxa={len(codon_out)}\n")
                if have_bits:
                    lf.write("colmask=bitstring\n")
                    lf.write(f"colmask_len={len(bits)}\n")
                    lf.write(f"ntmask_len={len(bits)*3}\n")
                else:
                    lf.write("colmask=absent\n")
                total_left=0
                for sp in headers:
                    used,total,left=stats[sp]; total_left+=left
                    lf.write(f"{sp}\tused_nt={used}\ttotal_nt={total}\tleftover_nt={left}\n")
                lf.write(f"leftover_nt_sum={total_left}\n")

            processed+=1
            if processed%200==0: log(f"进度：{processed}/{len(og_files)}")
        except Exception as e:
            err(f"OG={og} 失败：'{e}'")
            with open(oglog_path,"w",encoding="utf-8") as lf:
                lf.write(f"[ERROR] {e}\n")
            continue

    Path(SENTINEL_PATH).write_text("OK\n", encoding="utf-8")
    log(f"完成：成功处理 {processed}/{len(og_files)} 个 OG；哨兵写入 {SENTINEL_PATH}")

if __name__ == "__main__":
    main()

