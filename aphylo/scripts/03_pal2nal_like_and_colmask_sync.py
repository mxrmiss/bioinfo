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
RAW_MSA_SUFFIX = ".raw.fa"   # 与 10_publish_aphylo_ready.py 发布的 raw AA MSA 后缀一致
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


def _expand_publish_dir_placeholder(s: str, publish_dir: str) -> str:
    """替换 <publish_dir> 占位符（保持原逻辑不变：仅做字符串替换，不做扫描或猜测）"""
    if not isinstance(s, str):
        return s
    if "<publish_dir>" in s:
        return s.replace("<publish_dir>", publish_dir)
    return s

def _coalesce_from_inputs_block(cfg: dict):
    """
    兼容 aphylo 的统一配置结构：
      publish_dir: "../phylo/results/publish/aphylo_ready"
      inputs:
        sco_msa_dir: "<publish_dir>/strict_sco_msa"
        sco_msa_suffix: ".trim.faa"
        raw_msa_suffix: ".raw.fa"   # 新增
        colmask_dir: "<publish_dir>/colmask"
        colmask_suffix: ".colmask"
    说明：只做“显式键读取 + <publish_dir> 替换”，不做目录扫描/猜测。
    """
    global MSA_DIR, MSA_SUFFIX, RAW_MSA_SUFFIX, COLMASK_DIR, COLMASK_SUFFIX
    if not isinstance(cfg, dict):
        return
    pub = cfg.get("publish_dir", None)
    inputs = cfg.get("inputs", None)
    if not isinstance(inputs, dict):
        return
    if isinstance(pub, str) and pub.strip():
        pub = pub.strip()
        if isinstance(inputs.get("sco_msa_dir"), str):
            MSA_DIR = _expand_publish_dir_placeholder(inputs["sco_msa_dir"], pub)
        if isinstance(inputs.get("colmask_dir"), str):
            COLMASK_DIR = _expand_publish_dir_placeholder(inputs["colmask_dir"], pub)
    if isinstance(inputs.get("sco_msa_suffix"), str) and inputs.get("sco_msa_suffix").strip():
        MSA_SUFFIX = inputs["sco_msa_suffix"].strip()
    if isinstance(inputs.get("raw_msa_suffix"), str) and inputs.get("raw_msa_suffix").strip():
        RAW_MSA_SUFFIX = inputs["raw_msa_suffix"].strip()
    if isinstance(inputs.get("colmask_suffix"), str) and inputs.get("colmask_suffix").strip():
        COLMASK_SUFFIX = inputs["colmask_suffix"].strip()
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


def _apply_aa_mask(raw_aa: str, bits: str) -> str:
    """按 bits(0/1) 从 raw AA 对齐中取出保留列（'1'），返回 trim AA 字符串。"""
    if len(raw_aa) != len(bits):
        raise ValueError(f"[ERR] _apply_aa_mask 长度不一致：raw_aa={len(raw_aa)} vs bits={len(bits)}")
    out=[]
    for ch, b in zip(raw_aa, bits):
        if b == "1":
            out.append(ch)
    return "".join(out)

def _apply_codon_mask(raw_codon: str, bits: str) -> str:
    """
    raw_codon 为 pal2nal-like 输出的“raw 坐标”密码子对齐（长度 = raw_cols*3）。
    bits 长度 = raw_cols（raw AA 列数）。返回投影到 trim 坐标的密码子对齐（长度 = ones*3）。
    """
    raw_cols = len(bits)
    if len(raw_codon) != raw_cols * 3:
        raise ValueError(f"[ERR] _apply_codon_mask 长度不一致：raw_codon={len(raw_codon)} vs raw_cols*3={raw_cols*3}")
    out=[]
    for i, b in enumerate(bits):
        if b == "1":
            out.append(raw_codon[i*3:(i+1)*3])
    return "".join(out)

def _first_mismatch(a: str, b: str) -> int:
    """返回第一处不等的位置（0-based）；若相等返回 -1。"""
    n = min(len(a), len(b))
    for i in range(n):
        if a[i] != b[i]:
            return i
    if len(a) != len(b):
        return n
    return -1
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
        _coalesce_from_inputs_block(cfg)

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
            # 读取 trim AA（发布包中的 strict_sco_msa/*.trim.faa）
            aa_trim = read_fasta_dict(msa_path)
            headers=list(aa_trim.keys())
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

            # 新方案：使用 raw AA + raw→trim colmask 来生成“trim 坐标”的 codon 对齐
            #  - raw AA：strict_sco_msa/OGxxxx{RAW_MSA_SUFFIX}
            #  - colmask：publish/colmask/OGxxxx.colmask（长度=raw_cols；'1' 的数量=trim_cols）
            cm_path = Path(COLMASK_DIR) / f"{og}{COLMASK_SUFFIX}"
            raw_cols = None
            trim_cols = None
            bits = None
            if cm_path.exists():
                raw_path = msa_dir / f"{og}{RAW_MSA_SUFFIX}"
                if not raw_path.exists():
                    raise FileNotFoundError(f"[ERR] OG={og} 存在 colmask 但缺少 raw AA：{raw_path}")
                aa_raw = read_fasta_dict(raw_path)
                bad2=[h for h in aa_raw.keys() if "|" in h]
                if bad2:
                    raise ValueError(f"[契约违规] raw AA 表头必须仅为物种名；发现含 '|' 的 header：{bad2[:3]}")
                # 覆盖检查：raw/trim 必须有同一批物种
                if set(aa_raw.keys()) != set(headers):
                    miss = sorted(set(headers) - set(aa_raw.keys()))
                    extra = sorted(set(aa_raw.keys()) - set(headers))
                    raise ValueError(f"[ERR] raw/trim 物种集合不一致：missing={miss[:5]} extra={extra[:5]}")

                trim_cols = len(next(iter(aa_trim.values()))) if aa_trim else 0
                raw_cols  = len(next(iter(aa_raw.values())))  if aa_raw  else 0

                bits = load_colmask_bits(cm_path, raw_cols)
                ones = bits.count("1")
                if ones != trim_cols:
                    raise ValueError(f"[ERR] colmask 的 1 数量 != trim AA 列数：OG={og} ones={ones} trim_cols={trim_cols}")

                # 强校验：raw 上按 colmask 投影后必须逐物种严格等于 trim AA
                for sp in headers:
                    proj = _apply_aa_mask(aa_raw[sp], bits)
                    if proj != aa_trim[sp]:
                        mm = _first_mismatch(proj, aa_trim[sp])
                        raise ValueError(f"[ERR] raw→trim 投影不一致：OG={og} sp={sp} first_mismatch_col={mm}")

                # 以 raw AA 为坐标回译，再投影到 trim 坐标（避免“trim AA 与 CDS”在 gap 位置映射出错）
                from collections import OrderedDict
                aa_raw_ord = OrderedDict((sp, aa_raw[sp]) for sp in headers)
                codon_raw_out, stats = aa_to_codon_msa(aa_raw_ord, cds_map, strict=STRICT_CHECK)

                codon_out = OrderedDict()
                for sp, ntseq in codon_raw_out.items():
                    codon_out[sp] = _apply_codon_mask(ntseq, bits)

                # 输出的 codon 对齐列数（trim 坐标）
                n_cols = trim_cols

                # ntmask 仍需与 codon_out 一致：trim 坐标下全部为 1（因为已经投影/删列）
                ntmask = "1" * (n_cols * 3)
                (Path(NTMASK_DIR)/f"{og}.ntmask").write_text(ntmask+"\n", encoding="utf-8")
                have_bits=True
            else:
                # 兼容旧模式：没有 colmask 时，直接用 trim AA 回译（保持原逻辑）
                aa = aa_trim
                n_cols = len(next(iter(aa.values()))) if aa else 0
                codon_out, stats = aa_to_codon_msa(aa, cds_map, strict=STRICT_CHECK)
                warn(f"OG={og} 缺少发布的 colmask；使用旧模式直接回译 trim AA，跳过 ntmask 生成。")
                have_bits=False

            out_fa = Path(CODON_MSA_DIR) / f"{og}.codon.fna"
            write_fasta_dict(out_fa, codon_out, linewrap=LINE_WRAP)
            with open(oglog_path,"w",encoding="utf-8") as lf:
                lf.write(f"# OG={og}\n")
                lf.write(f"aa_cols={n_cols}\n")
                lf.write(f"n_taxa={len(codon_out)}\n")
                if have_bits:
                    lf.write("colmask=bitstring\n")
                    lf.write(f"colmask_len={len(bits)}\n")
                    lf.write(f"raw_aa_cols={raw_cols}\n")
                    lf.write(f"trim_aa_cols={n_cols}\n")
                    lf.write(f"colmask_ones={bits.count('1')}\n")
                    lf.write(f"colmask_zeros={len(bits)-bits.count('1')}\n")
                    lf.write(f"ntmask_len={n_cols*3}\n")
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


