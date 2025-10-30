#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, subprocess, gzip, re
from pathlib import Path
import yaml

def need(cfg, key):
    cur = cfg
    for k in key.split("."):
        if not isinstance(cur, dict) or k not in cur:
            sys.exit(f"[ERR] config 缺少 post.{key}")
        cur = cur[k]
    return cur

def normalize_gid(raw: str) -> str:
    name = raw.strip()
    m = re.search(r'cds_([A-Z]{2}_[0-9]+\.[0-9]+)', name)
    if m: return m.group(1)
    m = re.search(r'\b([NX][PM]_[0-9]+\.[0-9]+)\b', name)
    if m: return m.group(1)
    if '|' in name: name = name.split('|')[-1]
    name = name.split()[0]
    name = re.sub(r'[^A-Za-z0-9_.:-]', '_', name)
    return name

def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = ln[1:].strip(), []
            else:
                seq.append(ln.strip())
        if name is not None:
            yield name, "".join(seq).upper()

def fasta_len_no_gap(path):
    tot, names = 0, 0
    name, seq = None, []
    for ln in open(path, encoding="utf-8"):
        if ln.startswith(">"):
            if name is not None:
                tot += sum(1 for c in "".join(seq) if c not in "-?")
            name, seq = ln[1:].strip(), []
            names += 1
        else:
            seq.append(ln.strip())
    if name is not None:
        tot += sum(1 for c in "".join(seq) if c not in "-?")
    return names, tot

CFG_PATH = str(Path(__file__).resolve().parent.parent / "config.yaml")  # 只读 phylo-post/config.yaml
if not os.path.exists(CFG_PATH):
    sys.exit(f"[ERR] 缺少配置: {CFG_PATH}")
CFG  = yaml.safe_load(open(CFG_PATH, encoding="utf-8")) or {}
POST = CFG.get("post", {}) or {}

OUTDIR   = need(POST, "outdir")
PAL2NAL  = POST.get("env", {}).get("binaries", {}).get("pal2nal", None)
if not PAL2NAL:
    sys.exit("[ERR] post.env.binaries.pal2nal 未配置")

PROT_DIR = os.path.join(OUTDIR, "prot_msa_geneid")
CDS_DIR  = os.path.normpath(need(POST, "codon.cds_dir"))
TARGETS  = os.path.join(OUTDIR, "targets", "og_targets.list")
CODON_DIR= os.path.join(OUTDIR, "codon_sco")
RPT_DIR  = os.path.join(OUTDIR, "reports")

QC_MIN_TAXA  = int(need(POST, "codon.min_taxa"))
QC_MIN_NT    = int(need(POST, "codon.min_codon_len"))
QC_MAX_GAP   = float(need(POST, "codon.max_gap_frac"))
QC_DROP_STOP = bool(need(POST, "codon.drop_internal_stop"))

for p in [PROT_DIR, CDS_DIR, TARGETS]:
    if not os.path.exists(p):
        sys.exit(f"[ERR] 缺少输入: {p}")
Path(CODON_DIR).mkdir(parents=True, exist_ok=True)
Path(RPT_DIR).mkdir(parents=True, exist_ok=True)

gene2seq = {}
dup_warn = 0
file_cnt = 0
seq_cnt  = 0

ACC_VER    = re.compile(r'^([A-Za-z]{2}_[0-9]+)\.([0-9]+)$')
STRIP_VER  = lambda s: re.sub(r'\.[0-9]+$', '', s)

def paired_alias(core: str) -> str | None:
    if core.startswith("XP_"): return "XM_" + core[3:]
    if core.startswith("XM_"): return "XP_" + core[3:]
    if core.startswith("NP_"): return "NM_" + core[3:]
    if core.startswith("NM_"): return "NP_" + core[3:]
    return None

# 载入 CDS：索引原始号、去版本号、以及 XP<->XM / NP<->NM 配对别名
for fp in Path(CDS_DIR).iterdir():
    if not fp.is_file():
        continue
    file_cnt += 1
    for gid_raw, seq in read_fasta(fp):
        acc = normalize_gid(gid_raw).upper()
        if not acc:
            continue
        keys = set()
        keys.add(acc)
        core = STRIP_VER(acc)
        keys.add(core)
        alias = paired_alias(core)
        if alias:
            keys.add(alias)
        m = ACC_VER.match(acc)
        if m and alias:
            alias_with_ver = f"{alias}.{m.group(2)}"
            keys.add(alias_with_ver)
        for k in keys:
            if k in gene2seq:
                dup_warn += 1
                continue
            gene2seq[k] = seq
            seq_cnt += 1

sys.stderr.write(f"[INFO] 载入 CDS 数据：文件 {file_cnt} 个，索引键 {seq_cnt} 条（含去版本与成对别名），重复 {dup_warn} 次被忽略。\n")

qc_lines   = ["OG\tn_taxa\tnon_gap_len\tlen_nt\tpassed"]
drop_lines = ["OG\treason"]

def pal2nal_convert(og, aa_path, genes):
    tmp_cds = Path(CODON_DIR) / f"{og}.cds.tmp.fasta"
    with tmp_cds.open("w", encoding="utf-8") as w:
        for gid in genes:
            seq = gene2seq.get(gid, "")
            if not seq:
                g = gid
                if g.startswith(("XM_", "NM_")):
                    g_alt = ("XP_" + g[3:]) if g.startswith("XM_") else ("NP_" + g[3:])
                    seq = gene2seq.get(g_alt, "") or gene2seq.get(STRIP_VER(g_alt), "")
                elif g.startswith(("XP_", "NP_")):
                    g_alt = ("XM_" + g[3:]) if g.startswith("XP_") else ("NM_" + g[3:])
                    seq = gene2seq.get(g_alt, "") or gene2seq.get(STRIP_VER(g_alt), "")
            if not seq:
                tmp_cds.unlink(missing_ok=True)
                return None, f"missing_cds:{gid}"
            w.write(f">{gid}\n{seq}\n")

    out_codon = Path(CODON_DIR) / f"{og}.codon.fasta"
    cmd = [PAL2NAL, str(aa_path), str(tmp_cds), "-output", "paml", "-nogap"]
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        tmp_cds.unlink(missing_ok=True)
        return None, f"pal2nal_failed:{(e.stderr or e.stdout or '').strip()[:200]}"
    tmp_cds.unlink(missing_ok=True)

    txt = res.stdout.strip().splitlines()
    fasta, cur_id = [], None
    for ln in txt:
        if not ln.strip():
            continue
        if ln[0].isalnum() or ln[0] in "_.:-":
            parts = ln.strip().split()
            cur_id = parts[0]
            fasta.append(f">{cur_id}")
            if len(parts) > 1:
                fasta.append("".join(parts[1:]))
        else:
            if cur_id is None: continue
            fasta.append(ln.strip())
    with open(out_codon, "w", encoding="utf-8") as w:
        w.write("\n".join(fasta) + "\n")
    return str(out_codon), None

ok, tot = 0, 0

for og in [ln.strip() for ln in open(TARGETS, encoding="utf-8") if ln.strip()]:
    tot += 1
    aa_path = Path(PROT_DIR) / f"{og}.trim.faa"
    if not aa_path.exists():
        drop_lines.append(f"{og}\tmissing_aa"); continue

    gids = []
    with aa_path.open("r", encoding="utf-8", errors="replace") as r:
        for ln in r:
            if ln.startswith(">"):
                gids.append(ln[1:].strip().split()[0].upper())

    if len(gids) < QC_MIN_TAXA:
        drop_lines.append(f"{og}\tmin_taxa_fail"); continue

    out_codon, err = pal2nal_convert(og, aa_path, gids)
    # —— 此处为本次唯一修改（修复语法）——
    if err:
        drop_lines.append(f"{og}\t{err}")
        continue

    nseq, nongap = fasta_len_no_gap(out_codon)
    if nseq < QC_MIN_TAXA:
        drop_lines.append(f"{og}\tmin_taxa_fail"); Path(out_codon).unlink(missing_ok=True); continue

    first_seq = ""
    with open(out_codon, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith(">"):
                if first_seq: break
                continue
            if not first_seq:
                first_seq = ln.strip()
    # 保持原逻辑（不动判定/倍数）
    len_nt = len(first_seq.replace("-", "").replace("?", "")) * 3

    seqs, cur = [], []
    with open(out_codon, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else:
                cur.append(ln.strip())
        if cur: seqs.append("".join(cur))
    align_len = len(seqs[0]) if seqs else 0
    gaps = sum(s.count("-") + s.count("?") for s in seqs)
    gap_frac = gaps / max(1, align_len * len(seqs))
    if gap_frac > QC_MAX_GAP:
        drop_lines.append(f"{og}\tmax_gap_frac_fail"); Path(out_codon).unlink(missing_ok=True); continue

    if QC_DROP_STOP and any("*" in s[:-1] for s in seqs):
        drop_lines.append(f"{og}\tinternal_stop"); Path(out_codon).unlink(missing_ok=True); continue

    ok += 1
    qc_lines.append(f"{og}\t{nseq}\t{nongap}\t{len_nt}\tTRUE")

qc_path   = os.path.join(RPT_DIR, "codon_qc.tsv")
drop_path = os.path.join(RPT_DIR, "codon_drop.tsv")
with open(qc_path, "w", encoding="utf-8") as w:  w.write("\n".join(qc_lines) + "\n")
with open(drop_path, "w", encoding="utf-8") as w: w.write("\n".join(drop_lines) + "\n")
sys.stderr.write(f"[DONE] Codon 构建完成：通过 {ok} / 目标 {tot} → {qc_path}\n")

