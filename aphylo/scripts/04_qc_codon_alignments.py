#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
04_qc_codon_alignments.py — codon MSA 质量控制（按 config.codon 阈值）
"""
from __future__ import annotations
import sys, io, logging, subprocess, re
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable
import yaml

DEFAULT_CONFIG = "config.yaml"

def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(config_path)
    if not p.exists(): raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub: cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True); return p
def need_dir(p: Path, what: str):
    p = Path(p)
    if not p.is_dir(): raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p
def need_file(p: Path, what: str):
    p = Path(p)
    if not p.is_file(): raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name); lg.setLevel(level); lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt); fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout);     sh.setFormatter(fmt); sh.setLevel(level)
    lg.addHandler(fh); lg.addHandler(sh)
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout = _Flush(sys.stdout); sys.stderr = _Flush(sys.stderr)
    return lg

def banner(logger: logging.Logger, text: str):
    bar = "=" * max(10, len(text)+2); logger.info(bar); logger.info(f" {text} "); logger.info(bar)

def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    recs=[]; name=None; seq=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None: recs.append((name, "".join(seq)))
            name=line[1:].strip(); seq=[]
        else:
            seq.append(line.strip())
    if name is not None: recs.append((name, "".join(seq)))
    return recs

def write_fasta(path: Path, recs: Iterable[Tuple[str,str]]):
    with path.open("w", encoding="utf-8") as w:
        for h, s in recs:
            w.write(f">{h}\n")
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + "\n")

STD_CODE = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
"CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
"AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
"GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
"GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

def translate(nt: str) -> str:
    aa=[]; nt=nt.upper()
    for i in range(0, len(nt), 3):
        cod = nt[i:i+3]
        if len(cod)<3 or "-" in cod or "N" in cod: aa.append("X"); continue
        aa.append(STD_CODE.get(cod, "X"))
    return "".join(aa)

def avg_gap_ratio(seq: str) -> float:
    n = len(seq); gaps = seq.count("-")
    return gaps / max(1, n)

def internal_stop_count(nt: str) -> int:
    aa = translate(nt)
    if len(aa) == 0: return 0
    return aa[:-1].count("*")  # exclude terminal

def main():
    cfg = load_config()
    paths = cfg["paths"]; codon_cfg = cfg.get("codon", {})
    logs_dir = Path(paths["logs_dir"]); LOG_FILE = logs_dir / "04_qc_codon_alignments.log"
    log = get_logger("aphylo.04", LOG_FILE)
    banner(log, "APhylo 04 — codon 对齐质量控制")

    # 仅路径改动：bt_dir 固定为 results/03_codon
    bt_dir = need_dir(Path("results/03_codon"), "codon 输入目录")
    in_dir = bt_dir / "codon_msa"
    out_dir = ensure_dir(Path(paths["codon_qc_dir"]))

    min_taxa = int(codon_cfg.get("min_taxa", 0))
    min_len  = int(codon_cfg.get("min_codon_len", 150))
    max_gap  = float(codon_cfg.get("max_gap_frac", 0.5))
    stop_action = str(codon_cfg.get("stopcodon_action", "fail")).lower()

    codon_files = sorted(in_dir.glob("OG*.codon.fna"))
    if not codon_files: raise FileNotFoundError(f"[ERR] 未找到 codon MSA：{in_dir}")

    lines_qc=["OG\ttaxa\tnt_len\tinternal_stops\tavg_gap\tpass"]
    lines_drop=["OG\treason"]
    keep_ogs=[]
    for fa in codon_files:
        og = re.match(r"^(OG\d+)", fa.name).group(1)
        recs = read_fasta(fa)
        taxa = len(recs)
        seqlen = len(recs[0][1]) if recs else 0
        stops = sum(internal_stop_count(s) for _,s in recs)
        gaps  = sum(avg_gap_ratio(s) for _,s in recs)/max(1,taxa)
        stop_ok = (stops==0) if stop_action.startswith("fail") else True
        passed = (taxa>=min_taxa) and (seqlen//3 >= min_len) and stop_ok and (gaps <= max_gap)
        lines_qc.append("\t".join([og, str(taxa), str(seqlen), str(stops), f"{gaps:.3f}", str(int(passed))]))
        if passed: keep_ogs.append(og)
        else:
            reasons=[]
            if taxa < min_taxa: reasons.append("low_taxa")
            if seqlen//3 < min_len: reasons.append("short")
            if not stop_ok: reasons.append("internal_stop")
            if gaps > max_gap: reasons.append("high_gap")
            lines_drop.append("\t".join([og, ",".join(reasons)]))

    (out_dir/"codon_qc.tsv").write_text("\n".join(lines_qc) + "\n", encoding="utf-8")
    (out_dir/"codon_drop.tsv").write_text("\n".join(lines_drop) + "\n", encoding="utf-8")
    (out_dir/"keep_og.list").write_text("\n".join(keep_ogs) + "\n", encoding="utf-8")

    write_done(out_dir/".qc.done")
    log.info(f"[DONE] QC 完成：保留 {len(keep_ogs)} / 总计 {len(codon_files)}")

if __name__ == "__main__":
    try: main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)

