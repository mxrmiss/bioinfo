#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
11_joint_integration.py — 联合整合（输出 joint_dir）
"""
# ===== APhylo utils (config, logging, sentinels) =====
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
# ===== utils end =====

from pathlib import Path

def main():
    cfg = load_config()
    paths = cfg["paths"]
    logs_dir = Path(paths["logs_dir"]); LOG_FILE = logs_dir / "11_joint_integration.log"
    log = get_logger("aphylo.11", LOG_FILE)
    banner(log, "APhylo 11 — 联合整合")

    codeml = need_dir(Path(paths["codeml_agg_dir"]), "codeml 汇总目录")
    cafe   = need_dir(Path(paths["cafe_agg_dir"]), "cafe 汇总目录")
    out    = ensure_dir(Path(paths["joint_dir"]))

    D = need_file(codeml/"D_fdr_genes.tsv", "codeml D_fdr_genes.tsv")
    C = need_file(cafe/"cafe_significant_families.tsv", "cafe_significant_families.tsv")

    D_lines = [l for l in D.read_text(encoding="utf-8").splitlines() if l and not l.startswith("OG")]
    C_lines = [l for l in C.read_text(encoding="utf-8").splitlines() if l and not l.startswith("model")]

    (out/"integration_example.tsv").write_text(
        "note\tcount_codeml\tcount_cafe\n" + f"demo\t{len(D_lines)}\t{len(C_lines)}\n",
        encoding="utf-8"
    )

    write_done(out/".integration.done")
    log.info("[DONE] 联合整合完成")

if __name__ == "__main__":
    try: main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)
