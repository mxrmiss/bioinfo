#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys
from pathlib import Path
import yaml

def need(cfg, key):
    cur = cfg
    for k in key.split("."):
        if not isinstance(cur, dict) or k not in cur:
            sys.exit(f"[ERR] config 缺少 post.{key}")
        cur = cur[k]
    return cur

def find_latest_orthogroups(base_dir: Path) -> Path:
    pats = list(base_dir.glob("run_of_r*/Results_*/Orthogroups/Orthogroups.tsv"))
    if not pats: return None
    pats.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return pats[0]

def resolve_orthogroups(cfg_post: dict) -> Path:
    # 改：直接用配置提供的 results_txt（不再强拼 ../phylo）
    try:
        results_txt = os.path.normpath(need(cfg_post, "inputs.results_txt"))
        if os.path.exists(results_txt):
            resdir = Path(Path(results_txt).read_text(encoding="utf-8").strip())
            ogtsv = resdir / "Orthogroups" / "Orthogroups.tsv"
            if ogtsv.exists():
                sys.stderr.write(f"[INFO] 选用 results_txt 指向的 Orthogroups.tsv：{ogtsv}\n")
                return ogtsv
            else:
                sys.stderr.write(f"[WARN] results_txt 指向目录无 Orthogroups.tsv：{ogtsv}\n")
        else:
            sys.stderr.write(f"[WARN] results_txt 文件不存在：{results_txt}\n")
    except SystemExit:
        pass
    except Exception as e:
        sys.stderr.write(f"[WARN] 读取 results_txt 发生异常，将尝试自动扫描：{e}\n")

    # 兜底自动扫描（保留原逻辑，仍在 ../phylo 下）
    base_dir = Path("..") / "phylo" / "results" / "orthofinder"
    latest = find_latest_orthogroups(base_dir)
    if latest and latest.exists():
        sys.stderr.write(f"[OK] 自动选用最新 Orthogroups.tsv：{latest}\n")
        return latest

    sys.exit(f"[ERR] 未能定位 Orthogroups.tsv，请检查 results_txt 或 {base_dir}")

# 读取配置（改：本地）
CFG_PATH = str(Path(__file__).resolve().parent.parent / "config.yaml")
if not os.path.exists(CFG_PATH):
    sys.exit(f"[ERR] 未找到配置: {CFG_PATH}")
CFG = yaml.safe_load(open(CFG_PATH, encoding="utf-8")) or {}
POST = CFG.get("post", {}) or {}

OUTDIR = need(POST, "outdir")
AA_DIR = os.path.normpath(need(POST, "inputs.aa_align_dir"))  # 改：不再强拼 ../phylo
CODON_DIR = os.path.join(OUTDIR, "codon_sco")
TARGETS_DIR = os.path.join(OUTDIR, "targets")
SCO_KEEP = os.path.join("..", "phylo", "results", "orthogroups", "sco.keep.list")  # 保留
CAFE_SIG = POST.get("selection", {}).get("cafe_significant_tsv", None)

if not os.path.exists(AA_DIR):
    sys.exit(f"[ERR] 缺少输入路径: {AA_DIR}")

ogtsv = resolve_orthogroups(POST)

all_ogs = set()
with ogtsv.open("r", encoding="utf-8", errors="replace") as f:
    header = True
    for line in f:
        if header: header = False; continue
        if not line.strip(): continue
        all_ogs.add(line.split("\t", 1)[0].strip())

if not os.path.exists(SCO_KEEP):
    sys.exit(f"[ERR] 缺少严格SCO清单（由 phylo 流生成）: {SCO_KEEP}")
sco = {ln.strip() for ln in open(SCO_KEEP, encoding="utf-8") if ln.strip()}

cafe = set()
if CAFE_SIG and os.path.exists(CAFE_SIG):
    cafe = {ln.strip().split("\t")[0] for ln in open(CAFE_SIG, encoding="utf-8") if ln.strip()}
else:
    sys.stderr.write("[INFO] 未提供 CAFE 显著表：仅使用 SCO\n")

existing = set()
if os.path.isdir(CODON_DIR):
    for fp in Path(CODON_DIR).glob("*.codon.fasta"):
        og = fp.name.replace(".codon.fasta", "")
        existing.add(og)

cand = (sco | cafe) & all_ogs
final = sorted(cand - existing)

Path(TARGETS_DIR).mkdir(parents=True, exist_ok=True)
outf = os.path.join(TARGETS_DIR, "og_targets.list")
open(outf, "w", encoding="utf-8").write("\n".join(final) + "\n")
sys.stderr.write(f"[DONE] 目标OG={len(final)} (SCO={len(sco)}, CAFE={len(cafe)}, existing={len(existing)}) → {outf}\n")

