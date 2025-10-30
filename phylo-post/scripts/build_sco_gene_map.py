#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, csv
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
    try:
        results_txt = os.path.normpath(need(cfg_post,"inputs.results_txt"))  # 改
        if os.path.exists(results_txt):
            resdir = Path(Path(results_txt).read_text(encoding="utf-8").strip())
            ogtsv = resdir/"Orthogroups"/"Orthogroups.tsv"
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
        sys.stderr.write(f"[WARN] 读取 results_txt 异常，转入自动扫描：{e}\n")
    base_dir = Path("..")/"phylo"/"results"/"orthofinder"   # 保留兜底
    latest = find_latest_orthogroups(base_dir)
    if latest and latest.exists():
        sys.stderr.write(f"[OK] 自动选用最新 Orthogroups.tsv：{latest}\n")
        return latest
    sys.exit(f"[ERR] 未能定位 Orthogroups.tsv，请检查 results_txt 或 {base_dir}")

CFG_PATH = str(Path(__file__).resolve().parent.parent / "config.yaml")  # 改
if not os.path.exists(CFG_PATH): sys.exit(f"[ERR] 配置缺失: {CFG_PATH}")
CFG = yaml.safe_load(open(CFG_PATH, encoding="utf-8")) or {}
POST = CFG.get("post", {}) or {}

OUTDIR = need(POST, "outdir")
ogtsv = resolve_orthogroups(POST)

outp = Path(OUTDIR)/"codon_sco"/"og_species_to_gene.tsv"
outp.parent.mkdir(parents=True, exist_ok=True)

with ogtsv.open("r", encoding="utf-8", errors="replace", newline="") as f, \
     outp.open("w", encoding="utf-8") as w:
    rdr = csv.reader(f, delimiter="\t")
    header = next(rdr)
    species = header[1:]
    w.write("OG\tspecies\tgene_id\n")
    for row in rdr:
        if not row: continue
        og = row[0].strip()
        cells = row[1:]
        for sp, cell in zip(species, cells):
            cell = (cell or "").strip()
            if not cell: continue
            gid = cell.split(",")[0].strip()
            w.write(f"{og}\t{sp}\t{gid}\n")

sys.stderr.write(f"[DONE] 写出映射: {outp}\n")

