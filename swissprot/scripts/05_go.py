#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
05_go.py
================================================================================
[功能]
基于 Swiss-Prot besthit 注释结果（results/04_annot/{REL}/{species}/swissprot_annotation.tsv），
为每个物种追加 GO 注释（GO IDs），并输出“最终一张大表” annot.tsv + qc.tsv。

[核心策略（皇上已确认）]
1) 只允许一个 GO 缓存读取路径：
   data/db/go_cache/{REL}/acc2go.tsv.gz
2) 运行时先检查缓存：
   - 若缓存不存在 -> 自动构建缓存（只保留 Swiss-Prot reviewed accession）
   - 缓存构建完成后才进入逐物种 GO 注释
   - 注释阶段只从缓存读取；若缓存缺失/损坏 -> 直接报错
3) 不做旧目录兼容，不做白名单；config.yaml 的 species_prefixes 即为要做 GO 的物种列表
4) 产物命名短：
   - 最终表：annot.tsv（只有这一张表最常用）
   - QC：qc.tsv
   （可选：go.tsv，默认不输出）
5) 屏幕流式输出（阶段/进度可见）
6) 不接受命令行参数：参数集中在脚本顶部 + config.yaml

[输入依赖]
- results/04_annot/{REL}/{species}/swissprot_annotation.tsv
- data/db/raw/uniprot_sprot.fasta（由 01_db.py 解压生成；若缺失则尝试从 .gz 解压）
- data/db/raw/goa_uniprot_all.gaf.gz（你已用 Aspera 下载）

[输出]
- 每物种：
  results/05_go/{REL}/{species}/annot.tsv   （最终大表：Swiss-Prot 注释 + GO）
  results/05_go/{REL}/{species}/qc.tsv
  （可选）results/05_go/{REL}/{species}/go.tsv

- 缓存（唯一读取路径）：
  data/db/go_cache/{REL}/acc2go.tsv.gz
  data/db/go_cache/{REL}/manifest.tsv

================================================================================
"""

from __future__ import annotations

import os
import re
import gzip
import time
import shutil
from pathlib import Path
from typing import Dict, Any, Set, Tuple

# 让“项目根目录”进入 sys.path，保证可以 import scripts._common
import sys
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts._common import (
    find_project_root,
    read_yaml_config,
    parse_rel_from_01_db,
    get_species_list,
    get_only_species_from_env,
    ensure_dir,
)

# ----------------------------- 参数区（皇上可改：不接受命令行参数） -----------------------------

# 是否额外输出短表 go.tsv（默认 False：避免“多张表容易混淆”）
WRITE_GO_TSV = False

# 01_db.py 约定的 Swiss-Prot 解压文件名
SPROT_FASTA = "data/db/raw/uniprot_sprot.fasta"
SPROT_DL_GZ = "data/db/raw/uniprot_sprot.fasta.gz"     # 01_db.py 下载的 reviewed fasta.gz（同名）

# GOA GAF（默认路径，可在 config.yaml 的 go.gaf_gz 覆盖）
DEFAULT_GOA_GAF_GZ = "data/db/raw/goa_uniprot_all.gaf.gz"

# 缓存（唯一读取路径）
CACHE_REL_DIR = "data/db/go_cache/{REL}"
CACHE_FILE = "data/db/go_cache/{REL}/acc2go.tsv.gz"
CACHE_MANIFEST = "data/db/go_cache/{REL}/manifest.tsv"

# 05 输出目录
OUT_BASE = "results/05_go/{REL}/{SPECIES}"

# 最终大表文件名（皇上指定）
FINAL_TABLE_NAME = "annot.tsv"

# ----------------------------- 小工具 -----------------------------

def _cfg_get_go(cfg: Dict[str, Any]) -> Dict[str, Any]:
    d = cfg.get("go", {})
    return d if isinstance(d, dict) else {}

def _read_tsv_header_map(header_line: str) -> Dict[str, int]:
    cols = header_line.rstrip("\n").split("\t")
    return {c: i for i, c in enumerate(cols)}

def _stat_fingerprint(p: Path) -> Tuple[int, int]:
    """
    返回 (size, mtime_int)，用于快速一致性检查（避免每次算大文件 md5）
    """
    st = p.stat()
    return int(st.st_size), int(st.st_mtime)

def _gunzip_if_needed(src_gz: Path, dst: Path) -> None:
    """
    兜底：若 uniprot_sprot.fasta 不存在，则尝试从 gz 解压。
    """
    if dst.is_file() and dst.stat().st_size > 0:
        return
    if (not src_gz.is_file()) or src_gz.stat().st_size == 0:
        raise FileNotFoundError(f"Swiss-Prot FASTA missing: {dst} (and gz missing: {src_gz})")

    ensure_dir(dst.parent)
    print(f"[05_go] gunzip Swiss-Prot FASTA: {src_gz} -> {dst}", flush=True)

    with gzip.open(src_gz.as_posix(), "rb") as fin, dst.open("wb") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)

    if (not dst.is_file()) or dst.stat().st_size == 0:
        raise RuntimeError(f"Failed to gunzip: {dst}")

# ----------------------------- 缓存构建（只保留 Swiss-Prot reviewed accession） -----------------------------

def extract_sprot_accessions(fasta_path: Path) -> Set[str]:
    """
    从 uniprot_sprot.fasta 提取 accession 集合（Swiss-Prot reviewed）。
    典型 header：
      >sp|P27987|IP3KB_HUMAN ...
    也兼容：
      >P27987 ...
    """
    accs: Set[str] = set()
    n_hdr = 0

    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            n_hdr += 1
            s = line[1:].strip()
            if not s:
                continue

            if "|" in s:
                parts = s.split("|")
                if len(parts) >= 2 and parts[1]:
                    acc = parts[1]
                else:
                    acc = s.split()[0]
            else:
                acc = s.split()[0]

            acc = re.sub(r"\.[0-9]+$", "", acc)
            if acc:
                accs.add(acc)

            if n_hdr % 200000 == 0:
                print(f"[05_go] parsed {n_hdr} FASTA headers, unique acc={len(accs)}", flush=True)

    print(f"[05_go] Swiss-Prot accessions: {len(accs)}", flush=True)
    return accs

def build_acc2go_cache(rel: str, root: Path, cfg: Dict[str, Any]) -> None:
    """
    构建缓存：
      data/db/go_cache/{REL}/acc2go.tsv.gz

    缓存格式（每 accession 一行，紧凑好读）：
      swissprot_accession \t go_ids \t n_go
    """
    cfg_go = _cfg_get_go(cfg)

    cache_dir = root / CACHE_REL_DIR.format(REL=rel)
    cache_file = root / CACHE_FILE.format(REL=rel)
    manifest = root / CACHE_MANIFEST.format(REL=rel)

    ensure_dir(cache_dir)

    sprot_fa = root / SPROT_FASTA
    sprot_gz = root / SPROT_DL_GZ
    goa_gaf_gz = root / str(cfg_go.get("gaf_gz", DEFAULT_GOA_GAF_GZ))

    if (not sprot_fa.is_file()) or sprot_fa.stat().st_size == 0:
        _gunzip_if_needed(sprot_gz, sprot_fa)

    if (not goa_gaf_gz.is_file()) or goa_gaf_gz.stat().st_size == 0:
        raise FileNotFoundError(f"GOA GAF missing/empty: {goa_gaf_gz}")

    print(f"[05_go] build cache for REL={rel}", flush=True)
    print(f"[05_go] Swiss-Prot FASTA: {sprot_fa}", flush=True)
    print(f"[05_go] GOA GAF: {goa_gaf_gz}", flush=True)
    print(f"[05_go] CACHE_OUT: {cache_file}", flush=True)

    acc_set = extract_sprot_accessions(sprot_fa)

    acc2go: Dict[str, Set[str]] = {}
    n_total = 0
    n_keep = 0
    last_print = time.time()

    with gzip.open(goa_gaf_gz.as_posix(), "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("!"):
                continue
            n_total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            db = parts[0]
            acc = parts[1]
            go_id = parts[4]

            if db != "UniProtKB":
                continue

            acc = re.sub(r"\.[0-9]+$", "", acc)
            if acc not in acc_set:
                continue

            if not go_id or not go_id.startswith("GO:"):
                continue

            s = acc2go.get(acc)
            if s is None:
                s = set()
                acc2go[acc] = s
            s.add(go_id)
            n_keep += 1

            if time.time() - last_print >= 10:
                print(f"[05_go] scan GAF: lines={n_total:,} kept_records={n_keep:,} acc_with_go={len(acc2go):,}", flush=True)
                last_print = time.time()

    print(f"[05_go] scan GAF done: lines={n_total:,} kept_records={n_keep:,} acc_with_go={len(acc2go):,}", flush=True)

    tmp = cache_file.with_suffix(cache_file.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    with gzip.open(tmp.as_posix(), "wt", encoding="utf-8") as out:
        out.write("swissprot_accession\tgo_ids\tn_go\n")
        for acc in sorted(acc2go.keys()):
            go_ids = sorted(acc2go[acc])
            out.write(acc + "\t" + ";".join(go_ids) + "\t" + str(len(go_ids)) + "\n")

    tmp.replace(cache_file)

    if (not cache_file.is_file()) or cache_file.stat().st_size == 0:
        raise RuntimeError(f"cache build failed: {cache_file}")

    sprot_size, sprot_mtime = _stat_fingerprint(sprot_fa)
    goa_size, goa_mtime = _stat_fingerprint(goa_gaf_gz)
    cache_size, cache_mtime = _stat_fingerprint(cache_file)

    manifest_lines = [
        f"REL\t{rel}",
        f"SPROT_FASTA\t{sprot_fa.as_posix()}",
        f"SPROT_FASTA_SIZE\t{sprot_size}",
        f"SPROT_FASTA_MTIME\t{sprot_mtime}",
        f"GOA_GAF_GZ\t{goa_gaf_gz.as_posix()}",
        f"GOA_GAF_GZ_SIZE\t{goa_size}",
        f"GOA_GAF_GZ_MTIME\t{goa_mtime}",
        f"CACHE\t{cache_file.as_posix()}",
        f"CACHE_SIZE\t{cache_size}",
        f"CACHE_MTIME\t{cache_mtime}",
        f"BUILD_TIME\t{time.strftime('%Y-%m-%d %H:%M:%S')}",
    ]
    manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print(f"[05_go] cache built OK: {cache_file}", flush=True)
    print(f"[05_go] manifest: {manifest}", flush=True)

def check_or_build_cache(rel: str, root: Path, cfg: Dict[str, Any]) -> Path:
    """
    每次运行脚本都会做一次“轻量检查”：
    - cache 文件存在且非空
    - manifest 存在且基本完整
    - manifest 记录的 GOA 文件 size/mtime 与当前一致（快速检查）
    若 cache 缺失 -> 自动 build
    若 cache 存在但不一致/损坏 -> 报错（让你手动处理：删 cache 触发重建）
    """
    cache_file = root / CACHE_FILE.format(REL=rel)
    manifest = root / CACHE_MANIFEST.format(REL=rel)

    if (not cache_file.is_file()) or cache_file.stat().st_size == 0 or (not manifest.is_file()):
        print(f"[05_go] cache missing -> build", flush=True)
        build_acc2go_cache(rel, root, cfg)
        if (not cache_file.is_file()) or cache_file.stat().st_size == 0:
            raise RuntimeError(f"cache build failed: {cache_file}")
        if not manifest.is_file():
            raise RuntimeError(f"cache manifest missing after build: {manifest}")
        return cache_file

    kv: Dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        k, *rest = line.split("\t", 1)
        if rest:
            kv[k] = rest[0]

    goa_path = kv.get("GOA_GAF_GZ", "")
    goa_size = kv.get("GOA_GAF_GZ_SIZE", "")
    goa_mtime = kv.get("GOA_GAF_GZ_MTIME", "")

    if not goa_path or not goa_size or not goa_mtime:
        raise RuntimeError(f"cache manifest incomplete: {manifest}")

    goa_gaf = Path(goa_path)
    if (not goa_gaf.is_file()) or goa_gaf.stat().st_size == 0:
        raise RuntimeError(f"GOA file recorded in manifest missing/empty: {goa_gaf}")

    cur_size, cur_mtime = _stat_fingerprint(goa_gaf)
    if str(cur_size) != str(goa_size) or str(cur_mtime) != str(goa_mtime):
        raise RuntimeError(
            "GOA GAF changed since cache build. Please delete cache to rebuild:\n"
            f"  cache: {cache_file}\n"
            f"  manifest: {manifest}\n"
            f"  goa_gaf: {goa_gaf}\n"
        )

    print(f"[05_go] cache OK: {cache_file}", flush=True)
    return cache_file

def load_cache_to_memory(cache_file: Path) -> Dict[str, Tuple[str, int]]:
    """
    读取缓存到内存：acc -> (go_ids_str, n_go)
    """
    acc2: Dict[str, Tuple[str, int]] = {}
    n = 0
    with gzip.open(cache_file.as_posix(), "rt", encoding="utf-8", errors="replace") as f:
        header = next(f, None)
        if header is None:
            raise RuntimeError(f"cache empty: {cache_file}")
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc = parts[0]
            go_ids = parts[1]
            try:
                n_go = int(parts[2])
            except Exception:
                n_go = 0
            acc2[acc] = (go_ids, n_go)
            n += 1
            if n % 200000 == 0:
                print(f"[05_go] loaded cache rows: {n:,}", flush=True)
    print(f"[05_go] loaded cache rows total: {n:,}", flush=True)
    return acc2

# ----------------------------- 逐物种 GO 注释（直接输出最终大表 annot.tsv） -----------------------------

def annotate_one_species(root: Path, rel: str, species_id: str, acc2go: Dict[str, Tuple[str, int]]) -> None:
    in_annot = root / "results" / "04_annot" / rel / species_id / "swissprot_annotation.tsv"
    if (not in_annot.is_file()) or in_annot.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty input: {in_annot}")

    outdir = root / OUT_BASE.format(REL=rel, SPECIES=species_id)
    ensure_dir(outdir)

    out_final = outdir / FINAL_TABLE_NAME
    out_qc = outdir / "qc.tsv"
    out_go = outdir / "go.tsv"

    print(f"[05_go] species={species_id}", flush=True)
    print(f"[05_go] IN={in_annot}", flush=True)
    print(f"[05_go] OUT_FINAL={out_final}", flush=True)
    if WRITE_GO_TSV:
        print(f"[05_go] OUT_GO={out_go}", flush=True)

    n_total = 0
    n_has_go = 0
    n_no_go = 0

    with in_annot.open("r", encoding="utf-8", errors="replace") as fin, \
         out_final.open("w", encoding="utf-8", newline="") as ffinal:

        header = fin.readline()
        if not header:
            raise RuntimeError(f"empty file: {in_annot}")
        hmap = _read_tsv_header_map(header)

        if "protein_id" not in hmap or "swissprot_accession" not in hmap:
            raise RuntimeError(f"required columns missing in {in_annot}: protein_id / swissprot_accession")

        i_pid = hmap["protein_id"]
        i_acc = hmap["swissprot_accession"]

        # 写最终表表头：原表头 + go_ids + n_go
        orig_cols = header.rstrip("\n").split("\t")
        ffinal.write("\t".join(orig_cols + ["go_ids", "n_go"]) + "\n")

        # 可选：同时写短表 go.tsv（protein -> go）
        fgo = None
        if WRITE_GO_TSV:
            fgo = out_go.open("w", encoding="utf-8", newline="")
            fgo.write("protein_id\tswissprot_accession\tgo_ids\tn_go\n")

        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_pid, i_acc):
                continue

            pid = parts[i_pid]
            acc = parts[i_acc]
            acc = re.sub(r"\.[0-9]+$", "", acc)

            go_ids, n_go = ("", 0)
            hit = acc2go.get(acc)
            if hit is not None:
                go_ids, n_go = hit

            n_total += 1
            if n_go > 0:
                n_has_go += 1
            else:
                n_no_go += 1

            # 写最终表一行：原行 + GO
            ffinal.write("\t".join(parts + [go_ids, str(n_go)]) + "\n")

            # 可选写 go.tsv
            if fgo is not None:
                fgo.write(f"{pid}\t{acc}\t{go_ids}\t{n_go}\n")

            if n_total % 20000 == 0:
                print(f"[05_go] {species_id}: processed={n_total:,} has_go={n_has_go:,} no_go={n_no_go:,}", flush=True)

        if fgo is not None:
            fgo.close()

    out_qc.write_text(
        "species_id\tREL\tN_total\tN_has_go\tN_no_go\tHAS_GO_rate\n"
        f"{species_id}\t{rel}\t{n_total}\t{n_has_go}\t{n_no_go}\t"
        f"{(n_has_go / n_total if n_total else 0):.6f}\n",
        encoding="utf-8"
    )

    print(f"[05_go] DONE {species_id}: total={n_total:,} has_go={n_has_go:,} out={out_final}", flush=True)

# ----------------------------- main -----------------------------

def main() -> None:
    root = find_project_root()
    os.chdir(root)

    cfg = read_yaml_config(root)
    rel = parse_rel_from_01_db(root)

    # 物种列表：默认全跑；若设置环境变量 SPECIES_ID 则只跑一个（和 02/03/04 对齐）
    species_all = get_species_list(cfg)
    only = get_only_species_from_env()
    if only:
        species = [only]
        if only not in species_all:
            raise ValueError(f"SPECIES_ID={only} not in config species_prefixes")
    else:
        species = species_all

    print(f"[05_go] REL={rel}", flush=True)
    print(f"[05_go] N_species={len(species)}", flush=True)

    cache_file = check_or_build_cache(rel, root, cfg)

    if (not cache_file.is_file()) or cache_file.stat().st_size == 0:
        raise RuntimeError(f"cache missing/empty: {cache_file}")

    acc2go = load_cache_to_memory(cache_file)

    for sid in species:
        annotate_one_species(root, rel, sid, acc2go)

    print("[05_go] ALL DONE", flush=True)

if __name__ == "__main__":
    main()

