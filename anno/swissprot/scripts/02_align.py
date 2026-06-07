#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import re
import sys
import time
import subprocess
from pathlib import Path
from typing import Optional, List

# 让“项目根目录”进入 sys.path，保证可以 import scripts._common
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts._common import (
    find_project_root, read_yaml_config, parse_rel_from_01_db,
    get_species_list, get_only_species_from_env, resolve_db_dmnd,
    ensure_dir, count_fasta_seqs, md5sum_file, eprint
)

# ----------------------------- 流式输出降噪策略 -----------------------------
# 目标：
# 1) log 文件保留 diamond 原始完整输出（不丢任何行）
# 2) 屏幕端只显示关键信息，避免同一套 stage 日志反复刷屏

_DIAMOND_STAGE_PREFIXES = (
    "Building reference seed array",
    "Building query seed array",
    "Computing hash join",
    "Masking low complexity seeds",
    "Searching alignments",
    "Deallocating memory",
)

# 这些行非常重复，但“Processing ...”本身是有价值的进度信息，保留到屏幕
_DIAMOND_KEEP_REGEX = re.compile(r"^Processing query block\b")

# 常见错误/告警关键词：只要命中就强制显示到屏幕
_ERROR_REGEX = re.compile(r"(error|failed|fatal|terminate|exception)", re.IGNORECASE)

def _should_print_to_screen(line: str) -> bool:
    s = line.strip()
    if not s:
        return False
    if _ERROR_REGEX.search(s):
        return True
    if _DIAMOND_KEEP_REGEX.match(s):
        return True
    for p in _DIAMOND_STAGE_PREFIXES:
        if s.startswith(p):
            return False
    # 其余输出（例如 diamond 版本、参数、最后统计等）默认保留
    return True

def stream_command(cmd: List[str], log_path: Path, screen_rate_limit_sec: float = 0.2) -> int:
    """
    运行外部命令：stdout/stderr 合流，逐行写 log。
    屏幕端进行降噪：过滤 diamond 高频重复 stage 行；对进度行做轻度限频。
    """
    ensure_dir(log_path.parent)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("CMD: " + " ".join(cmd) + "\n")
        log.flush()

        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        assert p.stdout is not None
        last_print_t = 0.0
        last_processing = ""

        for line in p.stdout:
            # 永远完整写入 log（不做任何丢弃）
            log.write(line)
            log.flush()

            s = line.rstrip("\n")
            if not _should_print_to_screen(s):
                continue

            # 对 “Processing query block ...” 做轻度限频，避免 chunk 太多时仍然刷屏
            now = time.monotonic()
            if _DIAMOND_KEEP_REGEX.match(s):
                last_processing = s
                if now - last_print_t < screen_rate_limit_sec:
                    continue
                last_print_t = now
                print(s, flush=True)
                continue

            # 其它重要输出直接显示
            print(s, flush=True)

        rc = p.wait()

        # 如果最后一次 processing 被限频吞了，补打一行（让用户知道停在哪）
        if last_processing:
            # 只有在进程结束时，屏幕上最后一行不是这个 processing 才补打
            # 这里不追踪屏幕最后一行内容，用简单策略：结束后再提示一次
            print(last_processing, flush=True)

        return rc

# ----------------------------- 主逻辑 -----------------------------

def run_one(root: Path, cfg: dict, rel: str, species_id: str) -> None:
    threads = int(cfg.get("threads", 1))
    dcfg = cfg.get("diamond", {}) if isinstance(cfg.get("diamond", {}), dict) else {}
    sensitivity = str(dcfg.get("sensitivity", "sensitive")).strip()
    evalue = str(dcfg.get("evalue", "1e-5")).strip()
    max_target_seqs = int(dcfg.get("max_target_seqs", 50))

    query = root / "data" / "query" / f"{species_id}.faa"
    if not query.is_file() or query.stat().st_size == 0:
        raise FileNotFoundError(f"QUERY missing or empty: {query}")

    db_dmnd, db_prefix = resolve_db_dmnd(root, rel)
    if db_dmnd.stat().st_size == 0:
        raise FileNotFoundError(f"DB missing or empty: {db_dmnd}")

    outdir = root / "results" / "02_align" / rel / species_id
    ensure_dir(outdir)
    ensure_dir(root / "logs")

    out_hits = outdir / "hits.diamond.tsv"
    log_path = root / "logs" / f"02_align.{rel}.{species_id}.log"
    mani = outdir / "run_manifest.tsv"

    print(f"[02_align] species={species_id} REL={rel}", flush=True)
    print(f"[02_align] QUERY={query}", flush=True)
    print(f"[02_align] DB_DMND={db_dmnd}", flush=True)
    print(f"[02_align] OUT={out_hits}", flush=True)

    n_query = count_fasta_seqs(query)
    (outdir / "n_query_seqs.txt").write_text(str(n_query) + "\n", encoding="utf-8")
    (outdir / "md5_query.txt").write_text(f"{md5sum_file(query)}  {query.as_posix()}\n", encoding="utf-8")
    (outdir / "db_fileinfo.txt").write_text(f"{db_dmnd.as_posix()}\t{db_dmnd.stat().st_size}\n", encoding="utf-8")

    dv = subprocess.run(["diamond", "version"], capture_output=True, text=True)
    (outdir / "diamond_version.txt").write_text((dv.stdout + dv.stderr).strip() + "\n", encoding="utf-8")

    mani_lines = [
        f"REL\t{rel}",
        f"SPECIES_ID\t{species_id}",
        f"QUERY\t{query.as_posix()}",
        f"DB_DMND\t{db_dmnd.as_posix()}",
        f"DB_PREFIX\t{db_prefix}",
        f"THREADS\t{threads}",
        f"OUT\t{out_hits.as_posix()}",
        f"LOG\t{log_path.as_posix()}",
        f"SENSITIVITY\t{sensitivity}",
        f"EVALUE\t{evalue}",
        f"MAX_TARGET_SEQS\t{max_target_seqs}",
    ]
    mani.write_text("\n".join(mani_lines) + "\n", encoding="utf-8")

    cmd = [
        "diamond", "blastp",
        "--db", db_prefix,
        "--query", query.as_posix(),
        "--out", out_hits.as_posix(),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "stitle",
        "--threads", str(threads),
        "--evalue", evalue,
        "--max-target-seqs", str(max_target_seqs),
    ]
    if sensitivity == "sensitive":
        cmd.append("--sensitive")
    elif sensitivity == "very-sensitive":
        cmd.append("--very-sensitive")
    elif sensitivity == "fast":
        pass
    else:
        raise ValueError(f"Unsupported diamond.sensitivity: {sensitivity}")

    rc = stream_command(cmd, log_path, screen_rate_limit_sec=0.2)
    if rc != 0:
        raise RuntimeError(f"DIAMOND failed for {species_id} (exit={rc}). See {log_path}")

    if not out_hits.is_file() or out_hits.stat().st_size == 0:
        raise RuntimeError(f"HITS missing or empty: {out_hits}")

    print(f"[02_align] DONE species={species_id} hits={out_hits}", flush=True)

def main() -> None:
    root = find_project_root()
    os.chdir(root)

    cfg = read_yaml_config(root)
    rel = parse_rel_from_01_db(root)
    species_all = get_species_list(cfg)
    only = get_only_species_from_env()
    if only:
        species = [only]
        if only not in species_all:
            raise ValueError(f"SPECIES_ID={only} not in config species_prefixes")
    else:
        species = species_all

    for sid in species:
        run_one(root, cfg, rel, sid)

if __name__ == "__main__":
    main()

