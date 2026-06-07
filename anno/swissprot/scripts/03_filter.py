#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import csv
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import sys
from pathlib import Path

# 让“项目根目录”进入 sys.path，保证可以 import scripts._common
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


from scripts._common import (
    find_project_root, read_yaml_config, parse_rel_from_01_db,
    get_species_list, get_only_species_from_env, ensure_dir,
    clamp01, eprint
)

HEADER_OUT = ["qseqid","sseqid","evalue","bitscore","pident","length","qlen","slen","qcov","scov","stitle"]

def parse_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        # sometimes empty or weird; treat as NaN -> fail filter
        return float("nan")

def best_key(evalue: float, bitscore: float, qcov: float, pident: float) -> Tuple[float, float, float, float]:
    # smaller is better
    return (evalue, -bitscore, -qcov, -pident)

def run_one(root: Path, cfg: Dict[str, Any], rel: str, species_id: str) -> None:
    fcfg = cfg.get("filter", {}) if isinstance(cfg.get("filter", {}), dict) else {}
    EV_MAX = float(fcfg.get("EV_MAX", 1e-5))
    MIN_BITSCORE = float(fcfg.get("MIN_BITSCORE", 50))
    MIN_ALN_LEN = float(fcfg.get("MIN_ALN_LEN", 50))
    MIN_QCOV = float(fcfg.get("MIN_QCOV", 0.50))
    MIN_SCOV = float(fcfg.get("MIN_SCOV", 0.00))

    indir = root / "results" / "02_align" / rel / species_id
    in_hits = indir / "hits.diamond.tsv"
    if not in_hits.is_file() or in_hits.stat().st_size == 0:
        raise FileNotFoundError(f"Input hits missing/empty: {in_hits}")

    outdir = root / "results" / "03_filter" / rel / species_id
    ensure_dir(outdir)
    ensure_dir(root / "logs")

    out_filtered = outdir / "hits.filtered.tsv"
    out_besthit = outdir / "hits.besthit.tsv"
    out_params = outdir / "params.tsv"
    out_summary = outdir / "summary.tsv"

    print(f"[03_filter] species={species_id} REL={rel}", flush=True)
    print(f"[03_filter] IN={in_hits}", flush=True)
    print(f"[03_filter] OUT_FILTERED={out_filtered}", flush=True)
    print(f"[03_filter] OUT_BESTHIT={out_besthit}", flush=True)

    params_txt = (
        f"REL\t{rel}\n"
        f"SPECIES_ID\t{species_id}\n"
        f"IN\t{in_hits.as_posix()}\n"
        f"EV_MAX\t{EV_MAX}\n"
        f"MIN_BITSCORE\t{MIN_BITSCORE}\n"
        f"MIN_ALN_LEN\t{MIN_ALN_LEN}\n"
        f"MIN_QCOV\t{MIN_QCOV}\n"
        f"MIN_SCOV\t{MIN_SCOV}\n"
        f"QCOV_DEF\tqspan/qlen with clamp[0,1]\n"
        f"SCOV_DEF\tsspan/slen with clamp[0,1]\n"
    )
    out_params.write_text(params_txt, encoding="utf-8")

    n_total = 0
    n_pass = 0

    best: Dict[str, Tuple[Tuple[float,float,float,float], list]] = {}

    with in_hits.open("r", encoding="utf-8", errors="replace") as fin, \
         out_filtered.open("w", encoding="utf-8", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER_OUT)

        for row in reader:
            if not row:
                continue
            n_total += 1

            # Expected fixed columns:
            # 0 qseqid
            # 1 sseqid
            # 2 pident
            # 3 length
            # 4 mismatch
            # 5 gapopen
            # 6 qstart
            # 7 qend
            # 8 sstart
            # 9 send
            # 10 evalue
            # 11 bitscore
            # 12 qlen
            # 13 slen
            # 14 stitle (may spill to more columns, join them)
            if len(row) < 14:
                continue

            qseqid = row[0]
            sseqid = row[1]
            pident = parse_float(row[2])
            alen = parse_float(row[3])
            qstart = parse_float(row[6])
            qend = parse_float(row[7])
            sstart = parse_float(row[8])
            send = parse_float(row[9])
            evalue = parse_float(row[10])
            bitscore = parse_float(row[11])
            qlen = parse_float(row[12])
            slen = parse_float(row[13])
            stitle = row[14] if len(row) >= 15 else ""
            if len(row) > 15:
                stitle = " ".join([stitle] + row[15:])

            if qlen <= 0 or slen <= 0:
                continue

            # Coverage regularization: span-based, clamped
            qspan = abs(int(qend) - int(qstart)) + 1 if qstart > 0 and qend > 0 else 0
            sspan = abs(int(send) - int(sstart)) + 1 if sstart > 0 and send > 0 else 0
            qcov = clamp01((qspan / qlen) if qlen else 0.0)
            scov = clamp01((sspan / slen) if slen else 0.0)

            if not (evalue <= EV_MAX and bitscore >= MIN_BITSCORE and alen >= MIN_ALN_LEN and qcov >= MIN_QCOV and scov >= MIN_SCOV):
                continue

            n_pass += 1
            out_row = [qseqid, sseqid, f"{evalue:g}", f"{bitscore:g}", f"{pident:g}", f"{alen:g}", f"{int(qlen)}", f"{int(slen)}", f"{qcov:.6f}", f"{scov:.6f}", stitle]
            writer.writerow(out_row)

            k = best_key(evalue, bitscore, qcov, pident)
            cur = best.get(qseqid)
            if cur is None or k < cur[0]:
                best[qseqid] = (k, out_row)

    # Write besthit
    with out_besthit.open("w", encoding="utf-8", newline="") as fbh:
        writer = csv.writer(fbh, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER_OUT)
        for qseqid in sorted(best.keys()):
            writer.writerow(best[qseqid][1])

    n_best = len(best)

    # Query count from align
    n_query_file = root / "results" / "02_align" / rel / species_id / "n_query_seqs.txt"
    n_query = 0
    if n_query_file.is_file():
        try:
            n_query = int(n_query_file.read_text(encoding="utf-8").strip())
        except Exception:
            n_query = 0

    hit_rate = (n_best / n_query) if n_query > 0 else 0.0
    out_summary.write_text(
        f"REL\t{rel}\n"
        f"SPECIES_ID\t{species_id}\n"
        f"N_QUERY\t{n_query}\n"
        f"N_FILTERED_HITS\t{n_pass}\n"
        f"N_BESTHIT\t{n_best}\n"
        f"HIT_RATE\t{hit_rate:.6f}\n",
        encoding="utf-8"
    )

    if not out_besthit.is_file() or out_besthit.stat().st_size == 0:
        raise RuntimeError(f"besthit missing/empty: {out_besthit}")

    print(f"[03_filter] DONE species={species_id} filtered_hits={n_pass} besthit={n_best} hit_rate={hit_rate:.6f}", flush=True)

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

