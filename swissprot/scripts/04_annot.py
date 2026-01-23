#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import re
import csv
from pathlib import Path
from typing import Dict, Any, Optional
import sys
from pathlib import Path

# 让“项目根目录”进入 sys.path，保证可以 import scripts._common
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


from scripts._common import (
    find_project_root, read_yaml_config, parse_rel_from_01_db,
    get_species_list, get_only_species_from_env, ensure_dir
)

OUT_HEADER = ["protein_id","swissprot_accession","entry_name","protein_name","organism","gene_name","evalue","bitscore","pident","qcov","scov"]

def parse_accession_entry(sid: str) -> (str, str):
    acc = sid
    entry = ""
    if "|" in sid:
        parts = sid.split("|")
        if len(parts) >= 2:
            acc = parts[1]
        if len(parts) >= 3:
            entry = parts[2]
    acc = re.sub(r"\.[0-9]+$", "", acc)
    return acc, entry

def parse_organism(title: str) -> str:
    # Prefer OS=... OX=...
    m = re.search(r"OS=(.+?)\s+OX=", title)
    if m:
        return m.group(1).strip()
    # Fallback: OS= until GN/PE/SV/OX/end
    m = re.search(r"OS=(.+)$", title)
    if not m:
        return ""
    tmp = m.group(1)
    tmp = re.sub(r"\s+GN=.*$", "", tmp)
    tmp = re.sub(r"\s+PE=.*$", "", tmp)
    tmp = re.sub(r"\s+SV=.*$", "", tmp)
    tmp = re.sub(r"\s+OX=.*$", "", tmp)
    return tmp.strip()

def parse_gene_name(title: str) -> str:
    m = re.search(r"GN=([^\s]+)", title)
    return m.group(1).strip() if m else ""

def parse_protein_name(title: str) -> str:
    # take before OS=
    pname = re.sub(r"\s+OS=.*$", "", title).strip()
    # remove leading "sp|ACC|ENTRY " or "tr|ACC|ENTRY "
    pname = re.sub(r"^[a-z]{2}\|[^|]+\|[^ ]+\s+", "", pname)
    return pname.strip()

def run_one(root: Path, cfg: Dict[str, Any], rel: str, species_id: str) -> None:
    besthit = root / "results" / "03_filter" / rel / species_id / "hits.besthit.tsv"
    if not besthit.is_file() or besthit.stat().st_size == 0:
        raise FileNotFoundError(f"besthit missing/empty: {besthit}")

    outdir = root / "results" / "04_annot" / rel / species_id
    ensure_dir(outdir)

    out_annot = outdir / "swissprot_annotation.tsv"
    out_qc = outdir / "qc_missing_fields.tsv"
    out_md5 = outdir / "md5_swissprot_annotation.txt"

    # inherit params.tsv
    params_src = root / "results" / "03_filter" / rel / species_id / "params.tsv"
    if params_src.is_file():
        (outdir / "params.tsv").write_text(params_src.read_text(encoding="utf-8", errors="replace"), encoding="utf-8")

    print(f"[04_annot] species={species_id} REL={rel}", flush=True)
    print(f"[04_annot] IN={besthit}", flush=True)
    print(f"[04_annot] OUT={out_annot}", flush=True)

    n = 0
    no_os = 0
    no_gn = 0

    with besthit.open("r", encoding="utf-8", errors="replace") as fin, \
         out_annot.open("w", encoding="utf-8", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        header = next(reader, None)
        writer.writerow(OUT_HEADER)

        for row in reader:
            if not row:
                continue
            # besthit columns (from 03_filter): HEADER_OUT
            # qseqid, sseqid, evalue, bitscore, pident, length, qlen, slen, qcov, scov, stitle
            if len(row) < 11:
                continue
            q = row[0]
            sid = row[1]
            evalue = row[2]
            bitscore = row[3]
            pident = row[4]
            qcov = row[8]
            scov = row[9]
            title = row[10] if len(row) == 11 else " ".join(row[10:])

            acc, entry = parse_accession_entry(sid)
            org = parse_organism(title)
            gn = parse_gene_name(title)
            pname = parse_protein_name(title)

            n += 1
            if not org:
                no_os += 1
            if not gn:
                no_gn += 1

            writer.writerow([q, acc, entry, pname, org, gn, evalue, bitscore, pident, qcov, scov])

    if n == 0:
        raise RuntimeError(f"No annotations written: {out_annot}")

    no_os_rate = (no_os / n) if n else 0.0
    no_gn_rate = (no_gn / n) if n else 0.0
    out_qc.write_text(
        f"N_annot\t{n}\n"
        f"NO_OS\t{no_os}\n"
        f"NO_GN\t{no_gn}\n"
        f"NO_OS_rate\t{no_os_rate:.6f}\n"
        f"NO_GN_rate\t{no_gn_rate:.6f}\n",
        encoding="utf-8"
    )

    # md5
    import hashlib
    h = hashlib.md5()
    with out_annot.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    out_md5.write_text(f"{h.hexdigest()}  {out_annot.as_posix()}\n", encoding="utf-8")

    print(f"[04_annot] DONE species={species_id} N_annot={n} NO_GN_rate={no_gn_rate:.6f}", flush=True)

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

