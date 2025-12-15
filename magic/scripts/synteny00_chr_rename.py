#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny00_chr_rename.py
Chr length stats + rename table (Chr01..)

Contract:
- Read raw_data/synteny_species_meta.tsv (tab-separated; use first 2 columns: species_id, group).
- Auto-detect input files under raw_data/genomes/ and raw_data/gff/ (supports .gz).
- Output:
    output/synteny_00_chr_rename/chr_lengths_<species>.tsv
    output/synteny_00_chr_rename/chr_rename_<species>.tsv
    output/synteny_00_chr_rename/step00.summary.tsv
- No sentinel files.
- No command-line args. Edit parameters in the "USER PARAMETERS" section.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable


# ============================
# USER PARAMETERS (edit here)
# ============================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_DATA_DIR = PROJECT_ROOT / "raw_data"
META_TSV = RAW_DATA_DIR / "synteny_species_meta.tsv"
GENOMES_DIR = RAW_DATA_DIR / "genomes"
GFF_DIR = RAW_DATA_DIR / "gff"

OUTPUT_DIR = PROJECT_ROOT / "output" / "synteny_00_chr_rename"
LOG_DIR = PROJECT_ROOT / "logs"
CLEAN_OUTPUT = True

CHR_LENGTH_THRESHOLD_BP = 10_000_000
PREFER_GENOME_FASTA = True


# ============================
# Utilities
# ============================

def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


class Logger:
    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

    def _write(self, level: str, msg: str) -> None:
        line = f"{now_ts()} [{level}] {msg}"
        print(line, flush=True)
        with self.log_file.open("a", encoding="utf-8") as fw:
            fw.write(line + "\n")

    def info(self, msg: str) -> None:
        self._write("INFO", msg)

    def warn(self, msg: str) -> None:
        self._write("WARN", msg)

    def error(self, msg: str) -> None:
        self._write("ERROR", msg)


def abort(logger: Logger, msg: str, code: int = 1) -> None:
    logger.error(msg)
    raise SystemExit(code)


def clean_dir(path: Path, clean: bool, logger: Logger) -> None:
    if clean and path.exists():
        logger.info(f"Clean output dir: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def detect_first_existing(base_dir: Path, candidates: List[str]) -> Optional[Path]:
    for name in candidates:
        p = base_dir / name
        if p.exists():
            return p
    return None


def read_meta_species(meta_tsv: Path) -> Tuple[List[str], Dict[str, str]]:
    """
    Read meta: must have header starting with 'species_id' and 'group' (tab-separated).
    Extra columns (if any) are ignored; only first two are used.
    """
    species: List[str] = []
    group: Dict[str, str] = {}

    with open_text_auto(meta_tsv) as fr:
        header_line = fr.readline().rstrip("\n")
        if not header_line:
            raise ValueError("meta header is empty")
        header = header_line.split("\t")
        if len(header) < 2 or header[0] != "species_id" or header[1] != "group":
            raise ValueError(f"meta header must start with: species_id<TAB>group, got: {header}")

        for ln, line in enumerate(fr, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"meta line {ln} must have at least 2 columns: {line}")
            sid = parts[0].strip()
            grp = parts[1].strip()
            if not sid:
                raise ValueError(f"meta line {ln} species_id empty")
            species.append(sid)
            group[sid] = grp

    if len(species) < 2:
        raise ValueError("meta must contain at least 2 species")
    return species, group


def safe_int(x: str, default: int = 0) -> int:
    try:
        return int(x)
    except Exception:
        return default


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fw:
        w = csv.writer(fw, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def fasta_lengths(fa_path: Path) -> Dict[str, int]:
    lens: Dict[str, int] = {}
    cur_id: Optional[str] = None
    cur_len = 0

    with open_text_auto(fa_path) as fr:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    lens[cur_id] = cur_len
                cur_id = line[1:].split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id is not None:
            lens[cur_id] = cur_len

    return lens


def gff_gene_counts_and_maxend(gff_path: Path) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Return:
      n_genes per seqid (prefer gene; if no gene feature then use mRNA/transcript)
      max_end per seqid (as a rough length proxy)
    """
    n_gene: Dict[str, int] = {}
    n_mrna: Dict[str, int] = {}
    max_end: Dict[str, int] = {}
    has_gene = False

    with open_text_auto(gff_path) as fr:
        for line in fr:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid = parts[0]
            ftype = parts[2]
            end = parts[4]

            e = safe_int(end, 0)
            if e > max_end.get(seqid, 0):
                max_end[seqid] = e

            ft = ftype.lower()
            if ft == "gene":
                has_gene = True
                n_gene[seqid] = n_gene.get(seqid, 0) + 1
            elif ft in ("mrna", "transcript"):
                n_mrna[seqid] = n_mrna.get(seqid, 0) + 1

    return (n_gene if has_gene else n_mrna), max_end


def main() -> None:
    logger = Logger(LOG_DIR / "synteny00_chr_rename.log")
    t0 = time.time()

    logger.info("========== synteny00 — Chr rename table ==========")
    logger.info(f"PROJECT_ROOT={PROJECT_ROOT}")
    logger.info(f"META_TSV={META_TSV}")
    logger.info(f"OUTPUT_DIR={OUTPUT_DIR}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"CHR_LENGTH_THRESHOLD_BP={CHR_LENGTH_THRESHOLD_BP}")
    logger.info(f"PREFER_GENOME_FASTA={PREFER_GENOME_FASTA}")

    if not META_TSV.exists():
        abort(logger, f"Missing meta: {META_TSV}")

    try:
        species, _group = read_meta_species(META_TSV)
    except Exception as e:
        abort(logger, f"Meta parse error: {e}")

    logger.info(f"META species count={len(species)}")

    clean_dir(OUTPUT_DIR, CLEAN_OUTPUT, logger)

    summary_rows: List[List[str]] = []

    for sid in species:
        gff = detect_first_existing(GFF_DIR, [
            f"{sid}.gff3", f"{sid}.gff", f"{sid}.gff3.gz", f"{sid}.gff.gz"
        ])
        genome = detect_first_existing(GENOMES_DIR, [
            f"{sid}.fa", f"{sid}.fasta", f"{sid}.fna",
            f"{sid}.fa.gz", f"{sid}.fasta.gz", f"{sid}.fna.gz",
        ])

        if gff is None and genome is None:
            abort(logger, f"[{sid}] Missing both GFF and genome FASTA under raw_data/")

        n_genes: Dict[str, int] = {}
        gff_len_proxy: Dict[str, int] = {}
        if gff is not None:
            n_genes, gff_len_proxy = gff_gene_counts_and_maxend(gff)

        lens: Dict[str, int] = {}

        if PREFER_GENOME_FASTA and genome is not None:
            lens = fasta_lengths(genome)
            logger.info(f"[{sid}] length source=genome: {genome.name} (seqids={len(lens)})")
        else:
            if genome is not None and not PREFER_GENOME_FASTA:
                lens = fasta_lengths(genome)
                logger.info(f"[{sid}] length source=genome (forced): {genome.name} (seqids={len(lens)})")
            else:
                lens = {k: int(v) for k, v in gff_len_proxy.items()}
                logger.info(f"[{sid}] length source=gff (max_end proxy): {gff.name if gff else 'NA'} (seqids={len(lens)})")

        if not lens:
            abort(logger, f"[{sid}] No seqid length detected (empty lengths).")

        chr_lengths_path = OUTPUT_DIR / f"chr_lengths_{sid}.tsv"
        chr_len_rows: List[List[str]] = []
        for seqid_raw, L in sorted(lens.items(), key=lambda x: (-x[1], x[0])):
            chr_len_rows.append([seqid_raw, str(L), str(n_genes.get(seqid_raw, 0))])
        write_tsv(chr_lengths_path, ["seqid_raw", "length_bp", "n_genes"], chr_len_rows)

        kept = [(k, v) for k, v in lens.items() if v >= CHR_LENGTH_THRESHOLD_BP]
        kept_sorted = sorted(kept, key=lambda x: (-x[1], x[0]))

        if len(kept_sorted) == 0:
            abort(logger, f"[{sid}] n_kept_chromosomes==0 (threshold={CHR_LENGTH_THRESHOLD_BP}).")

        chr_rename_path = OUTPUT_DIR / f"chr_rename_{sid}.tsv"
        rename_rows: List[List[str]] = []
        for i, (seqid_raw, L) in enumerate(kept_sorted, start=1):
            chr_name = f"Chr{i:02d}" if i < 100 else f"Chr{i:03d}"
            rename_rows.append([sid, seqid_raw, chr_name, str(i), str(L), "yes"])

        kept_set = {k for k, _ in kept_sorted}
        for seqid_raw, L in sorted(lens.items(), key=lambda x: (-x[1], x[0])):
            if seqid_raw in kept_set:
                continue
            rename_rows.append([sid, seqid_raw, "NA", "0", str(L), "no"])

        write_tsv(
            chr_rename_path,
            ["species_id", "seqid_raw", "seqid_renamed", "rank", "length_bp", "is_chromosome"],
            rename_rows
        )

        total_kept = sum(L for _, L in kept_sorted)
        summary_rows.append([sid, str(len(lens)), str(len(kept_sorted)), str(CHR_LENGTH_THRESHOLD_BP), str(total_kept)])
        logger.info(
            f"[{sid}] n_raw_seqids={len(lens)} "
            f"n_kept_chromosomes={len(kept_sorted)} "
            f"total_kept_length_bp={total_kept}"
        )

    write_tsv(
        OUTPUT_DIR / "step00.summary.tsv",
        ["species_id", "n_raw_seqids", "n_kept_chromosomes", "min_chr_length_bp_threshold", "total_kept_length_bp"],
        summary_rows
    )

    logger.info(f"Done. runtime_sec={int(time.time() - t0)}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        logger = Logger(LOG_DIR / "synteny00_chr_rename.log")
        logger.error("Unhandled exception: " + repr(e))
        logger.error(traceback.format_exc())
        raise

