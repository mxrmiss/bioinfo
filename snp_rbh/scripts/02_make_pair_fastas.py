#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_make_pair_fastas.py
从 RBH pairs 表中抽取两物种蛋白+CDS，生成每对的小 fasta 文件。
- 输入 pairs 表：results/step2_rbh/rbh_pairs.tsv（第一行表头）
- 输出：
  results/step3_align_prot/pairs_faa/*.faa
  results/step4_codon/pairs_fna/*.fna
- 进度：每 100 对打印一次
"""

# ----------------------------- 用户自定义区（皇上可改） -----------------------------
PAIRS_TSV = "results/step2_rbh/rbh_pairs.tsv"

REF_PROTEOME = "data/proteomes/Sinonovacula_rivularis.faa"
REF_CDS      = "data/cds/Sinonovacula_rivularis.fna"

OTH_PROTEOME = "data/proteomes/Sinonovacula_constricta.faa"
OTH_CDS      = "data/cds/Sinonovacula_constricta.fna"

OUT_PROT_DIR = "results/step3_align_prot/pairs_faa"
OUT_CDS_DIR  = "results/step4_codon/pairs_fna"

PROGRESS_EVERY = 100
# ------------------------------------------------------------------------------

import os, sys

def project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def infer_code_from_path(path: str) -> str:
    bn = os.path.basename(path)
    for suf in [".gz", ".gff3", ".gff", ".faa", ".fna", ".fa", ".fasta"]:
        if bn.endswith(suf):
            bn = bn[: -len(suf)]
    parts = bn.split("_")
    if len(parts) < 2:
        raise ValueError(f"Cannot infer code from filename: {os.path.basename(path)} (need Genus_species.*)")
    genus, species = parts[0], parts[1]
    return (genus[0] + species[:2]).lower()

def read_pairs(pairs_tsv: str):
    pairs = []
    with open(pairs_tsv, "r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            a, b = line.split("\t")[:2]
            pairs.append((a, b))
    return pairs

def read_fasta_subset(path: str, wanted: set):
    seqs = {}
    cur_id = None
    cur_keep = False
    buf = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None and cur_keep:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]
                cur_keep = (cur_id in wanted)
                buf = []
            else:
                if cur_keep:
                    buf.append(line.strip())
        if cur_id is not None and cur_keep:
            seqs[cur_id] = "".join(buf)
    return seqs

def write_two_fasta(outpath: str, id1: str, seq1: str, id2: str, seq2: str):
    with open(outpath, "w", encoding="utf-8") as w:
        w.write(f">{id1}\n{seq1}\n>{id2}\n{seq2}\n")

def main():
    os.chdir(project_root())

    os.makedirs(OUT_PROT_DIR, exist_ok=True)
    os.makedirs(OUT_CDS_DIR, exist_ok=True)

    ref_code = infer_code_from_path(REF_PROTEOME)
    oth_code = infer_code_from_path(OTH_PROTEOME)

    pairs = read_pairs(PAIRS_TSV)
    n_pairs = len(pairs)
    sys.stderr.write(f"[02_make_pair_fastas] pairs={n_pairs} ref={ref_code} other={oth_code}\n")

    ref_ids = set(a for a, _ in pairs)
    oth_ids = set(b for _, b in pairs)

    sys.stderr.write("[02_make_pair_fastas] loading FASTA subsets...\n")
    ref_faa = read_fasta_subset(REF_PROTEOME, ref_ids)
    oth_faa = read_fasta_subset(OTH_PROTEOME, oth_ids)
    ref_fna = read_fasta_subset(REF_CDS, ref_ids)
    oth_fna = read_fasta_subset(OTH_CDS, oth_ids)

    n_ok = 0
    n_miss = 0
    for i, (rid, oid) in enumerate(pairs, start=1):
        ok = True
        if rid not in ref_faa or oid not in oth_faa:
            ok = False
        if rid not in ref_fna or oid not in oth_fna:
            ok = False

        if not ok:
            n_miss += 1
        else:
            tag = f"{rid}__{oid}"
            write_two_fasta(os.path.join(OUT_PROT_DIR, f"{tag}.faa"), rid, ref_faa[rid], oid, oth_faa[oid])
            write_two_fasta(os.path.join(OUT_CDS_DIR,  f"{tag}.fna"), rid, ref_fna[rid], oid, oth_fna[oid])
            n_ok += 1

        if (i % PROGRESS_EVERY) == 0 or i == n_pairs:
            sys.stderr.write(f"[02_make_pair_fastas] {i}/{n_pairs} done (ok={n_ok}, miss={n_miss})\n")

    sys.stderr.write(f"[02_make_pair_fastas] DONE ok={n_ok} miss={n_miss}\n")

if __name__ == "__main__":
    main()
