#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_matrix.py（极简）
- 读取拼接后的对齐（FASTA 或 .gz）
- 输出：species  total_sites  nongap_sites  occupancy_percent
用法：
  python scripts/report_matrix.py results/concat/supermatrix.aa.fas > results/reports/matrix.tsv
"""

import sys
import gzip

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode, encoding=None if "b" in mode else "utf-8", errors="ignore")
    return open(path, mode=mode, encoding=None if "b" in mode else "utf-8", errors="ignore")

def read_fasta(path):
    recs, hdr, seq = [], None, []
    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if hdr is not None:
            recs.append((hdr, "".join(seq)))
    return recs

def main():
    if len(sys.argv) != 2:
        print("Usage: report_matrix.py <supermatrix.fas>", file=sys.stderr)
        sys.exit(2)
    path = sys.argv[1]
    recs = read_fasta(path)
    print("species\ttotal_sites\tnongap_sites\toccupancy_percent")
    if not recs:
        sys.exit(0)
    total_len = len(recs[0][1])
    for hdr, seq in recs:
        sp = hdr.split()[0].split("|")[0]
        nongap = sum(1 for c in seq if c != "-" and c != "?")
        occ = (100.0 * nongap / total_len) if total_len > 0 else 0.0
        print(f"{sp}\t{total_len}\t{nongap}\t{occ:.2f}")

if __name__ == "__main__":
    main()
