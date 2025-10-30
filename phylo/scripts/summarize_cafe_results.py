#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
将 CAFÉ 结果标准化为一张 summary.tsv：
- 按 p_cutoff（默认与 config.cafe.p_cutoff 一致）标注 significant
- 合并 expanded/contracted（若存在）以补全报告
- 连接 family_id.map.tsv，保留 Orthogroup 原 ID 便于追溯
"""
import argparse, csv
from pathlib import Path

def read_tsv(p: Path):
    rows = []
    with p.open("r", newline='') as f:
        r = csv.reader(f, delimiter='\t')
        rows = [row for row in r]
    if not rows:
        return [], []
    return rows[0], rows[1:]

def write_tsv(p: Path, header, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(header)
        for row in rows:
            w.writerow(row)

def to_index_map(header):
    return {h: i for i, h in enumerate(header)}

def load_map(p: Path):
    if not p.exists(): return {}
    h, rows = read_tsv(p)
    idx = to_index_map(h)
    x = {}
    for r in rows:
        fid = r[idx.get("FamilyID", 0)]
        og  = r[idx.get("Orthogroup", 1)]
        x[fid] = og
    return x

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True)
    ap.add_argument("--expanded", required=False)
    ap.add_argument("--contracted", required=False)
    ap.add_argument("--fammap", required=True)
    ap.add_argument("--p-cutoff", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rep_h, rep = read_tsv(Path(args.report))
    idx = to_index_map(rep_h)

    # 尝试定位常用列：FamilyID、p值、物种/分支信息
    fid_k = "FamilyID" if "FamilyID" in idx else rep_h[0]
    pkeys = [k for k in rep_h if k.lower() in ["p", "pvalue", "p-value", "p_val", "pval"]]
    pkey  = pkeys[0] if pkeys else None

    # 附加 expanded / contracted 的并集信息（若存在）
    flags = {}
    def mark_flag(tsv, label):
        if not tsv: return
        h, rows = read_tsv(Path(tsv))
        if not rows: return
        i0 = 0  # FamilyID 默认在首列
        for r in rows:
            fid = r[i0]
            d = flags.setdefault(fid, set())
            d.add(label)

    mark_flag(args.expanded, "expanded")
    mark_flag(args.contracted, "contracted")

    fammap = load_map(Path(args.fammap))

    out_h = ["FamilyID", "Orthogroup", "flag", "p", "significant"] + [k for k in rep_h if k not in ("FamilyID", pkey)]
    out_rows = []
    for r in rep:
        fid = r[idx.get(fid_k, 0)]
        pval = r[idx[pkey]] if pkey else ""
        flag = ",".join(sorted(flags.get(fid, []))) if fid in flags else ""
        sig  = ""
        try:
            if pval != "":
                sig = "yes" if float(pval) <= args.p_cutoff else "no"
        except:
            sig = ""
        og = fammap.get(fid, fid)
        # 拼一行
        rest = [r[i] for i in range(len(rep_h)) if rep_h[i] not in ("FamilyID", pkey)]
        out_rows.append([fid, og, flag, pval, sig] + rest)

    write_tsv(Path(args.out), out_h, out_rows)

if __name__ == "__main__":
    main()

