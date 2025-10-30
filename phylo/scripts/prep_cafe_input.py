#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
prep_cafe_input.py —— 从 OrthoFinder GeneCount 生成“对齐且非零”的 CAFE5 输入

做什么：
1) 从 RESULTS_PATH.txt 找到 Orthogroups/Orthogroups.GeneCount.tsv
2) 读取等时树 tip，标准化命名（. - 空格→_，去不可见字符，确保唯一）
3) 用“标准化 + 去常见后缀 + 前缀最短匹配”把 GeneCount 的列名与树 tip 一一对齐
4) 构建 CAFE5 输入表：
   - 表头首列固定为 'family_id'（小写）
   - family_id 用 0,1,2,... 重新编号
   - 缺失补 0，非数转 0，负数截 0
5) 另写出：
   - 物种重命名表：species_rename.tsv（original→pretty）
   - family_id.map.tsv：FamilyID ↔ Orthogroup
"""

import argparse, csv, re, sys
from pathlib import Path
from collections import Counter

# ---------- 小工具 ----------

def norm_name(s: str) -> str:
    s = s.strip().strip("'\"")
    s = re.sub(r"[^\x20-\x7E]", "", s)
    s = re.sub(r"\|.*$", "", s)               # Species|SeqID -> Species
    s = s.replace(".", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9_]", "", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

_SUFFIX_PAT = re.compile(
    r"(?:_?proteome.*|_?protein.*|_?pep.*|_?aa.*|_?cds.*|_?faa.*|_?fa.*|_?seq.*|_?v\d+|_?ver\d+)$",
    re.IGNORECASE
)
def base_name(s: str) -> str:
    return _SUFFIX_PAT.sub("", s)

def ensure_unique(names):
    seen = Counter(); out = []
    for n in names:
        seen[n] += 1
        if seen[n] == 1:
            out.append(n)
        else:
            k = seen[n]; new = f"{n}_{k}"
            while new in seen: k += 1; new = f"{n}_{k}"
            seen[new] += 1; out.append(new)
    return out

def read_text_clean(p: Path) -> str:
    raw = p.read_bytes()
    if raw[:3] == b"\xef\xbb\xbf": raw = raw[3:]
    return raw.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-path", required=True, help="results/orthogroups/RESULTS_PATH.txt")
    ap.add_argument("--tree", required=True, help="等时树（已 pretty 过）")
    ap.add_argument("--out-tsv", required=True, help="输出：results/cafe/cafe_input.tsv")
    ap.add_argument("--rename-out", required=True, help="输出：results/cafe/species_rename.tsv")
    return ap.parse_args()

# ---------- 读取 OrthoFinder / 树 ----------

def gc_table_from_results(results_path: Path) -> Path:
    root = Path(results_path.read_text().strip())
    gc = root / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    if not gc.exists():
        sys.exit(f"[ERR] 缺少 GeneCount 表：{gc}")
    return gc

def read_tree_tips(tree_path: Path):
    txt = read_text_clean(tree_path)
    tips_raw = re.findall(r'[,(\s]([A-Za-z0-9._\-\|]+):\s*[0-9eE\+\-\.]', txt)
    if not tips_raw:
        sys.exit("[ERR] 无法从树解析到 tip 标签。")
    tips_norm = ensure_unique([norm_name(x) for x in tips_raw])
    return tips_norm

# ---------- GeneCount 列名 ↔ 树物种 对齐 ----------

def build_header_index(header_cells):
    exact, base, all_keys = {}, {}, set()
    for j, h in enumerate(header_cells):
        if j == 0:  # 跳过 'Orthogroup'
            continue
        nh = norm_name(h)
        nb = base_name(nh)
        if nh not in exact: exact[nh] = j
        if nb and nb not in base: base[nb] = j
        all_keys.add(nh); all_keys.add(nb)
    return exact, base, all_keys

def pick_column_for_species(sp, exact, base, all_keys):
    spn = norm_name(sp)
    spb = base_name(spn)
    if spn in exact: return exact[spn], spn
    if spn in base:  return base[spn],  spn
    if spb in exact: return exact[spb], spb
    if spb in base:  return base[spb],  spb
    cands = [k for k in all_keys if k and (k == spn or k.startswith(spn+"_") or k == spb or (spb and k.startswith(spb+"_")))]
    if cands:
        best = sorted(cands, key=len)[0]
        if best in exact: return exact[best], best
        if best in base:  return base[best],  best
    return None, None

def sanitize_int(x: str) -> int:
    s = re.sub(r"[^\x20-\x7E]", "", (x or "").strip())
    if not s or s.upper() == "NA": return 0
    try:
        v = int(float(s))
    except Exception:
        v = 0
    return 0 if v < 0 else v

# ---------- 主流程 ----------

def main():
    args = parse_args()
    results_path = Path(args.results_path)
    gc_path = gc_table_from_results(results_path)
    tips = read_tree_tips(Path(args.tree))

    # 读取 GeneCount
    with gc_path.open() as f:
        r = csv.reader(f, delimiter="\t")
        hdr = next(r)
        rows = list(r)

    exact, base, all_keys = build_header_index(hdr)

    # 为每个 tip 选择列号
    col_index = {}
    mapping_dbg = []
    for sp in tips:
        j, key = pick_column_for_species(sp, exact, base, all_keys)
        col_index[sp] = j
        mapping_dbg.append((sp, j, key))

    unmapped = [sp for sp, j, _ in mapping_dbg if j is None]
    if unmapped:
        print("[WARN] 以下物种在 GeneCount 表头无法匹配到列（将按 0 处理）：", ", ".join(unmapped))
    else:
        print("[INFO] GeneCount 与树的物种列已全部匹配。")
    print("[DEBUG] 物种→列 映射（前 10 个）：")
    for sp, j, key in mapping_dbg[:10]:
        print(f"    {sp:20s} <- col {j if j is not None else 'None'} ({key})")

    # family_id.map.tsv
    famap = Path(args.out_tsv).parent / "family_id.map.tsv"
    famap.parent.mkdir(parents=True, exist_ok=True)
    with famap.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["FamilyID", "Orthogroup"])
        for i, row in enumerate(rows):
            og = (row[0] if row else "").strip()
            w.writerow([str(i), og])

    # species_rename.tsv（记录“匹配到的原始列名”→“pretty 名”）
    rename_out = Path(args.rename_out)
    rename_out.parent.mkdir(parents=True, exist_ok=True)
    with rename_out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["original", "pretty"])
        for sp in tips:
            j, key = col_index[sp], None
            if j is not None and j < len(hdr): key = hdr[j]
            w.writerow([key if key is not None else "", sp])

    # 生成 CAFE 输入（首列 = family_id；0-based）
    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["family_id"] + tips)
        for i, row in enumerate(rows):  # 0-based
            vals = []
            for sp in tips:
                j = col_index[sp]; v = 0
                if j is not None and j < len(row):
                    v = sanitize_int(row[j])
                vals.append(str(v))
            w.writerow([str(i)] + vals)

    # 简要统计
    import itertools
    with out_tsv.open() as f:
        rr = csv.reader(f, delimiter="\t")
        _hdr = next(rr)
        n = z = 0
        for r in rr:
            n += 1
            if sum(int(x) for x in r[1:]) == 0: z += 1
    print(f"[OK] 写出 {out_tsv}；行数={n}；其中全 0 行={z} ({z/n:.2%})；非零行={n-z}")

if __name__ == "__main__":
    main()

