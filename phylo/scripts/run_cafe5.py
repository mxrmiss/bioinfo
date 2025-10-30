#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_cafe5.py —— 规范化树 + 强力对齐 GeneCount 列名与树 tip + 过滤全0家族 + 智能收集

关键增强：
- GeneCount 表头与树 tip 的对齐采用“标准化 + 剥常见后缀 + 前缀最短匹配”策略，
  解决 A_californica ↔ A_californica_protein 一类不一致导致的“全零”问题。
- 首列固定 family_id（0,1,2,...），剔除全零家族。
- 智能收集不同 CAFE5 变体的输出，统一到 report.tsv / expanded.tsv / contracted.tsv / lambda.txt
- 调试信息：输出哪些家族在某些物种中缺失，方便进一步定位。
"""

import argparse, csv, os, re, shlex, subprocess, sys, textwrap
from pathlib import Path
from collections import Counter

# ---------------- 工具函数 ------------------

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_tsv", required=True)
    ap.add_argument("--tree", dest="tree", required=True)
    ap.add_argument("--out", dest="outdir", required=True)
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--cafe5", default="cafe5")
    ap.add_argument("--p-cutoff", dest="pcut", type=float, default=0.05)
    return ap.parse_args()

def read_text_clean(p: Path) -> str:
    raw = p.read_bytes()
    if raw[:3] == b"\xef\xbb\xbf":
        raw = raw[3:]
    return raw.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")

def write_text(p: Path, txt: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(txt, encoding="utf-8", newline="\n")

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
    seen = Counter(); uniq = []
    for n in names:
        seen[n] += 1
        if seen[n] == 1:
            uniq.append(n)
        else:
            k = seen[n]
            new = f"{n}_{k}"
            while new in seen: k += 1; new = f"{n}_{k}"
            seen[new] += 1
            uniq.append(new)
    return uniq

TIP_PAT = re.compile(r"([\(,])\s*([A-Za-z0-9._\-\|]+)\s*(?=:\s*[0-9eE\+\-\.])")
def sanitize_tree_labels(text: str):
    tips_raw = []
    def _collect(m):
        tips_raw.append(m.group(2)); return m.group(0)
    TIP_PAT.sub(_collect, text)  # 先收集
    tips_norm = ensure_unique([norm_name(x) for x in tips_raw])
    mapping = {r: n for r, n in zip(tips_raw, tips_norm)}
    def _repl(m):
        return f"{m.group(1)}{mapping[m.group(2)]}"
    new_text = TIP_PAT.sub(_repl, text)
    if re.search(r":\s*[^0-9eE\+\.-]", new_text):
        raise SystemExit("ERROR: 树分支长度格式异常（冒号后不是数字）。")
    return new_text, tips_norm

# ---------- GeneCount 列名对齐 ----------
def build_header_index(header_cells):
    exact, base, all_keys = {}, {}, set()
    for j, h in enumerate(header_cells):
        if j == 0:  # 跳过 Orthogroup 列
            continue
        nh = norm_name(h); nb = base_name(nh)
        if nh not in exact: exact[nh] = j
        if nb and nb not in base: base[nb] = j
        all_keys.add(nh); all_keys.add(nb)
    return exact, base, all_keys

def pick_column_for_species(sp, exact, base, all_keys):
    spn = norm_name(sp); spb = base_name(spn)
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

def sanitize_int(cell: str) -> int:
    s = re.sub(r"[^\x20-\x7E]", "", (cell or "").strip())
    if not s or s.upper() == "NA": return 0
    try: v = int(float(s))
    except Exception: v = 0
    return 0 if v < 0 else v

# ---------- TSV 规范化 + 过滤全0家族，产出双格式 ----------
def prepare_two_inputs(tsv_text: str, tips_norm: list):
    lines = [ln for ln in tsv_text.split("\n") if ln != ""]
    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    if not rows:
        raise SystemExit("ERROR: cafe_input.tsv 为空。")
    hdr = [re.sub(r"[^\x20-\x7E]", "", c.strip()) for c in rows[0]]
    data = rows[1:]

    exact, base, all_keys = build_header_index(hdr)
    col_index = {}
    dbg = []
    for sp in tips_norm:
        j, key = pick_column_for_species(sp, exact, base, all_keys)
        col_index[sp] = j; dbg.append((sp, j, key))
    unmapped = [sp for sp,j,_ in dbg if j is None]
    if unmapped:
        print("[WARN] 以下物种在 GeneCount 中未匹配到列（按 0 处理）：", ", ".join(unmapped))
    else:
        print("[INFO] 所有物种都已成功映射到 GeneCount 列。")
    print("[DEBUG] 物种→列 映射示例（前 10 个）：")
    for sp,j,key in dbg[:10]:
        print(f"    {sp:20s} <- col {j if j is not None else 'None'} ({key})")

    # 聚合行：过滤全0家族
    kept = []
    drop0 = 0
    for row in data:
        vals = []
        for sp in tips_norm:
            j = col_index[sp]; v = 0
            if j is not None and j < len(row):
                v = sanitize_int(row[j])
            vals.append(v)
        if sum(vals) == 0:
            drop0 += 1
            continue
        kept.append(vals)
    if not kept:
        raise SystemExit(
            "ERROR: 规范化后所有家族对当前物种集合都是全0。\n"
            "请检查树 tip 与 GeneCount 列是否真的同一批物种（名字虽对，可能列全为0）。"
        )

    # A 格式（family_id + 0-based）
    famA = ["\t".join(["family_id"] + tips_norm)]
    for i, vals in enumerate(kept):
        famA.append("\t".join([str(i)] + [str(v) for v in vals]))
    tsv_A = "\n".join(famA) + "\n"

    # B 格式（desc + 1-based）
    famB = ["\t".join(["desc"] + tips_norm)]
    for i, vals in enumerate(kept, start=1):
        famB.append("\t".join([str(i)] + [str(v) for v in vals]))
    tsv_B = "\n".join(famB) + "\n"

    return tsv_A, tsv_B, drop0

# ---------- 文件收集与调试 ----------
def ls_tree(path: Path, max_lines=200):
    out = []
    for root, _, files in os.walk(path):
        r = Path(root)
        for f in sorted(files):
            p = (r / f)
            out.append(str(p.relative_to(path)))
            if len(out) >= max_lines:
                return "...\n" + "\n".join(out[-max_lines:])
    return "\n".join(out)

def _first(globs, base: Path):
    for pat in globs:
        for p in base.rglob(pat):
            if p.is_file() and p.stat().st_size > 0:
                return p
    return None

def _find_report_paths(outdir: Path):
    report = _first([ "*[Rr]eport*.tsv", "*summary*.tsv", "*results*.tsv",
                      "*[Rr]eport*.txt", "*[Rr]eport", "*[Rr]esults.txt" ], outdir)
    expanded = _first([ "*expand*.tsv", "*gain*.tsv", "*increase*.tsv" ], outdir)
    contracted = _first([ "*contract*.tsv", "*loss*.tsv", "*decrease*.tsv" ], outdir)
    lamb = _first([ "*lambda*.txt", "*lambda*.tsv", "*lambda" ], outdir)
    return report, expanded, contracted, lamb

def _copy_to(dst: Path, src: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(src.read_bytes())

def _derive_exp_con_from_report(rep: Path, out_exp: Path, out_con: Path):
    try:
        with rep.open() as f:
            hdr = f.readline().rstrip("\n").split("\t")
            rows = [ln.rstrip("\n").split("\t") for ln in f]
        hlow = [h.lower() for h in hdr]
        def col(*keys):
            for k in keys:
                if k in hlow: return hlow.index(k)
            return -1
        ci_type = col("type","status","event")
        ci_exp  = col("expanded","expansion","gain","increase")
        ci_con  = col("contracted","contraction","loss","decrease")
        if ci_type==-1 and ci_exp==-1 and ci_con==-1:
            return False
        exp_rows, con_rows = [], []
        for r in rows:
            flag_exp = False; flag_con = False
            if 0 <= ci_type < len(r):
                t = r[ci_type].lower()
                flag_exp |= ("expan" in t or "gain" in t or "incre" in t)
                flag_con |= ("contr" in t or "loss" in t or "decre" in t)
            if 0 <= ci_exp < len(r):
                flag_exp |= r[ci_exp].strip().lower() in ("1","true","yes","y")
            if 0 <= ci_con < len(r):
                flag_con |= r[ci_con].strip().lower() in ("1","true","yes","y")
            if flag_exp and not flag_con: exp_rows.append(r)
            if flag_con and not flag_exp: con_rows.append(r)
        if exp_rows:
            with out_exp.open("w", newline="") as f:
                f.write("\t".join(hdr) + "\n")
                for r in exp_rows: f.write("\t".join(r) + "\n")
        if con_rows:
            with out_con.open("w", newline="") as f:
                f.write("\t".join(hdr) + "\n")
                for r in con_rows: f.write("\t".join(r) + "\n")
        return bool(exp_rows or con_rows)
    except Exception:
        return False

def run_cafe(cmd, outdir: Path):
    env = os.environ.copy()
    for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        env[k] = env.get(k, str(cmd[-1])) if isinstance(cmd[-1], str) else "1"
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    (outdir / "cafe5.stdout.log").write_text(proc.stdout)
    return proc

def have_outputs(outdir: Path):
    std_rep = outdir / "report.tsv"
    std_exp = outdir / "expanded.tsv"
    std_con = outdir / "contracted.tsv"
    std_lmb = outdir / "lambda.txt"
    return all(p.exists() and p.stat().st_size > 0 for p in (std_rep, std_exp, std_con, std_lmb))

def collect_outputs(outdir: Path):
    std_rep = outdir / "report.tsv"
    std_exp = outdir / "expanded.tsv"
    std_con = outdir / "contracted.tsv"
    std_lmb = outdir / "lambda.txt"
    rep, exp, con, lamb = _find_report_paths(outdir)
    if rep and not std_rep.exists(): _copy_to(std_rep, rep)
    if lamb and not std_lmb.exists(): _copy_to(std_lmb, lamb)
    if exp and not std_exp.exists():  _copy_to(std_exp, exp)
    if con and not std_con.exists():  _copy_to(std_con, con)
    if std_rep.exists() and (not std_exp.exists() or not std_con.exists()):
        _derive_exp_con_from_report(std_rep, std_exp, std_con)
    return have_outputs(outdir)

# ---------------- 主流程 ----------------

def main():
    args = parse_args()
    in_tsv  = Path(args.in_tsv)
    in_tree = Path(args.tree)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sani = outdir / "_sanitized"; sani.mkdir(parents=True, exist_ok=True)

    # 1) 树规范化
    tree_txt = read_text_clean(in_tree)
    tree_sani_txt, tips = sanitize_tree_labels(tree_txt)
    write_text(sani / "tree.normalized.nwk", tree_sani_txt)

    # 2) TSV → 两种格式
    tsv_txt = read_text_clean(in_tsv)
    tsv_A, tsv_B, dropped = prepare_two_inputs(tsv_txt, tips)
    write_text(sani / "cafe_input.family_id.tsv", tsv_A)
    write_text(sani / "cafe_input.desc.tsv",      tsv_B)
    print("[DEBUG] sanitized TSV (family_id) head:\n" + "\n".join(tsv_A.splitlines()[:2]))
    if dropped > 0:
        print(f"[INFO] 过滤全0家族：{dropped} 条。")

    # 3) 尝试 A：family_id 版
    cmdA = [ args.cafe5, "-i", str(sani / "cafe_input.family_id.tsv"),
             "-t", str(sani / "tree.normalized.nwk"),
             "-o", str(outdir), "-k", str(args.threads) ]
    print("[INFO] Try(A):", " ".join(shlex.quote(c) for c in cmdA))
    procA = run_cafe(cmdA, outdir)
    _ = collect_outputs(outdir)

    need_fallback = (not have_outputs(outdir)) or ("not found in gene family" in procA.stdout.lower())
    if need_fallback:
        print("[WARN] A 格式未成功（或检测到 'not found in gene family'），尝试回退到 B 格式（desc + 1-based）。")
        # 清理可能的空文件，避免误判
        for name in ("report.tsv","expanded.tsv","contracted.tsv","lambda.txt"):
            p = outdir / name
            if p.exists() and p.stat().st_size == 0:
                try: p.unlink()
                except: pass
        cmdB = [ args.cafe5, "-i", str(sani / "cafe_input.desc.tsv"),
                 "-t", str(sani / "tree.normalized.nwk"),
                 "-o", str(outdir), "-k", str(args.threads) ]
        print("[INFO] Try(B):", " ".join(shlex.quote(c) for c in cmdB))
        procB = run_cafe(cmdB, outdir)
        _ = collect_outputs(outdir)

        if not have_outputs(outdir):
            print("\n[DEBUG] 输出目录结构（前 200 项）：\n" + ls_tree(outdir, 200))
            log_head = "\n".join((outdir / "cafe5.stdout.log").read_text().splitlines()[:120])
            print("\n[DEBUG] cafe5.stdout.log（前 120 行）：\n" + textwrap.dedent(log_head))
            raise FileNotFoundError("缺少 CAFÉ 关键产物：report.tsv, expanded.tsv, contracted.tsv, lambda.txt")
    else:
        print("[INFO] A 格式已成功生成产物。")

    print("[OK] CAFE5 finished and outputs collected:",
          ", ".join(x for x in ("report.tsv","expanded.tsv","contracted.tsv","lambda.txt")
                   if (outdir / x).exists() and (outdir / x).stat().st_size > 0))

if __name__ == "__main__":
    main()

