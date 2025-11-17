#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
09_publish_results.py — 发布工作簿（装订层，不改统计口径）
新增：强制把 BP/CC/MF & KEGG 的 by-term 宽表写入每本 <label>.xlsx（若文件存在）。
"""

import sys, re, yaml, pandas as pd
from pathlib import Path
from collections import defaultdict

DEFAULTS = {
    "paths": {
        "enrich_dir": "results/enrich",
        "quant_dir":  "results/quant",
        "deg_dir":    "results/deg",
        "logs_dir":   "logs",
        "publish_dir":"results/publish"
    }
}

def load_yaml(fp):
    with open(fp, "r", encoding="utf-8") as r:
        return yaml.safe_load(r) or {}

def deep_merge(a,b):
    out=dict(a)
    for k,v in (b or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k]=deep_merge(out[k],v)
        else:
            out[k]=v
    return out

def read_tsv(fp):  # 全字符串
    return pd.read_csv(fp, sep="\t", dtype=str).fillna("")

def read_tsv_num(fp):  # 数值列自动
    return pd.read_csv(fp, sep="\t")

def norm_subset(s: str) -> str:
    x = (s or "").strip().lower()
    if x == "all":  return "All"
    if x == "up":   return "Up"
    if x == "down": return "Down"
    return s

# --------- 收集富集表 ----------
def collect_enrich_tables(enrich_tables_dir: Path):
    m = defaultdict(lambda: {"GO":{}, "KEGG":{}})
    if not enrich_tables_dir.exists():
        raise FileNotFoundError(f"[ERR] 未找到富集表目录：{enrich_tables_dir}")
    files = list(enrich_tables_dir.glob("*.tsv"))
    if not files:
        raise FileNotFoundError(f"[ERR] 目录为空：{enrich_tables_dir}")

    pat = re.compile(r"^(go|kegg)_enrich_(.+)_(All|all|Up|up|Down|down)\.tsv$")
    matched = 0
    for fp in files:
        g = pat.match(fp.name)
        if not g: continue
        typ = "GO" if g.group(1).lower()=="go" else "KEGG"
        label = g.group(2)
        subset = norm_subset(g.group(3))
        df = read_tsv_num(fp)
        m[label][typ][subset] = df
        matched += 1
    if matched == 0:
        raise FileNotFoundError(f"[ERR] 未匹配到任何富集表：{enrich_tables_dir}")
    return m

# --------- 收集命中长表（可选） ----------
def collect_gene_hits(genes_dir: Path):
    hits = defaultdict(lambda: {"GO":{"All":{}, "Up":{}, "Down":{}},
                                "KEGG":{"All":{}, "Up":{}, "Down":{}}})
    if not genes_dir.exists():
        return hits
    pat = re.compile(r"^(GO|KEGG)_(all|up|down)__(hit_in_set|hit_in_sig)__(.+)\.tsv$", re.IGNORECASE)
    for fp in genes_dir.glob("*.tsv"):
        m = pat.match(fp.name)
        if not m: continue
        typ   = m.group(1).upper()
        subset= norm_subset(m.group(2))
        which = m.group(3)
        label = m.group(4)
        hits[label][typ][subset][which] = read_tsv(fp)
    return hits

# --------- 收集 by-term 宽表（强制装订，若存在就写） ----------
def collect_byterm(genes_dir: Path):
    """
    GO_BP/GO_CC/GO_MF_{all|up|down}_genes_by_term.tsv
    KEGG_pathway_{all|up|down}_genes_by_term.tsv
    与 label 无关（by-term 是跨 label 的视图），直接作为附表写入每本工作簿。
    """
    store = {"GO_BP":{}, "GO_CC":{}, "GO_MF":{}, "KEGG":{}}
    if not genes_dir.exists():
        return store
    for subset in ["All","Up","Down"]:
        # GO 三域
        for dom, prefix in [("GO_BP", "GO_BP"), ("GO_CC","GO_CC"), ("GO_MF","GO_MF")]:
            fp = genes_dir / f"{prefix}_{subset.lower()}_genes_by_term.tsv"
            if fp.exists(): store[dom][subset] = read_tsv(fp)
        # KEGG
        fp_k = genes_dir / f"KEGG_pathway_{subset.lower()}_genes_by_term.tsv"
        if fp_k.exists(): store["KEGG"][subset] = read_tsv(fp_k)
    return store

# --------- 写工作簿 ----------
def write_label_workbook(publish_dir: Path, label: str, bundle: dict, gene_hits: dict, byterm: dict):
    xlsx = publish_dir / f"{label}.xlsx"
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as w:
        # Summary
        summary_rows = []
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                df = bundle.get(typ, {}).get(subset)
                if df is None or len(df)==0: continue
                row = {
                    "type": typ, "subset": subset,
                    "bg_size": df.get("bg_size",[None]).iloc[0] if "bg_size" in df.columns else None,
                    "sig_size": df.get("sig_size",[None]).iloc[0] if "sig_size" in df.columns else None,
                    "n_terms": len(df)
                }
                summary_rows.append(row)
        if summary_rows:
            pd.DataFrame(summary_rows).to_excel(w, index=False, sheet_name="Summary")

        # 富集表
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                df = bundle.get(typ,{}).get(subset)
                if df is None or len(df)==0: continue
                df.to_excel(w, index=False, sheet_name=f"{typ}_{subset}")

        # 命中长表（若存在）
        gstore = gene_hits.get(label, {})
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                for which in ["hit_in_set","hit_in_sig"]:
                    df = (((gstore.get(typ, {})).get(subset, {})).get(which))
                    if df is None or len(df)==0: continue
                    sheet = f"{typ}_{subset}_{which}"
                    if len(sheet) > 31: sheet = sheet[:31]
                    df.to_excel(w, index=False, sheet_name=sheet)

        # by-term 宽表（强制写：若磁盘存在就写）
        for subset in ["All","Up","Down"]:
            for dom, sh in [("GO_BP","GO_BP"), ("GO_CC","GO_CC"), ("GO_MF","GO_MF")]:
                df = byterm.get(dom,{}).get(subset)
                if df is not None and len(df)>0:
                    name = f"{sh}_{subset}_by_term"
                    if len(name)>31: name = name[:31]
                    df.to_excel(w, index=False, sheet_name=name)
            dfk = byterm.get("KEGG",{}).get(subset)
            if dfk is not None and len(dfk)>0:
                name = f"KEGG_{subset}_by_term"
                if len(name)>31: name = name[:31]
                dfk.to_excel(w, index=False, sheet_name=name)

    return xlsx

# --------- 总览 ----------
def build_project_summary(publish_dir: Path, quant_dir: Path, collect: dict):
    xlsx = publish_dir / "project_summary.xlsx"
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as w:
        qsum = quant_dir / "summary.tsv"
        if qsum.exists():
            pd.read_csv(qsum, sep="\t").to_excel(w, index=False, sheet_name="Quant_QC")
        rows_go, rows_kegg = [], []
        for label, bundle in collect.items():
            for typ, store in bundle.items():
                for subset, df in store.items():
                    if df is None or len(df)==0: continue
                    tmp = df.copy()
                    tmp.insert(0, "label", label)
                    tmp.insert(1, "subset", subset)
                    (rows_go if typ=="GO" else rows_kegg).append(tmp)
        if rows_go:
            go_all = pd.concat(rows_go, ignore_index=True)
            if "padj" in go_all.columns and "pvalue" in go_all.columns:
                go_all = go_all.sort_values(["padj","pvalue"], na_position="last")
            go_all.to_excel(w, index=False, sheet_name="Top_GO")
        if rows_kegg:
            kegg_all = pd.concat(rows_kegg, ignore_index=True)
            if "padj" in kegg_all.columns and "pvalue" in kegg_all.columns:
                kegg_all = kegg_all.sort_values(["padj","pvalue"], na_position="last")
            kegg_all.to_excel(w, index=False, sheet_name="Top_KEGG")
    return xlsx

def main():
    proj = Path.cwd()
    cfg_fp = proj / "config.yaml"
    if not cfg_fp.exists():
        print("[ERR] 缺少 config.yaml（请在项目根目录运行）", file=sys.stderr); sys.exit(1)
    cfg = deep_merge(DEFAULTS, load_yaml(cfg_fp))
    paths = cfg.get("paths", {})

    enrich_dir        = proj / paths.get("enrich_dir", "results/enrich")
    enrich_tables_dir = enrich_dir / "tables"
    genes_dir         = enrich_dir / "genes"
    quant_dir         = proj / paths.get("quant_dir", "results/quant")
    publish_dir       = proj / paths.get("publish_dir","results/publish")
    publish_dir.mkdir(parents=True, exist_ok=True)

    if not enrich_dir.exists():
        print(f"[ERR] 未找到富集结果目录：{enrich_dir}（请先完成 08）", file=sys.stderr); sys.exit(2)
    if not enrich_tables_dir.exists():
        print(f"[ERR] 未找到富集表目录：{enrich_tables_dir}（请确认 08 是否成功）", file=sys.stderr); sys.exit(2)

    try:
        collect = collect_enrich_tables(enrich_tables_dir)
    except Exception as e:
        print(str(e), file=sys.stderr); sys.exit(2)

    # 可选：命中长表
    try:
        gene_hits = collect_gene_hits(genes_dir)
    except Exception as e:
        print(str(e), file=sys.stderr); sys.exit(3)

    # 强制：by-term 宽表（若存在就装订）
    byterm = collect_byterm(genes_dir)

    created = []
    for label, bundle in collect.items():
        xlsx = write_label_workbook(publish_dir, label, bundle, gene_hits, byterm)
        created.append(xlsx)

    try:
        summary_xlsx = build_project_summary(publish_dir, quant_dir, collect)
    except Exception as e:
        print(str(e), file=sys.stderr); sys.exit(4)

    print("[OK] 生成工作簿：")
    for x in created:
        print(" -", x)
    print(" -", summary_xlsx)

if __name__ == "__main__":
    main()

