#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
09_publish_results.py — 发布工作簿（仅用于查看，不参与下游）
输出：
  1) 每个 label（对比名或基因集名）一本 xlsx：GO/KEGG × All/Up/Down + 基因长/宽表 + Summary
  2) project_summary.xlsx：Quant_QC（Salmon 覆盖率）、Top_GO / Top_KEGG（跨 label 汇总）

约定：
  - 所有源表均来自 TSV（08 产物 & quant/summary.tsv）
  - 仅使用相对路径；config.yaml 覆盖脚本默认。
"""

import sys, os, re, yaml, pandas as pd
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

def read_tsv(fp):
    return pd.read_csv(fp, sep="\t", dtype=str).fillna("")

def read_tsv_num(fp):
    # 数值列自动类型；空值保持空串，避免科学计数法和 NA 污染
    df = pd.read_csv(fp, sep="\t")
    return df

def collect_enrich_tables(enrich_tables_dir: Path):
    """
    返回 dict[label][type]['All'/'Up'/'Down'] = DataFrame
      - type: 'GO' | 'KEGG'
      - label: go_enrich_<label>_<subset>.tsv / kegg_enrich_<label>_<subset>.tsv
    """
    m = defaultdict(lambda: {"GO":{}, "KEGG":{}})
    for fp in enrich_tables_dir.glob("*.tsv"):
        name = fp.name
        g = re.match(r"^(go|kegg)_enrich_(.+)_(All|Up|Down)\.tsv$", name, re.IGNORECASE)
        if not g: 
            continue
        typ = "GO" if g.group(1).lower()=="go" else "KEGG"
        label = g.group(2)
        subset = g.group(3)
        df = read_tsv_num(fp)
        m[label][typ][subset] = df
    return m

def collect_gene_tables(genes_dir: Path):
    """
    收集 08 生成的命中基因长/宽表（如果存在就写入 xlsx）
    约定文件名：
      go_genes_<Subset>_all_long.tsv / go_genes_<Subset>_all_byterm.tsv
      kegg_genes_<Subset>_all_long.tsv / ...
    这里与 label 无关（按 subset 聚合），因此仅作为补充 Sheet。
    """
    ret = {}
    for typ in ["go","kegg"]:
        for which in ["long","byterm"]:
            for subset in ["All","Up","Down"]:
                fp = genes_dir / f"{typ}_genes_{subset}_all_{which}.tsv"
                if fp.exists():
                    ret[f"{typ.upper()}_genes_{subset}_{which}"] = read_tsv(fp)
    return ret

def write_label_workbook(publish_dir: Path, label: str, bundle: dict, gene_sheets: dict):
    xlsx = publish_dir / f"{label}.xlsx"
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as w:
        # Summary（若能抽到 bg/sig 尺寸就放进去）
        summary_rows = []
        for typ in ["GO","KEGG"]:
            for subset, df in bundle.get(typ, {}).items():
                if df is None or df.empty: continue
                row = {
                    "type": typ, "subset": subset,
                    "bg_size": df.get("bg_size",[None]).iloc[0],
                    "sig_size": df.get("sig_size",[None]).iloc[0],
                    "n_terms": len(df)
                }
                summary_rows.append(row)
        if summary_rows:
            pd.DataFrame(summary_rows).to_excel(w, index=False, sheet_name="Summary")

        # 富集表
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                df = bundle.get(typ,{}).get(subset)
                if df is None or df.empty: continue
                df.to_excel(w, index=False, sheet_name=f"{typ}_{subset}")

        # 命中基因附表（如存在）
        for sheet, df in gene_sheets.items():
            if df is None or df.empty: continue
            df.to_excel(w, index=False, sheet_name=sheet[:31])  # Excel sheet名 ≤31
    return xlsx

def build_project_summary(publish_dir: Path, quant_dir: Path, collect: dict):
    xlsx = publish_dir / "project_summary.xlsx"
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as w:
        # Salmon 覆盖率
        qsum = quant_dir / "summary.tsv"
        if qsum.exists():
            pd.read_csv(qsum, sep="\t").to_excel(w, index=False, sheet_name="Quant_QC")

        # Top GO/KEGG（跨 label 合并，按 padj、pvalue 排序）
        rows_go, rows_kegg = [], []
        for label, bundle in collect.items():
            for typ, store in bundle.items():
                for subset, df in store.items():
                    if df is None or df.empty: continue
                    tmp = df.copy()
                    tmp.insert(0, "label", label)
                    tmp.insert(1, "subset", subset)
                    if typ=="GO":
                        rows_go.append(tmp)
                    else:
                        rows_kegg.append(tmp)
        if rows_go:
            go_all = pd.concat(rows_go, ignore_index=True)
            go_all = go_all.sort_values(["padj","pvalue"], na_position="last")
            go_all.to_excel(w, index=False, sheet_name="Top_GO")
        if rows_kegg:
            kegg_all = pd.concat(rows_kegg, ignore_index=True)
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

    enrich_dir   = proj / paths.get("enrich_dir", "results/enrich")
    enrich_tables_dir = enrich_dir / "tables"
    genes_dir    = enrich_dir / "genes"
    quant_dir    = proj / paths.get("quant_dir", "results/quant")
    publish_dir  = proj / paths.get("publish_dir","results/publish")
    publish_dir.mkdir(parents=True, exist_ok=True)

    collect = collect_enrich_tables(enrich_tables_dir)
    if not collect:
        print("[ERR] 未发现富集表（results/enrich/tables/*.tsv）", file=sys.stderr); sys.exit(2)

    # 附表：基因长/宽表（若存在）
    gene_sheets = collect_gene_tables(genes_dir)

    # 逐 label 写工作簿
    created = []
    for label, bundle in collect.items():
        xlsx = write_label_workbook(publish_dir, label, bundle, gene_sheets)
        created.append(xlsx)

    # 项目总览：含 Salmon 覆盖率
    summary_xlsx = build_project_summary(publish_dir, quant_dir, collect)

    print("[OK] 生成工作簿：")
    for x in created:
        print(" -", x)
    print(" -", summary_xlsx)

if __name__ == "__main__":
    main()

