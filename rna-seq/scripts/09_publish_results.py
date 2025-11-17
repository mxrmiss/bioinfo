#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
09_publish_results.py — 发布工作簿（仅用于查看，不参与下游）
输出：
  1) 每个 label（对比名或基因集名）一本 xlsx：GO/KEGG × All/Up/Down + 命中明细（若存在） + Summary
  2) project_summary.xlsx：Quant_QC（Salmon 覆盖率）、Top_GO / Top_KEGG（跨 label 汇总）

与 08 的产物完全对齐：
- tables/ 下文件名可大小写混用（all/All/up/Up/down/Down），此脚本统一规范为 All/Up/Down。
- genes/ 下命名采用：GO|KEGG_{all|up|down}__{hit_in_set|hit_in_sig}__<label>.tsv
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

# ---------- 小工具 ----------
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
    """所有列按字符串读入，避免 NA/科学计数污染"""
    return pd.read_csv(fp, sep="\t", dtype=str).fillna("")

def read_tsv_num(fp):
    """数值列保留原类型，用于 padj/pvalue 排序"""
    return pd.read_csv(fp, sep="\t")

def norm_subset(s: str) -> str:
    """把 all/up/down（任意大小写）统一为 All/Up/Down"""
    x = (s or "").strip().lower()
    if x == "all":  return "All"
    if x == "up":   return "Up"
    if x == "down": return "Down"
    return s  # 未识别时原样返回（后续一般会被忽略）

# ---------- 收集 08 的富集表 ----------
def collect_enrich_tables(enrich_tables_dir: Path):
    """
    返回 dict[label][type]['All'/'Up'/'Down'] = DataFrame
      - type: 'GO' | 'KEGG'
      - label: go_enrich_<label>_<subset>.tsv / kegg_enrich_<label>_<subset>.tsv
    """
    m = defaultdict(lambda: {"GO":{}, "KEGG":{}})

    if not enrich_tables_dir.exists():
        raise FileNotFoundError(f"[ERR] 未找到富集表目录：{enrich_tables_dir}")

    files = list(enrich_tables_dir.glob("*.tsv"))
    if not files:
        raise FileNotFoundError(f"[ERR] 目录为空：{enrich_tables_dir}（期待 go_enrich_*.tsv / kegg_enrich_*.tsv）")

    pat = re.compile(r"^(go|kegg)_enrich_(.+)_(All|Up|Down|all|up|down)\.tsv$", re.IGNORECASE)
    matched = 0
    for fp in files:
        g = pat.match(fp.name)
        if not g:
            continue
        typ = "GO" if g.group(1).lower()=="go" else "KEGG"
        label = g.group(2)
        subset = norm_subset(g.group(3))
        try:
            df = read_tsv_num(fp)
        except Exception as e:
            raise RuntimeError(f"[ERR] 读取富集表失败：{fp} ；原因：{e}")
        m[label][typ][subset] = df
        matched += 1

    if matched == 0:
        raise FileNotFoundError(f"[ERR] 未匹配到任何富集表文件，请检查命名：{enrich_tables_dir}")

    return m

# ---------- 收集 genes/ 命中明细（按当前命名规范） ----------
def collect_gene_hits(genes_dir: Path):
    """
    返回 dict[label][type]['All'/'Up'/'Down']['hit_in_set'|'hit_in_sig'] = DataFrame
    文件命名：GO|KEGG_{all|up|down}__{hit_in_set|hit_in_sig}__<label>.tsv
    """
    hits = defaultdict(lambda: {"GO": {"All":{}, "Up":{}, "Down":{}},
                               "KEGG":{"All":{}, "Up":{}, "Down":{}}})
    if not genes_dir.exists():
        # 没有 genes 目录不报错，视为可选
        return hits

    files = list(genes_dir.glob("*.tsv"))
    if not files:
        return hits

    pat = re.compile(r"^(GO|KEGG)_(all|up|down)__(hit_in_set|hit_in_sig)__(.+)\.tsv$", re.IGNORECASE)
    for fp in files:
        m = pat.match(fp.name)
        if not m:
            continue
        typ = m.group(1).upper()
        subset = norm_subset(m.group(2))
        which = m.group(3)  # 命名维持小写
        label = m.group(4)
        try:
            df = read_tsv(fp)
        except Exception as e:
            raise RuntimeError(f"[ERR] 读取命中明细失败：{fp} ；原因：{e}")
        hits[label][typ][subset][which] = df

    return hits

# ---------- 写单本工作簿 ----------
def write_label_workbook(publish_dir: Path, label: str, bundle: dict, gene_hits: dict):
    xlsx = publish_dir / f"{label}.xlsx"
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as w:
        # Summary
        summary_rows = []
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                df = bundle.get(typ, {}).get(subset)
                if df is None or len(df)==0:
                    continue
                row = {
                    "type": typ, "subset": subset,
                    "bg_size": df.get("bg_size",[None]).iloc[0] if "bg_size" in df.columns else None,
                    "sig_size": df.get("sig_size",[None]).iloc[0] if "sig_size" in df.columns else None,
                    "n_terms": len(df)
                }
                summary_rows.append(row)
        if summary_rows:
            pd.DataFrame(summary_rows).to_excel(w, index=False, sheet_name="Summary")

        # 富集表（GO/KEGG × All/Up/Down）
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                df = bundle.get(typ,{}).get(subset)
                if df is None or len(df)==0:
                    continue
                df.to_excel(w, index=False, sheet_name=f"{typ}_{subset}")

        # 命中明细（若存在则写入）
        gstore = gene_hits.get(label, {})
        for typ in ["GO","KEGG"]:
            for subset in ["All","Up","Down"]:
                for which in ["hit_in_set","hit_in_sig"]:
                    df = (((gstore.get(typ, {})).get(subset, {})).get(which))
                    if df is None or len(df)==0:
                        continue
                    # Sheet 名 ≤31 字符
                    sheet = f"{typ}_{subset}_{which}"
                    if len(sheet) > 31:
                        sheet = sheet[:31]
                    df.to_excel(w, index=False, sheet_name=sheet)

    return xlsx

# ---------- 项目总览 ----------
def build_project_summary(publish_dir: Path, quant_dir: Path, collect: dict):
    xlsx = publish_dir / "project_summary.xlsx"
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as w:
        # Salmon 覆盖率
        qsum = quant_dir / "summary.tsv"
        if qsum.exists():
            try:
                pd.read_csv(qsum, sep="\t").to_excel(w, index=False, sheet_name="Quant_QC")
            except Exception as e:
                raise RuntimeError(f"[ERR] 读取 quant/summary.tsv 失败：{qsum} ；原因：{e}")

        # Top GO/KEGG（跨 label 合并，按 padj、pvalue 排序）
        rows_go, rows_kegg = [], []
        for label, bundle in collect.items():
            for typ, store in bundle.items():
                for subset, df in store.items():
                    if df is None or len(df)==0:
                        continue
                    tmp = df.copy()
                    tmp.insert(0, "label", label)
                    tmp.insert(1, "subset", subset)
                    if typ=="GO":
                        rows_go.append(tmp)
                    else:
                        rows_kegg.append(tmp)
        if rows_go:
            go_all = pd.concat(rows_go, ignore_index=True)
            go_all = go_all.sort_values(["padj","pvalue"], na_position="last") if "padj" in go_all.columns and "pvalue" in go_all.columns else go_all
            go_all.to_excel(w, index=False, sheet_name="Top_GO")
        if rows_kegg:
            kegg_all = pd.concat(rows_kegg, ignore_index=True)
            kegg_all = kegg_all.sort_values(["padj","pvalue"], na_position="last") if "padj" in kegg_all.columns and "pvalue" in kegg_all.columns else kegg_all
            kegg_all.to_excel(w, index=False, sheet_name="Top_KEGG")
    return xlsx

# ---------- 主程 ----------
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

    # —— 输入检测（清晰报错）——
    if not enrich_dir.exists():
        print(f"[ERR] 未找到富集结果目录：{enrich_dir}（请先完成 08）", file=sys.stderr); sys.exit(2)
    if not enrich_tables_dir.exists():
        print(f"[ERR] 未找到富集表目录：{enrich_tables_dir}（请确认 08 是否成功）", file=sys.stderr); sys.exit(2)

    # 读取富集表（必须有）
    try:
        collect = collect_enrich_tables(enrich_tables_dir)
    except Exception as e:
        print(str(e), file=sys.stderr); sys.exit(2)

    # 读取命中明细（可选）
    try:
        gene_hits = collect_gene_hits(genes_dir)
    except Exception as e:
        print(str(e), file=sys.stderr); sys.exit(3)

    # 逐 label 写工作簿
    created = []
    for label, bundle in collect.items():
        xlsx = write_label_workbook(publish_dir, label, bundle, gene_hits)
        created.append(xlsx)

    # 项目总览
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

