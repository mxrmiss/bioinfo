#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py —— 整理 eggNOG-mapper 注释为“基因层映射”（手工前/后缀契约版）

契约要点：
  * 仅按 config.yaml 的“手工前/后缀开关 + 字面量列表”清理 emapper 的 query；不开关则不改动。
  * 不做任何“自动去数字/正则猜测/智能识别”。
  * 优先按（清理后）query == transcript_id 命中 → gene_id；若 allow_gene_space_match=true，则允许 query 直接等于 gene_id 命中。
  * 列名默认：query, GOs, KEGG_ko, KEGG_Pathway；可在 config.annotations.*_col 覆盖。
  * 输出：
      - results/07_annot/gene2go.tsv（gene_id, go_id）
      - results/07_annot/gene2ko.tsv（gene_id, ko_id）
      - results/07_annot/gene2pathway.tsv（gene_id, pathway_id）
      - results/07_annot/universe_coverage.tsv（覆盖率审计）
      - results/07_annot/unmapped_queries.list（未命中原始 query）
      - results/07_annot/id_cleanup.audit.tsv（raw→prefix→suffix 的命中率）
"""

from __future__ import annotations
import sys, csv, gzip, logging, re
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml

CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "reference": {"emapper": ""},
    "dirs": {"maps": "results/03_maps", "annotations": "results/07_annot"},
    "annotations": {
        "emapper_query_col": "",
        "emapper_gos_col": "",
        "emapper_ko_col": "",
        "emapper_pathway_col": "",
        "allow_gene_space_match": True,
        "id_cleanup": {
            "strip_prefix": False,
            "prefix": [],
            "strip_suffix": False,
            "suffix": [],
            "order": ["prefix","suffix"]  # 默认先前缀再后缀
        }
    },
    "logging": {"level": "INFO", "timestamp": True}
}

_SPLIT_RE = re.compile(r"[,\|\s;]+")

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def open_maybe_gz(p: Path):
    return gzip.open(p, "rt", encoding="utf-8", errors="ignore") if str(p).lower().endswith(".gz") \
           else open(p, "r", encoding="utf-8", errors="ignore")

def load_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        user = (yaml.safe_load(f) or {})
    def merge(u, d):
        out = dict(d)
        for k, v in u.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out
    return merge(user, DEFAULTS)

def clean_id_literal(x: str, prefixes: List[str], suffixes: List[str], order: List[str]) -> str:
    y = x
    for step in order:
        if step == "prefix":
            for pfx in prefixes:
                if y.startswith(pfx):
                    y = y[len(pfx):]
                    break
        elif step == "suffix":
            for sfx in suffixes:
                if y.endswith(sfx):
                    y = y[:-len(sfx)]
                    break
    return y

def normalize_ko(tok: str) -> Optional[str]:
    if not tok: return None
    s = tok.strip().upper()
    return s if re.fullmatch(r"K\d{5}", s) else None

def normalize_pathway(tok: str) -> Optional[str]:
    if not tok: return None
    s = tok.strip().replace("path:", "")
    if s.startswith("map"): s = "ko" + s[3:]
    return s if re.fullmatch(r"ko\d{5}", s) else None

def main() -> None:
    cfg = load_yaml(Path(CONFIG_PATH))
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True) else "[%(levelname)s] %(message)s",
    )

    logging.info("========== 07 — emapper 注释整理（手工前/后缀契约） ==========")

    # 路径
    emapper_fp = Path((cfg.get("reference", {}).get("emapper") or "").strip())
    maps_dir   = Path(cfg["dirs"]["maps"])
    out_dir    = Path(cfg["dirs"]["annotations"])
    mkdir_p(out_dir)
    if not emapper_fp.exists():
        raise FileNotFoundError(f"未找到 emapper 注释文件：{emapper_fp}")

    # 读取 tx2gene.clean.tsv
    tx2gene_fp = maps_dir / "tx2gene.clean.tsv"
    if not tx2gene_fp.exists():
        raise FileNotFoundError(f"未找到 tx2gene.clean.tsv：{tx2gene_fp}")
    t2g: Dict[str, str] = {}
    gene_set: set = set()
    with open(tx2gene_fp, "r", encoding="utf-8") as f:
        rdr = csv.DictReader(f, dialect=csv.excel_tab)
        need = ["transcript_id", "gene_id"]
        for c in need:
            if c not in (rdr.fieldnames or []):
                raise ValueError(f"tx2gene.clean.tsv 缺少必要列：{c}")
        for r in rdr:
            tid = (r["transcript_id"] or "").strip()
            gid = (r["gene_id"] or "").strip()
            if tid and gid:
                t2g[tid] = gid
                gene_set.add(gid)

    # 清理配置（仅按手工列表；不开关就不动）
    icfg = cfg["annotations"]["id_cleanup"]
    do_prefix = bool(icfg.get("strip_prefix", False))
    do_suffix = bool(icfg.get("strip_suffix", False))
    order    = list(icfg.get("order", ["prefix","suffix"]))
    pfx_list = list(icfg.get("prefix", [])) if do_prefix else []
    sfx_list = list(icfg.get("suffix", [])) if do_suffix else []
    allow_gene_space = bool(cfg["annotations"].get("allow_gene_space_match", True))
    logging.info("[Rule] 顺序=%s；strip_prefix=%s；strip_suffix=%s；allow_gene_space_match=%s",
                 order, do_prefix, do_suffix, allow_gene_space)
    if pfx_list: logging.info("[Rule] prefix=%s", pfx_list)
    if sfx_list: logging.info("[Rule] suffix=%s", sfx_list)

    # 解析 emapper 表头（仅认 '#query\t...' 行）
    header: Optional[List[str]] = None
    idx: Dict[str, int] = {}
    want = {
        "query": (cfg["annotations"].get("emapper_query_col") or "query"),
        "gos": (cfg["annotations"].get("emapper_gos_col") or "GOs"),
        "ko": (cfg["annotations"].get("emapper_ko_col") or "KEGG_ko"),
        "path": (cfg["annotations"].get("emapper_pathway_col") or "KEGG_Pathway"),
    }

    def locate(h: List[str]) -> Dict[str, int]:
        lo = [x.lower() for x in h]
        out = { "query": -1, "gos": -1, "ko": -1, "path": -1 }
        for k, nm in want.items():
            for i, v in enumerate(lo):
                if v == nm.lower():
                    out[k] = i; break
        return out

    # 先收集所有 query（便于审计 raw→prefix→suffix 的命中率）
    queries_raw: List[str] = []
    with open_maybe_gz(emapper_fp) as f:
        for line in f:
            if not line: continue
            s = line.rstrip("\n")
            if header is None:
                if s.startswith("#query\t") or s.startswith("#Query\t") or s.startswith("#QUERY\t"):
                    header = s[1:].split("\t")
                    idx = locate(header)
                    if idx["query"] < 0:
                        raise RuntimeError("未在 emapper 表头中找到 'query' 列")
                continue
            if s.startswith("#") or not s:
                continue
            parts = s.split("\t")
            queries_raw.append(parts[idx["query"]].strip())

    def apply_stage(x: str, stage: str) -> str:
        if stage == "prefix" and do_prefix:
            for pfx in pfx_list:
                if x.startswith(pfx):
                    return x[len(pfx):]
        if stage == "suffix" and do_suffix:
            for sfx in sfx_list:
                if x.endswith(sfx):
                    return x[:-len(sfx)]
        return x

    # 审计：raw→prefix→suffix 命中率（命中 = 命中任一集合 {transcript_id, gene_id(可选)}）
    audit_rows: List[List[str]] = [["stage","matched","total","match_rate"]]
    for stage in ["raw"] + order:
        matched = 0
        for q in queries_raw:
            qq = apply_stage(q, stage) if stage != "raw" else q
            if (qq in t2g) or (allow_gene_space and (qq in gene_set)):
                matched += 1
        total = len(queries_raw) if queries_raw else 1
        rate = matched / total
        audit_rows.append([stage, str(matched), str(total), f"{rate:.6f}"])

    # 正式产出
    def final_clean(x: str) -> str:
        return clean_id_literal(x, pfx_list if do_prefix else [], sfx_list if do_suffix else [], order)

    gene2go: Dict[str, set] = {}
    gene2ko: Dict[str, set] = {}
    gene2path: Dict[str, set] = {}
    mapped_tx = mapped_gene = unmapped = 0
    unmapped_examples: List[str] = []

    with open_maybe_gz(emapper_fp) as f:
        qcol = gos_col = ko_col = path_col = None
        for line in f:
            if not line: continue
            s = line.rstrip("\n")
            if qcol is None:
                if s.startswith("#query\t") or s.startswith("#Query\t") or s.startswith("#QUERY\t"):
                    hdr = s[1:].split("\t")
                    if header is None:
                        header = hdr
                    idx = locate(hdr)
                    qcol, gos_col, ko_col, path_col = idx["query"], idx["gos"], idx["ko"], idx["path"]
                    logging.info("[Info] 表头命中：query=%s, GOs=%s, KEGG_ko=%s, KEGG_Pathway=%s",
                                 hdr[qcol] if qcol>=0 else "NA",
                                 hdr[gos_col] if gos_col>=0 else "NA",
                                 hdr[ko_col]  if ko_col>=0  else "NA",
                                 hdr[path_col]if path_col>=0 else "NA")
                    continue
                else:
                    continue
            if s.startswith("#") or not s:
                continue

            parts = s.split("\t")
            q_raw = parts[qcol].strip()
            q = final_clean(q_raw)

            gid: Optional[str] = None
            if q in t2g:
                gid = t2g[q]; mapped_tx += 1
            elif allow_gene_space and (q in gene_set):
                gid = q; mapped_gene += 1
            else:
                unmapped += 1
                if len(unmapped_examples) < 50000:
                    unmapped_examples.append(q_raw)
                continue

            if gos_col is not None and gos_col >= 0 and gos_col < len(parts):
                s_g = parts[gos_col].strip()
                if s_g and s_g != "-":
                    for tok in _SPLIT_RE.split(s_g):
                        if tok and tok.startswith("GO:") and re.fullmatch(r"GO:\d{7}", tok):
                            gene2go.setdefault(gid, set()).add(tok)

            if ko_col is not None and ko_col >= 0 and ko_col < len(parts):
                s_k = parts[ko_col].strip()
                if s_k and s_k != "-":
                    for tok in _SPLIT_RE.split(s_k):
                        k = normalize_ko(tok)
                        if k:
                            gene2ko.setdefault(gid, set()).add(k)

            if path_col is not None and path_col >= 0 and path_col < len(parts):
                s_p = parts[path_col].strip()
                if s_p and s_p != "-":
                    for tok in _SPLIT_RE.split(s_p):
                        p = normalize_pathway(tok)
                        if p:
                            gene2path.setdefault(gid, set()).add(p)

    out_dir = Path(cfg["dirs"]["annotations"])
    mkdir_p(out_dir)

    def write_pairs(fp: Path, hdr: List[str], d: Dict[str, set]) -> None:
        with open(fp, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, dialect=csv.excel_tab)
            w.writerow(hdr)
            for k in sorted(d.keys()):
                for v in sorted(d[k]):
                    w.writerow([k, v])
        logging.info("[Out] %s", fp)

    g2go_fp   = out_dir / "gene2go.tsv"
    g2ko_fp   = out_dir / "gene2ko.tsv"
    g2path_fp = out_dir / "gene2pathway.tsv"
    write_pairs(g2go_fp, ["gene_id","go_id"], gene2go)
    write_pairs(g2ko_fp, ["gene_id","ko_id"], gene2ko)
    write_pairs(g2path_fp, ["gene_id","pathway_id"], gene2path)

    cov_fp   = out_dir / "universe_coverage.tsv"
    with open(cov_fp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["metric","value"])
        w.writerow(["emapper_data_rows", str(len(queries_raw))])
        w.writerow(["mapped_via_transcript", str(mapped_tx)])
        w.writerow(["mapped_via_gene", str(mapped_gene)])
        w.writerow(["unmapped", str(unmapped)])
        uniq_genes = len(set(list(gene2go.keys()) + list(gene2ko.keys()) + list(gene2path.keys())))
        w.writerow(["unique_genes_mapped", str(uniq_genes)])
    logging.info("[Out] %s", cov_fp)

    audit_fp = out_dir / "id_cleanup.audit.tsv"
    with open(audit_fp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerows(audit_rows)
    logging.info("[Out] %s", audit_fp)

    if unmapped_examples:
        with open(out_dir / "unmapped_queries.list", "w", encoding="utf-8") as f:
            f.write("\n".join(unmapped_examples))
        logging.info("[Out] %s", out_dir / "unmapped_queries.list")

    logging.info("========== 07 完成 ==========")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 07 执行失败：{e}", file=sys.stderr)
        sys.exit(1)

