#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py —— emapper 注释整理（契约：reference.emapper）
要点：
  * 只读 config.yaml；只认 reference.emapper（不找别名）
  * 用 tx2gene.clean.tsv 的 transcript_id ↔ emapper 的 query 映射
  * 输出：gene2go.tsv / gene2ko.tsv / gene2pathway.tsv / universe_coverage.tsv
  * 表头自愈（两处）：#query→query；KEGG_PathwayKEGG_Module → KEGG_Pathway + KEGG_Module
  * 可选的 query 前后缀清理（默认关闭，需在 config.emapper.query_id_clean.enable: true 才启用）
"""

from __future__ import annotations
import sys, csv, gzip, logging, re
from pathlib import Path
from typing import Dict, Any, List, Tuple

CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "dirs": {"maps": "results/03_maps", "annot": "results/07_annot", "logs": "logs"},
    "reference": {"emapper": "ref/annotations/annotations.tsv"},
    "emapper": {
        "query_id_clean": {
            "enable": False,
            "remove_prefixes": [],
            "remove_suffixes": [],
            "regex_suffixes": []  # 例：r"\.[0-9]+$", r"_t\d+$"
        }
    },
    "logging": {"level": "INFO", "timestamp": True}
}

def load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML：mamba/conda install pyyaml", file=sys.stderr); raise
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

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def norm_path(s: str) -> Path:
    val = (s or "").strip()
    return Path(val).expanduser() if val else Path("")

def show_resolved(p: Path) -> str:
    try:
        return str(p.resolve())
    except Exception:
        return str(p)

def open_auto(fp: Path):
    return gzip.open(fp, "rt") if str(fp).endswith(".gz") else open(fp, "r", encoding="utf-8")

def clean_query(q: str, cfg: Dict[str, Any]) -> str:
    qc = cfg.get("emapper", {}).get("query_id_clean", {})
    if not bool(qc.get("enable", False)):
        return q
    s = q
    for pre in qc.get("remove_prefixes", []) or []:
        if s.startswith(pre):
            s = s[len(pre):]
    for suf in qc.get("remove_suffixes", []) or []:
        if s.endswith(suf):
            s = s[: -len(suf)]
    for pat in qc.get("regex_suffixes", []) or []:
        try:
            s = re.sub(pat, "", s)
        except re.error:
            pass
    return s

def read_tx2gene(tx2: Path) -> Dict[str, str]:
    if not tx2.exists():
        raise FileNotFoundError(f"未找到 tx2gene.clean.tsv：{tx2}")
    m: Dict[str, str] = {}
    with open(tx2, "r", encoding="utf-8") as f:
        rdr = csv.DictReader(f, dialect=csv.excel_tab)
        need = ["transcript_id", "gene_id"]
        if rdr.fieldnames is None or any(c not in rdr.fieldnames for c in need):
            raise ValueError(f"tx2gene.clean.tsv 必须包含列：{', '.join(need)}（现有列：{rdr.fieldnames}）")
        for r in rdr:
            t = (r.get("transcript_id") or "").strip()
            g = (r.get("gene_id") or "").strip()
            if t:
                m[t] = g
    return m

def fix_header_line(line: str) -> Tuple[str, List[str]]:
    notes: List[str] = []
    s = line.rstrip("\n")
    if s.startswith("#query\t"):
        s = s[1:]; notes.append("header: #query -> query")
    if "KEGG_PathwayKEGG_Module" in s:
        s = s.replace("KEGG_PathwayKEGG_Module", "KEGG_Pathway\tKEGG_Module")
        notes.append("header: KEGG_PathwayKEGG_Module -> KEGG_Pathway + KEGG_Module")
    return s + "\n", notes

def main() -> None:
    # 读取契约配置（臣妾已阅读约定文件并按其键名实现）
    cfg = load_yaml(Path(CONFIG_PATH))
    # 日志
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True)
        else "[%(levelname)s] %(message)s")

    logging.info("========== 07 — emapper 注释整理（reference.emapper） ==========")

    emapper_fp = norm_path(cfg["reference"]["emapper"])
    maps_dir   = norm_path(cfg["dirs"]["maps"])
    out_dir    = norm_path(cfg["dirs"]["annot"])
    mkdir_p(out_dir); mkdir_p(out_dir / "_diag")

    logging.info(f"[Info] reference.emapper = {show_resolved(emapper_fp)}")
    logging.info(f"[Info] dirs.maps         = {show_resolved(maps_dir)}")
    logging.info(f"[Info] dirs.annot        = {show_resolved(out_dir)}")

    if not emapper_fp or not emapper_fp.exists():
        raise FileNotFoundError(f"未找到 emapper 注释文件：{emapper_fp}")

    tx2gene_fp = maps_dir / "tx2gene.clean.tsv"
    tx2map = read_tx2gene(tx2gene_fp)
    logging.info(f"[Info] tx2gene 映射条目数 = {len(tx2map)}")

    header_fixed_notes: List[str] = []
    header: List[str] = []
    data_lines: List[str] = []

    with open_auto(emapper_fp) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if not header:
                fixed, notes = fix_header_line(line)
                header_fixed_notes.extend(notes)
                header = fixed.rstrip("\n").split("\t")
                continue
            data_lines.append(line.rstrip("\n"))

    required = ["query", "GOs", "KEGG_ko", "KEGG_Pathway"]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"emapper 注释必须包含列：{', '.join(required)}（现有列：{header}）")

    if header_fixed_notes:
        with open(out_dir / "_diag" / "header_fix.note", "w", encoding="utf-8") as nf:
            nf.write("\n".join(header_fixed_notes) + "\n")
        for n in header_fixed_notes:
            logging.warning(f"[Fix] {n}")

    idx = {c: header.index(c) for c in required}

    go_fp = out_dir / "gene2go.tsv"
    ko_fp = out_dir / "gene2ko.tsv"
    pw_fp = out_dir / "gene2pathway.tsv"
    cov_fp = out_dir / "universe_coverage.tsv"

    go_w = csv.writer(open(go_fp, "w", encoding="utf-8", newline=""), dialect=csv.excel_tab); go_w.writerow(["gene_id","go_id"])
    ko_w = csv.writer(open(ko_fp, "w", encoding="utf-8", newline=""), dialect=csv.excel_tab); ko_w.writerow(["gene_id","ko"])
    pw_w = csv.writer(open(pw_fp, "w", encoding="utf-8", newline=""), dialect=csv.excel_tab); pw_w.writerow(["gene_id","pathway_id"])

    total = mapped = 0
    genes_seen = set()
    unmapped_examples: List[str] = []

    for ln in data_lines:
        if not ln or ln.startswith("#"): continue
        fields = ln.split("\t")
        if len(fields) < len(header): continue
        q_raw = fields[idx["query"]].strip()
        if not q_raw: continue
        total += 1

        q = clean_query(q_raw, cfg)  # 默认不改动；只有 enable:true 才生效
        gid = tx2map.get(q)
        if not gid:
            if len(unmapped_examples) < 10:
                unmapped_examples.append(q_raw)
            continue

        mapped += 1
        genes_seen.add(gid)

        gos = fields[idx["GOs"]].strip()
        if gos and gos != "-":
            for go in gos.split(","):
                go = go.strip()
                if go:
                    go_w.writerow([gid, go])

        kos = fields[idx["KEGG_ko"]].strip()
        if kos and kos != "-":
            for ko in kos.split(","):
                k = ko.strip()
                if k:
                    ko_w.writerow([gid, k])

        pws = fields[idx["KEGG_Pathway"]].strip()
        if pws and pws != "-":
            for pw in pws.split(","):
                p = pw.strip()
                if p:
                    pw_w.writerow([gid, p])

    with open(cov_fp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["total_queries","mapped_queries","unique_genes_mapped","unmapped_examples"])
        w.writerow([total, mapped, len(genes_seen), ",".join(unmapped_examples)])

    logging.info("========== 07 完成 ==========")
    logging.info(f"[Stat] emapper 总行={total}；mapped={mapped}；unique genes={len(genes_seen)}")
    logging.info(f"[Out] {go_fp}")
    logging.info(f"[Out] {ko_fp}")
    logging.info(f"[Out] {pw_fp}")
    logging.info(f"[Out] {cov_fp}")
    if unmapped_examples:
        logging.warning(f"[Warn] 未映射示例（≤10）：{unmapped_examples}")
        logging.warning("        若为蛋白/含前缀ID，请在 config.emapper.query_id_clean.* 中配置后再启用 enable:true")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 07 执行失败：{e}", file=sys.stderr); sys.exit(1)