#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py —— emapper 注释整理（严格契约 + 强鲁棒）
契约要点（与“约定文件”一致）：
  * 只读 config.yaml（项目根），不收命令行参数；
  * 关键输入：
      - annotations.emapper : ref/annotations/*.tsv   （emapper 的 TSV 导出）
      - dirs.maps           : results/03_maps         （内含 tx2gene.clean.tsv）
  * 关键输出（长表，英文列名）全部写入 dirs.annot（默认 results/07_annot）：
      - gene2go.tsv        ：gene_id, go_id
      - gene2ko.tsv        ：gene_id, ko
      - gene2pathway.tsv   ：gene_id, pathway_id
      - universe_coverage.tsv（覆盖率摘要）
      - _diag/header_fix.note（如触发表头自愈）
  * 严格字段要求：emapper 中必须可解析出列：query, GOs, KEGG_ko, KEGG_Pathway
    - 允许表头前导“#”自动去除（#query → query）
    - 允许修补黏连 “KEGG_PathwayKEGG_Module”（自动拆分为两列）
  * 映射策略：
      - 用 tx2gene.clean.tsv 的第一列 transcript_id 与 emapper 的 query 对齐；
      - 默认不改动 query；若 config.emapper.query_id_clean.enable=true，
        则按 remove_prefixes/remove_suffixes/regex_suffixes 依次清理后再匹配。
"""

from __future__ import annotations
import sys, os, io, csv, gzip, logging, re
from pathlib import Path
from typing import Dict, Any, List, Tuple

CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "dirs": {
        "maps": "results/03_maps",
        "annot": "results/07_annot",
        "logs": "logs",
    },
    "annotations": {
        "emapper": "ref/annotations/annotations.tsv"
    },
    "emapper": {
        "query_id_clean": {
            "enable": False,              # 默认关闭：不改动 query
            "remove_prefixes": [],        # 例：["Sinonovacula_constricta|", "sp|", "tr|"]
            "remove_suffixes": [],        # 例：[".t1", ".mRNA1"]
            "regex_suffixes": []          # 例：[r"\.[0-9]+$", r"_t\d+$"]
        }
    },
    "logging": {"level": "INFO", "timestamp": True}
}

# ------------------------ 小工具 ------------------------

def load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML，请先安装：mamba/conda install pyyaml", file=sys.stderr)
        raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}
    # 深合并：用户优先
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

def open_auto(fp: Path):
    return gzip.open(fp, "rt") if str(fp).endswith(".gz") else open(fp, "r", encoding="utf-8")

def clean_query(q: str, cfg: Dict[str, Any]) -> str:
    """按配置清理 query（默认不改动）。"""
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
    """修补一行表头文本并返回新的表头行 + 修补记录。"""
    notes: List[str] = []
    s = line.rstrip("\n")
    # 若以 '#query' 开头，去掉第一个 '#'
    if s.startswith("#query\t"):
        s = s[1:]
        notes.append("header: #query -> query")
    # 黏连修补：KEGG_PathwayKEGG_Module → KEGG_Pathway\tKEGG_Module
    if "KEGG_PathwayKEGG_Module" in s:
        s = s.replace("KEGG_PathwayKEGG_Module", "KEGG_Pathway\tKEGG_Module")
        notes.append("header: KEGG_PathwayKEGG_Module -> KEGG_Pathway + KEGG_Module")
    return s + "\n", notes

# ------------------------ 主流程 ------------------------

def main() -> None:
    # 读取契约配置
    cfg = load_yaml(Path(CONFIG_PATH))
    # 日志
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True)
        else "[%(levelname)s] %(message)s",
    )
    # 显式声明：臣妾已读取约定文件（config.yaml）
    logging.info("========== 07 — emapper 注释整理（严格契约 + 强正则） ==========")
    emapper_fp = Path(cfg["annotations"]["emapper"])
    maps_dir   = Path(cfg["dirs"]["maps"])
    out_dir    = Path(cfg["dirs"]["annot"])
    mkdir_p(out_dir)
    diag_dir = out_dir / "_diag"
    mkdir_p(diag_dir)

    logging.info(f"[Info] emapper  = {emapper_fp.resolve()}")
    logging.info(f"[Info] tx2gene  = {(maps_dir/'tx2gene.clean.tsv').resolve()}")
    logging.info(f"[Info] 输出目录  = {out_dir.resolve()}")

    # 读 tx2gene
    tx2gene_map = read_tx2gene(maps_dir / "tx2gene.clean.tsv")
    logging.info(f"[Info] tx2gene 映射条目数 = {len(tx2gene_map):d}")

    # 逐行读 emapper，定位表头并修补
    if not emapper_fp.exists():
        raise FileNotFoundError(f"未找到 emapper 注释文件：{emapper_fp}")

    header_fixed_notes: List[str] = []
    header: List[str] = []
    data_lines: List[str] = []

    with open_auto(emapper_fp) as f:
        for line in f:
            if line.startswith("##"):
                continue  # 注释
            if not header:
                fixed, notes = fix_header_line(line)
                header_fixed_notes.extend(notes)
                header = fixed.rstrip("\n").split("\t")
                continue
            data_lines.append(line.rstrip("\n"))

    # 若检测到列数与数据不一致，再尝试一次针对性修补
    if data_lines:
        cols = len(header)
        first_cols = len(data_lines[0].split("\t"))
        if "KEGG_PathwayKEGG_Module" in header and first_cols == cols + 1:
            # 罕见情况：header token 已经被 split，但仍是黏连字符串对象（防御）
            hi = header.index("KEGG_PathwayKEGG_Module")
            header = header[:hi] + ["KEGG_Pathway", "KEGG_Module"] + header[hi+1:]
            header_fixed_notes.append("header(tokens): KEGG_PathwayKEGG_Module -> KEGG_Pathway + KEGG_Module")

    # 确认必需列
    required = ["query", "GOs", "KEGG_ko", "KEGG_Pathway"]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"emapper 注释必须包含列：{', '.join(required)}（契约）")

    # 写入 header 修补说明（如有）
    if header_fixed_notes:
        with open(diag_dir / "header_fix.note", "w", encoding="utf-8") as nf:
            nf.write("\n".join(header_fixed_notes) + "\n")
        for n in header_fixed_notes:
            logging.warning(f"[Fix] {n}")

    # 建立列索引
    idx = {c: header.index(c) for c in required}

    # 产出文件准备
    go_out  = open(out_dir / "gene2go.tsv", "w", encoding="utf-8", newline="")
    ko_out  = open(out_dir / "gene2ko.tsv", "w", encoding="utf-8", newline="")
    pw_out  = open(out_dir / "gene2pathway.tsv", "w", encoding="utf-8", newline="")
    go_w  = csv.writer(go_out, dialect=csv.excel_tab);    go_w.writerow(["gene_id", "go_id"])
    ko_w  = csv.writer(ko_out, dialect=csv.excel_tab);    ko_w.writerow(["gene_id", "ko"])
    pw_w  = csv.writer(pw_out, dialect=csv.excel_tab);    pw_w.writerow(["gene_id", "pathway_id"])

    # 统计
    total = 0
    mapped = 0
    genes_seen = set()
    unmapped_examples: List[str] = []

    # 主循环
    for ln in data_lines:
        if not ln or ln.startswith("#"):
            continue
        fields = ln.split("\t")
        # 有些行可能列数不齐，防御跳过
        if len(fields) < len(header):
            continue
        q_raw = fields[idx["query"]].strip()
        if not q_raw:
            continue
        total += 1

        q = clean_query(q_raw, cfg)  # 按配置清理（默认不动）
        gene_id = tx2gene_map.get(q)
        if not gene_id:
            # 未映射，留少量样例便于排查
            if len(unmapped_examples) < 10:
                unmapped_examples.append(q_raw)
            continue

        mapped += 1
        genes_seen.add(gene_id)

        # 展开 GO（以逗号分隔；emapper 常见为 "GO:xxxx,GO:yyyy"）
        gos = fields[idx["GOs"]].strip()
        if gos and gos != "-":
            for go in gos.split(","):
                go = go.strip()
                if go:
                    go_w.writerow([gene_id, go])

        # 展开 KEGG_ko（形如 "ko:K12345,ko:K00001"）
        kos = fields[idx["KEGG_ko"]].strip()
        if kos and kos != "-":
            for ko in kos.split(","):
                k = ko.strip()
                if k:
                    ko_w.writerow([gene_id, k])

        # 展开 KEGG_Pathway（形如 "ko00010,map01100"）
        pws = fields[idx["KEGG_Pathway"]].strip()
        if pws and pws != "-":
            for pw in pws.split(","):
                p = pw.strip()
                if p:
                    pw_w.writerow([gene_id, p])

    go_out.close(); ko_out.close(); pw_out.close()

    # 覆盖率摘要
    cov_fp = out_dir / "universe_coverage.tsv"
    with open(cov_fp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["total_queries", "mapped_queries", "unique_genes_mapped", "unmapped_examples"])
        w.writerow([total, mapped, len(genes_seen), ",".join(unmapped_examples)])

    logging.info("========== 07 完成 ==========")
    logging.info(f"[Stat] emapper 总行={total}；mapped={mapped}；unique genes={len(genes_seen)}")
    logging.info(f"[Out] {out_dir/'gene2go.tsv'}")
    logging.info(f"[Out] {out_dir/'gene2ko.tsv'}")
    logging.info(f"[Out] {out_dir/'gene2pathway.tsv'}")
    logging.info(f"[Out] {cov_fp}")
    if unmapped_examples:
        logging.warning(f"[Warn] 有未映射 query（示例≤10）：{unmapped_examples}")
        logging.warning("        若是蛋白/含前缀的 ID，请在 config.emapper.query_id_clean 中开启并配置前/后缀。")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 07 执行失败：{e}", file=sys.stderr)
        sys.exit(1)