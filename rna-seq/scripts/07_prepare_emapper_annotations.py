#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py —— emapper 注释整理（vNext 版）

要点：
  * 只读 config.yaml；只认 reference.emapper（不找别名）
  * 用 tx2gene.clean.tsv 的 transcript_id ↔ emapper 的 query 映射
  * transcript_id / query 的 ID 清理规则统一走 annotations.id_cleanup
  * 输出：gene2go.tsv / gene2ko.tsv / gene2pathway.tsv / universe_coverage.tsv
  * 表头自愈（两处）：#query→query；KEGG_PathwayKEGG_Module → KEGG_Pathway + KEGG_Module
"""

from __future__ import annotations
import sys
import csv
import gzip
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable, Set

DEFAULT_CONFIG = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "reference": {
        "emapper": "ref/annotations/emapper.tsv",
    },
    "dirs": {
        "annot": "results/07_annot",
        "maps":  "results/03_maps",
    },
    "annotations": {
        "id_cleanup": {
            "strip_prefix": False,
            "prefix": [],
            "strip_suffix": False,
            "suffix": [],
            "order": ["prefix", "suffix"],
        }
    },
    "logging": {
        "level": "INFO",
        "timestamp": True,
    },
}


def load_config(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception as e:
        print("[ERR] 需要 PyYAML 支持：pip install pyyaml", file=sys.stderr)
        raise e

    if not path.exists():
        print(f"[ERR] 未找到配置文件：{path}", file=sys.stderr)
        sys.exit(1)

    with path.open("r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}

    def merge(base: Dict[str, Any], u: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(base)
        for k, v in u.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(out[k], v)
            else:
                out[k] = v
        return out

    return merge(DEFAULTS, user)


# ================= ID 清理（与 02/03 共用逻辑） =================

def apply_id_cleanup(raw: str, policy: Dict[str, Any]) -> str:
    s = raw
    order = policy.get("order") or ["prefix", "suffix"]
    strip_prefix = bool(policy.get("strip_prefix"))
    strip_suffix = bool(policy.get("strip_suffix"))
    prefixes: List[str] = policy.get("prefix") or []
    suffixes: List[str] = policy.get("suffix") or []

    for step in order:
        if step == "prefix" and strip_prefix:
            for p in prefixes:
                if p and s.startswith(p):
                    s = s[len(p):]
        if step == "suffix" and strip_suffix:
            for suf in suffixes:
                if suf and s.endswith(suf):
                    s = s[:-len(suf)]
    return s


# ================= 工具函数 =================

def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def normalize_header(header: List[str]) -> List[str]:
    """
    自愈 emapper header：
      - '#query' -> 'query'
      - KEGG_PathwayKEGG_Module -> KEGG_Pathway + KEGG_Module 拆列（这里仅改名，具体字段可按需要扩展）
    """
    out = []
    for h in header:
        if h.startswith("#"):
            h = h.lstrip("#")
        if h == "KEGG_PathwayKEGG_Module":
            # 兼容性：现阶段只用 KEGG_Pathway，对 Module 暂不展开
            h = "KEGG_Pathway"
        out.append(h)
    return out


def read_tx2gene(tx2gene_fp: Path, id_policy: Dict[str, Any], log: logging.Logger) -> Dict[str, str]:
    """
    读取 tx2gene.clean.tsv，并根据 id_policy 推导出：
      - transcript_id_clean -> gene_id 映射
    注意：03 已经应用过相同策略，这里只是假定“clean 后的 transcript_id 就是表里的样子”，
    为了稳妥起见，我们依然用 policy 再清一次，以防脚本版本不一致。
    """
    if not tx2gene_fp.exists():
        log.error("未找到 tx2gene.clean.tsv：%s", tx2gene_fp)
        sys.exit(1)

    mapping: Dict[str, str] = {}
    with tx2gene_fp.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            idx_tid = header.index("transcript_id")
            idx_gid = header.index("gene_id")
        except ValueError:
            log.error("tx2gene.clean.tsv 表头必须包含 transcript_id,gene_id")
            sys.exit(1)

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_tid, idx_gid):
                continue
            tid = parts[idx_tid]
            gid = parts[idx_gid]
            # 用相同策略再清一次，防止历史版本不一致
            tid_clean = apply_id_cleanup(tid, id_policy)
            mapping[tid_clean] = gid

    log.info("tx2gene 映射加载完成：transcript_id_clean 数量=%d", len(mapping))
    return mapping


# ================= 主流程 =================

def main() -> None:
    cfg = load_config(Path(DEFAULT_CONFIG))

    log_level = getattr(logging, str(cfg["logging"].get("level", "INFO")).upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [07_emapper] %(levelname)s: %(message)s" if cfg["logging"].get("timestamp") else "[07_emapper] %(levelname)s: %(message)s",
    )
    log = logging.getLogger("emapper")

    emapper_fp = Path(cfg["reference"]["emapper"])
    annot_dir = Path(cfg["dirs"]["annot"])
    maps_dir = Path(cfg["dirs"]["maps"])
    annot_dir.mkdir(parents=True, exist_ok=True)

    tx2gene_fp = maps_dir / "tx2gene.clean.tsv"
    gene2go_fp = annot_dir / "gene2go.tsv"
    gene2ko_fp = annot_dir / "gene2ko.tsv"
    gene2pathway_fp = annot_dir / "gene2pathway.tsv"
    universe_cov_fp = annot_dir / "universe_coverage.tsv"

    id_policy = cfg.get("annotations", {}).get("id_cleanup", {}) or DEFAULTS["annotations"]["id_cleanup"]

    log.info("emapper 输入：%s", emapper_fp)
    log.info("tx2gene.clean：%s", tx2gene_fp)
    log.info("输出目录：%s", annot_dir)

    if not emapper_fp.exists():
        log.error("emapper 注释文件不存在：%s", emapper_fp)
        sys.exit(1)

    tx2gene = read_tx2gene(tx2gene_fp, id_policy, log)

    # 读取 emapper（gz or plain）
    with open_maybe_gzip(emapper_fp) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            log.error("emapper 文件为空：%s", emapper_fp)
            sys.exit(1)
        header = normalize_header(header)
        col_index = {name: i for i, name in enumerate(header)}

        required_cols = ["query", "GOs", "KEGG_ko", "KEGG_Pathway"]
        for col in required_cols:
            if col not in col_index:
                log.warning("emapper 表头缺失列：%s（若不用对应注释可忽略）", col)

        idx_query = col_index.get("query")
        idx_go = col_index.get("GOs")
        idx_ko = col_index.get("KEGG_ko")
        idx_pathway = col_index.get("KEGG_Pathway")

        gene2go: Dict[str, Set[str]] = {}
        gene2ko: Dict[str, Set[str]] = {}
        gene2pathway: Dict[str, Set[str]] = {}

        n_lines = 0
        n_query_mapped = 0
        n_query_unmapped = 0

        for row in reader:
            n_lines += 1
            if idx_query is None or idx_query >= len(row):
                continue
            raw_query = row[idx_query].strip()
            if not raw_query:
                continue

            # 对 emapper 的 query 应用同一套 ID 清理规则，然后与 tx2gene 对齐
            query_clean = apply_id_cleanup(raw_query, id_policy)
            gid = tx2gene.get(query_clean)
            if gid is None:
                # 若严格模式可以在此处选择报错，这里采取记数 + 覆盖率报告
                n_query_unmapped += 1
                continue
            n_query_mapped += 1

            # GO
            if idx_go is not None and idx_go < len(row):
                gos_raw = row[idx_go].strip()
                if gos_raw and gos_raw != "-":
                    for go in gos_raw.replace(",", ";").split(";"):
                        go = go.strip()
                        if not go:
                            continue
                        gene2go.setdefault(gid, set()).add(go)

            # KO
            if idx_ko is not None and idx_ko < len(row):
                kos_raw = row[idx_ko].strip()
                if kos_raw and kos_raw != "-":
                    for ko in kos_raw.replace(",", ";").split(";"):
                        ko = ko.strip()
                        if not ko:
                            continue
                        gene2ko.setdefault(gid, set()).add(ko)

            # KEGG Pathway
            if idx_pathway is not None and idx_pathway < len(row):
                pw_raw = row[idx_pathway].strip()
                if pw_raw and pw_raw != "-":
                    for p in pw_raw.replace(",", ";").split(";"):
                        p = p.strip()
                        if not p:
                            continue
                        gene2pathway.setdefault(gid, set()).add(p)

    log.info("emapper 总记录数=%d，成功映射到 gene 的 query 数量=%d，未映射=%d",
             n_lines, n_query_mapped, n_query_unmapped)

    # 输出 gene2go.tsv
    with gene2go_fp.open("w", encoding="utf-8") as out:
        out.write("gene_id\tgo_id\n")
        for gid, gos in sorted(gene2go.items()):
            for go in sorted(gos):
                out.write(f"{gid}\t{go}\n")

    # 输出 gene2ko.tsv
    with gene2ko_fp.open("w", encoding="utf-8") as out:
        out.write("gene_id\tko_id\n")
        for gid, kos in sorted(gene2ko.items()):
            for ko in sorted(kos):
                out.write(f"{gid}\t{ko}\n")

    # 输出 gene2pathway.tsv（这里先用 ko->pathway 的字面映射，后续再由 kegg_local.sh 聚合）
    with gene2pathway_fp.open("w", encoding="utf-8") as out:
        out.write("gene_id\tpathway_id\n")
        for gid, pws in sorted(gene2pathway.items()):
            for pw in sorted(pws):
                out.write(f"{gid}\t{pw}\n")

    # 覆盖率报告
    universe_cov_fp.parent.mkdir(parents=True, exist_ok=True)
    with universe_cov_fp.open("w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"n_tx2gene_transcripts\t{len(tx2gene)}\n")
        out.write(f"n_emapper_records\t{n_lines}\n")
        out.write(f"n_query_mapped_to_gene\t{n_query_mapped}\n")
        out.write(f"n_query_unmapped\t{n_query_unmapped}\n")
        out.write(f"n_gene_with_GO\t{len(gene2go)}\n")
        out.write(f"n_gene_with_KO\t{len(gene2ko)}\n")
        out.write(f"n_gene_with_pathway\t{len(gene2pathway)}\n")

    logging.info("gene2go/gene2ko/gene2pathway/universe_coverage 已写出至：%s", annot_dir)


if __name__ == "__main__":
    main()