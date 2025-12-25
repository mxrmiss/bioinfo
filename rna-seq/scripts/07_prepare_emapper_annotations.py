#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py —— emapper 注释整理（vNext · 修复版）

要点：
  * 仅读 config.yaml；只认 reference.emapper（不找别名）
  * 用 tx2gene.clean.tsv 的 transcript_id ↔ emapper 的 query 映射
  * transcript_id / query 的 ID 清理规则统一走 annotations.id_cleanup
  * 自动跳过 emapper 文件开头的注释行，直到真正的表头行（#query ... GOs KEGG_ko KEGG_Pathway ...）
  * 输出：
      - results/07_annot/gene2go.tsv
      - results/07_annot/gene2ko.tsv
      - results/07_annot/gene2pathway.tsv
      - results/07_annot/universe_coverage.tsv
  * 表头自愈：
      - '#query' -> 'query'
      - 'KEGG_PathwayKEGG_Module' -> 'KEGG_Pathway'
"""

from __future__ import annotations
import sys
import csv
import gzip
import logging
from pathlib import Path
from typing import Dict, Any, List, Set
import itertools

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


# ================= 配置读取 =================

def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml，并与 DEFAULTS 递归合并。"""
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
    """按 annotations.id_cleanup 对 ID 进行前/后缀修剪（用于 transcript/query）。"""
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
    """支持 .gz 或普通文本两种 emapper 文件。"""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def normalize_header(header: List[str]) -> List[str]:
    """
    自愈 emapper header：
      - '#query' / 'query' 统一为 'query'
      - 'KEGG_PathwayKEGG_Module' -> 'KEGG_Pathway'
    """
    out = []
    for h in header:
        h = h.strip()
        if h.startswith("#"):
            h = h.lstrip("#")
        if h == "KEGG_PathwayKEGG_Module":
            h = "KEGG_Pathway"
        out.append(h)
    return out


def read_tx2gene(tx2gene_fp: Path,
                 id_policy: Dict[str, Any],
                 log: logging.Logger) -> Dict[str, str]:
    """
    读取 tx2gene.clean.tsv，并根据 id_policy 得到：
      transcript_id_clean -> gene_id 映射。
    说明：
      03 脚本已经用相同策略生成过一次，这里再清一遍只是保险，防止版本不一致。
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
            tid_clean = apply_id_cleanup(tid, id_policy)
            mapping[tid_clean] = gid

    log.info("tx2gene 映射加载完成：transcript_id_clean 数量=%d", len(mapping))
    return mapping


def find_emapper_header_and_reader(emapper_fp: Path,
                                   log: logging.Logger) -> (List[str], csv.reader):
    """
    从 emapper 文件中找到真正的 header 行（跳过前面的注释），
    并返回 (header_list, csv_reader)。
    header 判断规则：
      * 行中同时出现 'query' 和 'GOs'（大小写严格），且以 '#' 开头 或 从 'query' 起头。
    """
    fh = open_maybe_gzip(emapper_fp)

    header_line = None
    for line in fh:
        if not line.strip():
            continue
        # 典型 header 形态： "#query\tseed_ortholog\t...\tGOs\t...\tKEGG_ko\tKEGG_Pathway"
        stripped = line.lstrip()
        if (stripped.startswith("#query") or stripped.startswith("query")) and ("GOs" in line):
            header_line = line
            break

    if header_line is None:
        log.error("未在 emapper 文件中找到 header 行（包含 '#query' 和 'GOs'），请检查文件格式。")
        fh.close()
        sys.exit(1)

    # 将 header_line 和后续行组合成新的迭代器交给 csv.reader
    chained_iter = itertools.chain([header_line], fh)
    reader = csv.reader(chained_iter, delimiter="\t")
    header = next(reader, None)
    if header is None:
        log.error("emapper 文件 header 解析失败。")
        fh.close()
        sys.exit(1)

    return header, reader


# ================= 主流程 =================

def main() -> None:
    cfg = load_config(Path(DEFAULT_CONFIG))

    # 日志初始化
    log_level = getattr(logging,
                        str(cfg["logging"].get("level", "INFO")).upper(),
                        logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [07_emapper] %(levelname)s: %(message)s"
        if cfg["logging"].get("timestamp") else "[07_emapper] %(levelname)s: %(message)s",
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

    # 找到真正表头行并构造 csv.reader
    raw_header, reader = find_emapper_header_and_reader(emapper_fp, log)
    header = normalize_header(raw_header)
    col_index = {name: i for i, name in enumerate(header)}

    required_for_mapping = ["query"]
    for col in required_for_mapping:
        if col not in col_index:
            log.error("emapper header 缺少必需列：%s；实际表头为：%s", col, header)
            sys.exit(1)

    # 注释相关列可能缺失：给 WARNING，但不致命
    if "GOs" not in col_index:
        log.warning("emapper 表头缺失列：GOs（将无法生成 gene2go）")
    if "KEGG_ko" not in col_index:
        log.warning("emapper 表头缺失列：KEGG_ko（将无法生成 gene2ko）")
    if "KEGG_Pathway" not in col_index:
        log.warning("emapper 表头缺失列：KEGG_Pathway（将无法生成 gene2pathway）")

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

        # 对 emapper 的 query 应用统一 ID 清理规则，然后与 tx2gene 对齐
        query_clean = apply_id_cleanup(raw_query, id_policy)
        gid = tx2gene.get(query_clean)
        if gid is None:
            n_query_unmapped += 1
            continue
        n_query_mapped += 1

        # GO
        if idx_go is not None and idx_go < len(row):
            gos_raw = row[idx_go].strip()
            if gos_raw and gos_raw != "-":
                # emapper 里常见分隔符有 ',' 或 ';'
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

    log.info(
        "emapper 记录统计：总行数（含数据行）=%d，成功映射到 gene 的 query 数量=%d，未映射=%d",
        n_lines, n_query_mapped, n_query_unmapped
    )

    # ========= 写出 gene2go.tsv =========
    with gene2go_fp.open("w", encoding="utf-8") as out:
        out.write("gene_id\tgo_id\n")
        for gid, gos in sorted(gene2go.items()):
            for go in sorted(gos):
                out.write(f"{gid}\t{go}\n")

    # ========= 写出 gene2ko.tsv =========
    with gene2ko_fp.open("w", encoding="utf-8") as out:
        out.write("gene_id\tko_id\n")
        for gid, kos in sorted(gene2ko.items()):
            for ko in sorted(kos):
                out.write(f"{gid}\t{ko}\n")

    # ========= 写出 gene2pathway.tsv =========
    with gene2pathway_fp.open("w", encoding="utf-8") as out:
        out.write("gene_id\tpathway_id\n")
        for gid, pws in sorted(gene2pathway.items()):
            for pw in sorted(pws):
                out.write(f"{gid}\t{pw}\n")

    # ========= 覆盖率报告 =========
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

    log.info("gene2go/gene2ko/gene2pathway/universe_coverage 已写出至：%s", annot_dir)


if __name__ == "__main__":
    main()