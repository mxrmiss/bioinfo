#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_build_tx2gene_map.py —— 生成 tx2gene 映射表（vNext 版）

契约要点（与《转录组计划 vNext》一致）：
  1) 仅读取 config.yaml（不接命令行参数）。
  2) 输入：reference.ref_gtf（GFF3 或 GTF 均可）。
  3) 输出目录：dirs.maps
     - results/03_maps/tx2gene.clean.tsv        （表头：transcript_id, gene_id）
     - results/03_maps/tx2gene.blacklist.tsv    （表头：gene_id, n_tx；超阈值的一对多映射）
  4) ID 处理：
       - gene_id 作为 gene 层唯一主键，只做 strip 空白等轻度规范化，不做前/后缀修剪；
       - transcript_id 如需修剪，只能通过 annotations.id_cleanup 统一规则；
       - 03 与 02/04/07 共用同一套 ID 清理策略。
  5) 屏幕流式日志，关键统计与契约自检。
"""

from __future__ import annotations
import sys
import logging
from pathlib import Path
from typing import Dict, Any, Tuple, List
import re

DEFAULT_CONFIG = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "reference": {
        "ref_gtf": "",
    },
    "dirs": {
        "maps": "results/03_maps",
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
    "tx2gene": {
        "blacklist_tx_threshold": 10,  # 一对多映射超过该阈值的 gene 将进入 blacklist
    },
}


def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml，并与 DEFAULTS 合并。"""
    try:
        import yaml
    except Exception as e:
        print("[ERR] 需要 PyYAML 支持：pip install pyyaml", file=sys.stderr)
        raise e

    if not path.exists():
        print(f"[ERR] 未找到配置文件：{path}", file=sys.stderr)
        sys.exit(1)

    with open(path, "r", encoding="utf-8") as f:
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


# =============== ID 清理（仅用于 transcript_id） ====================

def apply_id_cleanup(raw: str, policy: Dict[str, Any]) -> str:
    """按 annotations.id_cleanup 对 transcript_id 进行前/后缀修剪。"""
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


def clean_gene_id(raw: str) -> str:
    """
    gene_id 的“轻度清理”：只去除首尾空白和明显奇怪的空格，不做任何版本号/前后缀剥离。
    符合“gene_id 永不动刀”的铁律。
    """
    if raw is None:
        return ""
    # 去掉两端空白和全角空格
    s = raw.strip().replace("\u3000", "")
    return s


# =============== GTF/GFF3 解析 ====================

def parse_attributes(attr: str) -> Dict[str, str]:
    """
    解析 GTF / GFF3 属性字段，兼容：
      - GTF: key " \"value\";" 风格
      - GFF3: key=value;key2=value2 风格
    """
    attr = attr.strip()
    if not attr:
        return {}
    d: Dict[str, str] = {}

    # 粗略判断 GTF vs GFF3
    if ";" in attr and "\"" in attr:
        # GTF 风格
        for part in attr.split(";"):
            part = part.strip()
            if not part:
                continue
            m = re.match(r'(\S+)\s+"(.+)"', part)
            if m:
                d[m.group(1)] = m.group(2)
    else:
        # GFF3 风格
        for part in attr.split(";"):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                k, v = part.split("=", 1)
                d[k] = v
    return d


def build_tx2gene_from_gtf(gtf: Path, log: logging.Logger) -> Dict[str, str]:
    """
    从 GTF/GFF3 中提取 transcript -> gene 映射关系。
    优先使用标准字段：
      - GTF: gene_id, transcript_id
      - GFF3: ID, Parent, 结合 feature type
    """
    tx2gene: Dict[str, str] = {}
    n_lines = 0
    n_tx = 0
    n_gene_level_only = 0

    with gtf.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr = fields
            attrs = parse_attributes(attr)

            # GTF 优先路径
            if "transcript_id" in attrs:
                tid = attrs.get("transcript_id", "").strip()
                gid = attrs.get("gene_id", "").strip()
                if not tid or not gid:
                    continue
                n_tx += 1
                tx2gene[tid] = gid
                continue

            # GFF3 风格：mRNA/ transcript 作为转录本，Parent 指向 gene
            if ftype.lower() in {"mrna", "transcript"} and "Parent" in attrs:
                tid = attrs.get("ID", "").strip()
                gid = attrs.get("Parent", "").strip()
                if "," in gid:
                    # 多 Parent 的情况：一般不太常见，这里只取第一个并报警
                    parents = [g.strip() for g in gid.split(",") if g.strip()]
                    if len(parents) > 1:
                        log.warning("mRNA/transcript 行有多个 Parent，仅保留第一个：%s", gid)
                    gid = parents[0] if parents else ""
                if not tid or not gid:
                    continue
                n_tx += 1
                tx2gene[tid] = gid
                continue

            # gene 行只统计，不直接输出
            if ftype.lower() == "gene":
                n_gene_level_only += 1

    log.info("GTF/GFF3 解析统计：总行数=%d，transcript 记录=%d，gene 级记录=%d", n_lines, n_tx, n_gene_level_only)
    if not tx2gene:
        log.error("未从 GTF/GFF3 中解析出任何 transcript->gene 映射，请检查文件格式。")
        sys.exit(1)
    return tx2gene


# =============== 主流程 ====================

def main() -> None:
    cfg = load_config(Path(DEFAULT_CONFIG))

    log_level = getattr(logging, str(cfg["logging"].get("level", "INFO")).upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [03_tx2gene] %(levelname)s: %(message)s" if cfg["logging"].get("timestamp") else "[03_tx2gene] %(levelname)s: %(message)s",
    )
    log = logging.getLogger("tx2gene")

    ref_gtf = Path(cfg["reference"]["ref_gtf"])
    maps_dir = Path(cfg["dirs"]["maps"])
    maps_dir.mkdir(parents=True, exist_ok=True)

    tx2gene_clean = maps_dir / "tx2gene.clean.tsv"
    tx2gene_blacklist = maps_dir / "tx2gene.blacklist.tsv"

    id_policy = cfg.get("annotations", {}).get("id_cleanup", {}) or DEFAULTS["annotations"]["id_cleanup"]
    blacklist_threshold = int(cfg.get("tx2gene", {}).get("blacklist_tx_threshold", 10))

    log.info("ref_gtf: %s", ref_gtf)
    log.info("输出目录：%s", maps_dir)
    log.info("blacklist_tx_threshold: %d", blacklist_threshold)

    if not ref_gtf.exists():
        log.error("参考注释文件不存在：%s", ref_gtf)
        sys.exit(1)

    # 1) 从 GTF/GFF3 中解析 transcript -> gene 映射
    raw_tx2gene = build_tx2gene_from_gtf(ref_gtf, log)

    # 2) 应用 ID 策略：
    #    - transcript_id 可按 policy 修剪
    #    - gene_id 只做 strip，不做前后缀修剪
    cleaned_tx2gene: Dict[str, str] = {}
    gene_tx_counts: Dict[str, int] = {}
    n_changed_tid = 0
    n_changed_gid = 0  # 理论上只是 strip，统计一下

    for tid_raw, gid_raw in raw_tx2gene.items():
        gid = clean_gene_id(gid_raw)
        if gid != gid_raw:
            n_changed_gid += 1
        tid_new = apply_id_cleanup(tid_raw, id_policy)
        if tid_new != tid_raw:
            n_changed_tid += 1

        if not tid_new or not gid:
            continue

        # 若同一个 transcript 被映射到多个 gene，以最后一次为准并给出警告
        if tid_new in cleaned_tx2gene and cleaned_tx2gene[tid_new] != gid:
            log.warning("transcript_id=%s 在 GTF 中被映射到多个 gene（%s, %s），以最后一次为准。",
                        tid_new, cleaned_tx2gene[tid_new], gid)
        cleaned_tx2gene[tid_new] = gid
        gene_tx_counts[gid] = gene_tx_counts.get(gid, 0) + 1

    log.info("应用 ID 策略后：transcript 记录数=%d，gene 数量=%d", len(cleaned_tx2gene), len(gene_tx_counts))
    log.info("被轻度规范化的 gene_id 数量（strip 空白等）=%d", n_changed_gid)
    log.info("被 ID 清理策略修改的 transcript_id 数量=%d", n_changed_tid)

    # 3) 输出 tx2gene.clean.tsv
    with tx2gene_clean.open("w", encoding="utf-8") as out:
        out.write("transcript_id\tgene_id\n")
        for tid, gid in sorted(cleaned_tx2gene.items()):
            out.write(f"{tid}\t{gid}\n")

    # 4) 输出 blacklist（超一对多阈值的 gene）
    with tx2gene_blacklist.open("w", encoding="utf-8") as out:
        out.write("gene_id\tn_tx\n")
        n_black = 0
        for gid, cnt in sorted(gene_tx_counts.items(), key=lambda kv: (-kv[1], kv[0])):
            if cnt >= blacklist_threshold:
                out.write(f"{gid}\t{cnt}\n")
                n_black += 1

    log.info("tx2gene.clean.tsv 写入：%s", tx2gene_clean)
    log.info("tx2gene.blacklist.tsv 写入：%s（n_blacklisted=%d）", tx2gene_blacklist, n_black)


if __name__ == "__main__":
    main()