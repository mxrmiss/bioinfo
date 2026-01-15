#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03a_build_tx2gene_map.py —— 生成 tx2gene 映射表（vNext 版·裸ID一致化）

契约要点（与《转录组计划 vNext》一致）：
  1) 仅读取 config.yaml（不接命令行参数）。
  2) 输入：reference.ref_gtf（GFF3 或 GTF 均可）。
  3) 输出目录：dirs.maps
     - results/03_maps/tx2gene.clean.tsv        （表头：transcript_id, gene_id）
     - results/03_maps/tx2gene.blacklist.tsv    （表头：gene_id, n_tx；超阈值的一对多映射）
  4) ID 处理（与 02/07 共用同一套规则）：
       - transcript_id 与 gene_id 均应用 annotations.id_cleanup（支持连环剥离）
       - 额外轻度规范化：strip 两端空白、去全角空格
       - 若清洗后出现 transcript->多个 gene 的冲突：直接报错退出（防止隐性错配）
  5) 屏幕流式日志，关键统计与契约自检。
"""

from __future__ import annotations
import sys
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
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
    """读取 config.yaml，并与 DEFAULTS 合并（浅层递归）。"""
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


# ========================= ID 清理（transcript + gene 共用） =========================

def _strip_spaces(raw: str) -> str:
    """轻度清理：去两端空白、去全角空格。"""
    if raw is None:
        return ""
    return raw.strip().replace("\u3000", "")


def apply_id_cleanup(raw: str, policy: Dict[str, Any]) -> str:
    """
    按 annotations.id_cleanup 对单个 ID 进行清理（支持连环剥离）。

    设计目标：
      - 与 02 的清理策略一致：prefix/suffix 可以连续剥离，不依赖列表顺序
      - 仅做“完全匹配”的前/后缀去除，不做任何模糊替换
    """
    s = raw
    order = policy.get("order") or ["prefix", "suffix"]
    strip_prefix = bool(policy.get("strip_prefix"))
    strip_suffix = bool(policy.get("strip_suffix"))
    prefixes: List[str] = policy.get("prefix") or []
    suffixes: List[str] = policy.get("suffix") or []

    prefixes = [p for p in prefixes if p]
    suffixes = [x for x in suffixes if x]

    def strip_prefixes_once(x: str) -> Tuple[str, bool]:
        changed = False
        for p in prefixes:
            if x.startswith(p):
                x = x[len(p):]
                changed = True
        return x, changed

    def strip_suffixes_once(x: str) -> Tuple[str, bool]:
        changed = False
        for suf in suffixes:
            if x.endswith(suf):
                x = x[:-len(suf)]
                changed = True
        return x, changed

    for step in order:
        if step == "prefix" and strip_prefix and prefixes:
            while True:
                s2, ch = strip_prefixes_once(s)
                s = s2
                if not ch:
                    break
        if step == "suffix" and strip_suffix and suffixes:
            while True:
                s2, ch = strip_suffixes_once(s)
                s = s2
                if not ch:
                    break

    return s


def normalize_id(raw: str, policy: Dict[str, Any]) -> str:
    """统一口径：先 strip 空白，再按 policy 去前/后缀。"""
    s0 = _strip_spaces(raw)
    return apply_id_cleanup(s0, policy)


# =============================== GTF/GFF3 解析 ===============================

def parse_attributes(attr: str) -> Tuple[Dict[str, str], str]:
    """
    解析 GTF / GFF3 属性字段。

    更稳判定逻辑：
      1) 先尝试按 GFF3: key=value;key2=value2 解析；若解析出 >=1 个键值对 => 认为是 GFF3
      2) 否则按 GTF: key "value"; key2 "value2"; 解析

    返回：(attrs_dict, style) 其中 style 为 "gff3" 或 "gtf" 或 "empty"
    """
    attr = attr.strip()
    if not attr:
        return {}, "empty"

    d_gff3: Dict[str, str] = {}
    for part in attr.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            k = k.strip()
            v = v.strip()
            if k:
                d_gff3[k] = v

    if d_gff3:
        return d_gff3, "gff3"

    # GTF 风格
    d_gtf: Dict[str, str] = {}
    for part in attr.split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r'(\S+)\s+"(.+)"', part)
        if m:
            d_gtf[m.group(1)] = m.group(2)

    if d_gtf:
        return d_gtf, "gtf"

    # 两种都没解析出来（极少见）
    return {}, "unknown"


def build_tx2gene_from_annot(annot: Path, log: logging.Logger) -> Dict[str, str]:
    """
    从 GTF/GFF3 中提取 transcript -> gene 映射关系。
    优先使用标准字段：
      - GTF: gene_id, transcript_id
      - GFF3: ID, Parent（在 mRNA/transcript 行）
    """
    tx2gene: Dict[str, str] = {}
    n_lines = 0
    n_tx = 0
    n_gene_level_only = 0
    n_gtf_hits = 0
    n_gff3_hits = 0

    with annot.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr = fields
            attrs, style = parse_attributes(attr)

            # GTF 路径：同时存在 transcript_id & gene_id
            if "transcript_id" in attrs and "gene_id" in attrs:
                tid = attrs.get("transcript_id", "").strip()
                gid = attrs.get("gene_id", "").strip()
                if not tid or not gid:
                    continue
                n_tx += 1
                n_gtf_hits += 1
                tx2gene[tid] = gid
                continue

            # GFF3 路径：mRNA / transcript 行，ID=tx, Parent=gene
            if ftype.lower() in {"mrna", "transcript"} and ("ID" in attrs) and ("Parent" in attrs):
                tid = attrs.get("ID", "").strip()
                gid = attrs.get("Parent", "").strip()

                # 多 Parent：只保留第一个并告警
                if "," in gid:
                    parents = [g.strip() for g in gid.split(",") if g.strip()]
                    if len(parents) > 1:
                        log.warning("mRNA/transcript 行有多个 Parent，仅保留第一个：%s", gid)
                    gid = parents[0] if parents else ""

                if not tid or not gid:
                    continue
                n_tx += 1
                n_gff3_hits += 1
                tx2gene[tid] = gid
                continue

            if ftype.lower() == "gene":
                n_gene_level_only += 1

    log.info("注释解析统计：总行数=%d，transcript 记录=%d，gene 级记录=%d", n_lines, n_tx, n_gene_level_only)
    log.info("命中路径统计：GTF命中=%d，GFF3命中=%d", n_gtf_hits, n_gff3_hits)

    if not tx2gene:
        log.error("未从注释中解析出任何 transcript->gene 映射，请检查 ref_gtf/ref_gff3 文件是否正确。")
        sys.exit(1)

    return tx2gene


# =============================== 主流程 ===============================

def main() -> None:
    cfg = load_config(Path(DEFAULT_CONFIG))

    log_level = getattr(logging, str(cfg["logging"].get("level", "INFO")).upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [03a_tx2gene] %(levelname)s: %(message)s" if cfg["logging"].get("timestamp") else "[03a_tx2gene] %(levelname)s: %(message)s",
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

    # 1) 从注释中解析 transcript -> gene 映射（未清理版本）
    raw_tx2gene = build_tx2gene_from_annot(ref_gtf, log)

    # 2) 应用统一 ID 策略（transcript + gene 都清理成“裸 ID”）
    cleaned_tx2gene: Dict[str, str] = {}
    gene_tx_counts: Dict[str, int] = {}

    # 统计与审计
    n_changed_tid = 0
    n_changed_gid = 0
    n_changed_tid_by_strip = 0
    n_changed_gid_by_strip = 0

    # 用于发现“过度剥离”的潜在风险
    tid_raw_to_new: Dict[str, str] = {}
    gid_raw_to_new: Dict[str, str] = {}

    # 冲突检测：同一 cleaned_tid 映射到多个 cleaned_gid -> 直接报错退出
    tid_to_gid_set: Dict[str, set] = {}

    # （仅告警）不同 raw_gid 清洗后合并到同一 gid 的情况
    gid_new_to_raw: Dict[str, set] = {}

    for tid_raw, gid_raw in raw_tx2gene.items():
        tid_raw0 = tid_raw
        gid_raw0 = gid_raw

        # 先 strip 空白，再按 policy 清理
        tid_strip = _strip_spaces(tid_raw0)
        gid_strip = _strip_spaces(gid_raw0)
        if tid_strip != tid_raw0:
            n_changed_tid_by_strip += 1
        if gid_strip != gid_raw0:
            n_changed_gid_by_strip += 1

        tid_new = apply_id_cleanup(tid_strip, id_policy)
        gid_new = apply_id_cleanup(gid_strip, id_policy)

        if tid_new != tid_raw0:
            n_changed_tid += 1
        if gid_new != gid_raw0:
            n_changed_gid += 1

        if not tid_new or not gid_new:
            continue

        tid_raw_to_new[tid_raw0] = tid_new
        gid_raw_to_new[gid_raw0] = gid_new

        tid_to_gid_set.setdefault(tid_new, set()).add(gid_new)
        gid_new_to_raw.setdefault(gid_new, set()).add(gid_raw0)

    # 冲突硬检查：cleaned_tid -> 多个 cleaned_gid
    n_conflict = 0
    for tid_new, gids in tid_to_gid_set.items():
        if len(gids) > 1:
            n_conflict += 1
            gids_show = ", ".join(sorted(list(gids))[:5])
            log.error("ID 清洗后发生 transcript->多个 gene 冲突：transcript_id=%s ; gene_id候选=%s", tid_new, gids_show)

    if n_conflict > 0:
        log.error("发现 %d 个 transcript->gene 冲突。请检查 annotations.id_cleanup.prefix 是否过度剥离。", n_conflict)
        sys.exit(1)

    # 通过检查后，写入最终映射与计数
    for tid_new, gids in tid_to_gid_set.items():
        gid_new = next(iter(gids))
        cleaned_tx2gene[tid_new] = gid_new
        gene_tx_counts[gid_new] = gene_tx_counts.get(gid_new, 0) + 1

    # gene 合并告警（不退出，只提醒）
    n_gid_merged = 0
    for gid_new, raws in gid_new_to_raw.items():
        if len(raws) > 1:
            n_gid_merged += 1
    if n_gid_merged > 0:
        # 只展示少量样例，避免刷屏
        shown = 0
        for gid_new, raws in sorted(gid_new_to_raw.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            if len(raws) <= 1:
                continue
            log.warning("gene_id 清洗后发生合并：gene_id=%s 由 %d 个原始ID合并（示例：%s）",
                        gid_new, len(raws), ", ".join(sorted(list(raws))[:3]))
            shown += 1
            if shown >= 10:
                break
        log.warning("gene_id 清洗合并总计：%d（仅展示前10条示例）", n_gid_merged)

    log.info("应用 ID 策略后：transcript 记录数=%d，gene 数量=%d", len(cleaned_tx2gene), len(gene_tx_counts))
    log.info("被 strip 轻度规范化影响的 transcript_id 数量=%d", n_changed_tid_by_strip)
    log.info("被 strip 轻度规范化影响的 gene_id 数量=%d", n_changed_gid_by_strip)
    log.info("被 ID 清理策略修改的 transcript_id 数量=%d", n_changed_tid)
    log.info("被 ID 清理策略修改的 gene_id 数量=%d", n_changed_gid)

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

