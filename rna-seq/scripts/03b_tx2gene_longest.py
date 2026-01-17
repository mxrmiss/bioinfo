#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03b_tx2gene_longest.py —— 用 AGAT 从注释中提取“每基因最长转录本”，为 tx2gene.clean.tsv 增补第三列（适配 03a 裸ID）

契约（按皇上最新要求 · 已修复 XM/XR 并存导致的冲突）：
1) 不改 03a_build_tx2gene_map.py，不覆盖 tx2gene.clean.tsv
2) 最长转录本提取：只用 AGAT（agat_sp_keep_longest_isoform.pl），不提供任何兜底策略
   - AGAT 不存在/运行失败/输出为空 -> 直接报错退出
3) 输出文件仅保留一个三列表：
   - tx2gene.longest.tsv
     transcript_id  gene_id  longest_transcript_id
     *前两列逐行完全照抄 tx2gene.clean.tsv（不做任何改动、顺序不变）*
     第三列来自注释中该 gene 的“代表最长转录本”（与 03a 同一套 id_cleanup 宇宙的“裸ID”）
4) 匹配全过程采用统一 ID 清洗策略（与 02/03a/07 共用）：
   - 对 keep_longest 注释中解析到的 Parent(gene) 与 ID(transcript) 都应用 annotations.id_cleanup（支持连环剥离）
   - tx2gene.clean.tsv 前两列原样输出；仅用于 lookup 的 gene_key 会按同一策略归一化
5) 自检策略（关键修复点）：
   - keep_longest 里同一 gene 可能同时出现 mRNA(XM_*) 与 transcript(XR_*) 两类代表，这是注释常态
   - 选择规则：优先选 mRNA；若该 gene 没有 mRNA 才使用 transcript
   - 仅在“同一优先级类别”内出现多个不同候选时才视为冲突并退出：
       * mRNA 候选 >1 -> 报错退出
       * transcript-only 候选 >1 -> 报错退出
   - 缺失 longest 的行允许存在（第三列填空），但会在 report 中统计并给出 gene 示例
6) 生成报告：
   - tx2gene.longest.report.tsv 记录关键统计、AGAT 命令、缺失情况、冲突情况等

输入：
- config.yaml
  - reference.ref_gtf         注释文件路径（GTF 或 GFF3；AGAT 均可接）
  - dirs.maps                 输出目录（默认 results/03_maps）
  - annotations.id_cleanup    统一 ID 清理策略（03a/02/07 共用）
- results/03_maps/tx2gene.clean.tsv（03a 产物）

输出（写入 dirs.maps）：
- tx2gene.longest.tsv
- tx2gene.longest.report.tsv
"""

from __future__ import annotations

import sys
import re
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set


# =============================================================================
# 参数区（按皇上习惯：不走命令行，只允许改这里 / config.yaml）
# =============================================================================

DEFAULT_CONFIG = "config.yaml"

# 找不到 longest 时第三列填什么（默认空）
MISSING_LONGEST_FILL = ""

# 临时目录名（创建在 dirs.maps 下）
TMP_SUBDIR = "_tmp"

# AGAT 程序名（契约：AGAT-only）
AGAT_KEEP_LONGEST = "agat_sp_keep_longest_isoform.pl"


# =============================================================================
# 默认配置（与 03/03a 同语义）
# =============================================================================

DEFAULTS: Dict[str, Any] = {
    "reference": {"ref_gtf": ""},
    "dirs": {"maps": "results/03_maps"},
    "annotations": {
        "id_cleanup": {
            "strip_prefix": False,
            "prefix": [],
            "strip_suffix": False,
            "suffix": [],
            "order": ["prefix", "suffix"],
        }
    },
    "logging": {"level": "INFO", "timestamp": True},
}


# =============================================================================
# 基础工具函数
# =============================================================================

def need(cmd: str) -> None:
    """检查依赖命令是否存在；不存在则直接退出（契约：AGAT-only）"""
    if subprocess.call(f"command -v {cmd} >/dev/null 2>&1", shell=True) != 0:
        print(f"[ERR] dependency not found: {cmd}", file=sys.stderr)
        sys.exit(1)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def run_capture(cmd: str, label: str = "CMD") -> Tuple[int, str, str]:
    """执行 shell 命令并捕获 stdout/stderr；失败直接退出"""
    print(f"[{label}] {cmd}")
    cp = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if cp.returncode != 0:
        print(f"[ERR] command failed ({label}), exit code: {cp.returncode}", file=sys.stderr)
        if cp.stderr:
            print(cp.stderr[:8000], file=sys.stderr)
        sys.exit(cp.returncode)
    return cp.returncode, (cp.stdout or ""), (cp.stderr or "")


def safe_preview(s: str, n: int = 300) -> str:
    """报告里预览用：把换行/制表符压成空格，截断"""
    if s is None:
        return ""
    x = s.replace("\r", " ").replace("\n", " ").replace("\t", " ")
    x = re.sub(r"\s+", " ", x).strip()
    return x[:n]


# =============================================================================
# 配置读取
# =============================================================================

def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml，并与 DEFAULTS 合并"""
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML 支持：pip install pyyaml", file=sys.stderr)
        sys.exit(1)

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


# =============================================================================
# ID 清理：与 02/03a 一致（连环剥离）
# =============================================================================

def _strip_spaces(raw: str) -> str:
    """轻度清理：去两端空白、去全角空格"""
    if raw is None:
        return ""
    return raw.strip().replace("\u3000", "")


def apply_id_cleanup(raw: str, policy: Dict[str, Any]) -> str:
    """
    按 annotations.id_cleanup 对单个 ID 进行清理（支持连环剥离）
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
    """统一口径：先 strip 空白，再按 policy 去前/后缀"""
    return apply_id_cleanup(_strip_spaces(raw), policy)


# =============================================================================
# GTF/GFF3 属性解析（供解析 keep_longest 输出）
# =============================================================================

def parse_attributes(attr: str) -> Tuple[Dict[str, str], str]:
    """
    更稳判定：
      1) 先按 GFF3 key=value 解析；若解析出 >=1 个键值对 -> gff3
      2) 否则按 GTF key "value" -> gtf
    返回：(attrs, style)
    """
    attr = (attr or "").strip()
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

    return {}, "unknown"


# =============================================================================
# 读取 tx2gene.clean.tsv（保持原顺序，前两列原样输出）
# =============================================================================

def read_tx2gene_clean(tx2gene_path: Path) -> Tuple[List[Tuple[str, str]], Set[str]]:
    rows: List[Tuple[str, str]] = []
    genes: Set[str] = set()

    with tx2gene_path.open("r", encoding="utf-8") as f:
        header = f.readline()
        if not header:
            print(f"[ERR] empty file: {tx2gene_path}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            tid_raw = parts[0]
            gid_raw = parts[1]
            rows.append((tid_raw, gid_raw))
            genes.add(gid_raw)

    return rows, genes


# =============================================================================
# AGAT keep_longest + 解析 gene->longest_transcript
# =============================================================================

def run_agat_keep_longest(ref_gtf: Path, out_gff: Path) -> Tuple[int, str, str]:
    """只用 AGAT：agat_sp_keep_longest_isoform.pl -g <ref> -o <out>"""
    need(AGAT_KEEP_LONGEST)
    ensure_dir(out_gff.parent)

    cmd = f"{AGAT_KEEP_LONGEST} -g '{ref_gtf}' -o '{out_gff}'"
    rc, out, err = run_capture(cmd, "AGAT:KEEP_LONGEST")

    if (not out_gff.exists()) or out_gff.stat().st_size == 0:
        print(f"[ERR] AGAT output missing or empty: {out_gff}", file=sys.stderr)
        sys.exit(1)

    return rc, out, err


def parse_keep_longest_gene2tx(
    keep_longest_gff: Path,
    id_cleanup_policy: Dict[str, Any],
    log: logging.Logger,
) -> Tuple[Dict[str, str], int, Dict[str, int]]:
    """
    从 keep_longest 注释解析 gene_id(cleaned) -> chosen_longest_transcript_id(cleaned)

    关键修复：
      - keep_longest 可能同时保留 mRNA 与 transcript 两类代表（常见 XM 与 XR 并存）
      - 选择规则：优先 mRNA；无 mRNA 才用 transcript
      - 冲突判定：仅在同一优先级类别内出现多个不同候选时才报错退出

    仅解析两类行：
      - ftype == mRNA
      - ftype == transcript

    期望属性（GFF3）：ID=tid ; Parent=gid
    Parent 多值：只取第一个，并计数告警

    返回：
      - gene2chosen_tx
      - multi_parent_warn_count
      - stats（用于 report）
    """
    gene2cands_mrna: Dict[str, Set[str]] = {}
    gene2cands_tx: Dict[str, Set[str]] = {}

    multi_parent_warn = 0

    # 统计用
    n_lines_mrna = 0
    n_lines_transcript = 0

    with keep_longest_gff.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            ftype = fields[2].strip().lower()
            if ftype not in {"mrna", "transcript"}:
                continue

            attrs, _style = parse_attributes(fields[8])
            tid_raw = (attrs.get("ID", "") or "").strip()
            gid_raw = (attrs.get("Parent", "") or "").strip()
            if not tid_raw or not gid_raw:
                continue

            if "," in gid_raw:
                parents = [g.strip() for g in gid_raw.split(",") if g.strip()]
                if len(parents) > 1:
                    multi_parent_warn += 1
                gid_raw = parents[0] if parents else ""

            if not gid_raw:
                continue

            gid_clean = normalize_id(gid_raw, id_cleanup_policy)
            tid_clean = normalize_id(tid_raw, id_cleanup_policy)

            if not gid_clean or not tid_clean:
                continue

            if ftype == "mrna":
                n_lines_mrna += 1
                gene2cands_mrna.setdefault(gid_clean, set()).add(tid_clean)
            else:
                n_lines_transcript += 1
                gene2cands_tx.setdefault(gid_clean, set()).add(tid_clean)

    # 汇总 gene 集合
    all_genes: Set[str] = set(gene2cands_mrna.keys()) | set(gene2cands_tx.keys())

    if not all_genes:
        print(f"[ERR] no gene candidates parsed from: {keep_longest_gff}", file=sys.stderr)
        sys.exit(1)

    # 选择 + 冲突检测
    gene2chosen: Dict[str, str] = {}

    conflict_mrna = 0
    conflict_transcript_only = 0

    # 为了避免日志刷屏，只预览前 N 个冲突
    conflict_preview_limit = 30
    conflict_preview_count = 0

    genes_with_mrna = 0
    genes_with_transcript = 0
    genes_with_both = 0
    genes_with_transcript_only = 0

    for g in sorted(all_genes):
        mrna_set = gene2cands_mrna.get(g, set())
        tx_set = gene2cands_tx.get(g, set())

        if mrna_set:
            genes_with_mrna += 1
        if tx_set:
            genes_with_transcript += 1
        if mrna_set and tx_set:
            genes_with_both += 1
        if (not mrna_set) and tx_set:
            genes_with_transcript_only += 1

        if mrna_set:
            if len(mrna_set) > 1:
                conflict_mrna += 1
                if conflict_preview_count < conflict_preview_limit:
                    log.error(
                        "同一 gene 在 mRNA 类别出现多个候选 longest_transcript（清洗后仍不同）：gene_id=%s ; candidates=%s",
                        g, ",".join(sorted(mrna_set))
                    )
                    conflict_preview_count += 1
                continue
            gene2chosen[g] = next(iter(mrna_set))
            continue

        if tx_set:
            if len(tx_set) > 1:
                conflict_transcript_only += 1
                if conflict_preview_count < conflict_preview_limit:
                    log.error(
                        "同一 gene 在 transcript 类别出现多个候选 longest_transcript（清洗后仍不同）：gene_id=%s ; candidates=%s",
                        g, ",".join(sorted(tx_set))
                    )
                    conflict_preview_count += 1
                continue
            gene2chosen[g] = next(iter(tx_set))
            continue

        # 理论上不会到这，因为 all_genes 来自两类 dict 的 key
        continue

    if conflict_mrna > 0 or conflict_transcript_only > 0:
        log.error(
            "keep_longest 解析冲突：mRNA 类冲突 gene=%d；transcript-only 类冲突 gene=%d。请检查注释或 ID 清洗规则。",
            conflict_mrna, conflict_transcript_only
        )
        sys.exit(1)

    if not gene2chosen:
        print(f"[ERR] no chosen gene->transcript after resolving: {keep_longest_gff}", file=sys.stderr)
        sys.exit(1)

    stats: Dict[str, int] = {
        "keep_longest_lines_mrna": n_lines_mrna,
        "keep_longest_lines_transcript": n_lines_transcript,
        "keep_longest_genes_total": len(all_genes),
        "keep_longest_genes_with_mrna": genes_with_mrna,
        "keep_longest_genes_with_transcript": genes_with_transcript,
        "keep_longest_genes_with_both": genes_with_both,
        "keep_longest_genes_with_transcript_only": genes_with_transcript_only,
        "keep_longest_conflict_mrna": conflict_mrna,
        "keep_longest_conflict_transcript_only": conflict_transcript_only,
    }

    return gene2chosen, multi_parent_warn, stats


# =============================================================================
# 主流程
# =============================================================================

def main() -> None:
    cfg = load_config(Path(DEFAULT_CONFIG))

    log_level = getattr(logging, str(cfg["logging"].get("level", "INFO")).upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [03b_longest] %(levelname)s: %(message)s" if cfg["logging"].get("timestamp") else "[03b_longest] %(levelname)s: %(message)s",
    )
    log = logging.getLogger("03b_longest")

    ref_gtf = Path(cfg["reference"]["ref_gtf"])
    maps_dir = Path(cfg["dirs"]["maps"])
    ensure_dir(maps_dir)

    tx2gene_clean = maps_dir / "tx2gene.clean.tsv"
    out_longest = maps_dir / "tx2gene.longest.tsv"
    out_report = maps_dir / "tx2gene.longest.report.tsv"

    tmp_dir = maps_dir / TMP_SUBDIR
    ensure_dir(tmp_dir)
    keep_longest_gff = tmp_dir / "ref.keep_longest.gff3"

    id_policy = cfg.get("annotations", {}).get("id_cleanup", {}) or DEFAULTS["annotations"]["id_cleanup"]

    log.info("ref_gtf: %s", ref_gtf)
    log.info("maps_dir: %s", maps_dir)
    log.info("tx2gene_clean: %s", tx2gene_clean)

    if not ref_gtf.exists():
        log.error("参考注释文件不存在：%s（请检查 config.yaml 的 reference.ref_gtf）", ref_gtf)
        sys.exit(1)

    if not tx2gene_clean.exists():
        log.error("未找到：%s（请先运行 03a_build_tx2gene_map.py）", tx2gene_clean)
        sys.exit(1)

    # 1) 读 tx2gene.clean.tsv（原样保留行框架；前两列将原样输出）
    rows_raw, genes_set_raw = read_tx2gene_clean(tx2gene_clean)
    log.info("读取 tx2gene.clean.tsv：rows=%d, genes=%d", len(rows_raw), len(genes_set_raw))

    # 2) 只用 AGAT 生成 keep_longest 注释
    rc, agat_stdout, agat_stderr = run_agat_keep_longest(ref_gtf, keep_longest_gff)

    # 3) 从 keep_longest 注释解析 gene(cleaned)->chosen_longest_tx(cleaned)
    gene2longest, multi_parent_warn, kl_stats = parse_keep_longest_gene2tx(keep_longest_gff, id_policy, log)
    log.info(
        "解析 keep_longest：chosen_genes=%d, parsed_genes_total=%d, multi_parent_warn=%d",
        len(gene2longest), kl_stats.get("keep_longest_genes_total", 0), multi_parent_warn
    )

    # 4) 写 tx2gene.longest.tsv（前两列原样照抄；第三列为 chosen_longest_transcript_id(cleaned)）
    miss_rows = 0
    miss_gene_examples: List[str] = []
    hit_genes: Set[str] = set()

    # 审计：用于 lookup 的 gene_key 是否与原 gene_id 不同（理论上 03a 已一致，这里仅统计）
    n_lookup_gene_changed = 0

    with out_longest.open("w", encoding="utf-8") as fo:
        fo.write("transcript_id\tgene_id\tlongest_transcript_id\n")

        for tid_raw, gid_raw in rows_raw:
            # 前两列：严格照抄
            tid_out = tid_raw
            gid_out = gid_raw

            # lookup 用：按同一策略归一化 gene
            gid_key = normalize_id(gid_raw, id_policy)
            if gid_key != gid_raw:
                n_lookup_gene_changed += 1

            longest_tid = gene2longest.get(gid_key, MISSING_LONGEST_FILL)

            if longest_tid:
                hit_genes.add(gid_raw)
            else:
                miss_rows += 1
                if len(miss_gene_examples) < 50 and gid_raw not in miss_gene_examples:
                    miss_gene_examples.append(gid_raw)

            fo.write(f"{tid_out}\t{gid_out}\t{longest_tid}\n")

    log.info("写入：%s", out_longest)

    # 5) 写 report
    agat_stdout_preview = safe_preview(agat_stdout, 500)
    agat_stderr_preview = safe_preview(agat_stderr, 500)

    # gene 命中率统计：以 tx2gene.clean.tsv 的 gene_id 原值为口径
    n_genes_total = len(genes_set_raw)
    n_genes_hit = len(hit_genes)
    hit_rate = (n_genes_hit / n_genes_total) if n_genes_total > 0 else 0.0

    with out_report.open("w", encoding="utf-8") as rep:
        rep.write("metric\tvalue\n")
        rep.write(f"ref_gtf\t{ref_gtf}\n")
        rep.write(f"maps_dir\t{maps_dir}\n")
        rep.write(f"tx2gene_clean\t{tx2gene_clean}\n")
        rep.write(f"out_longest\t{out_longest}\n")
        rep.write("strategy\tAGAT_only\n")
        rep.write("resolve_rule\tprefer_mRNA_else_transcript\n")
        rep.write(f"agat_exe\t{shutil.which(AGAT_KEEP_LONGEST) or ''}\n")
        rep.write(f"agat_returncode\t{rc}\n")
        rep.write(f"agat_stdout_preview\t{agat_stdout_preview}\n")
        rep.write(f"agat_stderr_preview\t{agat_stderr_preview}\n")

        rep.write(f"n_rows_tx2gene_clean\t{len(rows_raw)}\n")
        rep.write(f"n_genes_in_tx2gene_clean\t{n_genes_total}\n")
        rep.write(f"n_genes_hit_longest\t{n_genes_hit}\n")
        rep.write(f"gene_hit_rate\t{hit_rate:.6f}\n")

        rep.write(f"n_genes_in_keep_longest_chosen\t{len(gene2longest)}\n")
        rep.write(f"keep_longest_multi_parent_warn\t{multi_parent_warn}\n")

        # keep_longest 结构统计
        for k in [
            "keep_longest_lines_mrna",
            "keep_longest_lines_transcript",
            "keep_longest_genes_total",
            "keep_longest_genes_with_mrna",
            "keep_longest_genes_with_transcript",
            "keep_longest_genes_with_both",
            "keep_longest_genes_with_transcript_only",
            "keep_longest_conflict_mrna",
            "keep_longest_conflict_transcript_only",
        ]:
            rep.write(f"{k}\t{kl_stats.get(k, 0)}\n")

        rep.write(f"missing_longest_rows\t{miss_rows}\n")
        rep.write(f"missing_longest_unique_gene_preview_count\t{len(miss_gene_examples)}\n")
        rep.write(f"lookup_gene_changed_by_cleanup_count\t{n_lookup_gene_changed}\n")

        rep.write("\n")
        rep.write("missing_longest_gene_preview\t\n")
        for g in miss_gene_examples:
            rep.write(f"{g}\n")

    log.info("写入：%s", out_report)

    # 6) 友好告警（不退出）
    if miss_rows > 0:
        log.warning("存在缺失 longest 的行：%d（gene 示例已写入 report）", miss_rows)
    if n_lookup_gene_changed > 0:
        log.warning("注意：有 %d 行的 gene_id 在 lookup 时经过清洗发生变化（但输出前两列仍严格照抄）。", n_lookup_gene_changed)

    log.info("完成：rows=%d, genes=%d, gene_hit_rate=%.3f", len(rows_raw), n_genes_total, hit_rate)


if __name__ == "__main__":
    main()

