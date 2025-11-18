#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py —— 整理 eggNOG 注释为富集所需映射（严格契约 + 强化正则 + 硬报错）
契约与要点：
  * 仅读 config.yaml（不接命令行）
  * 固定读取：
      - reference.emapper                      -> eggNOG-mapper 注释表（TSV）
      - results/03_maps/tx2gene.clean.tsv     -> transcript_id → gene_id 映射
  * 固定输出（dirs.annot）：
      - gene2go.tsv（列：gene_id, go_id）
      - gene2ko.tsv（列：gene_id, ko_id）
      - gene2pathway.tsv（列：gene_id, pathway_id）
      - universe_coverage.tsv（覆盖率摘要）
      - unmapped_queries.tsv（无法映射的 query 列表，便于排查）
  * ID 清理（ids.cleanup_07）：
      - 先去前缀 → 后缀正则循环裁剪 → 可选 '|' 右侧截断
      - 若 hard_fail_on_unstripped_suffix = true 且侦测到疑似后缀未被去除 → 直接报错（提示把规则写入 config）
  * 仅使用契约列名：emapper 必须包含列 query, GOs, KEGG_ko, KEGG_Pathway；否则报错（不做“猜列名”）
"""

from __future__ import annotations
import sys, re, csv, logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable, Set

CONFIG_PATH = "config.yaml"

# ------------------------- 读取配置（带默认） -------------------------

DEFAULTS: Dict[str, Any] = {
    "reference": {
        "emapper": "ref/annotations/Sinonovacula_constricta_annotations.tsv"
    },
    "dirs": {
        "maps":  "results/03_maps",
        "annot": "results/07_annot"
    },
    "ids": {
        "cleanup_07": {
            "enable": True,
            "remove_prefixes": [],
            "remove_suffix_regex": [
                r"\.[0-9]+$",
                r"\.t[0-9]+$",
                r"_t[0-9]+$",
                r"-RA$", r"-RB$",
                r"-mRNA-[0-9]+$",
                r"\.p[0-9]+$"
            ],
            "strip_after_bar": True,
            "bar_chars": ["|"],
            "hard_fail_on_unstripped_suffix": True
        }
    },
    "logging": {"level": "INFO", "timestamp": True}
}

def load_cfg(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception:
        print("[ERR] 缺少 PyYAML，请安装：mamba/conda install pyyaml", file=sys.stderr)
        raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        u = yaml.safe_load(f) or {}

    def merge(user, defaults):
        out = dict(defaults)
        for k, v in user.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out
    return merge(u, DEFAULTS)

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def setup_logging(cfg: Dict[str, Any]) -> None:
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    fmt = "%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True) else "[%(levelname)s] %(message)s"
    logging.basicConfig(level=level, format=fmt)

# ------------------------- ID 清理器（先前缀 → 后缀 → 条形截断） -------------------------

class IDCleaner:
    def __init__(self, rules: Dict[str, Any]):
        self.enable = bool(rules.get("enable", True))
        self.remove_prefixes: List[str] = list(rules.get("remove_prefixes", []))
        self.remove_suffix_regex = [re.compile(p) for p in rules.get("remove_suffix_regex", [])]
        self.strip_after_bar = bool(rules.get("strip_after_bar", True))
        self.bar_chars: List[str] = list(rules.get("bar_chars", ["|"]))
        self.hard_fail = bool(rules.get("hard_fail_on_unstripped_suffix", True))
        # “可疑后缀”集合用于侦测：原始 ID 若被这些模式命中，但清理后仍命中 → 认为“未去除”
        self.suspect_patterns = [
            re.compile(r"\.[0-9]+$"),
            re.compile(r"\.t[0-9]+$"),
            re.compile(r"_t[0-9]+$"),
            re.compile(r"-R[A-Z]$"),
            re.compile(r"-mRNA-[0-9]+$"),
            re.compile(r"\.p[0-9]+$")
        ]

    def clean(self, s: str) -> Tuple[str, bool]:
        """
        返回：cleaned_id, had_unstripped_suspect_suffix
        had_unstripped_suspect_suffix = True 表示检测到疑似后缀但仍未被清理
        """
        if not self.enable:
            return s, False
        orig = s

        # 1) 去前缀（按列出的前缀逐一剥离，一个命中即剥离一次，不死循环）
        x = orig
        for pref in self.remove_prefixes:
            if x.startswith(pref):
                x = x[len(pref):]

        # 2) 后缀正则裁剪（循环，直到无命中）
        changed = True
        while changed:
            changed = False
            for pat in self.remove_suffix_regex:
                m = pat.search(x)
                if m:
                    x = x[:m.start()]
                    changed = True

        # 3) 条形截断（strip_after_bar）
        if self.strip_after_bar and any(ch in x for ch in self.bar_chars):
            idxs = [x.find(ch) for ch in self.bar_chars if ch in x]
            cut = min(idxs) if idxs else -1
            if cut >= 0:
                x = x[:cut]

        # 4) 可疑后缀未被去除的检测
        unstripped = False
        for sp in self.suspect_patterns:
            if sp.search(orig) and sp.search(x):
                unstripped = True
                break

        return x.strip(), unstripped

# ------------------------- 读入与解析 -------------------------

def read_tx2gene(tx2gene_fp: Path) -> Dict[str, str]:
    if not tx2gene_fp.exists():
        raise FileNotFoundError(f"未找到 tx2gene.clean.tsv：{tx2gene_fp}（请先完成 03）")
    m: Dict[str, str] = {}
    with open(tx2gene_fp, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        need = ["transcript_id", "gene_id"]
        if reader.fieldnames is None or any(c not in reader.fieldnames for c in need):
            raise ValueError(f"tx2gene.clean.tsv 必须包含列：{', '.join(need)}")
        for r in reader:
            tid = (r.get("transcript_id") or "").strip()
            gid = (r.get("gene_id") or "").strip()
            if tid and gid:
                m[tid] = gid
    return m

def read_emapper(emapper_fp: Path) -> Iterable[Dict[str, str]]:
    if not emapper_fp.exists():
        raise FileNotFoundError(f"未找到 emapper 注释：{emapper_fp}")
    with open(emapper_fp, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        need = ["query", "GOs", "KEGG_ko", "KEGG_Pathway"]
        if reader.fieldnames is None or any(c not in reader.fieldnames for c in need):
            raise ValueError(f"emapper 注释必须包含列：{', '.join(need)}（契约）")
        for r in reader:
            yield {k: (r.get(k) or "").strip() for k in need}

# ------------------------- 解析列内容（GO / KEGG） -------------------------

_GO_RE = re.compile(r"GO:\d+")
_KO_RE = re.compile(r"(?:ko:)?(K\d{5})")
_PATH_RE = re.compile(r"(?:ko|map)(\d{5})")

def split_terms(s: str) -> List[str]:
    # 兼容逗号/分号/空格分隔
    if not s:
        return []
    parts = re.split(r"[;,]\s*|\s+", s)
    return [p for p in parts if p]

def extract_go(s: str) -> List[str]:
    return _GO_RE.findall(s or "")

def extract_kos(s: str) -> List[str]:
    return [m.group(1) for m in _KO_RE.finditer(s or "")]

def extract_paths(s: str) -> List[str]:
    return [f"map{m.group(1)}" for m in _PATH_RE.finditer(s or "")]

# ------------------------- 主流程 -------------------------

def main() -> None:
    cfg = load_cfg(Path(CONFIG_PATH))
    setup_logging(cfg)

    emapper_fp = Path(cfg["reference"]["emapper"])
    maps_dir   = Path(cfg["dirs"]["maps"])
    annot_dir  = Path(cfg["dirs"]["annot"])
    mkdir_p(annot_dir)

    tx2gene_fp = maps_dir / "tx2gene.clean.tsv"

    logging.info("========== 07 — emapper 注释整理（严格契约 + 强正则） ==========")
    logging.info(f"[Info] emapper  = {emapper_fp.resolve()}")
    logging.info(f"[Info] tx2gene  = {tx2gene_fp.resolve()}")
    logging.info(f"[Info] 输出目录  = {annot_dir.resolve()}")

    cleaner = IDCleaner(cfg.get("ids", {}).get("cleanup_07", {}))

    # 读映射
    tx2gene = read_tx2gene(tx2gene_fp)
    logging.info(f"[Info] tx2gene 映射条目数 = {len(tx2gene)}")

    # 读 emapper & 构建映射
    gene2go: Set[Tuple[str, str]] = set()
    gene2ko: Set[Tuple[str, str]] = set()
    gene2path: Set[Tuple[str, str]] = set()

    n_rows = 0
    n_mapped = 0
    n_unmapped = 0
    unresolved_suffix_ids: List[str] = []
    unmapped_queries: List[str] = []

    for row in read_emapper(emapper_fp):
        n_rows += 1
        q = row["query"]
        q_clean, unstripped = cleaner.clean(q)
        if unstripped and cleaner.hard_fail:
            unresolved_suffix_ids.append(q)
            # 先收集，统一报错；避免中途退出看不到整体情况

        gid = tx2gene.get(q_clean, "")
        if not gid:
            n_unmapped += 1
            unmapped_queries.append(q)
            continue

        n_mapped += 1
        # GO
        for go in extract_go(row["GOs"]):
            gene2go.add((gid, go))
        # KO
        for ko in extract_kos(row["KEGG_ko"]):
            gene2ko.add((gid, ko))
        # Pathway
        for pw in extract_paths(row["KEGG_Pathway"]):
            gene2path.add((gid, pw))

    # 若有未去除的可疑后缀且强制报错
    if cleaner.hard_fail and unresolved_suffix_ids:
        examples = unresolved_suffix_ids[:20]
        msg = (
            "[ERR] 侦测到疑似后缀未被去除的 query ID（请在 config.ids.cleanup_07.remove_suffix_regex 中补充规则）：\n"
            "示例（最多20条）：\n  - " + "\n  - ".join(examples)
        )
        raise SystemExit(msg)

    # 写输出
    def write_pairs(p: Path, header: Tuple[str, str], items: Set[Tuple[str, str]]):
        with open(p, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, dialect=csv.excel_tab)
            w.writerow(list(header))
            for a, b in sorted(items):
                w.writerow([a, b])

    out_go   = annot_dir / "gene2go.tsv"
    out_ko   = annot_dir / "gene2ko.tsv"
    out_path = annot_dir / "gene2pathway.tsv"
    write_pairs(out_go,   ("gene_id", "go_id"),       gene2go)
    write_pairs(out_ko,   ("gene_id", "ko_id"),       gene2ko)
    write_pairs(out_path, ("gene_id", "pathway_id"),  gene2path)

    # 覆盖率摘要 + 未映射列表
    out_meta = annot_dir / "universe_coverage.tsv"
    with open(out_meta, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["n_emapper_rows", "mapped_queries", "unmapped_queries", "gene2go_pairs", "gene2ko_pairs", "gene2pathway_pairs"])
        w.writerow([n_rows, n_mapped, n_unmapped, len(gene2go), len(gene2ko), len(gene2path)])

    out_unmapped = annot_dir / "unmapped_queries.tsv"
    if unmapped_queries:
        with open(out_unmapped, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, dialect=csv.excel_tab)
            w.writerow(["query"])
            for q in unmapped_queries:
                w.writerow([q])

    logging.info("========== 07 完成 ==========")
    logging.info(f"[Stat] emapper 行数={n_rows}；mapped={n_mapped}；unmapped={n_unmapped}")
    logging.info(f"[Out]  {out_go}")
    logging.info(f"[Out]  {out_ko}")
    logging.info(f"[Out]  {out_path}")
    logging.info(f"[Out]  {out_meta}")
    if unmapped_queries:
        logging.info(f"[Out]  {out_unmapped}（未映射 query 列表）")

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        print(str(e), file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"[ERR] 07 执行失败：{e}", file=sys.stderr); sys.exit(1)