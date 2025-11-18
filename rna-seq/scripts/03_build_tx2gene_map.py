#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_build_tx2gene_map.py —— 生成 tx2gene 映射表（严格契约 · 仅清理 gene_id · transcript_id 原样）
契约要点（与约定文件一致）：
  1) 仅读取 config.yaml（不接命令行参数）。
  2) 输入：reference.ref_gtf（GFF3 或 GTF 均可）。
  3) 输出目录：dirs.maps
     - results/03_maps/tx2gene.clean.tsv        （表头：transcript_id, gene_id）
     - results/03_maps/tx2gene.blacklist.tsv    （表头：gene_id, n_tx；阈值见 tx2gene.blacklist_tx_threshold）
  4) ID 处理：只对 gene_id 应用“前缀/后缀规则”，transcript_id 保持原样（须与 quant.sf 的 Name 完全一致）。
  5) 屏幕流式日志，关键统计与契约自检。
"""

from __future__ import annotations
import sys, re, csv, logging, gzip
from pathlib import Path
from typing import Dict, Any, Iterable, Tuple, List

# ========================= 配置与默认值 =========================

CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "reference": {
        "ref_gtf": "ref/Sinonovacula_constricta_genome.gff3"
    },
    "dirs": {
        "maps": "results/03_maps"
    },
    "tx2gene": {
        # 黑名单阈值：一个 gene 对应的转录本数 >= 此阈值 → 列入 blacklist
        "blacklist_tx_threshold": 10,
        # 仅对 gene_id 生效的清理规则（transcript_id 一律不清理，必须保持与 quant.sf 一致）
        "cleanup": {
            "enable": True,
            # 先去前缀（按列出的前缀逐一剥离；不死循环）
            "remove_prefixes": [],
            # 再按末尾后缀正则循环裁剪（命中即剪，直到无命中）
            "remove_suffix_regex": [
                r"\.[0-9]+$",           # .1 / .2 / .10
                r"\.t[0-9]+$",          # .t1 / .t2
                r"_t[0-9]+$",           # _t1 / _t2
                r"-RA$", r"-RB$",       # -RA / -RB（Ensembl 风格）
                r"-mRNA-[0-9]+$",       # -mRNA-1
                r"\.p[0-9]+$",          # .p1（少见）
                r"_gene$"               # _gene（兜底）
            ],
            # 可选：若 gene_id 中存在 "X|Y"，仅保留左侧（如 UniProt|GeneID）
            "strip_after_bar": True,
            "bar_chars": ["|"],
            # 大小写不改动；如需统一可开启其一（默认 False）
            "upper": False,
            "lower": False
        }
    },
    "logging": {
        "level": "INFO",
        "timestamp": True
    }
}

# ========================= 工具函数 =========================

def load_cfg(path: Path) -> Dict[str, Any]:
    """读取 config.yaml 并与 DEFAULTS 递归合并（用户优先）。"""
    try:
        import yaml
    except Exception:
        print("[ERR] 缺少 PyYAML，请先安装：mamba/conda install pyyaml", file=sys.stderr)
        raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        user = (yaml.safe_load(f) or {})

    def merge(defaults, userdict):
        out = dict(defaults)
        for k, v in userdict.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(out[k], v)
            else:
                out[k] = v
        return out

    return merge(DEFAULTS, user)

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def setup_logging(cfg: Dict[str, Any]) -> None:
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    fmt = "%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True) else "[%(levelname)s] %(message)s"
    logging.basicConfig(level=level, format=fmt)

def open_any(path: Path):
    """支持 .gz 与普通文本的通用打开。"""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

# ------------------------- ID 清理（仅 gene_id 使用） -------------------------

class GeneIDCleaner:
    """仅清理 gene_id；transcript_id 不调用本清理器。"""
    def __init__(self, rules: Dict[str, Any]):
        self.enable = bool(rules.get("enable", True))
        self.remove_prefixes: List[str] = list(rules.get("remove_prefixes", []))
        self.remove_suffix_regex = [re.compile(p) for p in rules.get("remove_suffix_regex", [])]
        self.strip_after_bar = bool(rules.get("strip_after_bar", True))
        self.bar_chars: List[str] = list(rules.get("bar_chars", ["|"]))
        self.upper = bool(rules.get("upper", False))
        self.lower = bool(rules.get("lower", False))

    def clean(self, gid: str) -> str:
        if not self.enable:
            return gid
        x = gid

        # 1) 去前缀（列出即剥离一次）
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

        # 4) 大小写（可选）
        if self.upper and not self.lower:
            x = x.upper()
        elif self.lower and not self.upper:
            x = x.lower()

        return x.strip()

# ------------------------- GFF/GTF 解析 -------------------------

def parse_attrs_gff3(attr: str) -> Dict[str, str]:
    """解析 GFF3 第9列为字典：key=value;key2=value2"""
    out: Dict[str, str] = {}
    for kv in attr.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out

def parse_attrs_gtf(attr: str) -> Dict[str, str]:
    """解析 GTF 第9列为字典：key "value"; key2 "value2";"""
    out: Dict[str, str] = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr):
        out[m.group(1)] = m.group(2)
    return out

def sniff_format(p: Path) -> str:
    """根据文件扩展名嗅探格式：gtf/其它视为gff3。"""
    name = p.name.lower()
    if name.endswith(".gtf") or name.endswith(".gtf.gz"):
        return "gtf"
    return "gff3"

def iter_tx_gene_from_gff3(path: Path):
    """
    从 GFF3 中迭代 (transcript_id, gene_id)：
      - mRNA/transcript：用 ID 作为 transcript_id，Parent 的第一个值作为 gene_id
      - 少数注释把 transcript_id/gene_id 直接挂在exon/CDS上（兜底）
    """
    with open_any(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            _seqid, _src, ftype, _start, _end, _score, _strand, _phase, attr = parts
            a = parse_attrs_gff3(attr)
            if ftype in ("mRNA", "transcript"):
                tid = a.get("ID") or a.get("transcript_id")
                par = (a.get("Parent") or a.get("gene") or "").split(",")[0]
                if tid and par:
                    yield (tid, par)
            elif ("transcript_id" in a and "gene_id" in a):
                yield (a["transcript_id"], a["gene_id"])

def iter_tx_gene_from_gtf(path: Path):
    """
    从 GTF 中迭代 (transcript_id, gene_id)：
      - 大多数行（exon/CDS/transcript 等）都带 transcript_id 与 gene_id，任选读取。
    """
    with open_any(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attr = parts[8]
            a = parse_attrs_gtf(attr)
            tid = a.get("transcript_id")
            gid = a.get("gene_id")
            if tid and gid:
                yield (tid, gid)

# ========================= 主流程 =========================

def main() -> None:
    # 读取配置并初始化日志
    cfg = load_cfg(Path(CONFIG_PATH))
    setup_logging(cfg)

    ref_gtf = Path(cfg["reference"]["ref_gtf"])
    maps_dir = Path(cfg["dirs"]["maps"])
    mkdir_p(maps_dir)

    out_clean = maps_dir / "tx2gene.clean.tsv"
    out_black = maps_dir / "tx2gene.blacklist.tsv"

    threshold = int(cfg["tx2gene"]["blacklist_tx_threshold"])
    cleaner = GeneIDCleaner(cfg["tx2gene"].get("cleanup", {}))

    logging.info("========== 03 — 构建 tx2gene（仅清理 gene_id；transcript_id 原样） ==========")
    logging.info(f"[Info] reference.ref_gtf = {ref_gtf.resolve()}")
    logging.info(f"[Info] dirs.maps         = {maps_dir.resolve()}")
    logging.info(f"[Info] blacklist 阈值    = {threshold}")
    logging.info(f"[Info] gene_id 清理启用  = {cleaner.enable}")

    if not ref_gtf.exists():
        raise FileNotFoundError(f"未找到注释文件：{ref_gtf}")

    fmt = sniff_format(ref_gtf)
    logging.info(f"[Info] 识别为 {fmt.upper()} 格式")

    # 收集原始 (tid, gid) 配对；transcript_id 原样保留，gene_id 走清理器
    tid_to_gid: Dict[str, str] = {}
    conflicts = 0
    n_pairs = 0

    it = iter_tx_gene_from_gtf(ref_gtf) if fmt == "gtf" else iter_tx_gene_from_gff3(ref_gtf)
    for tid, gid in it:
        n_pairs += 1
        tid_clean = tid.strip()              # transcript_id 原样（严格契约）
        gid_clean = cleaner.clean(gid.strip())
        if not tid_clean or not gid_clean:
            continue
        # 若同一 tid 指向多个 gid（异常），保留首次并计冲突
        if tid_clean in tid_to_gid and tid_to_gid[tid_clean] != gid_clean:
            conflicts += 1
            continue
        tid_to_gid[tid_clean] = gid_clean

    if not tid_to_gid:
        raise RuntimeError("未从注释文件中解析到 transcript_id/gene_id 配对，请检查注释格式与属性字段。")

    # 统计 gene → n_tx
    gene2count: Dict[str, int] = {}
    for _tid, g in tid_to_gid.items():
        gene2count[g] = gene2count.get(g, 0) + 1

    # 写 clean 映射（按 transcript_id 排序，保证可复现）
    with open(out_clean, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["transcript_id", "gene_id"])
        for tid in sorted(tid_to_gid.keys()):
            w.writerow([tid, tid_to_gid[tid]])

    # 写 blacklist
    with open(out_black, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["gene_id", "n_tx"])
        for g, c in sorted(gene2count.items(), key=lambda x: -x[1]):
            if c >= threshold:
                w.writerow([g, c])

    logging.info("========== 03 完成 ==========")
    logging.info(f"[Stat] 原始配对行数≈{n_pairs}；有效配对={len(tid_to_gid)}；gene 数={len(gene2count)}；tid→多 gid 冲突={conflicts}")
    logging.info(f"[Out]  {out_clean}")
    logging.info(f"[Out]  {out_black}")
    logging.info("[Hint] transcript_id 与 quant.sf 的 Name 必须逐字一致；若 05 报不匹配，请核对 02/04 的 transcripts.fa 头与本表第一列。")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 03 执行失败：{e}", file=sys.stderr)
        sys.exit(1)