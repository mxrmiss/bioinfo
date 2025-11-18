#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_build_tx2gene_map.py —— 生成 tx2gene 映射（transcript_id → gene_id）
契约：
  - 仅读取 config.yaml；
  - 输入：reference.ref_gtf（GFF3/GTF 皆可）；
  - 输出：results/03_maps/tx2gene.clean.tsv（表头：transcript_id, gene_id）
          results/03_maps/tx2gene.blacklist.tsv（列：gene_id, n_tx）
  - 对 ID 进行“末尾后缀剪切”（回滚启用），剪切规则来自 config，可覆盖默认正则；
  - 屏幕流式输出 + 关键统计；不产生非契约副产物。
"""

from __future__ import annotations
import sys, re, csv, logging
from pathlib import Path
from typing import Dict, Any, Iterable, Tuple

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
        # ID 清理策略（回滚为剪切后缀 + 强化正则）
        "cleanup": {
            "enable": True,
            # 先切掉“|”右侧（如有）
            "strip_after_bar": True,
            "bar_chars": ["|"],
            # 可选：去除手动声明的前缀（按出现顺序逐一剥离；不启用则留空）
            "remove_prefixes": [],
            # 末尾后缀剪切（默认强正则，覆盖常见注释后缀）
            # 说明：按列表顺序依次应用；命中即剪，直到均不命中
            "remove_suffix_regex": [
                r"\.[0-9]+$",           # .1 / .2 / .10
                r"\.t[0-9]+$",          # .t1 / .t2
                r"_t[0-9]+$",           # _t1 / _t2
                r"-RA$", r"-RB$",       # -RA / -RB（Ensembl 风格）
                r"-mRNA-[0-9]+$",       # -mRNA-1
                r"\.p[0-9]+$",          # .p1（少见）
                r"_gene$"               # _gene（保守兜底）
            ],
            # 大小写不改动；如需统一可在此开启（默认 False）
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
    try:
        import yaml
    except Exception:
        print("[ERR] 缺少 PyYAML，请先安装：mamba/conda install pyyaml", file=sys.stderr)
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

def stream_log_setup(cfg: Dict[str, Any]) -> None:
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    fmt = "%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True) else "[%(levelname)s] %(message)s"
    logging.basicConfig(level=level, format=fmt)

# ------------------------- ID 清理 -------------------------

class IDCleaner:
    def __init__(self, rules: Dict[str, Any]):
        self.enable = bool(rules.get("enable", True))
        self.strip_after_bar = bool(rules.get("strip_after_bar", True))
        self.bar_chars = list(rules.get("bar_chars", ["|"]))
        self.remove_prefixes = list(rules.get("remove_prefixes", []))
        self.remove_suffix_regex = [re.compile(p) for p in rules.get("remove_suffix_regex", [])]
        self.upper = bool(rules.get("upper", False))
        self.lower = bool(rules.get("lower", False))

    def clean(self, s: str) -> str:
        if not self.enable:
            return s
        x = s
        # 1) bar 切割（如 UniProt|GeneID 等）
        if self.strip_after_bar and any(ch in x for ch in self.bar_chars):
            idx = min([i for i in [x.find(ch) for ch in self.bar_chars] if i != -1], default=-1)
            if idx >= 0:
                x = x[:idx]
        # 2) 去前缀（手动列举）
        for pref in self.remove_prefixes:
            if x.startswith(pref):
                x = x[len(pref):]
        # 3) 去后缀（强化正则，循环直到不命中）
        changed = True
        while changed:
            changed = False
            for pat in self.remove_suffix_regex:
                m = pat.search(x)
                if m:
                    x = x[:m.start()]
                    changed = True
        # 4) 大小写（可选）
        if self.upper and not self.lower:
            x = x.upper()
        elif self.lower and not self.upper:
            x = x.lower()
        return x.strip()

# ------------------------- GFF/GTF 解析 -------------------------

def parse_attrs_gff3(attr: str) -> Dict[str, str]:
    out = {}
    for kv in attr.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out

def parse_attrs_gtf(attr: str) -> Dict[str, str]:
    out = {}
    # 形如：key "value"; key2 "value2";
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr):
        out[m.group(1)] = m.group(2)
    return out

def iter_tx_gene_from_gff(path: Path) -> Iterable[Tuple[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            _seqid, _src, ftype, _start, _end, _score, _strand, _phase, attr = parts
            a = parse_attrs_gff3(attr)
            # 常见组合：mRNA/transcript 的 ID + Parent（gene）
            if ftype in ("mRNA", "transcript"):
                tid = a.get("ID") or a.get("transcript_id")
                # Parent 可能逗号分隔多值，取第一个
                par = (a.get("Parent") or a.get("gene") or "").split(",")[0]
                if tid and par:
                    yield (tid, par)
            # 某些注释把 transcript 直接用 transcript_id/gene_id 施加在 exon/CDS 上
            elif ("transcript_id" in a and "gene_id" in a):
                yield (a["transcript_id"], a["gene_id"])

def iter_tx_gene_from_gtf(path: Path) -> Iterable[Tuple[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
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

def sniff_format(p: Path) -> str:
    name = p.name.lower()
    if name.endswith(".gtf"):
        return "gtf"
    return "gff3"

# ========================= 主流程 =========================

def main() -> None:
    cfg = load_cfg(Path(CONFIG_PATH))
    stream_log_setup(cfg)

    ref_gtf = Path(cfg["reference"]["ref_gtf"])
    maps_dir = Path(cfg["dirs"]["maps"])
    mkdir_p(maps_dir)

    out_clean = maps_dir / "tx2gene.clean.tsv"
    out_black = maps_dir / "tx2gene.blacklist.tsv"

    threshold = int(cfg["tx2gene"]["blacklist_tx_threshold"])
    cleaner = IDCleaner(cfg["tx2gene"].get("cleanup", {}))

    logging.info("========== 03 — 构建 tx2gene（剪切后缀·强化正则） ==========")
    logging.info(f"[Info] reference.ref_gtf = {ref_gtf.resolve()}")
    logging.info(f"[Info] dirs.maps         = {maps_dir.resolve()}")
    logging.info(f"[Info] blacklist阈值     = {threshold}")
    logging.info(f"[Info] cleanup规则启用   = {cleaner.enable}")

    if not ref_gtf.exists():
        raise FileNotFoundError(f"未找到注释文件：{ref_gtf}")

    fmt = sniff_format(ref_gtf)
    logging.info(f"[Info] 识别为 {fmt.upper()} 格式")

    # 收集原始配对
    raw_pairs: Dict[str, str] = {}         # tid -> gid
    n_lines = 0
    it = iter_tx_gene_from_gtf(ref_gtf) if fmt == "gtf" else iter_tx_gene_from_gff(ref_gtf)
    for tid, gid in it:
        n_lines += 1
        # 清理 ID（按强正则剪掉后缀；strip_after_bar等）
        tid_c = cleaner.clean(tid)
        gid_c = cleaner.clean(gid)
        if tid_c and gid_c:
            raw_pairs[tid_c] = gid_c

    if not raw_pairs:
        raise RuntimeError("未从注释文件中解析到 transcript_id/gene_id 配对，请检查注释格式与属性字段。")

    # 统计 gene → n_tx
    gene2count: Dict[str, int] = {}
    for _tid, g in raw_pairs.items():
        gene2count[g] = gene2count.get(g, 0) + 1

    # 写 clean 映射
    with open(out_clean, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["transcript_id", "gene_id"])
        for tid, gid in raw_pairs.items():
            w.writerow([tid, gid])

    # 写 blacklist
    with open(out_black, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["gene_id", "n_tx"])
        for g, c in sorted(gene2count.items(), key=lambda x: -x[1]):
            if c >= threshold:
                w.writerow([g, c])

    n_pairs = len(raw_pairs)
    n_genes = len(gene2count)
    n_black = sum(1 for _g, c in gene2count.items() if c >= threshold)
    logging.info("========== 03 完成 ==========")
    logging.info(f"[Stat] 解析行数≈{n_lines}；得到配对 {n_pairs}；基因数 {n_genes}；黑名单基因 {n_black}")
    logging.info(f"[Out]  {out_clean}")
    logging.info(f"[Out]  {out_black}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 03 执行失败：{e}", file=sys.stderr)
        sys.exit(1)