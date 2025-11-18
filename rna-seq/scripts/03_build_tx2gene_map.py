#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_build_tx2gene_map.py —— 从 GFF3/GTF 生成 tx2gene（严格契约版）

契约要点（严格遵守《转录组计划1.md》）：
  • 仅读取项目根 config.yaml 中的 reference.ref_gtf（必须存在），不做别名兼容。
  • 固定输出目录：results/03_maps/
  • 固定主表表头：transcript_id, gene_id
  • 不进行任何 ID 前/后缀或版本号清理（原样保留）。
  • 记录来源指纹：TX2GENE.SOURCE（绝对路径与 MD5），便于下游审计。
  • 黑名单：gene 关联的转录本数 ≥ 阈值（tx2gene.blacklist_threshold，默认 10）。

输出文件：
  - results/03_maps/tx2gene.clean.tsv
  - results/03_maps/tx2gene.blacklist.tsv
  - results/03_maps/TX2GENE.SOURCE
"""

from __future__ import annotations
import sys, csv, hashlib
from pathlib import Path
from typing import Dict, Tuple, Iterable, List
from datetime import datetime

# 依赖：PyYAML
try:
    import yaml
except Exception:
    sys.stderr.write("[ERR] 缺少 PyYAML，请先安装：mamba/conda install pyyaml\n")
    sys.exit(1)

# ===================== 顶部参数区（所有参数集中于此） =====================
CONFIG_PATH = "config.yaml"   # 只读，不接受命令行参数

DEFAULTS: Dict[str, dict] = {
    "dirs": {
        "maps": "results/03_maps"
    },
    "tx2gene": {
        "blacklist_threshold": 10
    },
    "logging": {
        "timestamp": True  # 屏幕输出是否带时间戳
    }
}
# ========================================================================

_TS = True  # 是否打印时间戳，由配置控制

def log(msg: str) -> None:
    """标准化日志输出；是否带时间戳由配置决定。"""
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {msg}" if _TS else msg)

def load_cfg(p: Path) -> Dict:
    """读取 YAML 并与默认值浅合并（用户优先）；严格要求 reference.ref_gtf 存在。"""
    if not p.exists():
        sys.stderr.write(f"[ERR] 未找到配置文件：{p}\n")
        sys.exit(1)
    with open(p, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}

    # 浅合并（用户优先）
    out = dict(DEFAULTS)
    for k, v in user.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = {**out[k], **v}
        else:
            out[k] = v

    # 严格取 reference.ref_gtf（契约规定的唯一位置）
    ref_block = out.get("reference", {})
    if not isinstance(ref_block, dict):
        sys.stderr.write("[ERR] 配置缺少块：reference\n")
        sys.exit(1)
    ref_gtf = ref_block.get("ref_gtf")
    if not isinstance(ref_gtf, str) or not ref_gtf.strip():
        sys.stderr.write("[ERR] 配置缺少键：reference.ref_gtf（必须）\n")
        sys.exit(1)
    out["reference"]["ref_gtf"] = ref_gtf.strip()

    return out

def md5sum(path: Path) -> str:
    """计算文件 MD5（写入 SOURCE 用）。"""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(8192), b""):
            h.update(ch)
    return h.hexdigest()

# -------------------- 属性解析 --------------------
def parse_attr_gff3(s: str) -> Dict[str, str]:
    """解析 GFF3 属性字段：key=value;key2=value2;..."""
    out: Dict[str, str] = {}
    for part in s.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k.strip()] = v.strip()
    return out

def parse_attr_gtf(s: str) -> Dict[str, str]:
    """解析 GTF 属性字段：key "value"; key2 "value2"; ..."""
    out: Dict[str, str] = {}
    for seg in [x.strip() for x in s.rstrip(";").split(";") if x.strip()]:
        sp = seg.split()
        if len(sp) >= 2:
            out[sp[0].strip()] = " ".join(sp[1:]).strip().strip('"')
    return out

# -------------------- 提取 (transcript_id, gene_id) --------------------
def iter_tx_gene_pairs(gtf_path: Path) -> Iterable[Tuple[str, str]]:
    """
    仅产出“转录本→基因”的二元组 (transcript_id, gene_id)。
    - 对 GFF3：feature=mRNA 或 transcript 的行（属性 ID=<tx> 与 Parent=<gene>）
    - 对 GTF ：属性同时包含 transcript_id 与 gene_id
    """
    is_gff3 = gtf_path.suffix.lower() in (".gff3", ".gff")
    with open(gtf_path, "r", encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = (parts[2] or "").lower()
            attrs = parts[8]

            if is_gff3:
                if feature not in ("mrna", "transcript"):
                    continue
                a = parse_attr_gff3(attrs)
                tx = a.get("ID") or a.get("transcript_id")
                parent = a.get("Parent") or a.get("gene_id")
                if not tx or not parent:
                    continue
                # Parent 可能多值（逗号分隔），常规取第一个 gene
                gene = parent.split(",")[0].strip()
                yield (tx.strip(), gene)
            else:
                a = parse_attr_gtf(attrs)
                tx = (a.get("transcript_id") or "").strip()
                gene = (a.get("gene_id") or "").strip()
                if tx and gene:
                    yield (tx, gene)

# -------------------- 主流程 --------------------
def main() -> None:
    cfg = load_cfg(Path(CONFIG_PATH))
    global _TS
    _TS = bool(cfg.get("logging", {}).get("timestamp", True))

    maps_dir = Path(cfg["dirs"]["maps"])
    maps_dir.mkdir(parents=True, exist_ok=True)

    out_clean = maps_dir / "tx2gene.clean.tsv"
    out_black = maps_dir / "tx2gene.blacklist.tsv"
    out_src   = maps_dir / "TX2GENE.SOURCE"

    ref_gtf = Path(cfg["reference"]["ref_gtf"])
    if not ref_gtf.exists():
        sys.stderr.write(f"[ERR] 未找到注释文件：{ref_gtf}\n")
        sys.exit(1)

    log("========== 03 — 构建 tx2gene 映射 ==========")
    log(f"[Info] reference.ref_gtf = {ref_gtf.resolve()}")

    # 采集 (tx, gene)
    pairs: List[Tuple[str, str]] = []
    for tx, gene in iter_tx_gene_pairs(ref_gtf):
        if tx and gene:
            pairs.append((tx, gene))
    if not pairs:
        sys.stderr.write("[ERR] 未从注释中解析到任何 (transcript_id, gene_id)；请检查 GFF3/GTF 的 feature/属性是否标准。\n")
        sys.exit(1)

    # 去重
    seen = set()
    uniq: List[Tuple[str, str]] = []
    for tx, gene in pairs:
        k = (tx, gene)
        if k not in seen:
            uniq.append(k); seen.add(k)

    # 写 clean 主表
    with open(out_clean, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["transcript_id", "gene_id"])
        w.writerows(uniq)
    log(f"[Out] {out_clean}  行数：{len(uniq)}")

    # 黑名单（gene→转录本数 ≥ 阈值）
    thr = int(cfg.get("tx2gene", {}).get("blacklist_threshold", DEFAULTS["tx2gene"]["blacklist_threshold"]))
    g2n: Dict[str, int] = {}
    for _, g in uniq:
        g2n[g] = g2n.get(g, 0) + 1
    bad = [(g, n) for g, n in g2n.items() if n >= thr]
    with open(out_black, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["gene_id", "n_transcripts", "threshold"])
        for g, n in sorted(bad, key=lambda x: (-x[1], x[0])):
            w.writerow([g, n, thr])
    log(f"[Out] {out_black}  触发阈值基因数：{len(bad)}（阈值={thr}）")

    # 来源指纹（供 05 自检同源，不影响结果）
    try:
        with open(out_src, "w", encoding="utf-8") as f:
            f.write(f"ref_gtf={ref_gtf.resolve()}\n")
            f.write(f"md5={md5sum(ref_gtf)}\n")
    except Exception as e:
        log(f"[Warn] 写 SOURCE 指纹失败（可忽略）：{e}")

    log("========== 03 完成 ✅ ==========")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[ERR] 03 执行失败：{e}\n")
        sys.exit(1)

