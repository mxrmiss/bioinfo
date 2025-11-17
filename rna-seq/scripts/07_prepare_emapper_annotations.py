#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_prepare_emapper_annotations.py — 规范化 eggnog-mapper 注释为下游富集“原材料”
特性：
  1) 稳健读取 emapper *.annotations.tsv / *.tsv.gz：
     - 先嗅探以 '#' 开头的表头行（如 '#query ...'）
     - 再以 comment='#' 只读数据段，显式设置列名
  2) 列名别名对齐 + 缺列补空，统一输出规范化总表 emapper.parsed.tsv
  3) 将 emapper 的 query 字段规整为 gene_id（去物种前缀与转录本后缀），与 05/06 的基因矩阵天然对齐
  4) 生成三张“长表”原材料（TSV）：
       - gene2go.tsv        （gene_id, GO）
       - gene2ko.tsv        （gene_id, KO）
       - gene2pathway.tsv   （gene_id, pathway_id, pathway_name）
     其中 KEGG_Pathway 支持 koXXXXX / mapXXXXX / 'koXXXXX:Name' 混排并尽量解析名称
  5) 路径参数：默认优先读取项目根 config.yaml；若缺键则回退为脚本内默认。无绘图、无 Excel。
"""

from __future__ import annotations
import sys, io, os, gzip, re
from pathlib import Path
from typing import List, Dict, Any
import yaml
import pandas as pd

# ======================== 默认参数（可被 config.yaml 覆盖） ========================
LOCAL = {
    "config_yaml": "config.yaml",
    "tables": {
        # 皇上指定：使用 *.tsv，不用 xlsx；这里也会强校验
        "emapper": "ref/annotations/Sinonovacula_constricta_annotations.tsv"
    },
    "paths": {
        "anno_dir": "ref/annotations",
        "logs_dir": "logs"
    },
    "emapper": {
        # 列名别名映射（不同 emapper 版本可能有 query_name / target_name / kegg_ko 等）
        "col_alias": {
            "query_name": "query",
            "target_name": "seed_ortholog",
            "og": "eggNOG_OGs",
            "kegg_ko": "KEGG_ko",
            "kegg_pathway": "KEGG_Pathway",
            "kegg_module": "KEGG_Module",
            "kegg_reaction": "KEGG_Reaction",
            "brite": "BRITE",
            "cazy": "CAZy",
            "pfams": "PFAMs",
        },
        # 我们关心的标准列集合（缺失将补空列，以便下游稳定）
        "wanted_cols": [
            "query", "seed_ortholog", "evalue", "score",
            "eggNOG_OGs", "max_annot_lvl", "COG_category",
            "Description", "Preferred_name",
            "GOs", "EC",
            "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
            "KEGG_Reaction", "KEGG_rclass",
            "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction",
            "PFAMs"
        ],
        # 多值字段分隔符（eggnog-mapper 常用逗号）
        "multi_sep": {
            "GOs": ",",
            "KEGG_ko": ",",
            "KEGG_Pathway": ","
        }
    }
}

# ================================ 工具函数 ================================
def to_abs(p: str | Path) -> Path:
    """将相对路径转为基于当前工作目录的绝对路径（不修改原结构）"""
    p = Path(p)
    return p if p.is_absolute() else (Path.cwd() / p)

def load_yaml(path: str) -> Dict[str, Any]:
    """读取 YAML（缺失返回空 dict）"""
    p = to_abs(path)
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def merge_cfg(a: dict, b: dict) -> dict:
    """浅递归合并 dict：b 覆盖 a，同层优先 b"""
    if not isinstance(a, dict): return b
    if not isinstance(b, dict): return a
    out = dict(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge_cfg(out[k], v)
        else:
            out[k] = v
    return out

def ensure_dir(p: str | Path) -> Path:
    """确保目录存在"""
    p = to_abs(p)
    p.mkdir(parents=True, exist_ok=True)
    return p

def open_text_auto(path: Path):
    """自动识别 .gz 或普通文本"""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")

def sniff_emapper_header(path: Path) -> List[str]:
    """
    提取以 '#' 开头的表头行（如 '#query\t...'），去掉 '#' 后按制表符切列名。
    若未找到，返回空列表（后续用 wanted_cols 兜底）。
    """
    with open_text_auto(path) as fin:
        for line in fin:
            if not line:
                continue
            if line.startswith("#"):
                head = line.lstrip("#").rstrip("\n")
                cols = head.split("\t")
                if len(cols) >= 3:
                    return cols
            else:
                break
    return []

def read_emapper_tsv(path: Path, wanted_cols: List[str], alias_map: Dict[str, str]) -> pd.DataFrame:
    """
    稳健读取 emapper 的 TSV：
      1) 先 sniff '#'-header 作为列名；若无则用 wanted_cols 兜底
      2) 以 comment='#' 读取数据段；显式设置 names=列名；dtype=str
      3) 应用别名映射；缺列补空；仅保留“关心列”的超集（关心列 + 其余列）
    """
    cols = sniff_emapper_header(path)
    if not cols:
        cols = wanted_cols[:]  # 兜底列集合
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=cols,
        dtype=str
    ).fillna("")
    # 别名标准化
    col_map = {}
    for c in df.columns:
        key = c.strip()
        std = alias_map.get(key, alias_map.get(key.lower(), key))
        col_map[c] = std
    df = df.rename(columns=col_map)
    # 缺列补空
    for c in wanted_cols:
        if c not in df.columns:
            df[c] = ""
    keep = [c for c in wanted_cols if c in df.columns]
    others = [c for c in df.columns if c not in keep]
    df = df[keep + others]
    return df

def split_multi(s: str, sep: str) -> list[str]:
    """将多值字段按分隔符拆分为去空白的条目列表"""
    s = (s or "").strip()
    if not s:
        return []
    return [x.strip() for x in s.split(sep) if x.strip()]

# ---------- 基因 ID 规整：去物种前缀与转录本/版本后缀 ----------
_Q_SPLIT = re.compile(r"\|")
_TX_DOT_VER = re.compile(r"\.[0-9]+$")             # 例如 .1 /.2
_TX_SUFFIX = re.compile(r"[-_](?:RA|RB|RC|PA|PB)$", re.IGNORECASE)  # 常见转录本/蛋白后缀
def normalize_gene_id(query: str) -> str:
    """
    将 emapper 的 query 规整为基因级 ID：
      - 若含 'prefix|ID'，取最后一个 '|' 右侧的 ID
      - 去掉转录本/版本后缀（.数字 / -RA/-PA 等）
    """
    q = (query or "").strip()
    if not q:
        return q
    # 取最后一个 |
    parts = _Q_SPLIT.split(q)
    core = parts[-1] if parts else q
    core = _TX_DOT_VER.sub("", core)
    core = _TX_SUFFIX.sub("", core)
    return core

# ---------- Pathway 标准化 ----------
_PW_RE = re.compile(r"^(?:(?:ko|map))?(\d{5})(?::\s*(.+))?$", re.IGNORECASE)
def parse_pathway_token(token: str) -> tuple[str, str]:
    """
    将 Pathway 的条目解析为 (pathway_id, pathway_name)：
      - 支持 'ko04151' / 'map04151' / '04151' / 'ko04151:Apoptosis' 等
      - 未匹配到则返回 (原条目去空白, "")
    """
    t = (token or "").strip()
    if not t:
        return "", ""
    m = _PW_RE.match(t)
    if m:
        pid = f"ko{m.group(1)}"
        pname = (m.group(2) or "").strip()
        return pid, pname
    return t, ""

# ================================ 主流程 ================================
def main():
    # 1) 合并配置：config.yaml 覆盖 LOCAL
    cfg = merge_cfg(LOCAL, load_yaml(LOCAL["config_yaml"]))

    emap_fp  = to_abs(cfg.get("tables", {}).get("emapper", LOCAL["tables"]["emapper"]))
    anno_dir = ensure_dir(cfg.get("paths", {}).get("anno_dir", LOCAL["paths"]["anno_dir"]))
    ensure_dir(cfg.get("paths", {}).get("logs_dir", LOCAL["paths"]["logs_dir"]))

    # 2) 基础校验
    if not emap_fp.exists():
        print(f"[ERR] 缺少 emapper 注释文件：{emap_fp}", file=sys.stderr)
        sys.exit(1)
    if not (str(emap_fp).endswith(".tsv") or str(emap_fp).endswith(".tsv.gz")):
        print(f"[ERR] 输入不是 TSV/TSV.GZ：{emap_fp}\n请提供 eggnog-mapper 的 *.annotations.tsv（或 .tsv.gz）。", file=sys.stderr)
        sys.exit(1)

    # 3) 读取 emapper 表
    alias  = cfg["emapper"]["col_alias"]
    wanted = cfg["emapper"]["wanted_cols"]
    df = read_emapper_tsv(emap_fp, wanted_cols=wanted, alias_map=alias)

    # 4) 生成规范化总表 + gene_id 列
    if "query" not in df.columns:
        print("[ERR] emapper 表缺少 'query' 列，无法继续。", file=sys.stderr)
        sys.exit(1)

    df["gene_id"] = [normalize_gene_id(q) for q in df["query"]]

    # 将 Pathway 字段中显式 mapXXXXX 统一为 koXXXXX:Name 形式在长表解析时处理
    parsed_out = anno_dir / "emapper.parsed.tsv"
    df.to_csv(parsed_out, sep="\t", index=False)
    print(f"[OK] 规范化表写出：{parsed_out}（{len(df)} 行）")

    # 5) 产出三张“长表”
    gene_col = "gene_id"

    # 5.1 GO
    go_sep = cfg["emapper"]["multi_sep"].get("GOs", ",")
    rows_go: list[tuple[str, str]] = []
    gos = df["GOs"] if "GOs" in df.columns else []
    for g, s in zip(df[gene_col], gos):
        for go in split_multi(s, go_sep):
            rows_go.append((g, go))
    g2go_out = anno_dir / "gene2go.tsv"
    pd.DataFrame(rows_go, columns=["gene_id", "GO"]).to_csv(g2go_out, sep="\t", index=False)
    print(f"[OK] gene2go 写出：{g2go_out}（{len(rows_go)} 条）")

    # 5.2 KO
    ko_sep = cfg["emapper"]["multi_sep"].get("KEGG_ko", ",")
    rows_ko: list[tuple[str, str]] = []
    kos = df["KEGG_ko"] if "KEGG_ko" in df.columns else []
    for g, s in zip(df[gene_col], kos):
        for ko in split_multi(s, ko_sep):
            # emapper 有时会带 'ko:' 前缀，这里统一去掉 'ko:' 保留 K 号或直接保留原样都可
            ko_clean = ko.replace("ko:", "").strip()
            rows_ko.append((g, ko_clean))
    g2ko_out = anno_dir / "gene2ko.tsv"
    pd.DataFrame(rows_ko, columns=["gene_id", "KO"]).to_csv(g2ko_out, sep="\t", index=False)
    print(f"[OK] gene2ko 写出：{g2ko_out}（{len(rows_ko)} 条）")

    # 5.3 Pathway
    pw_sep = cfg["emapper"]["multi_sep"].get("KEGG_Pathway", ",")
    rows_pw: list[tuple[str, str, str]] = []
    pws = df["KEGG_Pathway"] if "KEGG_Pathway" in df.columns else []
    for g, s in zip(df[gene_col], pws):
        for it in split_multi(s, pw_sep):
            pid, pname = parse_pathway_token(it)
            if pid:  # 过滤空白
                rows_pw.append((g, pid, pname))
    g2pw_out = anno_dir / "gene2pathway.tsv"
    pd.DataFrame(rows_pw, columns=["gene_id", "pathway_id", "pathway_name"]).to_csv(g2pw_out, sep="\t", index=False)
    print(f"[OK] gene2pathway 写出：{g2pw_out}（{len(rows_pw)} 条）")

    print("完成 ✅")

if __name__ == "__main__":
    main()

