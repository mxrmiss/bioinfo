#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_build_tx2gene_map.py — 生成 tx2gene 映射表（transcript_id → gene_id）
契约要点：
  - 读取 GFF3/GTF 文件，提取 gene_id 和 transcript_id。
  - 输出文件：results/03_maps/tx2gene.clean.tsv
  - 不剪切 gene_id 后缀，保持原始的 gene_id。
"""

from __future__ import annotations
import sys
import csv
import logging
from pathlib import Path
from typing import Dict, Any, List
import gffutils

# 配置文件路径
CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "data": {
        "gff_file": "ref/Sinonovacula_constricta_genome.gff3",
        "tx2gene_out": "results/03_maps/tx2gene.clean.tsv"
    },
    "logging": {
        "level": "INFO",
        "timestamp": True,
        "file": "logs/03_build_tx2gene_map.log"
    }
}

def load_yaml(path: Path) -> Dict[str, Any]:
    """加载配置文件"""
    try:
        import yaml
    except ImportError:
        print("[ERR] 需要 PyYAML，请安装：mamba/conda install pyyaml", file=sys.stderr)
        raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        cfg = (yaml.safe_load(f) or {})
    return cfg

def mkdir_p(p: Path) -> None:
    """创建目录"""
    p.mkdir(parents=True, exist_ok=True)

def read_gff(gff_file: str) -> gffutils.FeatureDB:
    """读取 GFF 文件并创建数据库"""
    if not Path(gff_file).exists():
        raise FileNotFoundError(f"GFF 文件未找到：{gff_file}")
    return gffutils.create_db(gff_file, dbfn=gff_file + ".db", force=True, keep_order=True)

def extract_tx2gene(gff_db: gffutils.FeatureDB) -> Dict[str, str]:
    """从 GFF 文件中提取 transcript_id 和 gene_id"""
    tx2gene = {}
    for feature in gff_db.features_of_type("mRNA"):
        transcript_id = feature.attributes.get("ID", [None])[0]
        gene_id = feature.attributes.get("Parent", [None])[0]
        if transcript_id and gene_id:
            tx2gene[transcript_id] = gene_id
    return tx2gene

def write_tx2gene(tx2gene: Dict[str, str], output_file: Path) -> None:
    """将 tx2gene 映射写入文件"""
    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, dialect=csv.excel_tab)
        writer.writerow(["transcript_id", "gene_id"])  # 表头：按契约要求
        for transcript_id, gene_id in tx2gene.items():
            writer.writerow([transcript_id, gene_id])

def main() -> None:
    # 加载配置
    cfg = load_yaml(Path(CONFIG_PATH))
    gff_file = cfg["data"]["gff_file"]
    tx2gene_out = Path(cfg["data"]["tx2gene_out"])

    # 创建输出目录
    mkdir_p(tx2gene_out.parent)

    # 读取 GFF 文件
    gff_db = read_gff(gff_file)

    # 提取 tx2gene 映射
    tx2gene = extract_tx2gene(gff_db)

    # 写入输出文件
    write_tx2gene(tx2gene, tx2gene_out)

    logging.info(f"成功生成 tx2gene 映射：{tx2gene_out}")

if __name__ == "__main__":
    try:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
        main()
    except Exception as e:
        print(f"[ERR] 03 执行失败：{e}", file=sys.stderr)
        sys.exit(1)