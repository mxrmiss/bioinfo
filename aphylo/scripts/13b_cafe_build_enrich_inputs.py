#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13b_cafe_build_enrich_inputs.py —— CAFE5 → 富集分析桥接（构建 ORA 输入）

输入：
  1) config.yaml
       - paths.cafe_agg_dir
       - inputs.all_members_tsv = all_pep2cds_resolved.tsv
       - species.alias_map
       - cafe5.enrich_bridge:
           enable: true/false
           outputs_dir: "results/07_cafeagg/enrich_inputs"
           species_sets:
             <tag>: [Species1, Species2, ...]
           annotations:   # 可选
             gene2go: <path OR null>
             gene2pathway: <path OR null>

  2) 来自 13 步的聚合表：
       paths.cafe_agg_dir/cafe_significant_families_no_highfail.tsv
         必须含列：model, subset, family, P, Q, significant, source

  3) all_pep2cds_resolved.tsv（来自 phylo 发布包）
       列：OG  Species  protein_id  cds_id

输出结构：
  paths.cafe_agg_dir/enrich_inputs/
    README.txt
    <tag>/
      up.list          # 显著 family 的成员基因
      down.list        # 当前版本与 up.list 相同（代表“快速变化家族”全集）
      background.list  # 所有可评估 family 的成员基因
      meta.tsv         # 按 vNext 外部 tag 规范写 meta 行

注意：
  - 本版本不区分 expansion / contraction，仅提供“FDR 显著家族”的总集合；
  - 将来若基于 Base_asr.tre 能安全解析出方向，可以扩展为真正的 up/down。
"""

from __future__ import annotations

# ============================ 顶部可调区 ============================
CONFIG_PATH: str = "config.yaml"
LOG_LEVEL: str = "INFO"
LOG_FILE_BASENAME: str = "13b_cafe_build_enrich_inputs.log"
# ====================================================================

import sys
import yaml
import logging
import traceback
from pathlib import Path
from typing import Dict, List, Set
from datetime import datetime

# ---------------- 基础工具 ----------------

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def read_lines(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\r\n") for ln in f]

def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")

def setup_logging(logs_dir: Path) -> None:
    mkdir_p(logs_dir)
    log_file = logs_dir / LOG_FILE_BASENAME
    fmt = "[%(asctime)s] [%(levelname)s] %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format=fmt,
        datefmt=datefmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, encoding="utf-8"),
        ],
    )
    logging.info(f"[init] 使用配置文件：{CONFIG_PATH}")
    logging.info(f"[init] 日志写入：{log_file}")

# ---------------- 读取 CAFE 显著家族表 ----------------

def load_cafe_significant_nohighfail(path: Path, alpha: float) -> Dict[str, Set[str]]:
    """
    读取 cafe_significant_families_no_highfail.tsv，按 model 汇总：
      返回：model -> { sig_families }
    判定显著：Q <= alpha（忽略原 significant 字段，避免历史错误）
    """
    lines = read_lines(path)
    if not lines:
        raise RuntimeError(f"空文件：{path}")

    header = lines[0].split("\t")
    try:
        i_model = header.index("model")
        i_family = header.index("family")
        i_q = header.index("Q")
    except ValueError as e:
        raise RuntimeError(f"表头缺失必要列（model, family, Q）：{header}") from e

    model2sig: Dict[str, Set[str]] = {}
    n_all = 0
    n_sig = 0
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) <= max(i_model, i_family, i_q):
            continue
        model = parts[i_model]
        fam = parts[i_family]
        try:
            q = float(parts[i_q])
        except Exception:
            continue
        n_all += 1
        if q <= alpha:
            model2sig.setdefault(model, set()).add(fam)
            n_sig += 1
    logging.info(f"[cafe] 总 family 条目数：{n_all}，其中 Q<=α 的显著 family 数：{n_sig}")
    return model2sig

def load_cafe_family_universe(path: Path) -> Dict[str, Set[str]]:
    """
    读取 cafe_family_universe.tsv，返回 model -> {family}
    若找不到该文件，则退化为使用 cafe_significant_families_no_highfail.tsv 的全部 family。
    """
    if not path.exists():
        return {}
    lines = read_lines(path)
    if not lines:
        return {}
    header = lines[0].split("\t")
    try:
        i_model = header.index("model")
        i_family = header.index("family")
    except ValueError:
        return {}

    model2all: Dict[str, Set[str]] = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) <= max(i_model, i_family):
            continue
        model = parts[i_model]
        fam = parts[i_family]
        model2all.setdefault(model, set()).add(fam)
    return model2all

# ---------------- 读取 all_pep2cds_resolved.tsv ----------------

def load_all_members(path: Path, alias_map: Dict[str, str]) -> Dict[tuple, Set[str]]:
    """
    读取 all_pep2cds_resolved.tsv：
      OG  Species  protein_id  cds_id

    返回：
      members[(OG, species_norm)] = { cds_id1, cds_id2, ... }
    其中 species_norm 通过 alias_map 映射为规范名（若无映射则原样保留）。
    """
    lines = read_lines(path)
    if not lines:
        raise RuntimeError(f"空文件：{path}")
    header = lines[0].split("\t")
    if len(header) < 4 or header[0] != "OG":
        raise RuntimeError(f"all_members_tsv 表头不符合预期（首列应为 OG）：{header}")

    members: Dict[tuple, Set[str]] = {}
    n = 0
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 4:
            continue
        og = parts[0]
        sp = parts[1]
        cds = parts[3]
        sp_norm = alias_map.get(sp, sp)
        key = (og, sp_norm)
        members.setdefault(key, set()).add(cds)
        n += 1
    logging.info(f"[members] 从 {path} 读取到 {n} 条 (OG,Species,cds) 映射，键数：{len(members)}")
    return members

# ---------------- 读取注释字典（可选） ----------------

def load_annotation_universe(gene2go: Path | None, gene2path: Path | None) -> Set[str]:
    """
    读取 gene2go / gene2pathway，返回至少在任一字典中出现过的 gene_id 集合。
    若路径为 None 或不存在，则返回空集合。
    """
    ann_genes: Set[str] = set()

    def _load_two_col(path: Path, value_idx: int = 1) -> None:
        if not path or not path.exists():
            return
        lines = read_lines(path)
        if not lines:
            return
        # 跳首行表头
        for ln in lines[1:]:
            if not ln.strip():
                continue
            parts = ln.split("\t")
            if len(parts) <= value_idx:
                continue
            g = parts[value_idx]
            if g:
                ann_genes.add(g)

    if gene2go is not None:
        _load_two_col(gene2go, 1)
    if gene2path is not None:
        _load_two_col(gene2path, 1)

    logging.info(f"[annot] 注释字典中共有 {len(ann_genes)} 个带注释的基因")
    return ann_genes

# ---------------- 主流程 ----------------

def main() -> None:
    cfg = load_yaml(Path(CONFIG_PATH))
    paths_cfg = cfg.get("paths", {})
    cafe_agg_dir = Path(paths_cfg.get("cafe_agg_dir", "results/07_cafeagg")).resolve()
    logs_dir = Path(paths_cfg.get("logs_dir", "logs")).resolve()
    setup_logging(logs_dir)

    mkdir_p(cafe_agg_dir)
    logging.info(f"[init] cafe_agg_dir = {cafe_agg_dir}")

    inputs_cfg = cfg.get("inputs", {})
    all_members_path = inputs_cfg.get("all_members_tsv")
    if not all_members_path:
        logging.error("[ERR] config.inputs.all_members_tsv 未配置（应指向 all_pep2cds_resolved.tsv）")
        sys.exit(2)
    all_members_path = Path(all_members_path.replace("<publish_dir>", cfg.get("publish_dir", ""))).resolve()
    if not all_members_path.exists():
        logging.error(f"[ERR] all_members_tsv 不存在：{all_members_path}")
        sys.exit(2)

    species_cfg = cfg.get("species", {})
    alias_map: Dict[str, str] = species_cfg.get("alias_map", {}) or {}

    report_cfg = cfg.get("report", {})
    alpha = float(report_cfg.get("fdr_alpha", 0.05))
    logging.info(f"[init] FDR α = {alpha}")

    # enrich_bridge 配置
    cafe5_cfg = cfg.get("cafe5", {})
    bridge_cfg = cafe5_cfg.get("enrich_bridge", {}) or {}

    if not bridge_cfg.get("enable", False):
        logging.info("[init] cafe5.enrich_bridge.enable = false，本步不执行任何操作，直接退出。")
        return

    outputs_dir = Path(bridge_cfg.get("outputs_dir", cafe_agg_dir / "enrich_inputs")).resolve()
    mkdir_p(outputs_dir)
    logging.info(f"[init] enrich_outputs_dir = {outputs_dir}")

    species_sets_cfg = bridge_cfg.get("species_sets", {}) or {}
    if not species_sets_cfg:
        logging.error("[ERR] cafe5.enrich_bridge.species_sets 为空，至少需要定义一个 <tag> → [species] 映射")
        sys.exit(2)

    ann_cfg = bridge_cfg.get("annotations", {}) or {}
    gene2go_path = ann_cfg.get("gene2go")
    gene2path_path = ann_cfg.get("gene2pathway")

    gene2go = Path(gene2go_path).resolve() if gene2go_path else None
    gene2pathway = Path(gene2path_path).resolve() if gene2path_path else None

    # 1) 读取 CAFE 显著家族表
    sig_nofail_path = cafe_agg_dir / "cafe_significant_families_no_highfail.tsv"
    if not sig_nofail_path.exists():
        logging.error(f"[ERR] 缺少 13 步产物：{sig_nofail_path}")
        sys.exit(2)

    model2sig = load_cafe_significant_nohighfail(sig_nofail_path, alpha)

    # 2) 读取 family 宇宙（若有）
    fam_univ_path = cafe_agg_dir / "cafe_family_universe.tsv"
    model2all = load_cafe_family_universe(fam_univ_path)
    if not model2all:
        # 若缺少 family_universe.tsv，则退化使用 no_highfail 表中的全部 family
        logging.warning(f"[WARN] 未能加载 family 宇宙表 {fam_univ_path}，退化使用 no_highfail 表中的 family 集合作为背景 OG")
        model2all = {}
        lines = read_lines(sig_nofail_path)
        header = lines[0].split("\t")
        i_model = header.index("model")
        i_family = header.index("family")
        for ln in lines[1:]:
            if not ln.strip():
                continue
            parts = ln.split("\t")
            if len(parts) <= max(i_model, i_family):
                continue
            model = parts[i_model]
            fam = parts[i_family]
            model2all.setdefault(model, set()).add(fam)

    # 3) 读取 all_members_tsv
    members = load_all_members(all_members_path, alias_map)

    # 4) 注释宇宙（可选）
    annot_universe: Set[str] = set()
    if gene2go or gene2pathway:
        annot_universe = load_annotation_universe(gene2go, gene2pathway)

    # 5) 构建 enrich_inputs/README.txt（若不存在）
    readme_path = outputs_dir / "README.txt"
    if not readme_path.exists():
        with readme_path.open("w", encoding="utf-8") as f:
            f.write(
                "CAFE5 enrichment bridge\n"
                "=======================\n\n"
                "本目录由 13b_cafe_build_enrich_inputs.py 自动生成，用于把基因家族扩缩结果\n"
                "桥接到任意 rnaseq 富集模块。\n\n"
                "使用方法（示例）：\n"
                "  1) 在 aphylo 中运行 11/12/13/13b 步，得到本目录及子目录 <tag>/...\n"
                "  2) 选择一个 <tag>（例如 cafe_sc），将整个子目录拷贝到某个 RNA 工程的：\n"
                "       results/08_enrich/inputs/<tag>/\n"
                "  3) 在该 RNA 工程运行 08_enrich，对 up/down.list 和 background.list 进行 GO/KEGG 富集。\n\n"
                "当前版本说明：\n"
                "  - up.list / down.list 都表示 CAFE5 FDR 显著的基因家族成员基因集合\n"
                "    （尚未区分 expansion / contraction）；\n"
                "  - background.list 表示该物种集合在 CAFE 中成功评估的所有 family 的成员基因；\n"
                "  - meta.tsv 表头与 vNext 外部 tag 规范兼容：\n"
                "      label  n_detectable  n_annot_mapped  universe_size  detectable_rule  samples_used\n"
            )
        logging.info(f"[OK] 写出 README：{readme_path}")

    # 6) 针对每个 <tag> 生成 up/down/background/meta
    #    这里默认使用第一个 model（通常为 'global'）
    cafe_models = sorted(model2all.keys())
    if not cafe_models:
        logging.error("[ERR] 在 CAFE family 宇宙中没有找到任何 model，请检查 13 步产物。")
        sys.exit(2)
    model_used = cafe_models[0]
    logging.info(f"[init] enrich_inputs 将使用 model='{model_used}' 的 family 集合作为 OG 宇宙")

    fam_all = model2all[model_used]
    fam_sig = model2sig.get(model_used, set())
    logging.info(f"[cafe] model={model_used}: 背景 OG 数={len(fam_all)}，显著 OG 数={len(fam_sig)}")

    for tag, sp_list in species_sets_cfg.items():
        species_set = [alias_map.get(sp, sp) for sp in sp_list]
        logging.info(f"[tag={tag}] 物种集合：{species_set}")

        # 背景基因：所有背景 OG × species_set 的成员基因并集
        bg_genes: Set[str] = set()
        for og in fam_all:
            for sp in species_set:
                key = (og, sp)
                if key in members:
                    bg_genes.update(members[key])

        # 显著基因：显著 OG × species_set 的成员基因并集
        sig_genes: Set[str] = set()
        for og in fam_sig:
            for sp in species_set:
                key = (og, sp)
                if key in members:
                    sig_genes.update(members[key])

        logging.info(
            f"[tag={tag}] 背景 OG={len(fam_all)}, 背景基因={len(bg_genes)}, "
            f"显著 OG={len(fam_sig)}, 显著基因={len(sig_genes)}"
        )

        tag_dir = outputs_dir / tag
        mkdir_p(tag_dir)

        # 写 up.list / down.list / background.list
        up_path = tag_dir / "up.list"
        down_path = tag_dir / "down.list"
        bg_path = tag_dir / "background.list"

        with up_path.open("w", encoding="utf-8") as f:
            for g in sorted(sig_genes):
                f.write(g + "\n")

        # 当前版本 down.list 与 up.list 相同（代表“快速变化家族”全集）
        with down_path.open("w", encoding="utf-8") as f:
            for g in sorted(sig_genes):
                f.write(g + "\n")

        with bg_path.open("w", encoding="utf-8") as f:
            for g in sorted(bg_genes):
                f.write(g + "\n")

        # 构建 meta.tsv
        n_detectable = len(bg_genes)
        universe_size = n_detectable  # 外部 tag 背景 = CAFE 实际可评估宇宙
        if annot_universe:
            n_annot_mapped = len(bg_genes & annot_universe)
        else:
            n_annot_mapped = "NA"

        detectable_rule = "CAFE5_evaluable_gene_members_of_tested_families_any_direction"
        samples_used = ";".join(species_set)

        meta_path = tag_dir / "meta.tsv"
        write_tsv(
            meta_path,
            ["label", "n_detectable", "n_annot_mapped", "universe_size", "detectable_rule", "samples_used"],
            [[tag, n_detectable, n_annot_mapped, universe_size, detectable_rule, samples_used]],
        )

        logging.info(f"[tag={tag}] 写出 up/down/background/meta 到：{tag_dir}")

    logging.info("========== APhylo 13b — 完成 ==========")

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)