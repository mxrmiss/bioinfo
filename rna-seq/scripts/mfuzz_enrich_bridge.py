#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mfuzz_enrich_bridge.py

功能：
1. 读取 mfuzz_run.R 的输出结果
2. 选择 all / core / strict 三类 cluster 基因集之一作为富集输入
3. 复用现有 RNA-seq 流水线的 07 注释资源（gene2go / gene2pathway / GO / KEGG 词典）
4. 通过临时 runtime config 调用现有 08_g_enrich.R，避免重跑整个 08a 流程
5. 输出每个 cluster 的 GO / KEGG 富集结果，并生成汇总表

设计原则：
- 不修改原始 results/08_enrich
- 不重跑 RNA 对比富集
- 不重复实现一套新的 GO/KEGG 富集统计逻辑
- 尽量复用现有 07/08 资源和 08_g_enrich.R

运行方式：
    python3 scripts/mfuzz_enrich_bridge.py
"""

from __future__ import annotations

import csv
import os
import re
import sys
import shutil
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Any

try:
    import yaml
except Exception as e:
    print("[ERR] 需要 PyYAML，请先安装：pip install pyyaml", file=sys.stderr)
    raise e


###############################################################################
# 一、用户参数区（请按需修改）
###############################################################################

PROJECT_ROOT = "."
# 中文说明：项目根目录，通常保持为当前目录即可

MFUZZ_DIR = "results/11_mfuzz"
# 中文说明：mfuzz_run.R 的主输出目录

OUTDIR = "results/12_mfuzz_enrich"
# 中文说明：Mfuzz 富集桥接输出目录

GENESET_TYPE = "core"
# 中文说明：使用哪一类 cluster 基因集做富集
# 可选：all / core / strict
# - all：每簇全部基因
# - core：membership >= MEMBERSHIP_CORE
# - strict：membership >= MEMBERSHIP_STRICT

MIN_GENES_PER_CLUSTER = 10
# 中文说明：进入富集分析的最小基因数阈值
# 若某个 cluster 的基因数少于该值，则跳过该 cluster

RUN_ENRICH = True
# 中文说明：是否在准备好任务后，直接调用现有 08_g_enrich.R 执行富集
# True：准备任务后自动执行
# False：只生成任务与输入文件，不实际执行富集

RSCRIPT_BIN = None
# 中文说明：Rscript 路径
# 若设为 None，则优先从原 config.yaml 中读取 binaries.Rscript；若没有，则默认使用 "Rscript"

ORIGINAL_CONFIG = "config.yaml"
# 中文说明：原项目 config.yaml 路径

R_ENRICH_SCRIPT = "scripts/08_g_enrich.R"
# 中文说明：现有 RNA-seq 流水线中的 08_g_enrich.R 路径

GENEID_COL = "gene_id"
# 中文说明：单列基因列表文件中的列名

TASK_TEST_SET = "test"
# 中文说明：传给 08_g_enrich.R 的 test_set 名称
# 对 Mfuzz cluster 富集，统一使用 test

TASK_TYPE = "external"
# 中文说明：任务类型，沿用现有 08 流程中的 external 通道

LABEL_PREFIX = "mfuzz"
# 中文说明：输出 label 的前缀

USE_ABSOLUTE_PATHS_IN_RUNTIME_CONFIG = True
# 中文说明：临时 runtime config 中是否写入绝对路径
# 建议保持 True，更稳

TOP_N_TERMS_FOR_SUMMARY = 5
# 中文说明：汇总报告中，每个 cluster 保留前多少个 GO/KEGG term

VERBOSE = True
# 中文说明：是否输出详细日志


###############################################################################
# 二、基础函数区
###############################################################################

def log_msg(level: str, msg: str) -> None:
    if VERBOSE:
        ts = Path
        print(f"[{level}] {msg}")


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"找不到配置文件：{path}")
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def write_yaml(obj: Dict[str, Any], path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(obj, f, sort_keys=False, allow_unicode=True)


def read_single_column_list(path: Path, header_name: str = "gene_id") -> List[str]:
    """
    读取单列表文件：
    - 若有表头且包含 gene_id，则读取该列
    - 否则默认读取第一列
    """
    if not path.exists():
        return []

    genes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return []

        has_header = False
        idx = 0
        if header_name in first:
            has_header = True
            idx = first.index(header_name)
        elif len(first) == 1 and first[0].strip() == header_name:
            has_header = True
            idx = 0

        if not has_header:
            if len(first) > 0 and first[0].strip():
                genes.append(first[0].strip())

        for row in reader:
            if not row:
                continue
            if idx >= len(row):
                continue
            g = row[idx].strip()
            if g:
                genes.append(g)

    return sorted(set(genes))


def write_single_column_list(path: Path, items: List[str], header: str = "gene_id") -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{header}\n")
        for x in items:
            f.write(f"{x}\n")


def write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv_dict(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def get_nested(d: Dict[str, Any], keys: List[str], default=None):
    cur: Any = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def read_annotation_gene_universe(annot_dir: Path) -> Set[str]:
    """
    从 gene2go.tsv 和 gene2pathway.tsv 中读取注释基因宇宙
    """
    genes: Set[str] = set()

    for fp in [annot_dir / "gene2go.tsv", annot_dir / "gene2pathway.tsv"]:
        if not fp.exists():
            continue
        with fp.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gid = (row.get("gene_id") or "").strip()
                if gid:
                    genes.add(gid)
    return genes


def detect_cluster_files(cluster_dir: Path, geneset_type: str) -> List[Tuple[int, Path]]:
    """
    自动寻找 cluster 文件，并按 cluster 编号排序
    例如：
      cluster01_core.tsv
      cluster02_core.tsv
    """
    suffix_map = {
        "all": "_all.tsv",
        "core": "_core.tsv",
        "strict": "_strict.tsv",
    }
    if geneset_type not in suffix_map:
        raise ValueError("GENESET_TYPE 仅支持 all / core / strict")

    suffix = suffix_map[geneset_type]
    pat = re.compile(r"cluster(\d+)" + re.escape(suffix) + r"$")

    hits: List[Tuple[int, Path]] = []
    for fp in sorted(cluster_dir.glob(f"cluster*{suffix}")):
        m = pat.search(fp.name)
        if not m:
            continue
        cid = int(m.group(1))
        hits.append((cid, fp))

    hits.sort(key=lambda x: x[0])
    return hits


def summarize_cluster_name(cid: int, geneset_type: str) -> str:
    return f"{LABEL_PREFIX}_{geneset_type}_cluster{cid:02d}"


def safe_abs(path_str: str, project_root: Path) -> str:
    p = Path(path_str)
    if p.is_absolute():
        return str(p)
    return str((project_root / p).resolve())


def build_runtime_config(
    original_cfg: Dict[str, Any],
    runtime_enrich_dir: Path,
    annot_dir: Path,
    maps_dir: Path,
    project_root: Path,
) -> Dict[str, Any]:
    """
    生成给 08_g_enrich.R 使用的临时 config.yaml
    核心思路：
    - dirs.enrich 指向本次 Mfuzz 富集运行目录
    - dirs.annot / dirs.maps 复用原项目 07 资源
    - GSEA 强制关闭
    - 其他 ORA 相关参数尽量沿用原 config
    """
    cfg_new: Dict[str, Any] = {}

    # 目录
    if USE_ABSOLUTE_PATHS_IN_RUNTIME_CONFIG:
        cfg_new["dirs"] = {
            "enrich": str(runtime_enrich_dir.resolve()),
            "annot": str(annot_dir.resolve()),
            "maps": str(maps_dir.resolve()),
        }
    else:
        cfg_new["dirs"] = {
            "enrich": os.path.relpath(runtime_enrich_dir.resolve(), project_root.resolve()),
            "annot": os.path.relpath(annot_dir.resolve(), project_root.resolve()),
            "maps": os.path.relpath(maps_dir.resolve(), project_root.resolve()),
        }

    # 复制原 config 中相关参数
    cfg_new["go"] = original_cfg.get("go", {})
    cfg_new["kegg"] = original_cfg.get("kegg", {})
    cfg_new["enrich"] = original_cfg.get("enrich", {})

    # GSEA 对 Mfuzz cluster 不启用
    cfg_new["gsea"] = original_cfg.get("gsea", {}).copy() if isinstance(original_cfg.get("gsea", {}), dict) else {}
    cfg_new["gsea"]["enable"] = False

    return cfg_new


def parse_float(x: str, default: float = 1.0) -> float:
    try:
        return float(x)
    except Exception:
        return default


###############################################################################
# 三、主流程
###############################################################################

def main() -> None:
    project_root = Path(PROJECT_ROOT).resolve()
    mfuzz_dir = (project_root / MFUZZ_DIR).resolve()
    outdir = (project_root / OUTDIR).resolve()
    original_config_fp = (project_root / ORIGINAL_CONFIG).resolve()
    r_enrich_script_fp = (project_root / R_ENRICH_SCRIPT).resolve()

    if not mfuzz_dir.exists():
        raise FileNotFoundError(f"找不到 MFUZZ_DIR：{mfuzz_dir}")
    if not original_config_fp.exists():
        raise FileNotFoundError(f"找不到 ORIGINAL_CONFIG：{original_config_fp}")
    if not r_enrich_script_fp.exists():
        raise FileNotFoundError(f"找不到 R_ENRICH_SCRIPT：{r_enrich_script_fp}")

    ensure_dir(outdir)

    inputs_dir = outdir / "inputs"
    runtime_dir = outdir / "_runtime_08g"
    report_dir = outdir / "06_report"
    ensure_dir(inputs_dir)
    ensure_dir(runtime_dir)
    ensure_dir(report_dir)

    log_msg("INFO", f"读取原始配置：{original_config_fp}")
    original_cfg = read_yaml(original_config_fp)

    # 从原 config 推断注释资源路径
    annot_dir_str = get_nested(original_cfg, ["dirs", "annot"], "results/07_annot")
    maps_dir_str = get_nested(original_cfg, ["dirs", "maps"], "results/03_maps")

    annot_dir = Path(safe_abs(annot_dir_str, project_root))
    maps_dir = Path(safe_abs(maps_dir_str, project_root))

    if not annot_dir.exists():
        raise FileNotFoundError(f"找不到注释目录 dirs.annot：{annot_dir}")
    if not maps_dir.exists():
        log_msg("WARNING", f"maps 目录不存在：{maps_dir}；gene_info.tsv 将可能不可用")

    # Rscript 路径
    rscript_bin = RSCRIPT_BIN
    if rscript_bin is None:
        rscript_bin = get_nested(original_cfg, ["binaries", "Rscript"], "Rscript")

    # 读取 mfuzz universe
    mfuzz_universe_fp = mfuzz_dir / "10_for_enrichment" / "mfuzz_universe.tsv"
    if not mfuzz_universe_fp.exists():
        raise FileNotFoundError(f"找不到 mfuzz universe：{mfuzz_universe_fp}")

    universe_genes = read_single_column_list(mfuzz_universe_fp, header_name=GENEID_COL)
    universe_set = set(universe_genes)
    if len(universe_set) == 0:
        raise RuntimeError("mfuzz_universe.tsv 为空，无法继续")

    log_msg("INFO", f"Mfuzz universe 基因数：{len(universe_set)}")

    # 读取注释宇宙，用于 coverage 统计
    annot_gene_universe = read_annotation_gene_universe(annot_dir)
    log_msg("INFO", f"注释基因宇宙大小：{len(annot_gene_universe)}")

    # 读取 cluster gene lists
    cluster_lists_dir = mfuzz_dir / "08_cluster_gene_lists"
    if not cluster_lists_dir.exists():
        raise FileNotFoundError(f"找不到 cluster gene list 目录：{cluster_lists_dir}")

    cluster_files = detect_cluster_files(cluster_lists_dir, GENESET_TYPE)
    if len(cluster_files) == 0:
        raise RuntimeError(f"在 {cluster_lists_dir} 中未找到 GENESET_TYPE={GENESET_TYPE} 的 cluster 文件")

    log_msg("INFO", f"检测到 {len(cluster_files)} 个 cluster 基因集文件，类型={GENESET_TYPE}")

    # 先读取 cluster_summary，后面做摘要用
    cluster_summary_fp = mfuzz_dir / "07_final_clusters" / "cluster_summary.tsv"
    cluster_summary_rows = read_tsv_dict(cluster_summary_fp)
    cluster_summary_map: Dict[int, Dict[str, str]] = {}
    for row in cluster_summary_rows:
        try:
            cid = int(row.get("cluster", ""))
            cluster_summary_map[cid] = row
        except Exception:
            continue

    task_rows: List[Dict[str, Any]] = []
    coverage_rows: List[Dict[str, Any]] = []
    task_input_summary_rows: List[Dict[str, Any]] = []

    # 为每个 cluster 生成独立输入
    for cid, fp in cluster_files:
        label = summarize_cluster_name(cid, GENESET_TYPE)
        cluster_input_dir = inputs_dir / label
        cluster_outdir = outdir / label
        ensure_dir(cluster_input_dir)
        ensure_dir(cluster_outdir)

        genes_raw = read_single_column_list(fp, header_name=GENEID_COL)
        genes_raw_set = set(genes_raw)

        # 与 universe 对齐
        genes_in_universe = sorted(genes_raw_set & universe_set)

        # 注释覆盖统计
        n_with_go_or_kegg = len(set(genes_in_universe) & annot_gene_universe)

        task_input_summary_rows.append({
            "cluster": cid,
            "label": label,
            "geneset_type": GENESET_TYPE,
            "n_raw": len(genes_raw_set),
            "n_in_universe": len(genes_in_universe),
            "n_annotated": n_with_go_or_kegg,
            "source_file": str(fp),
        })

        # 太小的簇直接跳过
        if len(genes_in_universe) < MIN_GENES_PER_CLUSTER:
            log_msg(
                "WARNING",
                f"{label} 基因数过少（n_in_universe={len(genes_in_universe)} < {MIN_GENES_PER_CLUSTER}），跳过富集"
            )
            coverage_rows.append({
                "cluster": cid,
                "label": label,
                "n_input": len(genes_raw_set),
                "n_in_universe": len(genes_in_universe),
                "n_with_go_or_kegg": n_with_go_or_kegg,
                "status": "skipped_too_small",
            })
            continue

        # 写 test.list 和 background.list
        test_fp = cluster_input_dir / "test.list"
        bg_fp = cluster_input_dir / "background.list"
        meta_fp = cluster_input_dir / "meta.tsv"

        write_single_column_list(test_fp, genes_in_universe, header=GENEID_COL)
        write_single_column_list(bg_fp, sorted(universe_set), header=GENEID_COL)

        mapped_bg = len(universe_set & annot_gene_universe)
        coverage = (mapped_bg / len(universe_set)) if len(universe_set) > 0 else 0.0

        meta_rows = [{
            "label": label,
            "n_detectable": len(universe_set),
            "n_annot_mapped": mapped_bg,
            "universe_size": mapped_bg,
            "coverage": f"{coverage:.4f}",
            "detectable_rule": "mfuzz_filtered_universe",
            "samples_used": "NA",
        }]
        write_tsv(
            meta_fp,
            meta_rows,
            fieldnames=[
                "label",
                "n_detectable",
                "n_annot_mapped",
                "universe_size",
                "coverage",
                "detectable_rule",
                "samples_used",
            ],
        )

        task_rows.append({
            "label": label,
            "type": TASK_TYPE,
            "test_set": TASK_TEST_SET,
            "test_file": str(test_fp),
            "background_file": str(bg_fp),
            "universe_meta": str(meta_fp),
            "outdir": str(cluster_outdir),
            "n_deg_all": "NA",
        })

        coverage_rows.append({
            "cluster": cid,
            "label": label,
            "n_input": len(genes_raw_set),
            "n_in_universe": len(genes_in_universe),
            "n_with_go_or_kegg": n_with_go_or_kegg,
            "status": "ready",
        })

    if len(task_rows) == 0:
        raise RuntimeError("没有任何 cluster 通过最小基因数过滤，无法执行富集")

    # 写基础汇总
    write_tsv(
        outdir / "00_cluster_input_summary.tsv",
        task_input_summary_rows,
        fieldnames=["cluster", "label", "geneset_type", "n_raw", "n_in_universe", "n_annotated", "source_file"],
    )
    write_tsv(
        outdir / "01_annotation_coverage.tsv",
        coverage_rows,
        fieldnames=["cluster", "label", "n_input", "n_in_universe", "n_with_go_or_kegg", "status"],
    )

    # 写 tasks.tsv
    tasks_fp = outdir / "tasks.tsv"
    write_tsv(
        tasks_fp,
        task_rows,
        fieldnames=[
            "label",
            "type",
            "test_set",
            "test_file",
            "background_file",
            "universe_meta",
            "outdir",
            "n_deg_all",
        ],
    )
    log_msg("INFO", f"已写出 Mfuzz 富集 tasks.tsv：{tasks_fp}")

    # 生成 runtime config
    runtime_cfg = build_runtime_config(
        original_cfg=original_cfg,
        runtime_enrich_dir=outdir,
        annot_dir=annot_dir,
        maps_dir=maps_dir,
        project_root=project_root,
    )
    runtime_cfg_fp = runtime_dir / "config.yaml"
    write_yaml(runtime_cfg, runtime_cfg_fp)
    log_msg("INFO", f"已写出 runtime config：{runtime_cfg_fp}")

    # 记录运行参数
    run_params_rows = [{
        "parameter": "PROJECT_ROOT", "value": str(project_root)
    }, {
        "parameter": "MFUZZ_DIR", "value": str(mfuzz_dir)
    }, {
        "parameter": "OUTDIR", "value": str(outdir)
    }, {
        "parameter": "GENESET_TYPE", "value": str(GENESET_TYPE)
    }, {
        "parameter": "MIN_GENES_PER_CLUSTER", "value": str(MIN_GENES_PER_CLUSTER)
    }, {
        "parameter": "RUN_ENRICH", "value": str(RUN_ENRICH)
    }, {
        "parameter": "ORIGINAL_CONFIG", "value": str(original_config_fp)
    }, {
        "parameter": "R_ENRICH_SCRIPT", "value": str(r_enrich_script_fp)
    }, {
        "parameter": "RSCRIPT_BIN", "value": str(rscript_bin)
    }]
    write_tsv(outdir / "run_params.used.tsv", run_params_rows, fieldnames=["parameter", "value"])

    # 执行富集
    if RUN_ENRICH:
        log_msg("INFO", "开始调用现有 08_g_enrich.R 执行 Mfuzz cluster 富集")

        cmd = [str(rscript_bin), str(r_enrich_script_fp)]
        try:
            subprocess.run(
                cmd,
                cwd=str(runtime_dir),
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"08_g_enrich.R 执行失败，返回码={e.returncode}") from e

        log_msg("INFO", "08_g_enrich.R 执行完成")
    else:
        log_msg("INFO", "RUN_ENRICH=False，本次仅准备输入文件，不执行富集")

    # 生成汇总报告
    enrich_fdr = parse_float(str(get_nested(original_cfg, ["enrich", "fdr"], 0.05)), 0.05)

    report_rows: List[Dict[str, Any]] = []
    all_go_rows: List[Dict[str, Any]] = []
    all_kegg_rows: List[Dict[str, Any]] = []

    for task in task_rows:
        label = task["label"]
        out_label_dir = Path(task["outdir"])
        cid_match = re.search(r"cluster(\d+)", label)
        cid = int(cid_match.group(1)) if cid_match else None

        go_files = [
            out_label_dir / "GO_BP_by_term_test.tsv",
            out_label_dir / "GO_CC_by_term_test.tsv",
            out_label_dir / "GO_MF_by_term_test.tsv",
        ]
        kegg_file = out_label_dir / "KEGG_by_term_test.tsv"

        go_sig_terms: List[str] = []
        kegg_sig_terms: List[str] = []

        # 读取 GO
        for gfp in go_files:
            if not gfp.exists():
                continue
            rows = read_tsv_dict(gfp)
            for row in rows:
                row["_source_label"] = label
                all_go_rows.append(row)

            sig_rows = []
            for row in rows:
                try:
                    padj = float(row.get("p_adjust", "1"))
                except Exception:
                    padj = 1.0
                if padj <= enrich_fdr:
                    sig_rows.append(row)

            sig_rows.sort(key=lambda r: (
                float(r.get("p_adjust", "1")) if str(r.get("p_adjust", "")).strip() not in {"", "NA"} else 1.0,
                -float(r.get("gene_count", "0")) if str(r.get("gene_count", "")).strip() not in {"", "NA"} else 0.0
            ))

            for row in sig_rows[:TOP_N_TERMS_FOR_SUMMARY]:
                term_name = row.get("term_name", "")
                ontology = row.get("ontology", "")
                if term_name:
                    go_sig_terms.append(f"{ontology}:{term_name}")

        # 读取 KEGG
        if kegg_file.exists():
            rows = read_tsv_dict(kegg_file)
            for row in rows:
                row["_source_label"] = label
                all_kegg_rows.append(row)

            sig_rows = []
            for row in rows:
                try:
                    padj = float(row.get("p_adjust", "1"))
                except Exception:
                    padj = 1.0
                if padj <= enrich_fdr:
                    sig_rows.append(row)

            sig_rows.sort(key=lambda r: (
                float(r.get("p_adjust", "1")) if str(r.get("p_adjust", "")).strip() not in {"", "NA"} else 1.0,
                -float(r.get("gene_count", "0")) if str(r.get("gene_count", "")).strip() not in {"", "NA"} else 0.0
            ))

            for row in sig_rows[:TOP_N_TERMS_FOR_SUMMARY]:
                term_name = row.get("term_name", "")
                if term_name:
                    kegg_sig_terms.append(term_name)

        peak_stage = ""
        n_core_genes = ""
        if cid is not None and cid in cluster_summary_map:
            peak_stage = cluster_summary_map[cid].get("stage_peak", "")
            n_core_genes = cluster_summary_map[cid].get("n_core_ge_0.7", "")

        report_rows.append({
            "cluster": cid if cid is not None else "NA",
            "label": label,
            "geneset_type": GENESET_TYPE,
            "peak_stage": peak_stage,
            "n_core_genes": n_core_genes,
            "top_go_terms": " | ".join(go_sig_terms) if go_sig_terms else "NA",
            "top_kegg_terms": " | ".join(kegg_sig_terms) if kegg_sig_terms else "NA",
        })

    write_tsv(
        report_dir / "cluster_function_summary.tsv",
        report_rows,
        fieldnames=["cluster", "label", "geneset_type", "peak_stage", "n_core_genes", "top_go_terms", "top_kegg_terms"],
    )

    if len(all_go_rows) > 0:
        go_fields = list(all_go_rows[0].keys())
        write_tsv(outdir / "05_summary" / "all_clusters_go.tsv", all_go_rows, fieldnames=go_fields)
    if len(all_kegg_rows) > 0:
        kegg_fields = list(all_kegg_rows[0].keys())
        write_tsv(outdir / "05_summary" / "all_clusters_kegg.tsv", all_kegg_rows, fieldnames=kegg_fields)

    log_msg("INFO", f"Mfuzz 富集桥接完成，主输出目录：{outdir}")


if __name__ == "__main__":
    main()
