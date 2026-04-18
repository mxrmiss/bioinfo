#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import re
import sys
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    import yaml
except Exception as e:
    print("[ERR] 缺少 PyYAML，请先安装：pip install pyyaml", file=sys.stderr)
    raise e

###############################################################################
# 用户参数区
###############################################################################

CONFIG_YAML = "config.yaml"
# 中文说明：统一配置文件路径。脚本其余参数尽量从 config.yaml 读取，不重复造轮子。

VERBOSE = True
# 中文说明：是否打印详细日志。

###############################################################################
# 基础函数区
###############################################################################


def log_msg(level: str, msg: str) -> None:
    if VERBOSE:
        print(f"[{level}] {msg}")



def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)



def read_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"找不到配置文件：{path}")
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}



def write_yaml(obj: Dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(obj, f, sort_keys=False, allow_unicode=True)



def get_nested(d: Dict[str, Any], keys: List[str], default=None):
    cur: Any = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur



def to_bool(x: Any, default: bool = False) -> bool:
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in {"1", "true", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "no", "n", "off"}:
        return False
    return default



def to_int(x: Any, default: int) -> int:
    try:
        return int(x)
    except Exception:
        return default



def to_float(x: Any, default: float) -> float:
    try:
        return float(x)
    except Exception:
        return default



def as_str_list(x: Any) -> List[str]:
    if x is None:
        return []
    if isinstance(x, list):
        return [str(i).strip() for i in x if str(i).strip()]
    s = str(x).strip()
    return [s] if s else []



def safe_abs(path_str: str, project_root: Path) -> Path:
    p = Path(path_str)
    return p if p.is_absolute() else (project_root / p).resolve()



def sanitize_label(x: str) -> str:
    s = str(x).strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_\-.]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")



def read_tsv_dict(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))



def write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)



def read_single_column_list(path: Path, header_name: str = "gene_id") -> List[str]:
    if not path.exists():
        return []

    genes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return []

        idx = 0
        has_header = False
        if header_name in first:
            idx = first.index(header_name)
            has_header = True
        elif len(first) == 1 and first[0].strip() == header_name:
            has_header = True

        if not has_header and len(first) > 0:
            g = first[0].strip()
            if g:
                genes.append(g)

        for row in reader:
            if not row or idx >= len(row):
                continue
            g = row[idx].strip()
            if g:
                genes.append(g)

    return list(dict.fromkeys(genes))



def write_single_column_list(path: Path, items: List[str], header: str = "gene_id") -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{header}\n")
        for x in items:
            f.write(f"{x}\n")



def read_assignment_universe(assign_fp: Path) -> Set[str]:
    rows = read_tsv_dict(assign_fp)
    genes: Set[str] = set()
    for row in rows:
        gid = (row.get("gene_id") or "").strip()
        if gid:
            genes.add(gid)
    return genes



def read_annotation_gene_universe(annot_dir: Path) -> Set[str]:
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



def detect_module_list_files(list_dir: Path) -> List[Tuple[str, Path]]:
    hits: List[Tuple[str, Path]] = []
    for fp in sorted(list_dir.glob("module_*.list")):
        m = re.match(r"module_(.+)\.list$", fp.name)
        if m:
            hits.append((m.group(1), fp))
    return hits



def read_hub_top_genes(hub_fp: Path, top_n: int) -> Dict[str, List[str]]:
    rows = read_tsv_dict(hub_fp)
    out: Dict[str, List[str]] = {}
    for row in rows:
        mod = (row.get("module_color") or "").strip()
        gid = (row.get("gene_id") or "").strip()
        if not mod or not gid:
            continue
        out.setdefault(mod, [])
        if len(out[mod]) < top_n:
            out[mod].append(gid)
    for mod in list(out.keys()):
        out[mod] = list(dict.fromkeys(out[mod]))
    return out



def read_module_trait_matrix(path: Path) -> Tuple[List[str], Dict[str, Dict[str, float]]]:
    if not path.exists():
        raise FileNotFoundError(f"文件不存在：{path}")

    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise RuntimeError(f"矩阵表头异常：{path}")

        traits = header[1:]
        out: Dict[str, Dict[str, float]] = {}
        for row in reader:
            if not row:
                continue
            raw_name = row[0].strip()
            if not raw_name:
                continue
            module_name = re.sub(r"^ME", "", raw_name)
            val_map: Dict[str, float] = {}
            for trait, val in zip(traits, row[1:]):
                try:
                    val_map[trait] = float(val)
                except Exception:
                    val_map[trait] = float("nan")
            out[module_name] = val_map
    return traits, out



def module_pass_trait_filter(
    module_name: str,
    target_trait: str,
    cor_map: Dict[str, Dict[str, float]],
    fdr_map: Dict[str, Dict[str, float]],
    trait_direction: str,
    trait_cor_min: float,
    trait_fdr_max: float,
) -> Tuple[bool, float, float, str]:
    if module_name not in cor_map:
        return False, float("nan"), float("nan"), "missing_module_in_cor"
    if module_name not in fdr_map:
        return False, float("nan"), float("nan"), "missing_module_in_fdr"
    if target_trait not in cor_map[module_name]:
        return False, float("nan"), float("nan"), "missing_trait_in_cor"
    if target_trait not in fdr_map[module_name]:
        return False, float("nan"), float("nan"), "missing_trait_in_fdr"

    r = cor_map[module_name][target_trait]
    fdr = fdr_map[module_name][target_trait]

    if r != r or fdr != fdr:
        return False, r, fdr, "nan_value"
    if fdr > trait_fdr_max:
        return False, r, fdr, "fdr_too_large"

    direction = trait_direction.lower().strip()
    if direction == "positive":
        return (r >= trait_cor_min), r, fdr, "pass" if r >= trait_cor_min else "positive_cor_too_small"
    if direction == "negative":
        return (r <= -trait_cor_min), r, fdr, "pass" if r <= -trait_cor_min else "negative_cor_too_small"
    return (abs(r) >= trait_cor_min), r, fdr, "pass" if abs(r) >= trait_cor_min else "abs_cor_too_small"



def choose_module_candidates(
    wgcna_dir: Path,
    geneset_mode: str,
    geneid_col: str,
    hub_top_n: int,
    modules_include: List[str],
    modules_exclude: List[str],
) -> List[Tuple[str, List[str]]]:
    include_set = set(modules_include) if modules_include else None
    exclude_set = set(modules_exclude)
    module_rows: List[Tuple[str, List[str]]] = []

    if geneset_mode == "module_lists":
        list_dir = wgcna_dir / "lists"
        if not list_dir.exists():
            raise FileNotFoundError(f"找不到模块列表目录：{list_dir}")
        for mod, fp in detect_module_list_files(list_dir):
            genes = read_single_column_list(fp, header_name=geneid_col)
            module_rows.append((mod, genes))
    elif geneset_mode == "hub_top":
        hub_fp = wgcna_dir / "hub" / "hub_genes_by_module.tsv"
        if not hub_fp.exists():
            raise FileNotFoundError(f"找不到 hub 文件：{hub_fp}")
        mod2genes = read_hub_top_genes(hub_fp, hub_top_n)
        for mod in sorted(mod2genes.keys()):
            module_rows.append((mod, mod2genes[mod]))
    else:
        raise ValueError("geneset_mode 仅支持 module_lists / hub_top")

    filtered: List[Tuple[str, List[str]]] = []
    for mod, genes in module_rows:
        if include_set is not None and mod not in include_set:
            continue
        if mod in exclude_set:
            continue
        filtered.append((mod, genes))
    return filtered



def build_runtime_config(
    original_cfg: Dict[str, Any],
    runtime_enrich_dir: Path,
    annot_dir: Path,
    maps_dir: Path,
) -> Dict[str, Any]:
    cfg_new = dict(original_cfg)
    cfg_new.setdefault("dirs", {})
    cfg_new["dirs"]["enrich"] = str(runtime_enrich_dir.resolve())
    cfg_new["dirs"]["annot"] = str(annot_dir.resolve())
    cfg_new["dirs"]["maps"] = str(maps_dir.resolve())

    # WGCNA 桥接只跑 ORA，不跑 GSEA
    cfg_new.setdefault("gsea", {})
    cfg_new["gsea"]["enable"] = False
    return cfg_new


###############################################################################
# 主流程
###############################################################################


def main() -> None:
    project_root = Path(".").resolve()
    config_fp = (project_root / CONFIG_YAML).resolve()
    cfg = read_yaml(config_fp)

    wgcna_cfg = cfg.get("wgcna", {})
    bridge_cfg = cfg.get("wgcna_enrich", {})

    dirs_cfg = cfg.get("dirs", {})
    binaries_cfg = cfg.get("binaries", {})

    wgcna_dir = safe_abs(str(dirs_cfg.get("wgcna", "results/10_wgcna")), project_root)
    annot_dir = safe_abs(str(dirs_cfg.get("annot", "results/07_annot")), project_root)
    maps_dir = safe_abs(str(dirs_cfg.get("maps", "results/03_maps")), project_root)
    outdir = safe_abs(str(bridge_cfg.get("outdir", "results/12_wgcna_enrich")), project_root)
    r_enrich_script = safe_abs(str(bridge_cfg.get("r_enrich_script", "scripts/08_g_enrich.R")), project_root)
    rscript_bin = str(binaries_cfg.get("Rscript", "Rscript"))

    geneset_mode = str(bridge_cfg.get("geneset_mode", "module_lists"))
    geneid_col = str(bridge_cfg.get("geneid_col", "gene_id"))
    min_genes_per_set = to_int(bridge_cfg.get("min_genes_per_set", 10), 10)
    background_mode = str(bridge_cfg.get("background_mode", "wgcna_universe"))
    run_enrich = to_bool(bridge_cfg.get("run_enrich", True), True)
    modules_include = as_str_list(bridge_cfg.get("modules_include", []))
    modules_exclude = as_str_list(bridge_cfg.get("modules_exclude", ["grey"])) or ["grey"]

    hub_top_n = to_int(
        bridge_cfg.get(
            "hub_top_n",
            get_nested(wgcna_cfg, ["export", "hub_genes", "top_n_per_module"], 100),
        ),
        100,
    )

    trait_filter_cfg = bridge_cfg.get("trait_filter", {}) if isinstance(bridge_cfg.get("trait_filter", {}), dict) else {}
    trait_filter_enable = to_bool(trait_filter_cfg.get("enable", False), False)
    target_trait = str(
        trait_filter_cfg.get(
            "target_trait",
            get_nested(wgcna_cfg, ["export", "hub_genes", "primary_trait"], ""),
        )
    ).strip()
    trait_direction = str(trait_filter_cfg.get("direction", "both"))
    trait_cor_min = to_float(trait_filter_cfg.get("cor_min", 0.30), 0.30)
    trait_fdr_max = to_float(trait_filter_cfg.get("fdr_max", get_nested(cfg, ["enrich", "fdr"], 0.05)), 0.05)

    if not wgcna_dir.exists():
        raise FileNotFoundError(f"找不到 dirs.wgcna：{wgcna_dir}")
    if not annot_dir.exists():
        raise FileNotFoundError(f"找不到 dirs.annot：{annot_dir}")
    if not r_enrich_script.exists():
        raise FileNotFoundError(f"找不到 08 富集脚本：{r_enrich_script}")

    ensure_dir(outdir)
    inputs_dir = outdir / "inputs"
    runtime_dir = outdir / "_runtime_08g"
    ensure_dir(inputs_dir)
    ensure_dir(runtime_dir)

    assign_fp = wgcna_dir / "modules" / "gene_module_assignments.tsv"
    if not assign_fp.exists():
        raise FileNotFoundError(f"找不到模块分配表：{assign_fp}")

    wgcna_universe = read_assignment_universe(assign_fp)
    if not wgcna_universe:
        raise RuntimeError("gene_module_assignments.tsv 中没有读到任何 gene_id")
    annot_universe = read_annotation_gene_universe(annot_dir)
    if not annot_universe:
        log_msg("WARNING", f"在 {annot_dir} 下没有读到 GO/KEGG 注释基因，后续富集可能为空")

    if background_mode == "wgcna_universe":
        background_set = set(wgcna_universe)
    elif background_mode == "annotated_wgcna_universe":
        background_set = set(wgcna_universe) & set(annot_universe)
    else:
        raise ValueError("background_mode 仅支持 wgcna_universe / annotated_wgcna_universe")

    if not background_set:
        raise RuntimeError("背景基因集合为空，请检查 wgcna_enrich.background_mode")

    log_msg("INFO", f"WGCNA 基因宇宙大小：{len(wgcna_universe)}")
    log_msg("INFO", f"注释基因宇宙大小：{len(annot_universe)}")
    log_msg("INFO", f"最终背景基因数：{len(background_set)}")

    module_candidates = choose_module_candidates(
        wgcna_dir=wgcna_dir,
        geneset_mode=geneset_mode,
        geneid_col=geneid_col,
        hub_top_n=hub_top_n,
        modules_include=modules_include,
        modules_exclude=modules_exclude,
    )
    if not module_candidates:
        raise RuntimeError("未找到任何可用模块基因集")
    log_msg("INFO", f"检测到 {len(module_candidates)} 个候选模块基因集")

    cor_map: Dict[str, Dict[str, float]] = {}
    fdr_map: Dict[str, Dict[str, float]] = {}
    if trait_filter_enable:
        if not target_trait:
            raise ValueError("开启 trait_filter 时，target_trait 不能为空")
        cor_fp = wgcna_dir / "assoc" / "module_trait_cor.tsv"
        fdr_fp = wgcna_dir / "assoc" / "module_trait_padj.tsv"
        cor_traits, cor_map = read_module_trait_matrix(cor_fp)
        fdr_traits, fdr_map = read_module_trait_matrix(fdr_fp)
        if target_trait not in cor_traits or target_trait not in fdr_traits:
            raise ValueError(f"target_trait 不在模块相关矩阵列名中：{target_trait}")
        log_msg(
            "INFO",
            f"启用 trait 过滤：target_trait={target_trait}, direction={trait_direction}, cor_min={trait_cor_min}, fdr_max={trait_fdr_max}",
        )

    task_rows: List[Dict[str, Any]] = []
    skipped = 0

    for mod, raw_genes in module_candidates:
        genes_in_bg = sorted(set(raw_genes) & background_set)

        if trait_filter_enable:
            passed, r_val, fdr_val, reason = module_pass_trait_filter(
                module_name=mod,
                target_trait=target_trait,
                cor_map=cor_map,
                fdr_map=fdr_map,
                trait_direction=trait_direction,
                trait_cor_min=trait_cor_min,
                trait_fdr_max=trait_fdr_max,
            )
            if not passed:
                skipped += 1
                log_msg("INFO", f"跳过模块 {mod}：{reason}")
                continue
            trait_r = f"{r_val:.6g}"
            trait_fdr = f"{fdr_val:.6g}"
            trait_sign = "positive" if r_val > 0 else ("negative" if r_val < 0 else "zero")
        else:
            trait_r = ""
            trait_fdr = ""
            trait_sign = ""

        if len(genes_in_bg) < min_genes_per_set:
            skipped += 1
            log_msg("INFO", f"跳过模块 {mod}：背景内基因数 {len(genes_in_bg)} < min_genes_per_set {min_genes_per_set}")
            continue

        if trait_filter_enable:
            label = f"wgcna_{sanitize_label(target_trait)}_{trait_sign}_{sanitize_label(mod)}"
        else:
            label = f"wgcna_{sanitize_label(mod)}"

        module_input_dir = inputs_dir / label
        module_outdir = outdir / label
        ensure_dir(module_input_dir)
        ensure_dir(module_outdir)

        test_fp = module_input_dir / "test.list"
        bg_fp = module_input_dir / "background.list"
        meta_fp = module_input_dir / "meta.tsv"

        write_single_column_list(test_fp, genes_in_bg, header=geneid_col)
        write_single_column_list(bg_fp, sorted(background_set), header=geneid_col)

        mapped_bg = len(set(background_set) & annot_universe)
        coverage = (mapped_bg / len(background_set)) if background_set else 0.0
        write_tsv(
            meta_fp,
            [{
                "label": label,
                "module": mod,
                "geneset_mode": geneset_mode,
                "n_detectable": len(background_set),
                "n_annot_mapped": mapped_bg,
                "universe_size": mapped_bg,
                "coverage": f"{coverage:.4f}",
                "detectable_rule": background_mode,
                "trait": target_trait,
                "trait_r": trait_r,
                "trait_fdr": trait_fdr,
            }],
            fieldnames=[
                "label", "module", "geneset_mode", "n_detectable", "n_annot_mapped",
                "universe_size", "coverage", "detectable_rule", "trait", "trait_r", "trait_fdr",
            ],
        )

        task_rows.append({
            "label": label,
            "type": "external",
            "test_set": "test",
            "test_file": str(test_fp),
            "background_file": str(bg_fp),
            "universe_meta": str(meta_fp),
            "outdir": str(module_outdir),
            "n_deg_all": "NA",
        })

    if not task_rows:
        raise RuntimeError("没有任何模块通过筛选，无法执行富集")

    tasks_fp = outdir / "tasks.tsv"
    write_tsv(
        tasks_fp,
        task_rows,
        fieldnames=["label", "type", "test_set", "test_file", "background_file", "universe_meta", "outdir", "n_deg_all"],
    )
    log_msg("INFO", f"已写出 tasks.tsv：{tasks_fp}")
    log_msg("INFO", f"通过模块数：{len(task_rows)}；跳过模块数：{skipped}")

    runtime_cfg = build_runtime_config(
        original_cfg=cfg,
        runtime_enrich_dir=outdir,
        annot_dir=annot_dir,
        maps_dir=maps_dir,
    )
    runtime_cfg_fp = runtime_dir / "config.yaml"
    write_yaml(runtime_cfg, runtime_cfg_fp)
    log_msg("INFO", f"已写出 runtime config：{runtime_cfg_fp}")

    if run_enrich:
        cmd = [rscript_bin, str(r_enrich_script)]
        log_msg("INFO", f"开始调用 08 富集脚本：{' '.join(cmd)}")
        try:
            subprocess.run(cmd, cwd=str(runtime_dir), check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"08_g_enrich.R 执行失败，返回码={e.returncode}") from e
        log_msg("INFO", "08_g_enrich.R 执行完成")
    else:
        log_msg("INFO", "run_enrich=False，本次仅准备桥接输入文件")


if __name__ == "__main__":
    main()

