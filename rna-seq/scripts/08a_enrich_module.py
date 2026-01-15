#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
08a_enrich_module.py —— ORA + GSEA 调度模块（08a：原 08_enrich_module.py 改名版）

严格遵守皇上需求清单（关键点）：
- 参数优先级：优先从 config.yaml 读取；仅当 config.yaml 缺失相关配置时才使用脚本内部默认值。
- 增量兼容：不改变原有目录结构、表头约定、产物命名；原本能跑通的任务与输出保持不变。
- 启动即确保外包目录存在：results/08_enrich/inputs/ 不存在则创建；存在则不覆盖、不清空。
- 若 config.yaml 启用了 one_vs_rest（08b）功能：在写完 base tasks.tsv 后触发 08b；
  并确保 08_g_enrich.R 读取到的是“最终合并后的 tasks.tsv”（而不是 base 版）。
- 禁止生成任何 tasks.tsv.bak_* 之类备份文件（本脚本不生成；08b 也不生成）。

注意：
- 全程使用相对路径（默认从项目根目录运行：./scripts/08a_enrich_module.py）
- gene_id 为主键：本脚本不做任何 ID 改写。
"""

from __future__ import annotations

import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set, Optional

import yaml
import csv



# ====================== 用户自定义参数区（仅在 config 缺失相关配置时启用） ======================
# 皇上注意：这里是“兜底默认值”，只有在 config.yaml 缺失对应字段时才会被使用。
# 正常情况下（你的项目已提供完善 config.yaml），这些值不会生效。

USER_DEFAULTS = {
    # 二进制
    "binaries": {
        "Rscript": "Rscript",
    },
    # dirs
    "dirs": {
        "matrix": "results/05_matrix",
        "deg": "results/06_deg",
        "annot": "results/07_annot",
        "enrich": "results/08_enrich",
    },
    # data
    "data": {
        "samples_tsv": "data/samples.tsv",
        "contrasts_tsv": "data/contrasts.tsv",
    },
    # background
    "background": {
        "ora_inputs_dir": "results/08_enrich/inputs",
        "detectable": "TPM>0_or_count>0_in>=1_sample",
    },
    # gsea
    "gsea": {
        "enable": False,
        "score_from": "stat",
    },
    # one_vs_rest（08b）开关：缺失时默认关闭
    "one_vs_rest": {
        "enable": False,
        "targets": [],
    },
}

# ====================== 基础工具函数 ======================

def load_config(config_path: Path) -> Dict[str, Any]:
    """读取 config.yaml，并返回 dict。"""
    if not config_path.exists():
        print(f"[ERR] 找不到配置文件：{config_path}", file=sys.stderr)
        sys.exit(1)
    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg


def get_nested(cfg: Dict[str, Any], keys: List[str], default=None):
    """从多层 dict 安全取值。"""
    cur: Any = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def setup_logging(cfg: Dict[str, Any]) -> None:
    """根据 config.logging 初始化日志。"""
    log_cfg = cfg.get("logging", {}) or {}
    level_str = str(log_cfg.get("level", "INFO")).upper()
    level = getattr(logging, level_str, logging.INFO)
    ts = bool(log_cfg.get("timestamp", True))
    fmt = "%(asctime)s [08_enrich] %(levelname)s: %(message)s" if ts else "[08_enrich] %(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=fmt)


def read_single_column_list(path: Path) -> List[str]:
    """
    读取一列基因列表文件：
      - 若存在表头 gene_id，则取该列；
      - 否则取第一列。
    """
    genes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return []
        has_header = ("gene_id" in first) or (len(first) == 1 and first[0] == "gene_id")
        if has_header:
            header = first
            try:
                idx = header.index("gene_id")
            except ValueError:
                idx = 0
        else:
            idx = 0
            if first and first[0]:
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


def write_single_column_list(path: Path, items: List[str], header: str | None = None) -> None:
    """写出单列表文件，可选写 header。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        if header:
            f.write(header + "\n")
        for it in items:
            f.write(str(it) + "\n")


def read_gene_universe_from_annotations(annot_dir: Path) -> Set[str]:
    """
    从 gene2go.tsv 与 gene2pathway.tsv 中提取注释基因宇宙：
      gene_universe = union(gene2go.gene_id, gene2pathway.gene_id)
    """
    import csv as _csv
    gene2go_fp = annot_dir / "gene2go.tsv"
    gene2path_fp = annot_dir / "gene2pathway.tsv"

    genes: Set[str] = set()
    for fp in [gene2go_fp, gene2path_fp]:
        if fp.exists():
            with fp.open("r", encoding="utf-8") as f:
                reader = _csv.DictReader(f, delimiter="\t")
                for row in reader:
                    gid = (row.get("gene_id") or "").strip()
                    if gid:
                        genes.add(gid)
    return genes


def read_gene_matrix(path: Path) -> Tuple[List[str], List[str], Dict[str, List[float]]]:
    """
    读取 gene_counts 或 gene_tpm：
      返回：（gene_ids, sample_ids, matrix[gene_id] = [v1, v2, ...]）
    """
    if not path.exists():
        raise FileNotFoundError(str(path))
    matrix: Dict[str, List[float]] = {}
    gene_ids: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, [])
        if not header:
            raise RuntimeError(f"{path} 表头为空")
        sample_ids = header[1:]
        for row in reader:
            if not row:
                continue
            gid = row[0].strip()
            if not gid:
                continue
            gene_ids.append(gid)
            vals: List[float] = []
            for x in row[1:]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(0.0)
            matrix[gid] = vals
    return gene_ids, sample_ids, matrix


def align_matrix_to_samples(
    base_gene_ids: List[str],
    base_samples: List[str],
    mat_other: Dict[str, List[float]],
    other_samples: List[str],
) -> Dict[str, List[float]]:
    """
    将另一个矩阵（如 TPM）对齐到 base_samples 的顺序。
    - gene 以 base_gene_ids 为准：缺失的 gene/样本补 0。
    - 避免“counts 与 TPM 列顺序不同”造成 detectable 判定逻辑漏洞。
    """
    other_index = {s: i for i, s in enumerate(other_samples)}
    out: Dict[str, List[float]] = {}
    for gid in base_gene_ids:
        src = mat_other.get(gid, [])
        aligned: List[float] = []
        for s in base_samples:
            j = other_index.get(s, None)
            if j is None:
                aligned.append(0.0)
            else:
                aligned.append(src[j] if j < len(src) else 0.0)
        out[gid] = aligned
    return out


def filter_detectable_genes(
    gene_ids: List[str],
    samples: List[str],
    counts: Dict[str, List[float]],
    tpms: Optional[Dict[str, List[float]]],
    sample_subset: List[str],
) -> Set[str]:
    """
    根据“TPM>0 或 count>0（任一对比样本）”规则，返回可检测基因集合。
    - 若 tpms 为 None：退化为 counts>0 规则（仍可跑，但与 config 的 detectable 语义略有差异）。
    """
    sample_index = {s: i for i, s in enumerate(samples)}
    idxs = [sample_index[s] for s in sample_subset if s in sample_index]
    if not idxs:
        return set()

    detectable: Set[str] = set()
    for gid in gene_ids:
        cvals = counts.get(gid, [])
        tvals = tpms.get(gid, []) if tpms is not None else []
        flag = False
        for j in idxs:
            cv = cvals[j] if j < len(cvals) else 0.0
            tv = tvals[j] if j < len(tvals) else 0.0
            if cv > 0 or tv > 0:
                flag = True
                break
        if flag:
            detectable.add(gid)
    return detectable


def ensure_dir_exists(p: Path) -> None:
    """确保目录存在：不存在则创建；存在则不覆盖、不清空。"""
    p.mkdir(parents=True, exist_ok=True)


# ====================== 主流程 ======================

def main() -> None:
    config_path = Path("config.yaml")
    cfg = load_config(config_path)
    setup_logging(cfg)
    log = logging.getLogger("enrich")

    # 目录与路径（全部保持相对路径）
    dirs_cfg = cfg.get("dirs", USER_DEFAULTS["dirs"]) or USER_DEFAULTS["dirs"]
    data_cfg = cfg.get("data", USER_DEFAULTS["data"]) or USER_DEFAULTS["data"]
    background_cfg = cfg.get("background", USER_DEFAULTS["background"]) or USER_DEFAULTS["background"]
    gsea_cfg = cfg.get("gsea", USER_DEFAULTS["gsea"]) or USER_DEFAULTS["gsea"]

    gsea_enable = bool(gsea_cfg.get("enable", USER_DEFAULTS["gsea"]["enable"]))
    gsea_score_from = str(gsea_cfg.get("score_from", USER_DEFAULTS["gsea"]["score_from"])).lower()

    project_root = Path(".")  # 保持相对
    dir_matrix = project_root / dirs_cfg.get("matrix", USER_DEFAULTS["dirs"]["matrix"])
    dir_deg = project_root / dirs_cfg.get("deg", USER_DEFAULTS["dirs"]["deg"])
    dir_annot = project_root / dirs_cfg.get("annot", USER_DEFAULTS["dirs"]["annot"])
    dir_enrich = project_root / dirs_cfg.get("enrich", USER_DEFAULTS["dirs"]["enrich"])

    dir_bg = dir_enrich / "background"

    # 外包 inputs 目录：启动即确保存在（皇上硬要求）
    dir_inputs = project_root / background_cfg.get("ora_inputs_dir", USER_DEFAULTS["background"]["ora_inputs_dir"])
    ensure_dir_exists(dir_inputs)
    log.info("已检查/创建外包 inputs 目录（存在不覆盖）：%s", dir_inputs)

    dir_matrix_counts = dir_matrix / "counts"
    dir_matrix_tpms = dir_matrix / "tpms"

    counts_fp = dir_matrix_counts / "gene_counts.tsv"
    tpm_fp = dir_matrix_tpms / "gene_tpm.tsv"
    samples_tsv = project_root / data_cfg.get("samples_tsv", USER_DEFAULTS["data"]["samples_tsv"])
    contrasts_tsv = project_root / data_cfg.get("contrasts_tsv", USER_DEFAULTS["data"]["contrasts_tsv"])

    # 注释词典路径检查（存在性后交给 R 再精细报错）
    gene2go_fp = dir_annot / "gene2go.tsv"
    gene2path_fp = dir_annot / "gene2pathway.tsv"
    go_term2name_fp = dir_annot / "go" / "term2name.tsv"
    go_obsolete_fp = dir_annot / "go" / "obsolete_map.tsv"
    kegg_term2gene_fp = dir_annot / "kegg" / "term2gene.tsv"
    kegg_term2name_fp = dir_annot / "kegg" / "term2name.tsv"

    for fp in [gene2go_fp, gene2path_fp, go_term2name_fp, go_obsolete_fp, kegg_term2gene_fp, kegg_term2name_fp]:
        if not fp.exists():
            log.warning("注释文件缺失（可在 08_g_enrich.R 中进一步检查）：%s", fp)

    # 读取 gene 矩阵（若失败则仅使用外部 tag 通道）
    matrix_available = True
    try:
        gene_ids_counts, samples_counts, counts_mat = read_gene_matrix(counts_fp)
    except Exception as e:
        log.warning("读取 counts 矩阵失败，将跳过 RNA 对比，仅使用外部 tag；原因：%s", e)
        matrix_available = False
        gene_ids_counts, samples_counts, counts_mat = [], [], {}

    # TPM 可选：若失败不致命（detectable 退化）
    tpms_mat_aligned: Optional[Dict[str, List[float]]] = None
    if matrix_available:
        try:
            gene_ids_tpm, samples_tpm, tpm_mat = read_gene_matrix(tpm_fp)
            if samples_tpm != samples_counts:
                log.info("检测到 TPM 列顺序与 counts 不一致：将对齐 TPM 到 counts 的 samples 顺序（避免 detectable 逻辑漏洞）。")
            # gene 以 counts 为准；TPM 对齐到 counts 的样本顺序
            tpms_mat_aligned = align_matrix_to_samples(
                base_gene_ids=gene_ids_counts,
                base_samples=samples_counts,
                mat_other=tpm_mat,
                other_samples=samples_tpm,
            )
        except Exception as e:
            log.warning("读取 TPM 矩阵失败：detectable 将退化为 counts>0；原因：%s", e)
            tpms_mat_aligned = None

    gene_ids = gene_ids_counts
    samples = samples_counts

    # 样本与对比信息（若缺失则仅使用外部 tag 通道）
    sample_group: Dict[str, str] = {}
    contrasts: List[Tuple[str, str, str]] = []

    if matrix_available:
        if not samples_tsv.exists() or not contrasts_tsv.exists():
            log.warning(
                "samples.tsv 或 contrasts.tsv 不存在，将跳过 RNA 对比，仅处理外部 tag：%s, %s",
                samples_tsv,
                contrasts_tsv,
            )
            matrix_available = False
        else:
            with samples_tsv.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if "sample" not in (reader.fieldnames or []) or "group" not in (reader.fieldnames or []):
                    log.warning("samples.tsv 必须包含列 sample, group，将跳过 RNA 对比，仅处理外部 tag")
                    matrix_available = False
                else:
                    for row in reader:
                        sid = (row.get("sample") or "").strip()
                        grp = (row.get("group") or "").strip()
                        if sid:
                            sample_group[sid] = grp

            if matrix_available:
                with contrasts_tsv.open("r", encoding="utf-8") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    required_cols = ["contrast", "case", "control"]
                    miss = [c for c in required_cols if c not in (reader.fieldnames or [])]
                    if miss:
                        log.warning(
                            "contrasts.tsv 需要包含列：%s，将跳过 RNA 对比，仅处理外部 tag",
                            ",".join(required_cols),
                        )
                        matrix_available = False
                    else:
                        for row in reader:
                            label = (row.get("contrast") or "").strip()
                            case = (row.get("case") or "").strip()
                            ctrl = (row.get("control") or "").strip()
                            if not label:
                                continue
                            contrasts.append((label, case, ctrl))

    if matrix_available:
        log.info("检测到 RNA 对比数量：%d", len(contrasts))
    else:
        log.info("RNA 对比模块已关闭，本次仅使用外部 tag 的 ORA 背景（如存在）。")

    # 注释基因宇宙
    annot_gene_universe = read_gene_universe_from_annotations(dir_annot)
    log.info("注释基因宇宙大小（gene2go ∪ gene2pathway）：%d", len(annot_gene_universe))

    detectable_rule = background_cfg.get("detectable", USER_DEFAULTS["background"]["detectable"])

    ensure_dir_exists(dir_bg)

    tasks: List[Dict[str, str]] = []
    label_deg_all_fp: Dict[str, Path] = {}

    # ---------- RNA 对比（ORA + GSEA 排名输入） ----------
    if matrix_available:
        for label, g_case, g_ctrl in contrasts:
            log.info("处理 RNA 对比：%s (%s vs %s)", label, g_case, g_ctrl)
            deg_dir_label = dir_deg / label
            deg_all_fp = deg_dir_label / "DEG_all.tsv"
            deg_up_fp = deg_dir_label / "DEG_up.tsv"
            deg_down_fp = deg_dir_label / "DEG_down.tsv"

            if not deg_all_fp.exists():
                log.warning("DEG_all.tsv 不存在，跳过对比：%s", label)
                continue

            # 统计 n_deg_all
            with deg_all_fp.open("r", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                _header = next(reader, None)
                n_deg_all = sum(1 for _ in reader)

            label_deg_all_fp[label] = deg_all_fp

            # 确定对比使用的样本（按 group 匹配 case/control）
            used_samples = [s for s, g in sample_group.items() if g in {g_case, g_ctrl}]
            used_samples = [s for s in used_samples if s in samples]
            if not used_samples:
                log.warning("对比 %s 没有匹配到任何样本（根据 group），背景计算可能不可靠", label)

            # 可检测基因
            detectable = filter_detectable_genes(
                gene_ids=gene_ids,
                samples=samples,
                counts=counts_mat,
                tpms=tpms_mat_aligned,
                sample_subset=used_samples,
            )
            n_detectable = len(detectable)
            if n_detectable == 0:
                log.warning("对比 %s 的可检测基因数为 0，请检查表达矩阵与样本分组", label)

            # 背景路径
            bg_list_fp = dir_bg / f"{label}.list"
            bg_meta_fp = dir_bg / f"{label}.meta.tsv"

            if bg_list_fp.exists() and bg_meta_fp.exists():
                log.info("对比 %s 已存在背景文件，将直接使用现有 background.list/meta.tsv", label)
            else:
                universe = sorted(detectable & annot_gene_universe)
                universe_size = len(universe)
                coverage = (universe_size / float(n_detectable)) if n_detectable > 0 else 0.0

                write_single_column_list(bg_list_fp, universe, header="gene_id")

                bg_meta_fp.parent.mkdir(parents=True, exist_ok=True)
                with bg_meta_fp.open("w", encoding="utf-8") as meta_f:
                    meta_f.write(
                        "label\tn_detectable\tn_annot_mapped\tuniverse_size\tcoverage\tdetectable_rule\tsamples_used\n"
                    )
                    meta_f.write(
                        "\t".join(
                            [
                                label,
                                str(n_detectable),
                                str(universe_size),
                                str(universe_size),
                                f"{coverage:.4f}",
                                detectable_rule,
                                ",".join(sorted(used_samples)) if used_samples else "NA",
                            ]
                        )
                        + "\n"
                    )
                log.info(
                    "已生成对比 %s 的背景：universe_size=%d, coverage=%.3f",
                    label,
                    universe_size,
                    coverage,
                )

            outdir_label = dir_enrich / label
            outdir_label.mkdir(parents=True, exist_ok=True)

            # 为 ORA 构建 test.list（all/up/down）
            def build_test_list_from_deg(deg_fp: Path, test_set: str) -> Optional[Path]:
                if not deg_fp.exists():
                    log.info("对比 %s 缺少 %s，跳过该 test_set", label, deg_fp.name)
                    return None
                gene_ids_local: List[str] = []
                with deg_fp.open("r", encoding="utf-8") as f_in:
                    reader = csv.DictReader(f_in, delimiter="\t")
                    if "gene_id" not in (reader.fieldnames or []):
                        log.warning("DEG 文件缺少 gene_id 列：%s", deg_fp)
                        return None
                    for row in reader:
                        gid = (row.get("gene_id") or "").strip()
                        if gid:
                            gene_ids_local.append(gid)
                gene_ids_local = sorted(set(gene_ids_local))
                if not gene_ids_local:
                    log.info("对比 %s 的 %s 集合为空，将在富集中产生空结果表", label, test_set)
                test_fp = outdir_label / f"{test_set}.list"
                write_single_column_list(test_fp, gene_ids_local, header="gene_id")
                return test_fp

            for test_set, fp in [("all", deg_all_fp), ("up", deg_up_fp), ("down", deg_down_fp)]:
                test_fp = build_test_list_from_deg(fp, test_set)
                if test_fp is None:
                    continue
                tasks.append(
                    {
                        "label": label,
                        "type": "rna",
                        "test_set": test_set,
                        "test_file": str(test_fp),
                        "background_file": str(bg_list_fp),
                        "universe_meta": str(bg_meta_fp),
                        "outdir": str(outdir_label),
                        "n_deg_all": str(n_deg_all),
                    }
                )

    # ---------- 外部 tag 通道（通用 ORA 引擎） ----------
    # 注意：dir_inputs 由启动阶段保证存在；这里仅枚举其中的 tag 文件夹（如存在）。
    for tag_dir in sorted([p for p in dir_inputs.iterdir() if p.is_dir()]):
        tag = tag_dir.name
        log.info("处理外部 ORA tag：%s", tag)
        bg_list_fp = tag_dir / "background.list"
        if not bg_list_fp.exists():
            log.warning("tag=%s 缺少 background.list，跳过", tag)
            continue

        bg_meta_fp = tag_dir / "meta.tsv"
        if not bg_meta_fp.exists():
            bg_genes = read_single_column_list(bg_list_fp)
            n_bg = len(bg_genes)
            mapped = len(set(bg_genes) & annot_gene_universe)
            coverage = (mapped / float(n_bg)) if n_bg > 0 else 0.0
            bg_meta_fp.parent.mkdir(parents=True, exist_ok=True)
            with bg_meta_fp.open("w", encoding="utf-8") as meta_f:
                meta_f.write(
                    "label\tn_detectable\tn_annot_mapped\tuniverse_size\tcoverage\tdetectable_rule\tsamples_used\n"
                )
                meta_f.write(
                    "\t".join(
                        [
                            tag,
                            "NA",
                            str(mapped),
                            str(mapped),
                            f"{coverage:.4f}",
                            "external_provided_background",
                            "NA",
                        ]
                    )
                    + "\n"
                )
            log.info("已为 tag=%s 生成 meta.tsv：universe_size=%d, coverage=%.3f", tag, mapped, coverage)

        outdir_label = dir_enrich / tag
        outdir_label.mkdir(parents=True, exist_ok=True)

        test_main_fp = tag_dir / "test.list"
        up_fp = tag_dir / "up.list"
        down_fp = tag_dir / "down.list"

        if test_main_fp.exists():
            tasks.append(
                {
                    "label": tag,
                    "type": "external",
                    "test_set": "test",
                    "test_file": str(test_main_fp),
                    "background_file": str(bg_list_fp),
                    "universe_meta": str(bg_meta_fp),
                    "outdir": str(outdir_label),
                    "n_deg_all": "NA",
                }
            )
        else:
            if up_fp.exists():
                tasks.append(
                    {
                        "label": tag,
                        "type": "external",
                        "test_set": "up",
                        "test_file": str(up_fp),
                        "background_file": str(bg_list_fp),
                        "universe_meta": str(bg_meta_fp),
                        "outdir": str(outdir_label),
                        "n_deg_all": "NA",
                    }
                )
            if down_fp.exists():
                tasks.append(
                    {
                        "label": tag,
                        "type": "external",
                        "test_set": "down",
                        "test_file": str(down_fp),
                        "background_file": str(bg_list_fp),
                        "universe_meta": str(bg_meta_fp),
                        "outdir": str(outdir_label),
                        "n_deg_all": "NA",
                    }
                )

    # ---------- GSEA 排名文件（仅 RNA 对比，按 config.gsea 决定） ----------
    if gsea_enable and matrix_available:
        log.info("GSEA 功能已启用，将为 RNA 对比生成 gsea_ranks.tsv（score_from=%s）", gsea_score_from)
        for label, deg_all_fp in label_deg_all_fp.items():
            outdir_label = dir_enrich / label
            outdir_label.mkdir(parents=True, exist_ok=True)
            ranks_fp = outdir_label / "gsea_ranks.tsv"

            with deg_all_fp.open("r", encoding="utf-8") as f_in:
                reader = csv.DictReader(f_in, delimiter="\t")
                fieldnames = reader.fieldnames or []
                score_col = None
                if gsea_score_from == "stat" and "stat" in fieldnames:
                    score_col = "stat"
                elif gsea_score_from in ("log2fc", "log2foldchange", "log2_fold_change"):
                    if "log2fc" in fieldnames:
                        score_col = "log2fc"
                    elif "log2FoldChange" in fieldnames:
                        score_col = "log2FoldChange"
                elif "stat" in fieldnames:
                    score_col = "stat"

                if score_col is None:
                    log.warning(
                        "对比 %s 的 DEG_all.tsv 中找不到适合作为 GSEA 排名的列（stat/log2FC），跳过 GSEA 排名文件。",
                        label,
                    )
                    continue

                rows: List[Tuple[str, float]] = []
                for row in reader:
                    gid = (row.get("gene_id") or "").strip()
                    if not gid:
                        continue
                    val_str = (row.get(score_col) or "").strip()
                    if val_str == "":
                        continue
                    try:
                        val = float(val_str)
                    except ValueError:
                        continue
                    rows.append((gid, val))

            if not rows:
                log.warning("对比 %s 的 GSEA 排名为空，跳过 gsea_ranks.tsv 写出。", label)
                continue

            # 去重：若同一 gene 多次出现，保留第一次（与原脚本一致）
            seen: Set[str] = set()
            unique_rows: List[Tuple[str, float]] = []
            for gid, val in rows:
                if gid in seen:
                    continue
                seen.add(gid)
                unique_rows.append((gid, val))

            unique_rows.sort(key=lambda x: x[1], reverse=True)

            with ranks_fp.open("w", encoding="utf-8") as f_out:
                f_out.write("gene_id\tscore\trank_source\n")
                source = "DESeq2.stat" if gsea_score_from == "stat" else "DESeq2.log2fc"
                for gid, val in unique_rows:
                    f_out.write(f"{gid}\t{val}\t{source}\n")

            log.info("已为对比 %s 写出 GSEA 排名文件：%s", label, ranks_fp)

    # ---------- 写 base 任务清单（先写 base，再按需触发 08b 合并） ----------
    if not tasks:
        log.error("未检测到任何 ORA 任务（RNA 对比与外部 tag 均为空），08 模块无需执行。")
        sys.exit(1)

    tasks_fp = dir_enrich / "tasks.tsv"
    tasks_fp.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "label",
        "type",
        "test_set",
        "test_file",
        "background_file",
        "universe_meta",
        "outdir",
        "n_deg_all",
    ]
    with tasks_fp.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for t in tasks:
            writer.writerow(t)

    log.info("已写出 base 任务清单：%s，任务数量=%d", tasks_fp, len(tasks))

    # ---------- 按 config 开关触发 08b（one-vs-rest） ----------
    one_vs_rest_cfg = cfg.get("one_vs_rest", USER_DEFAULTS["one_vs_rest"]) or USER_DEFAULTS["one_vs_rest"]
    one_vs_rest_enable = bool(one_vs_rest_cfg.get("enable", False))
    targets = one_vs_rest_cfg.get("targets", None)
    # 若 targets 缺失/为空：按“启用但无 targets”视作不运行（避免看似成功但无任务）
    targets_ok = isinstance(targets, list) and any(str(x).strip() for x in targets)

    if one_vs_rest_enable and targets_ok:
        rscript_bin = get_nested(cfg, ["binaries", "Rscript"], default="Rscript")
        r_08b = project_root / "scripts" / "08b_any_vs_other.R"
        if not r_08b.exists():
            log.error("config.yaml 启用了 one_vs_rest，但未找到脚本：%s", r_08b)
            sys.exit(1)

        cmd_08b = [rscript_bin, str(r_08b)]
        log.info("one_vs_rest 已启用：触发 08b 合并 tasks.tsv：%s", " ".join(cmd_08b))
        try:
            subprocess.run(cmd_08b, check=True)
        except subprocess.CalledProcessError as e:
            log.error("08b_any_vs_other.R 执行失败，返回码=%d", e.returncode)
            sys.exit(e.returncode)

        # 轻量复核：确保 tasks.tsv 仍可读且非空（不做任何备份/覆盖）
        try:
            with tasks_fp.open("r", encoding="utf-8") as f:
                n_lines = sum(1 for _ in f)
            if n_lines <= 1:
                log.error("08b 合并后 tasks.tsv 为空（只有表头），请检查 08b 输出。")
                sys.exit(1)
            log.info("08b 合并完成：最终 tasks.tsv 行数（含表头）=%d", n_lines)
        except Exception as e:
            log.error("无法读取 08b 合并后的 tasks.tsv：%s", e)
            sys.exit(1)
    else:
        if one_vs_rest_enable and not targets_ok:
            log.warning("one_vs_rest.enable=true 但 targets 为空/缺失：本次不触发 08b。")
        else:
            log.info("one_vs_rest 未启用：本次不触发 08b。")

    # ---------- 调用 R 富集脚本（必须读最终 tasks.tsv） ----------
    rscript_bin = get_nested(cfg, ["binaries", "Rscript"], default="Rscript")
    r_script = project_root / "scripts" / "08_g_enrich.R"
    if not r_script.exists():
        log.error("未找到 R 脚本：%s", r_script)
        sys.exit(1)

    cmd = [rscript_bin, str(r_script)]
    log.info("调用 Rscript 执行富集分析（读取最终 tasks.tsv）：%s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        log.error("R 富集脚本执行失败，返回码=%d", e.returncode)
        sys.exit(e.returncode)

    log.info("===== 08a_enrich_module.py 执行完成 =====")


if __name__ == "__main__":
    main()

