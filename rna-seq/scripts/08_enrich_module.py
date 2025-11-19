#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
08_enrich_module.py —— ORA + GSEA 调度模块（vNext + GSEA 最终版）

职责：
  1) 读取 config.yaml，锁定目录与参数；
  2) 识别富集任务：
       - RNA 对比：results/06_deg/<contrast>/
       - 外部 tag：results/08_enrich/inputs/<tag>/
  3) 为每个 RNA 对比构建背景集合：
       - 可检测基因 = TPM>0 或 count>0（任一对比样本）
       - 背景 = 可检测 ∩ 注释宇宙（gene2go ∪ gene2pathway）
       - 写出 background/<contrast>.list + background/<contrast>.meta.tsv
  4) 为外部 tag 检查 background.list，并补齐 meta.tsv（如缺失）；
  5) 为 ORA 生成 test.list（RNA 由 DEG_* 转换而来，外部 tag 直接使用用户提供）；
  6) 若 config.gsea.enable = true，为每个 RNA 对比生成 GSEA 排名文件：
       - results/08_enrich/<contrast>/gsea_ranks.tsv
  7) 生成统一任务清单 tasks.tsv 交给 08_g_enrich.R；
  8) 调用 Rscript scripts/08_g_enrich.R 执行 ORA + GSEA。

注意：
  - 所有参数只从 config.yaml 读取；
  - 不改变现有目录结构与表头约定；
  - gene_id 是全流程主键，不在本脚本中做任何 ID 改写。
"""

from __future__ import annotations

import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set

import yaml
import csv


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
            header = None
            if first[0]:
                genes.append(first[0].strip())
            idx = 0

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
    import csv as _csv  # 局部导入避免覆盖
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


def filter_detectable_genes(
    gene_ids: List[str],
    samples: List[str],
    counts: Dict[str, List[float]],
    tpms: Dict[str, List[float]],
    sample_subset: List[str],
) -> Set[str]:
    """
    根据“TPM>0 或 count>0（任一对比样本）”规则，返回可检测基因集合。
    """
    sample_index = {s: i for i, s in enumerate(samples)}
    idxs = [sample_index[s] for s in sample_subset if s in sample_index]
    if not idxs:
        return set()

    detectable: Set[str] = set()
    for gid in gene_ids:
        cvals = counts.get(gid, [])
        tvals = tpms.get(gid, [])
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


# ====================== 主流程 ======================

def main() -> None:
    config_path = Path("config.yaml")
    cfg = load_config(config_path)
    setup_logging(cfg)
    log = logging.getLogger("enrich")

    # 目录与路径
    dirs_cfg = cfg.get("dirs", {}) or {}
    data_cfg = cfg.get("data", {}) or {}
    background_cfg = cfg.get("background", {}) or {}
    gsea_cfg = cfg.get("gsea", {}) or {}

    gsea_enable = bool(gsea_cfg.get("enable", False))
    gsea_score_from = str(gsea_cfg.get("score_from", "stat")).lower()

    project_root = Path(".").resolve()
    dir_matrix = project_root / dirs_cfg.get("matrix", "results/05_matrix")
    dir_deg = project_root / dirs_cfg.get("deg", "results/06_deg")
    dir_annot = project_root / dirs_cfg.get("annot", "results/07_annot")
    dir_enrich = project_root / dirs_cfg.get("enrich", "results/08_enrich")

    dir_bg = dir_enrich / "background"
    dir_inputs = project_root / background_cfg.get("ora_inputs_dir", "results/08_enrich/inputs")

    dir_matrix_counts = dir_matrix / "counts"
    dir_matrix_tpms = dir_matrix / "tpms"

    counts_fp = dir_matrix_counts / "gene_counts.tsv"
    tpm_fp = dir_matrix_tpms / "gene_tpm.tsv"
    samples_tsv = project_root / data_cfg.get("samples_tsv", "data/samples.tsv")
    contrasts_tsv = project_root / data_cfg.get("contrasts_tsv", "data/contrasts.tsv")

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

    # 读取 gene 矩阵
    try:
        gene_ids_counts, samples_counts, counts_mat = read_gene_matrix(counts_fp)
        gene_ids_tpm, samples_tpm, tpm_mat = read_gene_matrix(tpm_fp)
    except Exception as e:
        log.error("读取表达矩阵失败：%s", e)
        sys.exit(1)

    if gene_ids_counts != gene_ids_tpm or samples_counts != samples_tpm:
        log.warning("counts 与 TPM 的行/列不完全一致，将以 counts 为准对齐 TPM")

    gene_ids = gene_ids_counts
    samples = samples_counts

    # 样本与对比信息
    if not samples_tsv.exists() or not contrasts_tsv.exists():
        log.error("samples.tsv 或 contrasts.tsv 不存在：%s, %s", samples_tsv, contrasts_tsv)
        sys.exit(1)

    sample_group: Dict[str, str] = {}
    with samples_tsv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "sample" not in reader.fieldnames or "group" not in reader.fieldnames:
            log.error("samples.tsv 必须包含列：sample, group")
            sys.exit(1)
        for row in reader:
            sid = (row.get("sample") or "").strip()
            grp = (row.get("group") or "").strip()
            if sid:
                sample_group[sid] = grp

    contrasts: List[Tuple[str, str, str]] = []
    with contrasts_tsv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required_cols = ["contrast", "case", "control"]
        miss = [c for c in required_cols if c not in (reader.fieldnames or [])]
        if miss:
            log.error("contrasts.tsv 需要包含列：%s", ",".join(required_cols))
            sys.exit(1)
        for row in reader:
            label = (row.get("contrast") or "").strip()
            case = (row.get("case") or "").strip()
            ctrl = (row.get("control") or "").strip()
            if not label:
                continue
            contrasts.append((label, case, ctrl))

    log.info("检测到 RNA 对比数量：%d", len(contrasts))

    # 注释基因宇宙
    annot_gene_universe = read_gene_universe_from_annotations(dir_annot)
    log.info("注释基因宇宙大小（gene2go ∪ gene2pathway）：%d", len(annot_gene_universe))

    detectable_rule = background_cfg.get("detectable", "TPM>0_or_count>0_in>=1_sample")

    dir_bg.mkdir(parents=True, exist_ok=True)

    tasks: List[Dict[str, str]] = []
    label_deg_all: Dict[str, int] = {}
    label_deg_all_fp: Dict[str, Path] = {}

    # ---------- RNA 对比（ORA + GSEA 排名输入） ----------
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
        label_deg_all[label] = n_deg_all
        label_deg_all_fp[label] = deg_all_fp

        # 确定对比使用的样本
        used_samples = [s for s, g in sample_group.items() if g in {g_case, g_ctrl}]
        used_samples = [s for s in used_samples if s in samples]
        if not used_samples:
            log.warning("对比 %s 没有匹配到任何样本（根据 group），背景计算可能不可靠", label)

        # 可检测基因
        detectable = filter_detectable_genes(
            gene_ids=gene_ids,
            samples=samples,
            counts=counts_mat,
            tpms=tpm_mat,
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
        def build_test_list_from_deg(deg_fp: Path, test_set: str) -> Path | None:
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
                    "n_deg_all": str(label_deg_all.get(label, 0)),
                }
            )

    # ---------- 外部 tag 通道（通用 ORA 引擎） ----------
    if dir_inputs.exists():
        for tag_dir in sorted(p for p in dir_inputs.iterdir() if p.is_dir()):
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
                log.info(
                    "已为 tag=%s 生成 meta.tsv：universe_size=%d, coverage=%.3f",
                    tag,
                    mapped,
                    coverage,
                )

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
    if gsea_enable:
        log.info("GSEA 功能已启用，将为 RNA 对比生成 gsea_ranks.tsv（score_from=%s）", gsea_score_from)
        for label, deg_all_fp in label_deg_all_fp.items():
            outdir_label = dir_enrich / label
            outdir_label.mkdir(parents=True, exist_ok=True)
            ranks_fp = outdir_label / "gsea_ranks.tsv"

            # 读取 DEG_all.tsv
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
                    if not val_str:
                        continue
                    try:
                        val = float(val_str)
                    except ValueError:
                        continue
                    # 允许 score 为 0，但 GSEA 时会自动处理
                    rows.append((gid, val))

            if not rows:
                log.warning("对比 %s 的 GSEA 排名为空，跳过 gsea_ranks.tsv 写出。", label)
                continue

            # 去重：若同一 gene 多次出现，保留第一次
            seen: Set[str] = set()
            unique_rows: List[Tuple[str, float]] = []
            for gid, val in rows:
                if gid in seen:
                    continue
                seen.add(gid)
                unique_rows.append((gid, val))

            # 按 score 降序排序
            unique_rows.sort(key=lambda x: x[1], reverse=True)

            # 写出 gsea_ranks.tsv
            ranks_fp.parent.mkdir(parents=True, exist_ok=True)
            with ranks_fp.open("w", encoding="utf-8") as f_out:
                f_out.write("gene_id\tscore\trank_source\n")
                source = "DESeq2.stat" if gsea_score_from == "stat" else "DESeq2.log2fc"
                for gid, val in unique_rows:
                    f_out.write(f"{gid}\t{val}\t{source}\n")

            log.info("已为对比 %s 写出 GSEA 排名文件：%s", label, ranks_fp)

    # ---------- 写任务清单 ----------
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

    log.info("已写出任务清单：%s，任务数量=%d", tasks_fp, len(tasks))

    # ---------- 调用 R 富集脚本 ----------
    rscript_bin = get_nested(cfg, ["binaries", "Rscript"], default="Rscript")
    r_script = project_root / "scripts" / "08_g_enrich.R"
    if not r_script.exists():
        log.error("未找到 R 脚本：%s", r_script)
        sys.exit(1)

    cmd = [rscript_bin, str(r_script)]
    log.info("调用 Rscript 执行富集分析：%s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        log.error("R 富集脚本执行失败，返回码=%d", e.returncode)
        sys.exit(e.returncode)

    log.info("===== 08_enrich_module.py 执行完成 =====")


if __name__ == "__main__":
    main()