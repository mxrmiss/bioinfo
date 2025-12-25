#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13b_cafe_build_enrich_inputs.py —— 为 CAFE5 扩张/收缩构建发布包（叶节点富集包 + 全节点OG包）

【皇上最新规则（本脚本已落实）】
  A）不设置任何“节点开关”：所有节点（祖先节点 + 叶节点）必须出现在“该出现的表格”中；
  B）不再读取 config.yaml 中的 cafe5.enrich_bridge：
        - 叶节点富集发布包：对“所有叶节点物种”逐物种输出；
        - 全节点OG发布包：对“所有节点（祖先+叶）”逐节点输出，便于共扩张/共收缩交集分析；
  C）输出位置与命名：
        1）叶节点富集发布包：写到
              <cafe_agg_dir>/enrich_inputs/<S_constricta_cafe>/
           文件固定四件套，且“表头/内容格式”保持原样：
              background.list / up.list / down.list / meta.tsv
        2）全节点OG发布包：写到
              <cafe_agg_dir>/enrich_inputs/node/<node_pkg_name>/
           其中：
              - 祖先节点：{node_id}node_cafe，例如 40node_cafe
              - 叶节点：S_constricta_node_cafe
           OG发布包为新增产物，表头由本脚本自定义（不受“已有表头不准改”的限制）

【与 13 步对齐（避免衔接误差）】
  - 13 输出：cafe_significant_families_no_highfail.tsv
      表头包含：model / family / Q
    本脚本严格按列名读取，不使用列号猜测。

【显著扩/缩判定口径（与旧13b一致）】
  认定“某个 family 在某个 node 上显著扩/缩”需满足：
    1) family 的 Q_family <= report.fdr_alpha（全局FDR显著）
    2) 若存在 Base_branch_probabilities.tab：
         p_branch <= cafe5.branch_p_alpha
       否则在 ALLOW_NO_BRANCH_PROB=True 时退化为仅 Q + Δ
    3) Base_change.tab 的 Δ != 0
  并据 Δ 正负标记 status：expanded / contracted

【表格输出（旧表头不改）】
  - cafe_family_branch_status.tsv：表头不变，仅增加祖先节点行
  - cafe_significant_clade_summary.tsv：★加回，生成内容应与旧13b一致（仅叶节点/物种）
  - cafe_node_significant_clade_summary.tsv：★新增，补齐祖先节点信息（只给 family 计数，不给基因计数）

注意：
  - 祖先节点没有真实“基因ID”，因此本脚本不为祖先节点输出基因级富集包；
    祖先节点只输出OG级发布包（用于共扩张/共收缩交集与后续统计）。
"""

from __future__ import annotations

import sys
import logging
import traceback
import math
import re
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Set, Any

import yaml

# -------------------------- 顶部参数区 --------------------------

CONFIG_PATH = "config.yaml"
LOG_FILE_BASENAME = "13b_cafe_build_enrich_inputs.log"
LOG_LEVEL = "INFO"

# 当 Base_branch_probabilities.tab 缺失时，是否允许退化为仅使用 Q + Δ 判定分支显著性
ALLOW_NO_BRANCH_PROB = True


# ----------------------- 基础工具 & 配置读取 -----------------------


def mkdir_p(p: Path) -> Path:
    """类似 mkdir -p 功能，若目录存在则不报错。"""
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_lines(path: Path) -> List[str]:
    """读取文本文件全部行，去掉行尾换行符。"""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\n\r") for ln in f]


def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    """写 TSV 文件，None 写为空字符串。"""
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")


def _expand_publish_placeholders(obj: Any, publish_dir: str) -> Any:
    """将对象中的 <publish_dir> 占位符展开。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj


def load_config(config_path: str) -> Dict[str, Any]:
    """
    读取 config.yaml，并对 inputs 段落中的 <publish_dir> 占位符进行展开：
      - publish_dir: "../phylo/results/publish/aphylo_ready"
      - inputs.*: "<publish_dir>/xxxx"
    """
    p = Path(config_path).resolve()
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg


@dataclass
class NodeInfo:
    node_id: int
    taxon_label: str
    species_name: Optional[str]
    is_leaf: bool


# ----------------------- CAFE 结果解析工具 -----------------------


def parse_node_label_to_id(label: str) -> Optional[int]:
    """
    从列名或 Taxon_ID 中解析 node_id：
      - "<22>"           → 22
      - "Name<10>"       → 10
    """
    label = label.strip()
    m_int = re.fullmatch(r"<(\d+)>", label)
    if m_int:
        return int(m_int.group(1))
    m_leaf = re.fullmatch(r".+<(\d+)>", label)
    if m_leaf:
        return int(m_leaf.group(1))
    return None


def parse_clade_results(path: Path) -> Dict[int, NodeInfo]:
    """
    解析 Base_clade_results.txt，构建 node_id -> NodeInfo 映射。
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_clade_results.txt 为空：{path}")
    node_map: Dict[int, NodeInfo] = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        taxon = parts[0].strip()

        m_int = re.fullmatch(r"<(\d+)>", taxon)
        if m_int:
            nid = int(m_int.group(1))
            node_map[nid] = NodeInfo(
                node_id=nid,
                taxon_label=taxon,
                species_name=None,
                is_leaf=False,
            )
            continue

        m_leaf = re.fullmatch(r"(.+)<(\d+)>", taxon)
        if m_leaf:
            name = m_leaf.group(1).strip()
            nid = int(m_leaf.group(2))
            node_map[nid] = NodeInfo(
                node_id=nid,
                taxon_label=taxon,
                species_name=name,
                is_leaf=True,
            )
            continue

        logging.warning(f"[WARN] 未能解析 Taxon_ID 字段：{taxon}（{path}）")

    if not node_map:
        raise ValueError(f"[ERR] 未从 {path} 解析出任何节点信息")
    return node_map


def parse_change_table(path: Path) -> List[Tuple[str, int, float]]:
    """
    解析 Base_change.tab，提取 (family, node_id, delta) 列表。
    兼容矩阵格式与旧长表格式。
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_change.tab 为空：{path}")
    header = re.split(r"\s+", lines[0].strip())
    data_lines = lines[1:]

    is_matrix = False
    if header and header[0].lower().startswith("family"):
        for col in header[1:]:
            if parse_node_label_to_id(col) is not None:
                is_matrix = True
                break

    out: List[Tuple[str, int, float]] = []

    if is_matrix:
        node_ids_by_col: Dict[int, int] = {}
        for idx, col in enumerate(header[1:], start=1):
            nid = parse_node_label_to_id(col)
            if nid is not None:
                node_ids_by_col[idx] = nid

        for ln in data_lines:
            if not ln.strip():
                continue
            parts = re.split(r"\s+", ln.strip())
            if len(parts) < 2:
                continue
            fam = parts[0].strip()
            if not fam:
                continue
            for idx, node_id in node_ids_by_col.items():
                if idx >= len(parts):
                    continue
                delta_s = parts[idx].strip()
                if not delta_s or delta_s.upper() in ("NA", "N/A", "."):
                    continue
                try:
                    delta = float(delta_s)
                except Exception:
                    continue
                if math.isnan(delta):
                    continue
                out.append((fam, node_id, delta))
    else:
        data_lines2 = lines
        first = lines[0].strip()
        if first.startswith("FamilyID") or first.startswith("#FamilyID"):
            data_lines2 = lines[1:]
        for ln in data_lines2:
            if not ln.strip():
                continue
            parts = re.split(r"\s+", ln.strip())
            if len(parts) < 5:
                continue
            fam = parts[0].strip()
            node_s = parts[2].strip()
            delta_s = parts[4].strip()
            try:
                node_id = int(node_s)
                delta = float(delta_s)
            except Exception:
                continue
            if math.isnan(delta):
                continue
            out.append((fam, node_id, delta))

    return out


def parse_branch_probabilities(path: Path) -> Dict[Tuple[str, int], float]:
    """
    解析 Base_branch_probabilities.tab，提取 (family, node_id) -> p_branch 映射。
    兼容矩阵格式与旧长表格式。
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_branch_probabilities.tab 为空：{path}")
    header = re.split(r"\s+", lines[0].strip())
    hlow = [h.lower() for h in header]
    data_lines = lines[1:]

    out: Dict[Tuple[str, int], float] = {}

    is_matrix = False
    if header and header[0].lower().startswith("family"):
        for col in header[1:]:
            if parse_node_label_to_id(col) is not None:
                is_matrix = True
                break

    if is_matrix:
        node_ids_by_col: Dict[int, int] = {}
        for idx, col in enumerate(header[1:], start=1):
            nid = parse_node_label_to_id(col)
            if nid is not None:
                node_ids_by_col[idx] = nid

        for ln in data_lines:
            if not ln.strip():
                continue
            parts = re.split(r"\s+", ln.strip())
            if len(parts) < 2:
                continue
            fam = parts[0].strip()
            if not fam:
                continue
            for idx, node_id in node_ids_by_col.items():
                if idx >= len(parts):
                    continue
                prob_s = parts[idx].strip()
                if not prob_s or prob_s.upper() in ("NA", "N/A", "."):
                    continue
                try:
                    p = float(prob_s)
                except Exception:
                    continue
                if math.isnan(p):
                    continue
                out[(fam, node_id)] = p
        return out

    family_idx = None
    node_idx = None
    prob_idx = None

    for i, h in enumerate(hlow):
        if family_idx is None and "famil" in h:
            family_idx = i
        if node_idx is None and ("node" in h or "clade" in h or "taxon" in h):
            node_idx = i
        if prob_idx is None and ("prob" in h or "pvalue" in h or h in ("p", "p.", "pval")):
            prob_idx = i

    if family_idx is None or node_idx is None or prob_idx is None:
        raise ValueError(f"[ERR] 无法从表头识别 family/node/prob 列：{header}")

    for ln in data_lines:
        if not ln.strip():
            continue
        parts = re.split(r"\s+", ln.strip())
        if len(parts) <= max(family_idx, node_idx, prob_idx):
            continue
        fam = parts[family_idx].strip()
        node_s = parts[node_idx].strip()
        prob_s = parts[prob_idx].strip()
        if prob_s.upper() in ("NA", "N/A", "."):
            continue
        try:
            node_id = int(node_s)
            p = float(prob_s)
        except Exception:
            continue
        if math.isnan(p):
            continue
        out[(fam, node_id)] = p

    return out


def species_to_leaf_pkg_prefix(species: str) -> str:
    """Sinonovacula_constricta -> S_constricta"""
    parts = species.strip().split("_")
    if len(parts) >= 2 and parts[0] and parts[1]:
        return f"{parts[0][0].upper()}_{parts[1]}"
    return species.strip()


def write_list_plain(path: Path, items: Set[str]) -> None:
    """写 list 文件：一行一个字符串，无表头。"""
    with path.open("w", encoding="utf-8") as f:
        for x in sorted(items):
            f.write(x + "\n")


# ----------------------- 主逻辑：生成发布包 -----------------------


def main() -> None:
    # 1) 读取配置
    cfg_path = Path(CONFIG_PATH).resolve()
    cfg = load_config(str(cfg_path))

    paths = cfg.get("paths", {})
    cafe_run_dir = Path(paths.get("cafe_run_dir", "results/06_cafe")).resolve()
    cafe_agg_dir = Path(paths.get("cafe_agg_dir", "results/07_cafeagg")).resolve()
    logs_dir = Path(paths.get("logs_dir", "logs")).resolve()

    mkdir_p(logs_dir)
    mkdir_p(cafe_agg_dir)

    # 2) 日志初始化
    logfile = logs_dir / LOG_FILE_BASENAME
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(logfile, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("========== APhylo 13b — CAFE5 发布包（叶节点富集 + 全节点OG） ==========")
    logging.info(f"[init] 使用配置文件：{cfg_path}")
    logging.info(f"[init] cafe_run_dir      = {cafe_run_dir}")
    logging.info(f"[init] cafe_agg_dir      = {cafe_agg_dir}")
    logging.info(f"[init] logs_dir          = {logs_dir}")

    # 3) 前置检查：13 步产物
    sig_nohf_path = cafe_agg_dir / "cafe_significant_families_no_highfail.tsv"
    if not sig_nohf_path.exists():
        logging.error(f"[ERR] 未找到 13 步产物：{sig_nohf_path}")
        sys.exit(2)

    sentinel = cafe_run_dir / ".cafe.done"
    if not sentinel.exists():
        logging.warning(f"[WARN] 未发现 12 步完成哨兵：{sentinel}，请确认 CAFE 已完整运行。")

    # 4) 阈值
    cafe_cfg = cfg.get("cafe5", {})
    report_cfg = cfg.get("report", {})

    fdr_alpha = float(report_cfg.get("fdr_alpha", 0.05))
    if not (0.0 < fdr_alpha <= 1.0):
        logging.warning(f"[WARN] report.fdr_alpha 非法：{fdr_alpha}，重置为 0.05")
        fdr_alpha = 0.05
    logging.info(f"[init] FDR α = {fdr_alpha}")

    branch_p_alpha = float(cafe_cfg.get("branch_p_alpha", 0.05))
    if not (0.0 < branch_p_alpha <= 1.0):
        logging.warning(f"[WARN] cafe5.branch_p_alpha 非法：{branch_p_alpha}，重置为 0.05")
        branch_p_alpha = 0.05
    logging.info(f"[init] 分支显著性阈值 branch_p_alpha = {branch_p_alpha}")

    # 5) 读取 13 的 FDR 表：构建 (model,family)->Q；以及 fam_all / fam_sig
    fam2q: Dict[Tuple[str, str], float] = {}
    fam_all_by_model: Dict[str, Set[str]] = {}
    fam_sig_by_model: Dict[str, Set[str]] = {}

    lines = read_lines(sig_nohf_path)
    if not lines:
        logging.error(f"[ERR] {sig_nohf_path} 为空")
        sys.exit(2)
    header = lines[0].split("\t")
    col_map = {name: i for i, name in enumerate(header)}
    required_cols = ["model", "family", "Q"]
    for c in required_cols:
        if c not in col_map:
            logging.error(f"[ERR] cafe_significant_families_no_highfail.tsv 缺少列：{c}")
            sys.exit(2)

    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) <= max(col_map.values()):
            continue
        model = parts[col_map["model"]].strip()
        fam = parts[col_map["family"]].strip()
        q_s = parts[col_map["Q"]].strip()
        if not model or not fam:
            continue
        try:
            q = float(q_s)
        except Exception:
            continue
        if math.isnan(q):
            continue
        fam2q[(model, fam)] = q
        fam_all_by_model.setdefault(model, set()).add(fam)
        if q <= fdr_alpha:
            fam_sig_by_model.setdefault(model, set()).add(fam)

    if not fam_all_by_model:
        logging.error("[ERR] 未从 FDR 表中解析到任何 family 记录")
        sys.exit(2)

    total_fam = sum(len(v) for v in fam_all_by_model.values())
    total_sig = sum(len(v) for v in fam_sig_by_model.values())
    logging.info(f"[cafe] 总 family 条目数：{total_fam}，其中 Q<=α 的显著 family 数：{total_sig}")

    # 6) 选择模型：优先 report.selected_cafe_model，否则 cafe5.models[0]
    selected_model = report_cfg.get("selected_cafe_model")
    cafe_models = cafe_cfg.get("models", ["global"])
    if not selected_model:
        selected_model = cafe_models[0]
        logging.warning(f"[WARN] report.selected_cafe_model 未设置，将使用 cafe5.models[0]：{selected_model}")
    else:
        logging.info(f"[init] selected_cafe_model = {selected_model}")

    models = [selected_model]
    models_dir = cafe_run_dir / "models"

    # 7) 解析 CAFE 输出：构建“全节点”显著扩/缩 family 集合 + “叶节点（物种）”集合
    expanded_fams_node: Dict[Tuple[str, int], Set[str]] = {}
    contracted_fams_node: Dict[Tuple[str, int], Set[str]] = {}

    expanded_fams_leaf: Dict[Tuple[str, str], Set[str]] = {}
    contracted_fams_leaf: Dict[Tuple[str, str], Set[str]] = {}

    node_map_by_model: Dict[str, Dict[int, NodeInfo]] = {}
    leaf_species_by_model: Dict[str, Set[str]] = {}

    branch_status_rows: List[List[object]] = []

    for model in models:
        mdir = models_dir / model
        primary_dir = mdir / "primary_global"
        if not primary_dir.exists():
            logging.error(f"[ERR] 缺少 primary_global 目录：{primary_dir}")
            sys.exit(2)

        clade_path = primary_dir / "Base_clade_results.txt"
        change_path = primary_dir / "Base_change.tab"
        branch_prob_path = primary_dir / "Base_branch_probabilities.tab"

        if not clade_path.exists():
            logging.error(f"[ERR] 缺少 Base_clade_results.txt：{clade_path}")
            sys.exit(2)
        if not change_path.exists():
            logging.error(f"[ERR] 缺少 Base_change.tab：{change_path}")
            sys.exit(2)

        node_map = parse_clade_results(clade_path)
        node_map_by_model[model] = node_map

        leaf_species: Set[str] = set()
        for info in node_map.values():
            if info.is_leaf and info.species_name:
                leaf_species.add(info.species_name)
        leaf_species_by_model[model] = leaf_species
        logging.info(f"[model={model}] 叶节点物种数：{len(leaf_species)}；全节点数：{len(node_map)}")

        change_records = parse_change_table(change_path)
        logging.info(f"[model={model}] 从 {change_path} 解析到 {len(change_records)} 条 Δ 记录")

        branch_prob: Dict[Tuple[str, int], float] = {}
        has_branch_prob = False
        if branch_prob_path.exists():
            try:
                branch_prob = parse_branch_probabilities(branch_prob_path)
                has_branch_prob = True
                logging.info(f"[model={model}] 从 {branch_prob_path} 解析到 {len(branch_prob)} 条分支概率记录")
            except Exception as e:
                logging.warning(f"[WARN] 解析 Base_branch_probabilities.tab 失败，将忽略分支 p 值：{e}")
                if not ALLOW_NO_BRANCH_PROB:
                    logging.error("[ERR] 禁止在缺少分支概率的情况下继续执行（ALLOW_NO_BRANCH_PROB=false）")
                    sys.exit(2)
        else:
            logging.warning(f"[WARN] 未找到 Base_branch_probabilities.tab：{branch_prob_path}，仅使用 Q + Δ 判定分支显著性")
            if not ALLOW_NO_BRANCH_PROB:
                logging.error("[ERR] 禁止在缺少分支概率的情况下继续执行（ALLOW_NO_BRANCH_PROB=false）")
                sys.exit(2)

        fam_all = fam_all_by_model.get(model, set())
        logging.info(f"[model={model}] fam_all={len(fam_all)}，fam_sig={len(fam_sig_by_model.get(model, set()))}（Q<=α）")

        for fam, node_id, delta in change_records:
            if fam_all and fam not in fam_all:
                continue

            q = fam2q.get((model, fam))
            if q is None or q > fdr_alpha:
                continue

            node_info = node_map.get(node_id)
            if node_info is None:
                continue

            if delta == 0:
                continue

            p_branch = None
            if has_branch_prob:
                p_branch = branch_prob.get((fam, node_id))
                if p_branch is None:
                    continue
                if p_branch > branch_p_alpha:
                    continue

            status = "expanded" if delta > 0 else "contracted"

            # ★全节点 per-family 状态表（旧表头不变）
            branch_status_rows.append([
                model,
                fam,
                node_id,
                node_info.taxon_label,
                "" if node_info.species_name is None else node_info.species_name,
                status,
                f"{delta:.6g}",
                "" if p_branch is None else f"{p_branch:.6g}",
                f"{q:.6g}",
            ])

            key_node = (model, node_id)
            if status == "expanded":
                expanded_fams_node.setdefault(key_node, set()).add(fam)
            else:
                contracted_fams_node.setdefault(key_node, set()).add(fam)

            if node_info.is_leaf and node_info.species_name:
                key_leaf = (model, node_info.species_name)
                if status == "expanded":
                    expanded_fams_leaf.setdefault(key_leaf, set()).add(fam)
                else:
                    contracted_fams_leaf.setdefault(key_leaf, set()).add(fam)

    # 8) 写出 per-family 状态表（已有表格：表头不变，只增行）
    status_path = cafe_agg_dir / "cafe_family_branch_status.tsv"
    write_tsv(
        status_path,
        ["model", "family", "node_id", "taxon_label", "species", "status", "delta", "p_branch", "Q_family"],
        branch_status_rows,
    )
    logging.info(f"[OK] 写出 per-family 状态表：{status_path}（{len(branch_status_rows)} 行；含祖先+叶）")

    # 9) 读取成员映射：inputs.all_members_tsv（all_pep2cds_resolved.tsv）
    inputs_cfg = cfg.get("inputs", {})
    all_members_path = inputs_cfg.get("all_members_tsv")
    if not all_members_path:
        logging.error("[ERR] config.inputs.all_members_tsv 未设置，无法构建 OG->基因 映射")
        sys.exit(2)
    all_members_path = Path(all_members_path).resolve()
    if not all_members_path.exists():
        logging.error(f"[ERR] all_members_tsv 文件不存在：{all_members_path}")
        sys.exit(2)

    lines = read_lines(all_members_path)
    if not lines:
        logging.error(f"[ERR] all_members_tsv 文件为空：{all_members_path}")
        sys.exit(2)

    header = lines[0].split("\t")
    hmap = {name: i for i, name in enumerate(header)}
    for c in ("OG", "Species", "cds_id"):
        if c not in hmap:
            logging.error(f"[ERR] all_members_tsv 缺少列：{c}")
            sys.exit(2)

    member_map: Dict[Tuple[str, str], Set[str]] = {}
    total_records = 0
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) <= max(hmap.values()):
            continue
        og = parts[hmap["OG"]].strip()
        sp = parts[hmap["Species"]].strip()
        cds = parts[hmap["cds_id"]].strip()
        if not og or not sp or not cds:
            continue
        member_map.setdefault((og, sp), set()).add(cds)
        total_records += 1
    logging.info(f"[members] 读取到 {total_records} 条 (OG,Species,cds) 映射，键数：{len(member_map)}")

    # ==================== ★ 加回：cafe_significant_clade_summary.tsv（内容与旧13b一致口径） ====================
    # 说明：
    #   - 该表“只统计叶节点物种”，列：model, species, n_expanded_fams, n_contracted_fams, n_expanded_genes, n_contracted_genes
    #   - 统计方式：对每个物种，取 expanded_fams_leaf / contracted_fams_leaf，
    #              并把这些 family 在该物种上的成员基因并起来计数（基于 all_members_tsv）
    sig_summary_rows: List[List[object]] = []

    for model in models:
        # 该模型下所有叶节点物种名（来自 Base_clade_results）
        species_set: Set[str] = set(leaf_species_by_model.get(model, set()))
        # 把 expanded/contracted 里出现过的物种也补进去，防止极端缺失
        for (m, sp) in expanded_fams_leaf.keys():
            if m == model:
                species_set.add(sp)
        for (m, sp) in contracted_fams_leaf.keys():
            if m == model:
                species_set.add(sp)

        for sp in sorted(species_set):
            key = (model, sp)
            fam_exp = expanded_fams_leaf.get(key, set())
            fam_con = contracted_fams_leaf.get(key, set())

            exp_genes: Set[str] = set()
            for fam in fam_exp:
                genes = member_map.get((fam, sp))
                if genes:
                    exp_genes.update(genes)

            con_genes: Set[str] = set()
            for fam in fam_con:
                genes = member_map.get((fam, sp))
                if genes:
                    con_genes.update(genes)

            sig_summary_rows.append([
                model,
                sp,
                len(fam_exp),
                len(fam_con),
                len(exp_genes),
                len(con_genes),
            ])

    sig_summary_path = cafe_agg_dir / "cafe_significant_clade_summary.tsv"
    write_tsv(
        sig_summary_path,
        ["model", "species", "n_expanded_fams", "n_contracted_fams", "n_expanded_genes", "n_contracted_genes"],
        sig_summary_rows,
    )
    logging.info(
        f"[OK] 写出 cafe_significant_clade_summary.tsv：{sig_summary_path}（{len(sig_summary_rows)} 行）"
    )
    # ==================== ★ 加回结束 ====================

    # ==================== ★ 新增：cafe_node_significant_clade_summary.tsv（补祖先节点） ====================
    # 说明：
    #   - 该表统计“所有节点”（祖先+叶），但只给 family 计数（不做祖先节点基因计数）
    #   - 列自定义，不影响旧表头约束
    node_summary_rows: List[List[object]] = []
    node_summary_header = [
        "model",
        "node_id",
        "taxon_label",
        "is_leaf",
        "species",
        "n_expanded_fams",
        "n_contracted_fams",
    ]

    for model in models:
        node_map = node_map_by_model.get(model, {}) or {}
        for nid in sorted(node_map.keys()):
            info = node_map[nid]
            key_node = (model, nid)
            n_up = len(expanded_fams_node.get(key_node, set()))
            n_down = len(contracted_fams_node.get(key_node, set()))
            node_summary_rows.append([
                model,
                nid,
                info.taxon_label,
                "TRUE" if info.is_leaf else "FALSE",
                "" if info.species_name is None else info.species_name,
                n_up,
                n_down,
            ])

    node_summary_path = cafe_agg_dir / "cafe_node_significant_clade_summary.tsv"
    write_tsv(node_summary_path, node_summary_header, node_summary_rows)
    logging.info(
        f"[OK] 写出 cafe_node_significant_clade_summary.tsv：{node_summary_path}（{len(node_summary_rows)} 行；含祖先+叶）"
    )
    # ==================== ★ 新增结束 ====================

    # 10) 输出根目录（脚本控制）
    enrich_out_root = cafe_agg_dir / "enrich_inputs"
    node_out_root = enrich_out_root / "node"
    mkdir_p(enrich_out_root)
    mkdir_p(node_out_root)
    logging.info(f"[out] leaf_enrich_root = {enrich_out_root}")
    logging.info(f"[out] node_og_root     = {node_out_root}")

    # 11) 叶节点富集发布包（所有叶物种）
    model = models[0]
    fam_all = fam_all_by_model.get(model, set())
    if not fam_all:
        logging.error(f"[ERR] model={model} 在 FDR 表中无 family 记录")
        sys.exit(2)

    leaf_species = sorted(leaf_species_by_model.get(model, set()))
    if not leaf_species:
        logging.error(f"[ERR] model={model} 未解析到任何叶节点物种（Base_clade_results.txt）")
        sys.exit(2)

    for sp in leaf_species:
        prefix = species_to_leaf_pkg_prefix(sp)
        pkg_name = f"{prefix}_cafe"  # 叶节点富集包：S_constricta_cafe
        pkg_dir = enrich_out_root / pkg_name
        mkdir_p(pkg_dir)

        bg_genes: Set[str] = set()
        for fam in fam_all:
            genes = member_map.get((fam, sp))
            if genes:
                bg_genes.update(genes)

        key_leaf = (model, sp)
        exp_fams = expanded_fams_leaf.get(key_leaf, set())
        con_fams = contracted_fams_leaf.get(key_leaf, set())

        up_genes: Set[str] = set()
        for fam in exp_fams:
            genes = member_map.get((fam, sp))
            if genes:
                up_genes.update(genes)

        down_genes: Set[str] = set()
        for fam in con_fams:
            genes = member_map.get((fam, sp))
            if genes:
                down_genes.update(genes)

        write_list_plain(pkg_dir / "background.list", bg_genes)
        write_list_plain(pkg_dir / "up.list", up_genes)
        write_list_plain(pkg_dir / "down.list", down_genes)

        meta_path = pkg_dir / "meta.tsv"
        n_detectable = len(bg_genes)
        universe_size = n_detectable
        with meta_path.open("w", encoding="utf-8") as f:
            f.write("label\tn_detectable\tn_annot_mapped\tuniverse_size\tdetectable_rule\tsamples_used\n")
            f.write(
                "\t".join([
                    pkg_name,
                    str(n_detectable),
                    "NA",
                    str(universe_size),
                    "CAFE_evaluable_gene_members_of_OGs_in_species_set",
                    sp,
                ]) + "\n"
            )

        logging.info(
            f"[leaf_enrich:{sp}] bg_genes={len(bg_genes)} up_genes={len(up_genes)} down_genes={len(down_genes)} -> {pkg_dir}"
        )

    # 12) 全节点OG发布包（祖先+叶节点）
    node_map = node_map_by_model.get(model, {})
    if not node_map:
        logging.error(f"[ERR] model={model} 未获取到 node_map")
        sys.exit(2)

    def node_pkg_name(info: NodeInfo) -> str:
        if info.is_leaf and info.species_name:
            prefix2 = species_to_leaf_pkg_prefix(info.species_name)
            return f"{prefix2}_node_cafe"  # 叶节点OG包：S_constricta_node_cafe
        return f"{info.node_id}node_cafe"  # 祖先节点OG包：40node_cafe

    og_bg_header = ["model", "node_id", "taxon_label", "is_leaf", "species", "family"]
    og_updown_header = ["model", "node_id", "taxon_label", "is_leaf", "species", "family", "Q_family", "p_branch_rule", "delta_rule"]
    node_meta_header = ["model", "node_id", "taxon_label", "is_leaf", "species", "pkg_name", "n_bg_og", "n_up_og", "n_down_og"]

    for nid in sorted(node_map.keys()):
        info = node_map[nid]
        pkg = node_pkg_name(info)
        pkg_dir = node_out_root / pkg
        mkdir_p(pkg_dir)

        key_node = (model, nid)
        up_ogs = expanded_fams_node.get(key_node, set())
        down_ogs = contracted_fams_node.get(key_node, set())

        bg_rows = []
        for fam in sorted(fam_all):
            bg_rows.append([
                model,
                nid,
                info.taxon_label,
                "TRUE" if info.is_leaf else "FALSE",
                "" if info.species_name is None else info.species_name,
                fam,
            ])
        write_tsv(pkg_dir / "og_background.tsv", og_bg_header, bg_rows)

        up_rows = []
        for fam in sorted(up_ogs):
            q = fam2q.get((model, fam))
            up_rows.append([
                model,
                nid,
                info.taxon_label,
                "TRUE" if info.is_leaf else "FALSE",
                "" if info.species_name is None else info.species_name,
                fam,
                "" if q is None else f"{q:.6g}",
                f"p_branch<= {branch_p_alpha} (or NA if no file)",
                "delta>0 & delta!=0",
            ])
        write_tsv(pkg_dir / "og_up.tsv", og_updown_header, up_rows)

        down_rows = []
        for fam in sorted(down_ogs):
            q = fam2q.get((model, fam))
            down_rows.append([
                model,
                nid,
                info.taxon_label,
                "TRUE" if info.is_leaf else "FALSE",
                "" if info.species_name is None else info.species_name,
                fam,
                "" if q is None else f"{q:.6g}",
                f"p_branch<= {branch_p_alpha} (or NA if no file)",
                "delta<0 & delta!=0",
            ])
        write_tsv(pkg_dir / "og_down.tsv", og_updown_header, down_rows)

        meta_rows = [[
            model,
            nid,
            info.taxon_label,
            "TRUE" if info.is_leaf else "FALSE",
            "" if info.species_name is None else info.species_name,
            pkg,
            len(fam_all),
            len(up_ogs),
            len(down_ogs),
        ]]
        write_tsv(pkg_dir / "meta.tsv", node_meta_header, meta_rows)

        logging.info(
            f"[node_og:{info.taxon_label}] bg_og={len(fam_all)} up_og={len(up_ogs)} down_og={len(down_ogs)} -> {pkg_dir}"
        )

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

