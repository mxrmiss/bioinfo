#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13b_cafe_build_enrich_inputs.py —— 为 CAFE5 per-species 扩张/收缩构建富集输入（适配矩阵版 Base_* 文件）

核心设计（与皇上的 config.yaml + CAFE5 输出对齐）：
  1）基于 13 步的 FDR 表（cafe_significant_families_no_highfail.tsv）确定全局显著家族：
        - 使用 report.fdr_alpha 作为 Q 阈值；
        - 每个 model 拿到 fam_all（用于背景）与 fam_sig（用于候选扩缩）。
  2）从 primary_global 下读取：
        - Base_clade_results.txt           —— 提取 node_id ↔ 物种名 ↔ 是否叶节点；
        - Base_change.tab                  —— 每个 (family, node) 的 Δ（child-parent），方向信息；
        - Base_branch_probabilities.tab    —— 每个 (family, node) 的分支显著性概率/ p 值。
     ★ 适配两种格式：
        a) 传统“长表”格式：FamilyID  Node  Probability
        b) 新版“矩阵”格式：FamilyID  Branchiostrea<1>  Homo_sapiens<2>  …  <29>
  3）仅在满足以下条件时，认定“某个 family 在某个物种叶节点上发生了显著扩/缩”：
        a. family 的 Q_family <= report.fdr_alpha（全局 FDR 显著）；
        b. 若存在 Base_branch_probabilities.tab，则该分支 p_branch <= cafe5.branch_p_alpha（默认 0.05）；
           若 branch_probabilities 缺失且 ALLOW_NO_BRANCH_PROB=True，则仅依赖 Q + Δ；
        c. Base_change.tab 中该 (family, node) 的 Δ != 0。
     并据 Δ 的正负将 status 标记为 expanded / contracted。
  4）在 cafe_agg_dir 写出：
        - cafe_family_branch_status.tsv：记录所有叶节点的 per-family 扩缩状态。
        - cafe_significant_clade_summary.tsv：
              按物种汇总“显著扩张/收缩家族数 + 显著扩张/收缩基因数”（★ 给树图与结果段使用）。
  5）利用 inputs.all_members_tsv（当前为 "<publish_dir>/all_pep2cds_resolved.tsv"）：
        - 通过 (OG, Species) -> {cds_id} 映射，把 family 状态折叠成 per-species 基因集合：
            * 背景基因 = 该物种在 CAFE 可评估家族中的所有基因；
            * 扩张基因 = 所有 expanded family × 该物种 的基因并集；
            * 收缩基因 = 所有 contracted family × 该物种 的基因并集。
  6）结合 config.cafe5.enrich_bridge：
        cafe5:
          enrich_bridge:
            enable: true
            outputs_dir: "results/07_cafeagg/enrich_inputs"
            species_sets:
              cafe_sc:
                - "Sinonovacula_constricta"
              cafe_sr:
                - "Sinonovacula_rivularis"
        - 对每个 tag（例如 cafe_sc / cafe_sr），按 species 列表构建：
            * <outputs_dir>/<tag>/background.list
            * <outputs_dir>/<tag>/up.list     （扩张）
            * <outputs_dir>/<tag>/down.list   （收缩）
            * <outputs_dir>/<tag>/meta.tsv
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

# 配置文件路径（皇上统一在脚本顶部配置）
CONFIG_PATH = "config.yaml"

# 日志文件名（写在 paths.logs_dir 下）
LOG_FILE_BASENAME = "13b_cafe_build_enrich_inputs.log"

# 日志级别
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
    """
    将对象中的 <publish_dir> 占位符展开：
      - 字符串：直接替换
      - 列表/字典：递归展开
      - 其他类型：原样返回
    """
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
    """存放节点与物种的映射信息。"""
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
      - 其他             → None
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

    文件格式示例：
      #Taxon_ID       Increase        Decrease
      <22>    125     33
      Sinonovacula_constricta<10>     541     441
      Cyclina_sinensis<9>             1587    519
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
        # 内部节点：<22>
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
        # 叶节点：Name<10>
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

    适配两种格式：

    1）矩阵格式（当前皇上的 CAFE5 输出）：
        表头：
          FamilyID  Branchiostrea_floridae<1>  Homo_sapiens<2>  ...  <29>
        每一行：
          OG0000001  0.1  0.0  ...  NA

        行解析逻辑：
          - 第 1 列为 family；
          - 后续每一列名都可以用 parse_node_label_to_id() 解析成 node_id；
          - 单元格为 Δ（change），NA/N/A/./空 跳过。

    2）旧版“长表”格式（少数脚本使用）：
        FamilyID  Parent  Node  ParentCount  Change  ...
        - 第1列：family
        - 第3列：node_id
        - 第5列：Change
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_change.tab 为空：{path}")
    header = re.split(r"\s+", lines[0].strip())
    data_lines = lines[1:]

    # 判断是否为矩阵格式：首列 FamilyID，且第二列及之后的列名中有 "<id>" 结构
    is_matrix = False
    if header and header[0].lower().startswith("family"):
        for col in header[1:]:
            if parse_node_label_to_id(col) is not None:
                is_matrix = True
                break

    out: List[Tuple[str, int, float]] = []

    if is_matrix:
        # 矩阵模式
        # 先解析每一列对应的 node_id
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
        # 退化为旧版长表解析
        data_lines = lines
        first = lines[0].strip()
        if first.startswith("FamilyID") or first.startswith("#FamilyID"):
            data_lines = lines[1:]
        for ln in data_lines:
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

    适配两种格式：

    1）矩阵格式（当前皇上的 CAFE5 输出）：
        表头：
          FamilyID  Branchiostrea_floridae<1>  Homo_sapiens<2>  ...  <29>
        每一行：
          OG0000001  0.01  0.23  ...  NA

        行解析逻辑：
          - 第 1 列为 family；
          - 后续每一列名都可以用 parse_node_label_to_id() 解析成 node_id；
          - 单元格为 p_branch，NA/N/A/./空 跳过。

    2）旧版“长表”格式：
        FamilyID  Node  Probability
        - 自动根据列名查找 family/node/prob 列。
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_branch_probabilities.tab 为空：{path}")
    header = re.split(r"\s+", lines[0].strip())
    hlow = [h.lower() for h in header]
    data_lines = lines[1:]

    out: Dict[Tuple[str, int], float] = {}

    # 尝试判定是否矩阵格式
    is_matrix = False
    if header and header[0].lower().startswith("family"):
        for col in header[1:]:
            if parse_node_label_to_id(col) is not None:
                is_matrix = True
                break

    if is_matrix:
        # 矩阵模式
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

    # 否则走旧版长表格式：识别 family/node/prob 列
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
        raise ValueError(
            f"[ERR] 无法从表头识别 family/node/prob 列：{header}"
        )

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


# ----------------------- 主逻辑：生成富集输入 -----------------------


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
    logging.info("========== APhylo 13b — CAFE5 per-species 扩张/收缩富集输入 ==========")
    logging.info(f"[init] 使用配置文件：{cfg_path}")
    logging.info(f"[init] cafe_run_dir      = {cafe_run_dir}")
    logging.info(f"[init] cafe_agg_dir      = {cafe_agg_dir}")
    logging.info(f"[init] logs_dir          = {logs_dir}")

    # 3) 检查 enrich_bridge 是否启用
    cafe_cfg = cfg.get("cafe5", {})
    enrich_cfg = cafe_cfg.get("enrich_bridge", {}) or {}
    if not enrich_cfg.get("enable", False):
        logging.info("[init] cafe5.enrich_bridge.enable = false，本脚本不执行任何操作，直接退出。")
        return

    outputs_dir_cfg = enrich_cfg.get("outputs_dir", "results/07_cafeagg/enrich_inputs")
    enrich_out_root = Path(outputs_dir_cfg).resolve()
    mkdir_p(enrich_out_root)
    logging.info(f"[init] enrich_outputs_dir = {enrich_out_root}")

    species_sets_cfg = enrich_cfg.get("species_sets", {})
    if not species_sets_cfg:
        logging.error("[ERR] cafe5.enrich_bridge.species_sets 未配置，无法生成 CAFE 富集输入")
        sys.exit(2)

    # 4) 前置检查：13 步产物是否存在
    sig_nohf_path = cafe_agg_dir / "cafe_significant_families_no_highfail.tsv"
    if not sig_nohf_path.exists():
        logging.error(f"[ERR] 未找到 13 步产物：{sig_nohf_path}")
        sys.exit(2)

    sentinel = cafe_run_dir / ".cafe.done"
    if not sentinel.exists():
        logging.warning(f"[WARN] 未发现 12 步完成哨兵：{sentinel}，请确认 CAFE 已完整运行。")

    # 5) FDR 与分支 p 值阈值
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

    # 6) 读取 13 的 FDR 表，构建 family → Q 映射
    fam2q: Dict[Tuple[str, str], float] = {}      # (model, family) -> Q
    fam_all_by_model: Dict[str, Set[str]] = {}    # model -> 所有 family
    fam_sig_by_model: Dict[str, Set[str]] = {}    # model -> Q<=α family

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
        if not fam:
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

    total_fam = sum(len(v) for v in fam_all_by_model.values())
    total_sig = sum(len(v) for v in fam_sig_by_model.values())
    logging.info(f"[cafe] 总 family 条目数：{total_fam}，其中 Q<=α 的显著 family 数：{total_sig}")

    if not fam_all_by_model:
        logging.error("[ERR] 未从 FDR 表中解析到任何 family 记录")
        sys.exit(2)

    # 7) 决定要使用的 CAFE 模型：优先 report.selected_cafe_model，否则用 cafe5.models[0]
    selected_model = report_cfg.get("selected_cafe_model")
    cafe_models = cafe_cfg.get("models", ["global"])
    if not selected_model:
        selected_model = cafe_models[0]
        logging.warning(f"[WARN] report.selected_cafe_model 未设置，将使用 cafe5.models[0]：{selected_model}")
    else:
        logging.info(f"[init] selected_cafe_model = {selected_model}")
    models = [selected_model]

    models_dir = cafe_run_dir / "models"

    # per-branch 状态记录
    branch_status_rows: List[List[object]] = []
    # (model, species) -> 显著扩张 / 收缩 family 集合
    expanded_fams: Dict[Tuple[str, str], Set[str]] = {}
    contracted_fams: Dict[Tuple[str, str], Set[str]] = {}
    # 记录每个 model 下所有叶节点的物种名，用于补 0
    leaf_species_by_model: Dict[str, Set[str]] = {}

    # 8) 遍历模型，解析节点信息、Δ 与分支概率，填充 expanded_fams / contracted_fams
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
        change_records = parse_change_table(change_path)
        logging.info(f"[model={model}] 从 {change_path} 解析到 {len(change_records)} 条 Δ 记录")

        # 收集该 model 下所有叶节点物种名（用于后续显著家族 / 基因统计）
        leaf_species: Set[str] = set()
        for info in node_map.values():
            if info.is_leaf and info.species_name:
                leaf_species.add(info.species_name)
        leaf_species_by_model[model] = leaf_species
        logging.info(f"[model={model}] 叶节点物种数：{len(leaf_species)}")

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
            # fam 必须存在于 13 步 FDR 表（fam_all），否则跳过
            if fam_all and fam not in fam_all:
                continue
            # family 必须先在全局 FDR 上显著
            q = fam2q.get((model, fam))
            if q is None or q > fdr_alpha:
                continue

            node_info = node_map.get(node_id)
            if node_info is None:
                continue
            if not node_info.is_leaf:
                # 仅 per-species：内部节点忽略
                continue
            species = node_info.species_name
            if not species:
                continue

            if delta == 0:
                # Δ=0 无扩缩
                continue

            p_branch = None
            if has_branch_prob:
                p_branch = branch_prob.get((fam, node_id))
                if p_branch is None:
                    # 没有记录一般说明 parent=child，视为不显著
                    continue
                if p_branch > branch_p_alpha:
                    continue

            status = "expanded" if delta > 0 else "contracted"

            branch_status_rows.append([
                model,
                fam,
                node_id,
                node_info.taxon_label,
                species,
                status,
                f"{delta:.6g}",
                "" if p_branch is None else f"{p_branch:.6g}",
                f"{q:.6g}",
            ])

            key = (model, species)
            if status == "expanded":
                expanded_fams.setdefault(key, set()).add(fam)
            else:
                contracted_fams.setdefault(key, set()).add(fam)

    # 写出 per-branch 状态表
    status_path = cafe_agg_dir / "cafe_family_branch_status.tsv"
    write_tsv(
        status_path,
        ["model", "family", "node_id", "taxon_label", "species", "status", "delta", "p_branch", "Q_family"],
        branch_status_rows,
    )
    logging.info(f"[OK] 写出 per-branch 状态表：{status_path}（{len(branch_status_rows)} 行）")

    # 9) 成员映射：使用 inputs.all_members_tsv （all_pep2cds_resolved.tsv）
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

    # member_map[(OG, Species)] -> {cds_id}
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
        key = (og, sp)
        member_map.setdefault(key, set()).add(cds)
        total_records += 1

    logging.info(
        f"[members] 从 {all_members_path} 读取到 {total_records} 条 (OG,Species,cds) 映射，键数：{len(member_map)}"
    )

    # ==================== ★ 新增：显著扩/缩 家族数 + 基因数 汇总 ====================
    # 目标：按物种统计“显著扩张/收缩家族数 + 显著扩张/收缩基因数”，
    #       写出到 cafe_significant_clade_summary.tsv，供系统发育树和结果段使用。
    sig_summary_rows: List[List[object]] = []

    for model in models:
        # 该模型下所有叶节点物种名（来自 Base_clade_results）
        species_set: Set[str] = set(leaf_species_by_model.get(model, set()))
        # 把 expanded_fams / contracted_fams 里出现过的物种也补进去，防止极端缺失
        for (m, sp) in expanded_fams.keys():
            if m == model:
                species_set.add(sp)
        for (m, sp) in contracted_fams.keys():
            if m == model:
                species_set.add(sp)

        for sp in sorted(species_set):
            key = (model, sp)
            fam_exp = expanded_fams.get(key, set())
            fam_con = contracted_fams.get(key, set())

            # 显著扩张基因集合：所有显著扩张 family 在该物种上的成员基因并集
            exp_genes: Set[str] = set()
            for fam in fam_exp:
                genes = member_map.get((fam, sp))
                if genes:
                    exp_genes.update(genes)

            # 显著收缩基因集合
            con_genes: Set[str] = set()
            for fam in fam_con:
                genes = member_map.get((fam, sp))
                if genes:
                    con_genes.update(genes)

            n_exp_fams = len(fam_exp)
            n_con_fams = len(fam_con)
            n_exp_genes = len(exp_genes)
            n_con_genes = len(con_genes)

            sig_summary_rows.append([
                model,
                sp,
                n_exp_fams,
                n_con_fams,
                n_exp_genes,
                n_con_genes,
            ])

    sig_summary_path = cafe_agg_dir / "cafe_significant_clade_summary.tsv"
    write_tsv(
        sig_summary_path,
        ["model", "species", "n_expanded_fams", "n_contracted_fams", "n_expanded_genes", "n_contracted_genes"],
        sig_summary_rows,
    )
    logging.info(
        f"[OK] 写出显著扩/缩 家族 + 基因 汇总表：{sig_summary_path}（{len(sig_summary_rows)} 行；"
        f"列为 model, species, n_expanded_fams, n_contracted_fams, n_expanded_genes, n_contracted_genes）"
    )
    # ==================== ★ 新增部分结束 ====================

    # 10) species.alias_map：用于统一物种命名（对 enrich_inputs 的 tag 生效）
    alias_map = (cfg.get("species", {}) or {}).get("alias_map", {}) or {}

    # 11) 针对每个 tag 构建 up/down/background/meta（保持原有逻辑不变）
    model = models[0]
    fam_all = fam_all_by_model.get(model, set())
    if not fam_all:
        logging.error(f"[ERR] model={model} 在 FDR 表中无 family 记录")
        sys.exit(2)

    for tag, sp_list in species_sets_cfg.items():
        # config 中 species_sets 的值是物种列表
        if not isinstance(sp_list, list) or len(sp_list) == 0:
            logging.warning(f"[WARN] tag={tag} 的 species_sets 不是非空列表，跳过")
            continue

        sp_list_resolved: List[str] = []
        for sp in sp_list:
            sp_norm = alias_map.get(sp, sp)
            sp_list_resolved.append(sp_norm)

        logging.info(f"[tag={tag}] 物种集合：{sp_list_resolved}")

        # 背景基因：所有 fam_all × species_set 的并集
        bg_genes: Set[str] = set()
        for fam in fam_all:
            for sp in sp_list_resolved:
                key = (fam, sp)
                genes = member_map.get(key)
                if genes:
                    bg_genes.update(genes)

        # 扩张 / 收缩 family：根据 per-branch 状态表汇总
        exp_fams_tag: Set[str] = set()
        con_fams_tag: Set[str] = set()
        for sp in sp_list_resolved:
            key = (model, sp)
            exp_fams_tag.update(expanded_fams.get(key, set()))
            con_fams_tag.update(contracted_fams.get(key, set()))

        # 扩张基因集合
        up_genes: Set[str] = set()
        for fam in exp_fams_tag:
            for sp in sp_list_resolved:
                key = (fam, sp)
                genes = member_map.get(key)
                if genes:
                    up_genes.update(genes)

        # 收缩基因集合
        down_genes: Set[str] = set()
        for fam in con_fams_tag:
            for sp in sp_list_resolved:
                key = (fam, sp)
                genes = member_map.get(key)
                if genes:
                    down_genes.update(genes)

        logging.info(
            f"[tag={tag}] 背景 OG={len(fam_all)}, 背景基因={len(bg_genes)}, "
            f"扩张 OG={len(exp_fams_tag)}, 扩张基因={len(up_genes)}, "
            f"收缩 OG={len(con_fams_tag)}, 收缩基因={len(down_genes)}"
        )

        tag_dir = enrich_out_root / tag
        mkdir_p(tag_dir)

        def write_list(path: Path, genes: Set[str]) -> None:
            with path.open("w", encoding="utf-8") as f:
                for g in sorted(genes):
                    f.write(g + "\n")

        write_list(tag_dir / "background.list", bg_genes)
        write_list(tag_dir / "up.list", up_genes)
        write_list(tag_dir / "down.list", down_genes)

        # meta.tsv：按 vNext 约定
        meta_path = tag_dir / "meta.tsv"
        n_detectable = len(bg_genes)
        universe_size = n_detectable
        with meta_path.open("w", encoding="utf-8") as f:
            f.write("label\tn_detectable\tn_annot_mapped\tuniverse_size\tdetectable_rule\tsamples_used\n")
            f.write(
                "\t".join([
                    tag,
                    str(n_detectable),
                    "NA",
                    str(universe_size),
                    "CAFE_evaluable_gene_members_of_OGs_in_species_set",
                    ";".join(sp_list_resolved),
                ]) + "\n"
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

