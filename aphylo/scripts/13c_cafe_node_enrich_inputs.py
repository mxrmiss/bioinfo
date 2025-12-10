#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13c_cafe_node_enrich_inputs.py —— 按 CAFE 节点构建富集外包（节点版）

设计目标：
  - 皇上在 config.yaml 中指定一个或多个 CAFE 节点（node_id），
    以及每个节点对应的“代表物种”（用于取该物种的基因做富集）。
  - 本脚本读取：
       1) 13 步的 FDR 表：cafe_significant_families_no_highfail.tsv
       2) CAFE primary_global 下的：
            Base_clade_results.txt
            Base_change.tab
            Base_branch_probabilities.tab（可选）
       3) inputs.all_members_tsv（<publish_dir>/all_pep2cds_resolved.tsv）
  - 对每个节点（例如 node_id=30）：
       1) 识别在该节点上显著扩张 / 收缩的家族（OG）：
            - family 的 Q_family <= report.fdr_alpha
            - Base_change.tab 中该节点的 Δ != 0
            - 若存在 Base_branch_probabilities.tab，则该节点 p_branch <= cafe5.branch_p_alpha
       2) 使用 represent_species（例如 Sinonovacula_constricta）
          将 family 映射为基因集合：
            - 背景基因 = 所有 fam_all × represent_species 的基因并集
            - 扩张基因 = expanded_fams_node × represent_species 的基因并集
            - 收缩基因 = contracted_fams_node × represent_species 的基因并集
       3) 为每个节点输出一套富集外包：
            <outputs_dir>/<tag>/background.list
            <outputs_dir>/<tag>/up.list
            <outputs_dir>/<tag>/down.list
            <outputs_dir>/<tag>/meta.tsv

  - 额外输出：
       cafe_node_enrich_summary.tsv（写在 paths.cafe_agg_dir）：
         model  tag  node_id  taxon_label  represent_species
         n_expanded_fams  n_contracted_fams  n_expanded_genes  n_contracted_genes
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

# ======================== 顶部参数区 ========================

# 配置文件路径（在脚本顶部统一配置）
CONFIG_PATH = "config.yaml"

# 日志文件名（写在 paths.logs_dir 下）
LOG_FILE_BASENAME = "13c_cafe_node_enrich_inputs.log"

# 日志级别
LOG_LEVEL = "INFO"

# 当 Base_branch_probabilities.tab 缺失时，是否允许退化为仅使用 Q + Δ 判定分支显著性
ALLOW_NO_BRANCH_PROB = True


# ==================== 基础工具 & 配置读取 ====================

def mkdir_p(p: Path) -> Path:
    """类似 mkdir -p，若目录存在则不报错。"""
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


# ======================= CAFE 结果解析 =======================

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
      Sinonovacula_constricta<10>     541 441
      Cyclina_sinensis<9>             1587    519

    解析规则：
      - "<22>" 视为内部节点，species_name=None, is_leaf=False
      - "Name<10>" 视为叶节点，species_name="Name", is_leaf=True
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_clade_results.txt 为空：{path}")
    node_map: Dict[int, NodeInfo] = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if not parts:
            continue
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

    适配矩阵格式（当前皇上的 CAFE5 输出）：
      表头：
        FamilyID  Argopecten_irradians<1> ... Sinonovacula_rivularis<16> Lottia_gigantea<17> <18> ... <33>
      每一行：
        OG0000015  29  -4  ...  -7

      解析逻辑：
        - 第 1 列为 family；
        - 后续每一列名都可以用 parse_node_label_to_id() 解析成 node_id；
        - 单元格为 Δ（change），NA/N/A/./空 跳过。
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

    if not is_matrix:
        raise ValueError("[ERR] 暂未实现非矩阵格式的 Base_change.tab 解析，请检查文件格式")

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

    return out


def parse_branch_probabilities(path: Path) -> Dict[Tuple[str, int], float]:
    """
    解析 Base_branch_probabilities.tab，提取 (family, node_id) -> p_branch 映射。

    适配矩阵格式：
      表头：
        FamilyID  Argopecten_irradians<1> ... <33>
      每一行：
        OG0000015  0.1  0.03  ... NA

      解析逻辑：
        - 第 1 列为 family；
        - 后续列名中带 <id> 的解析为 node_id；
        - 单元格为 p_branch，NA/N/A/./空 跳过。
    """
    lines = read_lines(path)
    if not lines:
        raise ValueError(f"[ERR] Base_branch_probabilities.tab 为空：{path}")
    header = re.split(r"\s+", lines[0].strip())
    data_lines = lines[1:]

    out: Dict[Tuple[str, int], float] = {}

    # 判定是否矩阵格式
    is_matrix = False
    if header and header[0].lower().startswith("family"):
        for col in header[1:]:
            if parse_node_label_to_id(col) is not None:
                is_matrix = True
                break

    if not is_matrix:
        raise ValueError("[ERR] 暂未实现非矩阵格式的 Base_branch_probabilities.tab 解析，请检查文件格式")

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


# ===================== 主逻辑：节点富集外包 =====================

def main() -> None:
    # 1) 读取配置
    cfg_path = Path(CONFIG_PATH).resolve()
    cfg = load_config(str(cfg_path))

    paths = cfg.get("paths", {}) or {}
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
    logging.info("========== APhylo 13c — CAFE 节点富集外包（node_bridge） ==========")
    logging.info(f"[init] 使用配置文件：{cfg_path}")
    logging.info(f"[init] cafe_run_dir      = {cafe_run_dir}")
    logging.info(f"[init] cafe_agg_dir      = {cafe_agg_dir}")
    logging.info(f"[init] logs_dir          = {logs_dir}")

    # 3) 检查 node_bridge 是否启用
    cafe_cfg = cfg.get("cafe5", {}) or {}
    node_cfg = cafe_cfg.get("node_bridge", {}) or {}
    if not node_cfg.get("enable", False):
        logging.info("[init] cafe5.node_bridge.enable = false，本脚本不执行任何操作，直接退出。")
        return

    outputs_dir_cfg = node_cfg.get("outputs_dir", "results/07_cafeagg/enrich_inputs_nodes")
    node_out_root = Path(outputs_dir_cfg).resolve()
    mkdir_p(node_out_root)
    logging.info(f"[init] node_enrich_outputs_dir = {node_out_root}")

    sets_cfg = node_cfg.get("sets", {}) or {}
    if not isinstance(sets_cfg, dict) or not sets_cfg:
        logging.error("[ERR] cafe5.node_bridge.sets 未配置或为空，无法生成节点富集输入")
        sys.exit(2)

    # 4) 前置检查：13 步产物是否存在
    sig_nohf_path = cafe_agg_dir / "cafe_significant_families_no_highfail.tsv"
    if not sig_nohf_path.exists():
        logging.error(f"[ERR] 未找到 13 步产物：{sig_nohf_path}")
        sys.exit(2)

    # 5) FDR 与分支 p 值阈值
    report_cfg = cfg.get("report", {}) or {}
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

    total_fam = sum(len(v) for v in fam_all_by_model.values())
    logging.info(f"[fdr] 总 family 条目数：{total_fam}")

    if not fam_all_by_model:
        logging.error("[ERR] 未从 FDR 表中解析到任何 family 记录")
        sys.exit(2)

    # 7) 决定要使用的 CAFE 模型：优先 report.selected_cafe_model，否则用 cafe5.models[0]
    cafe_models = cafe_cfg.get("models", ["global"])
    selected_model = report_cfg.get("selected_cafe_model")
    if not selected_model:
        selected_model = cafe_models[0]
        logging.warning(f"[WARN] report.selected_cafe_model 未设置，将使用 cafe5.models[0]：{selected_model}")
    else:
        logging.info(f"[init] selected_cafe_model = {selected_model}")

    model = selected_model
    fam_all = fam_all_by_model.get(model, set())
    if not fam_all:
        logging.error(f"[ERR] model={model} 在 FDR 表中无 family 记录")
        sys.exit(2)

    models_dir = cafe_run_dir / "models"
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

    # 8) 成员映射：使用 inputs.all_members_tsv （all_pep2cds_resolved.tsv）
    inputs_cfg = cfg.get("inputs", {}) or {}
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

    # species.alias_map：用于统一物种命名（对 represent_species 生效）
    alias_map = (cfg.get("species", {}) or {}).get("alias_map", {}) or {}

    # 记录节点汇总结果
    summary_rows: List[List[object]] = []

    # 9) 遍历每个 node_bridge.sets，构建该节点的扩/缩外包
    for tag, info in sets_cfg.items():
        if not isinstance(info, dict):
            logging.warning(f"[WARN] node_bridge.sets.{tag} 不是字典，跳过")
            continue
        if "node_id" not in info or "represent_species" not in info:
            logging.warning(f"[WARN] node_bridge.sets.{tag} 缺少 node_id 或 represent_species，跳过")
            continue

        node_id_raw = info.get("node_id")
        try:
            node_id = int(node_id_raw)
        except Exception:
            logging.error(f"[ERR] tag={tag} 的 node_id 非法：{node_id_raw}")
            sys.exit(2)

        rep_sp_raw = str(info.get("represent_species")).strip()
        if not rep_sp_raw:
            logging.error(f"[ERR] tag={tag} 的 represent_species 为空")
            sys.exit(2)
        rep_sp = alias_map.get(rep_sp_raw, rep_sp_raw)

        node_info = node_map.get(node_id)
        if node_info is None:
            logging.error(f"[ERR] tag={tag} 指定的 node_id={node_id} 不存在于 Base_clade_results.txt 中")
            sys.exit(2)

        logging.info(
            f"[tag={tag}] 节点 node_id={node_id}, taxon_label={node_info.taxon_label}, "
            f"represent_species={rep_sp_raw} -> {rep_sp}"
        )

        # 该节点的显著扩张 / 收缩 family 集合
        expanded_fams: Set[str] = set()
        contracted_fams: Set[str] = set()

        # 遍历所有 change_records，挑出该 node_id 的条目
        # change_records: List[(family, node_id, delta)]
        for fam, nid, delta in change_records:
            if nid != node_id:
                continue
            # fam 必须存在于该 model 的 fam_all 中
            if fam_all and fam not in fam_all:
                continue
            # family 必须先在全局 FDR 上显著
            q = fam2q.get((model, fam))
            if q is None or q > fdr_alpha:
                continue
            if delta == 0:
                continue

            p_branch = None
            if has_branch_prob:
                p_branch = branch_prob.get((fam, node_id))
                # 没有记录一般说明 parent=child，视为不显著
                if p_branch is None or p_branch > branch_p_alpha:
                    continue

            if delta > 0:
                expanded_fams.add(fam)
            else:
                contracted_fams.add(fam)

        logging.info(
            f"[tag={tag}] 节点 {node_id} 上显著扩张 family={len(expanded_fams)}，显著收缩 family={len(contracted_fams)}"
        )

        # 背景基因：fam_all × represent_species
        bg_genes: Set[str] = set()
        for fam in fam_all:
            genes = member_map.get((fam, rep_sp))
            if genes:
                bg_genes.update(genes)

        # 扩张基因
        up_genes: Set[str] = set()
        for fam in expanded_fams:
            genes = member_map.get((fam, rep_sp))
            if genes:
                up_genes.update(genes)

        # 收缩基因
        down_genes: Set[str] = set()
        for fam in contracted_fams:
            genes = member_map.get((fam, rep_sp))
            if genes:
                down_genes.update(genes)

        logging.info(
            f"[tag={tag}] 背景基因={len(bg_genes)}, 扩张基因={len(up_genes)}, 收缩基因={len(down_genes)}"
        )

        # 写出到 <outputs_dir>/<tag>/
        tag_dir = node_out_root / tag
        mkdir_p(tag_dir)

        def write_list(path: Path, genes: Set[str]) -> None:
            with path.open("w", encoding="utf-8") as f:
                for g in sorted(genes):
                    f.write(g + "\n")

        write_list(tag_dir / "background.list", bg_genes)
        write_list(tag_dir / "up.list", up_genes)
        write_list(tag_dir / "down.list", down_genes)

        # meta.tsv：格式与 13b 保持一致
        meta_path = tag_dir / "meta.tsv"
        n_detectable = len(bg_genes)
        universe_size = n_detectable
        detectable_rule = f"CAFE_node_{node_id}_represent_species_{rep_sp}"
        with meta_path.open("w", encoding="utf-8") as f:
            f.write("label\tn_detectable\tn_annot_mapped\tuniverse_size\tdetectable_rule\tsamples_used\n")
            f.write(
                "\t".join([
                    tag,
                    str(n_detectable),
                    "NA",
                    str(universe_size),
                    detectable_rule,
                    rep_sp,
                ]) + "\n"
            )

        logging.info(f"[tag={tag}] 写出节点外包 up/down/background/meta 到：{tag_dir}")

        # 汇总表记录
        summary_rows.append([
            model,
            tag,
            node_id,
            node_info.taxon_label,
            rep_sp,
            len(expanded_fams),
            len(contracted_fams),
            len(up_genes),
            len(down_genes),
        ])

    # 10) 写出节点汇总表
    summary_path = cafe_agg_dir / "cafe_node_enrich_summary.tsv"
    write_tsv(
        summary_path,
        [
            "model",
            "tag",
            "node_id",
            "taxon_label",
            "represent_species",
            "n_expanded_fams",
            "n_contracted_fams",
            "n_expanded_genes",
            "n_contracted_genes",
        ],
        summary_rows,
    )
    logging.info(f"[OK] 写出节点汇总表：{summary_path}（{len(summary_rows)} 行）")

    logging.info("========== APhylo 13c — 完成 ==========")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

