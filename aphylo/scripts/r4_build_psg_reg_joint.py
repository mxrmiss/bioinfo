#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
r4_build_psg_reg_joint.py

功能：
1）读取现有 PSG 主线结果表
2）读取 r3_make_reg_table.py 生成的 REG_plot_input.tsv
3）额外读取 07b 生成的高可信 PSG 表（BEB>=0.95）
4）按 OG + foreground 合并，生成 PSG × REG 联合表
5）输出：
   results/05_sbranch_model/summary/
     - PSG_REG_joint.tsv
     - PSG_REG_counts.tsv
     - PSG_REG_labels.tsv
     - .r4.done

说明：
- 本脚本默认优先使用 PSG 的 Q 值来判定 gene-level 显著性
- 可通过 config.yaml: branch_model.psg_alpha_mode 改为用 P 值
- 最终 psg_sig = gene-level 显著 AND BEB>=0.95 支持
- joint_class 分四类：
    None
    PSG_only
    REG_only
    Both
- x_psg_p = -log10(psg_p)
- y_reg_fdr 直接取 REG_plot_input.tsv 中的 y_reg_fdr

注意：
- 本版严格保持原版 r4 输出表头不变，不新增任何列
- 只改变 psg_sig 和 joint_class 的判定逻辑

运行方式：
  python scripts/r4_build_psg_reg_joint.py
"""

from __future__ import annotations

import io
import sys
import math
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set


# ==============================
# 配置区（皇上主要看这里）
# ==============================

PROJECT_ROOT = Path(".").resolve()
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

REG_INPUT_FILENAME = "REG_plot_input.tsv"
OUT_JOINT = "PSG_REG_joint.tsv"
OUT_COUNTS = "PSG_REG_counts.tsv"
OUT_LABELS = "PSG_REG_labels.tsv"
DONE_FILENAME = ".r4.done"

# PSG 主表：必须继续使用全量 gene-level 结果表
DEFAULT_PSG_TABLE = PROJECT_ROOT / "results/05_cmlagg/D_fdr_genes.tsv"

# 07b 输出的高可信 PSG 表：仅用来判断某个 OG×FG 是否满足 BEB>=0.95
DEFAULT_BEB_SIG_TABLE = PROJECT_ROOT / "results/05_cmlagg/PSG_qsig_bebsig_OGxFG.tsv"


# ==============================
# 配置读取
# ==============================

def _expand_publish_placeholders(obj, publish_dir: str):
    """递归替换 <publish_dir> 占位符。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj


def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml。"""
    if not path.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{path}")

    try:
        import yaml
    except Exception as e:
        raise RuntimeError(f"[ERR] 读取 YAML 需要 pyyaml，但导入失败：{e}")

    cfg = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg


# ==============================
# 基础工具
# ==============================

def ensure_dir(p: Path) -> Path:
    """确保目录存在。"""
    p.mkdir(parents=True, exist_ok=True)
    return p


def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """创建日志器，同时输出到屏幕与文件。"""
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name)
    lg.setLevel(level)
    lg.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(level)

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(level)

    lg.addHandler(fh)
    lg.addHandler(sh)

    class _Flush(io.TextIOBase):
        def __init__(self, s):
            self.s = s

        def write(self, x):
            self.s.write(x)
            self.s.flush()
            return len(x)

    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)
    return lg


def banner(log: logging.Logger, text: str):
    """打印横幅。"""
    bar = "=" * max(10, len(text) + 2)
    log.info(bar)
    log.info(f" {text} ")
    log.info(bar)


def parse_bool(x: Any) -> bool:
    """把各种真假值转为布尔。"""
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    return s in {"1", "true", "t", "yes", "y"}


def safe_float(x: Any, default: float = float("nan")) -> float:
    """安全转 float。"""
    try:
        return float(x)
    except Exception:
        return default


def split_semicolon(s: str) -> List[str]:
    """按分号拆分字符串，去空。"""
    if s is None:
        return []
    return [x.strip() for x in str(s).split(";") if x.strip()]


def write_tsv(path: Path, header: List[str], rows: List[Dict[str, Any]]) -> None:
    """写 TSV。"""
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as w:
        w.write("\t".join(header) + "\n")
        for row in rows:
            w.write("\t".join(str(row.get(k, "")) for k in header) + "\n")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    """读取 TSV 到 dict 列表。"""
    if not path.is_file():
        raise FileNotFoundError(f"[ERR] 缺少输入表：{path}")

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines:
        raise ValueError(f"[ERR] 输入表为空：{path}")

    header = lines[0].rstrip("\n").split("\t")
    rows: List[Dict[str, str]] = []

    for raw in lines[1:]:
        parts = raw.rstrip("\n").split("\t")
        if len(parts) < len(header):
            parts.extend([""] * (len(header) - len(parts)))
        row = {header[i]: parts[i] for i in range(len(header))}
        rows.append(row)

    return rows


def find_col(cols: Dict[str, int], candidates: List[str]) -> int:
    """从表头字典中寻找列索引，支持模糊匹配。"""
    for c in candidates:
        if c in cols:
            return cols[c]
    for k, i in cols.items():
        for c in candidates:
            if c in k:
                return i
    return -1


def resolve_path_maybe_relative(p: str | Path) -> Path:
    """相对路径相对项目根目录解析。"""
    pp = Path(str(p)).expanduser()
    if pp.is_absolute():
        return pp
    return (PROJECT_ROOT / pp).resolve()


# ==============================
# 路径与输入
# ==============================

def get_summary_dir(cfg: Dict[str, Any]) -> Path:
    """获取 results/05_sbranch_model/summary。"""
    bm = cfg.get("branch_model") or {}
    run_dir = bm.get("run_dir", "results/05_sbranch_model")
    return resolve_path_maybe_relative(run_dir) / "summary"


def get_psg_table(cfg: Dict[str, Any]) -> Path:
    """获取 PSG 主表路径。"""
    bm = cfg.get("branch_model") or {}
    p = bm.get("psg_table", str(DEFAULT_PSG_TABLE))
    return resolve_path_maybe_relative(p)


def get_beb_sig_table(cfg: Dict[str, Any]) -> Path:
    """获取 07b 的高可信 PSG 表路径。"""
    bm = cfg.get("branch_model") or {}
    p = bm.get("beb_sig_table", str(DEFAULT_BEB_SIG_TABLE))
    return resolve_path_maybe_relative(p)


# ==============================
# PSG 表读取与标准化
# ==============================

def load_psg_rows(psg_table: Path) -> List[Dict[str, str]]:
    """
    读取现有 PSG 表。
    兼容列名：
    - OG / orthogroup
    - foreground
    - P / pvalue
    - Q / qvalue
    - foreground_gene_ids
    """
    if not psg_table.is_file():
        raise FileNotFoundError(f"[ERR] PSG 主表不存在：{psg_table}")

    lines = psg_table.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines:
        raise ValueError(f"[ERR] PSG 主表为空：{psg_table}")

    header = lines[0].rstrip("\n").split("\t")
    cols = {h.strip().lower(): i for i, h in enumerate(header)}

    i_og = find_col(cols, ["og", "orthogroup"])
    i_fg = find_col(cols, ["foreground", "fg"])
    i_p = find_col(cols, ["p", "pvalue"])
    i_q = find_col(cols, ["q", "qvalue", "fdr"])
    i_gid = find_col(cols, ["foreground_gene_ids", "gene_ids", "gene_id"])

    if i_og < 0 or i_fg < 0 or i_p < 0:
        raise ValueError(
            f"[ERR] PSG 主表表头无法识别所需列。至少需要 OG / foreground / P。文件：{psg_table}"
        )

    rows: List[Dict[str, str]] = []
    for raw in lines[1:]:
        if not raw.strip():
            continue
        parts = raw.rstrip("\n").split("\t")
        if len(parts) <= max(i_og, i_fg, i_p):
            continue

        row = {
            "OG": parts[i_og].strip(),
            "foreground": parts[i_fg].strip(),
            "psg_p": parts[i_p].strip(),
            "psg_q": parts[i_q].strip() if i_q >= 0 and i_q < len(parts) else "",
            "foreground_gene_ids": parts[i_gid].strip() if i_gid >= 0 and i_gid < len(parts) else "",
        }
        if row["OG"] and row["foreground"]:
            rows.append(row)

    return rows


def load_beb_supported_keys(beb_sig_table: Path) -> Set[Tuple[str, str]]:
    """
    读取 07b 输出表，只提取满足 BEB>=0.95 的 (OG, foreground) 键。
    不把任何新列写入联合表，避免影响 r5。
    """
    if not beb_sig_table.is_file():
        raise FileNotFoundError(f"[ERR] BEB 支持表不存在：{beb_sig_table}")

    lines = beb_sig_table.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines:
        raise ValueError(f"[ERR] BEB 支持表为空：{beb_sig_table}")

    header = lines[0].rstrip("\n").split("\t")
    cols = {h.strip().lower(): i for i, h in enumerate(header)}

    i_og = find_col(cols, ["og", "orthogroup"])
    i_fg = find_col(cols, ["foreground", "fg"])

    if i_og < 0 or i_fg < 0:
        raise ValueError(
            f"[ERR] BEB 支持表表头无法识别所需列。至少需要 OG / foreground。文件：{beb_sig_table}"
        )

    keys: Set[Tuple[str, str]] = set()
    for raw in lines[1:]:
        if not raw.strip():
            continue
        parts = raw.rstrip("\n").split("\t")
        if len(parts) <= max(i_og, i_fg):
            continue
        og = parts[i_og].strip()
        fg = parts[i_fg].strip()
        if og and fg:
            keys.add((og, fg))

    return keys


# ==============================
# gene label 与显著性判定
# ==============================

def choose_joint_label(psg_gene_ids: str, reg_label: str, og: str) -> str:
    """
    选择联合表的 gene_label。
    优先：
    1）REG_plot_input 中的 gene_label
    2）PSG 的 foreground_gene_ids 第一个
    3）OG
    """
    if str(reg_label).strip():
        return str(reg_label).strip()

    gids = split_semicolon(psg_gene_ids)
    if gids:
        return gids[0]

    return og


def calc_minus_log10(x: Any, zero_cap: float = 300.0) -> str:
    """
    计算 -log10(x)
    - 空值或非法值返回空
    - x<=0 时返回封顶值
    """
    v = safe_float(x)
    if math.isnan(v):
        return ""
    if v <= 0:
        return f"{zero_cap:.6f}"
    return f"{(-math.log10(v)):.6f}"


def psg_is_sig(psg_p: str, psg_q: str, mode: str, alpha: float) -> bool:
    """
    PSG gene-level 显著性判定。
    mode:
      - q：优先用 q 值
      - p：用 p 值
    """
    mode = str(mode).strip().lower()
    if mode == "q":
        q = safe_float(psg_q)
        if math.isnan(q):
            return False
        return q < alpha

    p = safe_float(psg_p)
    if math.isnan(p):
        return False
    return p < alpha


def reg_is_sig_from_plot(reg_sig: str) -> bool:
    """直接使用 r3 输出的 reg_sig。"""
    return parse_bool(reg_sig)


def joint_class_of(psg_sig: bool, reg_sig: bool) -> str:
    """四分类。"""
    if psg_sig and reg_sig:
        return "Both"
    if psg_sig and not reg_sig:
        return "PSG_only"
    if (not psg_sig) and reg_sig:
        return "REG_only"
    return "None"


# ==============================
# 主流程
# ==============================

def main() -> int:
    cfg = load_config(CONFIG_PATH)

    paths = cfg.get("paths") or {}
    logs_dir = resolve_path_maybe_relative(paths.get("logs_dir", "logs"))
    log_file = logs_dir / "r4_build_psg_reg_joint.log"
    log = get_logger("aphylo.r4", log_file)

    banner(log, "APhylo r4 — build PSG × REG joint table")

    bm = cfg.get("branch_model") or {}
    psg_alpha_mode = str(bm.get("psg_alpha_mode", "q")).strip().lower()
    psg_alpha = float(bm.get("psg_alpha", 0.05))
    label_top_n = int(bm.get("label_top_n", 20))

    summary_dir = get_summary_dir(cfg)
    reg_plot_tsv = summary_dir / REG_INPUT_FILENAME
    psg_table = get_psg_table(cfg)
    beb_sig_table = get_beb_sig_table(cfg)

    if psg_alpha_mode not in {"q", "p"}:
        raise ValueError(f"[ERR] branch_model.psg_alpha_mode 只能是 q 或 p，当前为：{psg_alpha_mode}")

    log.info(f"[INIT] PSG_TABLE={psg_table}")
    log.info(f"[INIT] BEB_SIG_TABLE={beb_sig_table}")
    log.info(f"[INIT] REG_PLOT={reg_plot_tsv}")
    log.info(f"[RULE] psg_alpha_mode={psg_alpha_mode}; psg_alpha={psg_alpha}; label_top_n={label_top_n}")
    log.info("[RULE] 最终 psg_sig = gene-level 显著 AND BEB-supported")

    psg_rows = load_psg_rows(psg_table)
    reg_rows = read_tsv(reg_plot_tsv)
    beb_supported_keys = load_beb_supported_keys(beb_sig_table)

    if not psg_rows:
        log.error("[ERR] PSG 表无有效数据")
        return 2
    if not reg_rows:
        log.error("[ERR] REG_plot_input.tsv 无有效数据")
        return 2

    # 建索引：key = (OG, foreground)
    psg_map: Dict[Tuple[str, str], Dict[str, str]] = {}
    reg_map: Dict[Tuple[str, str], Dict[str, str]] = {}

    for r in psg_rows:
        key = (r["OG"], r["foreground"])
        psg_map[key] = r

    for r in reg_rows:
        key = (r.get("OG", ""), r.get("foreground", ""))
        if key[0] and key[1]:
            reg_map[key] = r

    # 用并集，保证四类都能出现
    all_keys = sorted(set(psg_map.keys()) | set(reg_map.keys()))
    log.info(f"[INIT] PSG_keys={len(psg_map)}; REG_keys={len(reg_map)}; UNION_keys={len(all_keys)}")
    log.info(f"[INIT] BEB_supported_keys={len(beb_supported_keys)}")

    rows_joint: List[Dict[str, Any]] = []
    class_count = {"None": 0, "PSG_only": 0, "REG_only": 0, "Both": 0}

    for og, fg in all_keys:
        key = (og, fg)
        pr = psg_map.get(key, {})
        rr = reg_map.get(key, {})

        psg_p = pr.get("psg_p", "")
        psg_q = pr.get("psg_q", "")
        reg_p = rr.get("reg_p", "")
        reg_q = rr.get("reg_q", "")
        reg_sig = rr.get("reg_sig", "")

        omega_fg = rr.get("omega_foreground", "")
        omega_bg = rr.get("omega_background", "")
        delta = rr.get("delta_omega", "")
        y_reg_fdr = rr.get("y_reg_fdr", "")

        fg_gene_ids = pr.get("foreground_gene_ids", "")
        reg_label = rr.get("gene_label", "")

        gene_label = choose_joint_label(fg_gene_ids, reg_label, og)

        # gene-level 显著
        psg_sig_gene = psg_is_sig(psg_p, psg_q, psg_alpha_mode, psg_alpha)

        # 叠加 BEB 支持
        psg_sig = psg_sig_gene and (key in beb_supported_keys)

        reg_sig_bool = reg_is_sig_from_plot(reg_sig)
        joint_class = joint_class_of(psg_sig, reg_sig_bool)
        class_count[joint_class] += 1

        row = {
            "OG": og,
            "foreground": fg,
            "gene_label": gene_label,
            "foreground_gene_ids": fg_gene_ids,
            "psg_p": psg_p,
            "psg_q": psg_q,
            "reg_p": reg_p,
            "reg_q": reg_q,
            "omega_foreground": omega_fg,
            "omega_background": omega_bg,
            "delta_omega": delta,
            "psg_sig": "TRUE" if psg_sig else "FALSE",
            "reg_sig": "TRUE" if reg_sig_bool else "FALSE",
            "joint_class": joint_class,
            "x_psg_p": calc_minus_log10(psg_p),
            "y_reg_fdr": y_reg_fdr,
        }
        rows_joint.append(row)

    # 写联合表
    # 表头严格保持原版 r4 一模一样
    out_joint = summary_dir / OUT_JOINT
    head_joint = [
        "OG",
        "foreground",
        "gene_label",
        "foreground_gene_ids",
        "psg_p",
        "psg_q",
        "reg_p",
        "reg_q",
        "omega_foreground",
        "omega_background",
        "delta_omega",
        "psg_sig",
        "reg_sig",
        "joint_class",
        "x_psg_p",
        "y_reg_fdr",
    ]
    write_tsv(out_joint, head_joint, rows_joint)

    # 写四类计数表
    out_counts = summary_dir / OUT_COUNTS
    rows_counts = [
        {"class": "None", "n": class_count["None"]},
        {"class": "PSG_only", "n": class_count["PSG_only"]},
        {"class": "REG_only", "n": class_count["REG_only"]},
        {"class": "Both", "n": class_count["Both"]},
    ]
    write_tsv(out_counts, ["class", "n"], rows_counts)

    # 生成标签推荐表
    # 规则：
    # 1) Both 优先
    # 2) 其次 PSG_only / REG_only
    # 3) 再按 x_psg_p + y_reg_fdr 降序
    class_priority = {"Both": 0, "PSG_only": 1, "REG_only": 2, "None": 3}

    def label_score(row: Dict[str, Any]) -> float:
        x = safe_float(row.get("x_psg_p", ""), default=0.0)
        y = safe_float(row.get("y_reg_fdr", ""), default=0.0)
        if math.isnan(x):
            x = 0.0
        if math.isnan(y):
            y = 0.0
        return x + y

    rows_for_label = []
    for r in rows_joint:
        rr = dict(r)
        rr["_score"] = label_score(r)
        rows_for_label.append(rr)

    rows_for_label.sort(
        key=lambda r: (
            class_priority.get(r.get("joint_class", "None"), 9),
            -r.get("_score", 0.0),
            r.get("OG", ""),
            r.get("foreground", "")
        )
    )

    rows_labels: List[Dict[str, Any]] = []
    for idx, r in enumerate(rows_for_label[:max(0, label_top_n)], start=1):
        rows_labels.append({
            "OG": r.get("OG", ""),
            "foreground": r.get("foreground", ""),
            "gene_label": r.get("gene_label", ""),
            "psg_p": r.get("psg_p", ""),
            "psg_q": r.get("psg_q", ""),
            "reg_p": r.get("reg_p", ""),
            "reg_q": r.get("reg_q", ""),
            "x_psg_p": r.get("x_psg_p", ""),
            "y_reg_fdr": r.get("y_reg_fdr", ""),
            "joint_class": r.get("joint_class", ""),
            "label_rank": idx,
        })

    out_labels = summary_dir / OUT_LABELS
    head_labels = [
        "OG",
        "foreground",
        "gene_label",
        "psg_p",
        "psg_q",
        "reg_p",
        "reg_q",
        "x_psg_p",
        "y_reg_fdr",
        "joint_class",
        "label_rank",
    ]
    write_tsv(out_labels, head_labels, rows_labels)

    done_file = summary_dir / DONE_FILENAME
    done_file.touch()

    log.info(f"[DONE] 写出：{out_joint}")
    log.info(f"[DONE] 写出：{out_counts}")
    log.info(f"[DONE] 写出：{out_labels}")
    log.info(
        f"[STAT] None={class_count['None']} "
        f"PSG_only={class_count['PSG_only']} "
        f"REG_only={class_count['REG_only']} "
        f"Both={class_count['Both']}"
    )
    log.info(f"[DONE] sentinel={done_file}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断")
        sys.exit(130)
