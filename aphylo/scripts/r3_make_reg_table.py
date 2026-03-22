#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
r3_make_reg_table.py

功能：
1）读取 r2_branch_model_aggregate.py 生成的 branch_model_all.tsv
2）根据 config.yaml: branch_model 中的阈值规则，判定 rapidly evolving genes (REG)
3）输出：
   results/05_sbranch_model/summary/
     - REG_all.tsv
     - REG_qsig.tsv
     - REG_qsig_gene.tsv
     - REG_plot_input.tsv
     - .r3.done

默认判定规则：
- status == ok
- is_testable == TRUE
- Q < fdr_alpha
- omega_foreground > omega_background

可选附加规则：
- require_omega_fg_gt_1: true/false
- min_delta_omega: float

说明：
- 本脚本不重新计算统计量，只做筛选与整理
- gene_label 优先级：
  foreground_gene_names > foreground_gene_ids > OG
- 如果一个 OG 对应多个 foreground gene_id / gene_name，
  REG_qsig_gene.tsv 会尽量按分号拆开逐行展开；
  若 gene_name 数量与 gene_id 数量不一致，则按位置尽力匹配，不足处留空
"""

from __future__ import annotations

import io
import sys
import math
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple


# ==============================
# 配置区（皇上主要看这里）
# ==============================

PROJECT_ROOT = Path(".").resolve()
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

INPUT_FILENAME = "branch_model_all.tsv"
OUT_REG_ALL = "REG_all.tsv"
OUT_REG_QSIG = "REG_qsig.tsv"
OUT_REG_QSIG_GENE = "REG_qsig_gene.tsv"
OUT_REG_PLOT = "REG_plot_input.tsv"
DONE_FILENAME = ".r3.done"


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


def ts() -> str:
    """返回时间字符串。"""
    import datetime
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


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
    """把各种真假字符串转成 bool。"""
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
    """按分号拆分，去空。"""
    if s is None:
        return []
    return [x.strip() for x in str(s).split(";") if x.strip()]


def choose_gene_label(gene_names: str, gene_ids: str, og: str) -> str:
    """
    选择作图标签：
    优先 gene_name，其次 gene_id，最后 OG
    如果有多个，取第一个
    """
    gns = split_semicolon(gene_names)
    if gns:
        return gns[0]

    gids = split_semicolon(gene_ids)
    if gids:
        return gids[0]

    return og


def paired_gene_rows(gene_ids: str, gene_names: str) -> List[Tuple[str, str]]:
    """
    将 foreground_gene_ids / foreground_gene_names 尽量配对展开。
    规则：
    - 若两边长度都 >0，按位置配对
    - 较短的一侧不足时补空字符串
    - 若一侧为空，则另一侧正常展开
    """
    gids = split_semicolon(gene_ids)
    gns = split_semicolon(gene_names)

    if not gids and not gns:
        return []

    n = max(len(gids), len(gns))
    out = []
    for i in range(n):
        gid = gids[i] if i < len(gids) else ""
        gnm = gns[i] if i < len(gns) else ""
        out.append((gid, gnm))
    return out


def reg_class_of_row(
    row: Dict[str, str],
    fdr_alpha: float,
    require_omega_fg_gt_bg: bool,
    require_omega_fg_gt_1: bool,
    min_delta_omega: float
) -> Tuple[str, bool, str]:
    """
    根据规则判定当前行的 reg_class 与 is_reg。
    返回：
      (reg_class, is_reg, note)
    """

    status = str(row.get("status", "")).strip()
    is_testable = parse_bool(row.get("is_testable", "FALSE"))

    if status != "ok":
        return "failed", False, "status_not_ok"

    if not is_testable:
        return "failed", False, "not_testable"

    q = safe_float(row.get("Q", ""))
    p = safe_float(row.get("P", ""))
    omega_fg = safe_float(row.get("omega_foreground", ""))
    omega_bg = safe_float(row.get("omega_background", ""))
    delta = safe_float(row.get("delta_omega", ""))

    if math.isnan(q):
        return "failed", False, "missing_Q"
    if math.isnan(p):
        return "failed", False, "missing_P"
    if math.isnan(omega_fg):
        return "failed", False, "missing_omega_fg"
    if math.isnan(omega_bg):
        return "failed", False, "missing_omega_bg"
    if math.isnan(delta):
        delta = omega_fg - omega_bg

    if q >= fdr_alpha:
        return "not_sig", False, f"Q>={fdr_alpha}"

    if require_omega_fg_gt_bg and not (omega_fg > omega_bg):
        return "fg_not_higher", False, "omega_fg_not_gt_bg"

    if require_omega_fg_gt_1 and not (omega_fg > 1.0):
        return "fg_not_gt_1", False, "omega_fg_not_gt_1"

    if delta < min_delta_omega:
        return "delta_too_small", False, f"delta<{min_delta_omega}"

    return "REG", True, "pass"


def get_summary_dir(cfg: Dict[str, Any]) -> Path:
    """获取 results/05_sbranch_model/summary 目录。"""
    bm = cfg.get("branch_model") or {}
    run_dir = bm.get("run_dir", "results/05_sbranch_model")
    return PROJECT_ROOT / str(run_dir) / "summary"


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


def write_tsv(path: Path, header: List[str], rows: List[Dict[str, Any]]) -> None:
    """写 TSV。"""
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as w:
        w.write("\t".join(header) + "\n")
        for row in rows:
            w.write("\t".join(str(row.get(k, "")) for k in header) + "\n")


# ==============================
# 主流程
# ==============================

def main() -> int:
    cfg = load_config(CONFIG_PATH)

    paths = cfg.get("paths") or {}
    logs_dir = PROJECT_ROOT / str(paths.get("logs_dir", "logs"))
    log_file = logs_dir / "r3_make_reg_table.log"
    log = get_logger("aphylo.r3", log_file)

    banner(log, "APhylo r3 — make REG tables")

    bm = cfg.get("branch_model") or {}
    fdr_alpha = float(bm.get("fdr_alpha", 0.05))
    require_omega_fg_gt_bg = bool(bm.get("require_omega_fg_gt_bg", True))
    require_omega_fg_gt_1 = bool(bm.get("require_omega_fg_gt_1", False))
    min_delta_omega = float(bm.get("min_delta_omega", 0.0))

    summary_dir = get_summary_dir(cfg)
    input_tsv = summary_dir / INPUT_FILENAME

    rows = read_tsv(input_tsv)
    if not rows:
        log.error(f"[ERR] 输入表无数据：{input_tsv}")
        return 2

    log.info(f"[INIT] input={input_tsv}")
    log.info(f"[INIT] rows={len(rows)}")
    log.info(
        f"[RULE] fdr_alpha={fdr_alpha}; "
        f"require_omega_fg_gt_bg={require_omega_fg_gt_bg}; "
        f"require_omega_fg_gt_1={require_omega_fg_gt_1}; "
        f"min_delta_omega={min_delta_omega}"
    )

    rows_reg_all: List[Dict[str, Any]] = []
    rows_reg_qsig: List[Dict[str, Any]] = []
    rows_reg_qsig_gene: List[Dict[str, Any]] = []
    rows_reg_plot: List[Dict[str, Any]] = []

    stat_total = 0
    stat_testable = 0
    stat_reg = 0
    stat_failed = 0
    stat_not_sig = 0
    stat_fg_not_higher = 0
    stat_fg_not_gt_1 = 0
    stat_delta_small = 0

    for row in rows:
        stat_total += 1

        if parse_bool(row.get("is_testable", "FALSE")):
            stat_testable += 1

        og = row.get("OG", "")
        fg = row.get("foreground", "")
        gene_ids = row.get("foreground_gene_ids", "")
        gene_names = row.get("foreground_gene_names", "")

        reg_class, is_reg, note_r3 = reg_class_of_row(
            row=row,
            fdr_alpha=fdr_alpha,
            require_omega_fg_gt_bg=require_omega_fg_gt_bg,
            require_omega_fg_gt_1=require_omega_fg_gt_1,
            min_delta_omega=min_delta_omega
        )

        if reg_class == "REG":
            stat_reg += 1
        elif reg_class == "failed":
            stat_failed += 1
        elif reg_class == "not_sig":
            stat_not_sig += 1
        elif reg_class == "fg_not_higher":
            stat_fg_not_higher += 1
        elif reg_class == "fg_not_gt_1":
            stat_fg_not_gt_1 += 1
        elif reg_class == "delta_too_small":
            stat_delta_small += 1

        gene_label = choose_gene_label(gene_names, gene_ids, og)

        out_row = {
            "OG": og,
            "foreground": fg,
            "foreground_gene_ids": gene_ids,
            "foreground_gene_names": gene_names,
            "P": row.get("P", ""),
            "Q": row.get("Q", ""),
            "lnL_one": row.get("lnL_one", ""),
            "lnL_two": row.get("lnL_two", ""),
            "LRT": row.get("LRT", ""),
            "omega_one": row.get("omega_one", ""),
            "omega_foreground": row.get("omega_foreground", ""),
            "omega_background": row.get("omega_background", ""),
            "delta_omega": row.get("delta_omega", ""),
            "omega_fg_gt_bg": row.get("omega_fg_gt_bg", ""),
            "omega_fg_gt_1": row.get("omega_fg_gt_1", ""),
            "reg_class": reg_class,
            "is_reg": "TRUE" if is_reg else "FALSE",
            "note": note_r3 if note_r3 else row.get("note", ""),
        }
        rows_reg_all.append(out_row)

        if is_reg:
            qsig_row = {
                "OG": og,
                "foreground": fg,
                "foreground_gene_ids": gene_ids,
                "foreground_gene_names": gene_names,
                "Q": row.get("Q", ""),
                "P": row.get("P", ""),
                "omega_foreground": row.get("omega_foreground", ""),
                "omega_background": row.get("omega_background", ""),
                "delta_omega": row.get("delta_omega", ""),
            }
            rows_reg_qsig.append(qsig_row)

            gene_pairs = paired_gene_rows(gene_ids, gene_names)
            if gene_pairs:
                for gid, gnm in gene_pairs:
                    rows_reg_qsig_gene.append({
                        "OG": og,
                        "foreground": fg,
                        "gene_id": gid,
                        "gene_name": gnm,
                        "Q": row.get("Q", ""),
                        "P": row.get("P", ""),
                        "omega_foreground": row.get("omega_foreground", ""),
                        "omega_background": row.get("omega_background", ""),
                        "delta_omega": row.get("delta_omega", ""),
                    })
            else:
                rows_reg_qsig_gene.append({
                    "OG": og,
                    "foreground": fg,
                    "gene_id": "",
                    "gene_name": "",
                    "Q": row.get("Q", ""),
                    "P": row.get("P", ""),
                    "omega_foreground": row.get("omega_foreground", ""),
                    "omega_background": row.get("omega_background", ""),
                    "delta_omega": row.get("delta_omega", ""),
                })

        # REG_plot_input.tsv：全部都保留，后面 r4 合并更方便
        reg_q = safe_float(row.get("Q", ""))
        y_reg_fdr = ""
        if not math.isnan(reg_q):
            if reg_q > 0:
                y_reg_fdr = f"{(-math.log10(reg_q)):.6f}"
            else:
                # 避免 q=0 导致无穷，给一个很大的封顶值
                y_reg_fdr = "300.000000"

        rows_reg_plot.append({
            "OG": og,
            "foreground": fg,
            "gene_label": gene_label,
            "reg_p": row.get("P", ""),
            "reg_q": row.get("Q", ""),
            "omega_foreground": row.get("omega_foreground", ""),
            "omega_background": row.get("omega_background", ""),
            "delta_omega": row.get("delta_omega", ""),
            "reg_sig": "TRUE" if is_reg else "FALSE",
            "y_reg_fdr": y_reg_fdr,
        })

    # 排序
    def _sort_q_then_p(r: Dict[str, Any]):
        q = safe_float(r.get("Q", ""))
        p = safe_float(r.get("P", ""))
        if math.isnan(q):
            q = 999.0
        if math.isnan(p):
            p = 999.0
        return (q, p)

    rows_reg_qsig.sort(key=_sort_q_then_p)
    rows_reg_qsig_gene.sort(key=_sort_q_then_p)
    rows_reg_plot.sort(
        key=lambda r: (
            0 if r.get("reg_sig", "") == "TRUE" else 1,
            safe_float(r.get("reg_q", ""), default=999.0),
            safe_float(r.get("reg_p", ""), default=999.0)
        )
    )

    # 写出 4 张表
    out1 = summary_dir / OUT_REG_ALL
    head1 = [
        "OG",
        "foreground",
        "foreground_gene_ids",
        "foreground_gene_names",
        "P",
        "Q",
        "lnL_one",
        "lnL_two",
        "LRT",
        "omega_one",
        "omega_foreground",
        "omega_background",
        "delta_omega",
        "omega_fg_gt_bg",
        "omega_fg_gt_1",
        "reg_class",
        "is_reg",
        "note",
    ]
    write_tsv(out1, head1, rows_reg_all)

    out2 = summary_dir / OUT_REG_QSIG
    head2 = [
        "OG",
        "foreground",
        "foreground_gene_ids",
        "foreground_gene_names",
        "Q",
        "P",
        "omega_foreground",
        "omega_background",
        "delta_omega",
    ]
    write_tsv(out2, head2, rows_reg_qsig)

    out3 = summary_dir / OUT_REG_QSIG_GENE
    head3 = [
        "OG",
        "foreground",
        "gene_id",
        "gene_name",
        "Q",
        "P",
        "omega_foreground",
        "omega_background",
        "delta_omega",
    ]
    write_tsv(out3, head3, rows_reg_qsig_gene)

    out4 = summary_dir / OUT_REG_PLOT
    head4 = [
        "OG",
        "foreground",
        "gene_label",
        "reg_p",
        "reg_q",
        "omega_foreground",
        "omega_background",
        "delta_omega",
        "reg_sig",
        "y_reg_fdr",
    ]
    write_tsv(out4, head4, rows_reg_plot)

    done_file = summary_dir / DONE_FILENAME
    done_file.touch()

    log.info(f"[DONE] 写出：{out1}")
    log.info(f"[DONE] 写出：{out2}")
    log.info(f"[DONE] 写出：{out3}")
    log.info(f"[DONE] 写出：{out4}")
    log.info(f"[STAT] total={stat_total}")
    log.info(f"[STAT] testable={stat_testable}")
    log.info(f"[STAT] REG={stat_reg}")
    log.info(f"[STAT] failed={stat_failed}")
    log.info(f"[STAT] not_sig={stat_not_sig}")
    log.info(f"[STAT] fg_not_higher={stat_fg_not_higher}")
    log.info(f"[STAT] fg_not_gt_1={stat_fg_not_gt_1}")
    log.info(f"[STAT] delta_too_small={stat_delta_small}")
    log.info(f"[DONE] sentinel={done_file}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断")
        sys.exit(130)
