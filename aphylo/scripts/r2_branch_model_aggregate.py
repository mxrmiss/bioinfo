#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
r2_branch_model_aggregate.py

功能：
1）读取 r1_codeml_branch_batch.py 生成的 branch model 结果
2）对每个 OG × foreground，解析：
   - one-ratio 的 lnL 与 omega_one
   - two-ratio 的 lnL 与 omega_background / omega_foreground
3）计算：
   - LRT = 2 * (lnL_two - lnL_one)
   - P（卡方 df=1 上尾）
   - Q（BH-FDR）
4）输出：
   results/05_sbranch_model/summary/
     - branch_model_all.tsv
     - branch_model_failed.tsv
     - foreground_gene_map.tsv
     - .r2.done

设计原则：
- 路径与现有 aphylo 保持一致
- 尽量复用 07_codeml_aggregate.py 的 foreground_gene_ids 构造逻辑
- gene_name 若上游 pep2cds_map 无对应列，则输出为空
- 本脚本只做“聚合与统计”，不做最终 REG 判定（交给 r3）

运行方式：
  python scripts/r2_branch_model_aggregate.py
"""

from __future__ import annotations

import re
import io
import sys
import math
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

# ==============================
# 配置区（皇上主要看这里）
# ==============================

PROJECT_ROOT = Path(".").resolve()
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

SUMMARY_DIRNAME = "summary"
ALL_TSV = "branch_model_all.tsv"
FAILED_TSV = "branch_model_failed.tsv"
FG_MAP_TSV = "foreground_gene_map.tsv"
DONE_FILENAME = ".r2.done"

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
# 日志与基础工具
# ==============================

def ensure_dir(p: Path) -> Path:
    """确保目录存在。"""
    p.mkdir(parents=True, exist_ok=True)
    return p


def need_file(p: Path, what: str) -> Path:
    """要求文件必须存在。"""
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
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


def write_done(path: Path):
    """写哨兵文件。"""
    path.touch()


# ==============================
# 统计函数
# ==============================

def bh_fdr(pvals: List[float]) -> List[float]:
    """
    Benjamini-Hochberg FDR 校正。
    返回与输入同顺序的 q 值列表。
    """
    m = len(pvals)
    if m == 0:
        return []

    order = sorted(range(m), key=lambda i: pvals[i])  # 升序
    q = [1.0] * m
    prev = 1.0

    for rank_from_end, i in enumerate(reversed(order), start=1):
        rank = m - rank_from_end + 1
        val = pvals[i] * m / rank
        if val > 1.0:
            val = 1.0
        prev = min(prev, val)
        q[i] = prev

    return q


def chi2_df1_sf(x: float) -> float:
    """
    卡方分布 df=1 的上尾概率。
    对于 df=1，有：
      SF(x) = erfc(sqrt(x/2))
    """
    if x <= 0:
        return 1.0
    return math.erfc(math.sqrt(x / 2.0))


# ==============================
# Newick 与 tips 解析
# ==============================

_TOK_TIP = re.compile(r'(?<=\(|,)\s*([^()\s:,]+)\s*(?=[:),])')


def parse_tips(nwk_text: str) -> List[str]:
    """从 Newick 文本提取 tip。"""
    return _TOK_TIP.findall(nwk_text)


def strip_mark(name: str) -> str:
    """去掉 #1 标记。"""
    return name.replace("#1", "")


# ==============================
# pep2cds / order / foreground 映射
# ==============================

def normalize_key(x: str) -> str:
    """
    规范化 ID：
    - 去空白
    - 只取空白前 token
    - 如果有 |，只取最后一段
    """
    x = (x or "").strip()
    if not x:
        return ""
    x = x.split()[0]
    return x.split("|")[-1]


def find_col(cols: Dict[str, int], candidates: List[str]) -> int:
    """按候选列名寻找列索引，支持模糊匹配。"""
    for c in candidates:
        if c in cols:
            return cols[c]
    for k, i in cols.items():
        for c in candidates:
            if c in k:
                return i
    return -1


def load_pep2cds_info(map_path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """
    读取 pep2cds_resolved.tsv，返回：
      info[(OG, species)] = {
          "protein_id": ...,
          "cds_id": ...,
          "gene_id": ...,
          "gene_name": ...
      }

    注：
    - gene_id / gene_name 列可能不存在；不存在则留空
    """
    need_file(map_path, "pep2cds_map")
    lines = map_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines:
        raise ValueError(f"[ERR] 映射表为空：{map_path}")

    header = lines[0].rstrip("\n").split("\t")
    cols = {h.strip().lower(): i for i, h in enumerate(header)}

    i_og = find_col(cols, ["og", "orthogroup"])
    i_sp = find_col(cols, ["species", "sp"])
    i_pid = find_col(cols, ["protein_id", "pep_id", "peptide_id", "protein", "pep", "prot", "peptide"])
    i_cds = find_col(cols, ["cds_id", "cds", "transcript_id", "mrna_id"])
    i_gid = find_col(cols, ["gene_id", "gene"])
    i_gn = find_col(cols, ["gene_name", "symbol", "gene_symbol", "name"])

    if i_og < 0 or i_sp < 0:
        raise ValueError("[ERR] pep2cds_map 缺少 OG / species 列。")

    info: Dict[Tuple[str, str], Dict[str, str]] = {}

    for raw in lines[1:]:
        if not raw.strip():
            continue
        parts = raw.rstrip("\n").split("\t")
        if len(parts) <= max(i_og, i_sp):
            continue

        og = parts[i_og].strip()
        sp = parts[i_sp].strip()
        if not og or not sp:
            continue

        pid = parts[i_pid].strip() if i_pid >= 0 and i_pid < len(parts) else ""
        cds = parts[i_cds].strip() if i_cds >= 0 and i_cds < len(parts) else ""
        gid = parts[i_gid].strip() if i_gid >= 0 and i_gid < len(parts) else ""
        gnm = parts[i_gn].strip() if i_gn >= 0 and i_gn < len(parts) else ""

        info[(og, sp)] = {
            "protein_id": normalize_key(pid),
            "cds_id": normalize_key(cds),
            "gene_id": gid,
            "gene_name": gnm,
        }

    return info


def get_fg_tips(sets_dir: Path, fg: str, cache: Dict[str, set], log: logging.Logger) -> set:
    """读取 foreground 的物种列表。"""
    if fg in cache:
        return cache[fg]

    fp = sets_dir / f"{fg}.list"
    if not fp.is_file():
        log.warning(f"[WARN] 找不到前景列表文件：{fp}")
        cache[fg] = set()
        return cache[fg]

    tips = [ln.strip() for ln in fp.read_text(encoding="utf-8").splitlines() if ln.strip()]
    cache[fg] = set(tips)
    return cache[fg]


def get_og_order_rows(order_dir: Path, og: str, cache: Dict[str, List[Tuple[str, str]]], log: logging.Logger) -> List[Tuple[str, str]]:
    """
    读取 02 产生的 order/OG.order.tsv
    返回：
      [(Species, cds_id), ...]
    """
    if og in cache:
        return cache[og]

    tsv = order_dir / f"{og}.order.tsv"
    if not tsv.is_file():
        log.warning(f"[WARN] 找不到 OG={og} 的 order.tsv：{tsv}")
        cache[og] = []
        return cache[og]

    rows: List[Tuple[str, str]] = []
    for raw in tsv.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split("\t")
        if parts[0] == "OG":
            continue
        if len(parts) < 4:
            continue
        species = parts[1].strip()
        cds_id = normalize_key(parts[3])
        rows.append((species, cds_id))

    cache[og] = rows
    return cache[og]


def build_foreground_gene_bundle(
    og: str,
    fg: str,
    sets_dir: Path,
    order_dir: Path,
    pep2cds_info: Dict[Tuple[str, str], Dict[str, str]],
    fg_tips_cache: Dict[str, set],
    og_order_cache: Dict[str, List[Tuple[str, str]]],
    log: logging.Logger
) -> Dict[str, str]:
    """
    构建前景基因信息：
    - foreground_species
    - foreground_gene_ids
    - foreground_gene_names
    """
    tips = get_fg_tips(sets_dir, fg, fg_tips_cache, log)
    rows = get_og_order_rows(order_dir, og, og_order_cache, log)

    fg_species = sorted(tips)
    fg_gene_ids: List[str] = []
    fg_gene_names: List[str] = []

    if rows:
        cds_by_sp = {sp: cds for sp, cds in rows}
        for sp in fg_species:
            cds_id = cds_by_sp.get(sp, "")
            meta = pep2cds_info.get((og, sp), {})

            gene_id = meta.get("gene_id", "").strip()
            gene_name = meta.get("gene_name", "").strip()

            # 优先 gene_id；若没有则退回 cds_id
            if gene_id:
                fg_gene_ids.append(gene_id)
            elif cds_id:
                fg_gene_ids.append(cds_id)

            if gene_name:
                fg_gene_names.append(gene_name)

    # 去重并保持排序
    fg_gene_ids = sorted(set(x for x in fg_gene_ids if x))
    fg_gene_names = sorted(set(x for x in fg_gene_names if x))

    return {
        "foreground_species": ";".join(fg_species),
        "foreground_gene_ids": ";".join(fg_gene_ids),
        "foreground_gene_names": ";".join(fg_gene_names),
    }


# ==============================
# mlc 解析
# ==============================

_RE_LNL = re.compile(r"lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.Ee+]+)")
_RE_OMEGA_DNDS = re.compile(r"omega\s*\(dN/dS\)\s*=\s*([-\d.Ee+]+)", re.I)
_RE_OMEGA_EQ = re.compile(r"(?<!fix_)omega\s*=\s*([-\d.Ee+]+)", re.I)

# 例如：
# omega (0) = 0.12345
# omega (1) = 0.45678
_RE_OMEGA_LABEL = re.compile(r"omega\s*\((\d+)\)\s*=\s*([-\d.Ee+]+)", re.I)

# 例如：
# w (dN/dS) for branches: 0.12345 0.45678
_RE_W_BRANCHES = re.compile(r"w\s*\(dN/dS\)\s*for\s*branches\s*:\s*(.+)$", re.I)


def safe_float(x: str) -> Optional[float]:
    """安全转 float。"""
    try:
        return float(x)
    except Exception:
        return None


def last_lnl(txt: str) -> Optional[float]:
    """提取最后一个 lnL。"""
    vals = [safe_float(m.group(1)) for m in _RE_LNL.finditer(txt)]
    vals = [x for x in vals if x is not None]
    if not vals:
        return None
    return vals[-1]


def parse_omega_one(txt: str) -> Optional[float]:
    """
    解析 one-ratio 的 omega。
    优先：
    1）omega (dN/dS) = ...
    2）普通 omega = ...
    """
    m = _RE_OMEGA_DNDS.search(txt)
    if m:
        return safe_float(m.group(1))

    candidates = []
    for mm in _RE_OMEGA_EQ.finditer(txt):
        val = safe_float(mm.group(1))
        if val is not None:
            candidates.append(val)

    if candidates:
        return candidates[-1]

    return None


def parse_omega_bg_fg(txt: str) -> Tuple[Optional[float], Optional[float], str]:
    """
    解析 two-ratio 的 background / foreground omega。

    解析优先级：
    1）显式 omega (0), omega (1)
       - 认为 0=background, 1=foreground
    2）w (dN/dS) for branches: x y
       - 当恰好有两个值时，默认第 1 个是 background，第 2 个是 foreground

    返回：
      (omega_bg, omega_fg, parse_method)
    """
    # 方法1：显式 omega(0), omega(1)
    found = {}
    for mm in _RE_OMEGA_LABEL.finditer(txt):
        idx = mm.group(1)
        val = safe_float(mm.group(2))
        if val is not None:
            found[idx] = val

    if "0" in found and "1" in found:
        return found["0"], found["1"], "omega(label)"

    # 方法2：w (dN/dS) for branches: ...
    for line in txt.splitlines():
        m = _RE_W_BRANCHES.search(line.strip())
        if m:
            nums = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", m.group(1))
            vals = [safe_float(x) for x in nums]
            vals = [x for x in vals if x is not None]
            if len(vals) == 2:
                return vals[0], vals[1], "w_for_branches"

    return None, None, ""


# ==============================
# 结果遍历
# ==============================

def get_run_root(cfg: Dict[str, Any]) -> Path:
    """获取 results/05_sbranch_model 根目录。"""
    bm = cfg.get("branch_model") or {}
    run_dir = bm.get("run_dir", "results/05_sbranch_model")
    return PROJECT_ROOT / str(run_dir)


def collect_all_units(run_root: Path) -> List[Tuple[str, str, Path, Path]]:
    """
    收集所有 OG × FG 组合，返回：
      [(og, fg, one_mlc, two_mlc), ...]
    """
    out: List[Tuple[str, str, Path, Path]] = []
    runs_root = run_root / "runs"
    if not runs_root.is_dir():
        return out

    for og_dir in sorted(runs_root.iterdir()):
        if not og_dir.is_dir():
            continue
        og = og_dir.name
        for fg_dir in sorted(og_dir.iterdir()):
            if not fg_dir.is_dir():
                continue
            fg = fg_dir.name
            one_mlc = fg_dir / "one" / "mlc.txt"
            two_mlc = fg_dir / "two" / "mlc.txt"
            out.append((og, fg, one_mlc, two_mlc))

    return out


# ==============================
# 主流程
# ==============================

def main() -> int:
    cfg = load_config(CONFIG_PATH)

    paths = cfg.get("paths") or {}
    logs_dir = PROJECT_ROOT / str(paths.get("logs_dir", "logs"))
    log_file = logs_dir / "r2_branch_model_aggregate.log"
    log = get_logger("aphylo.r2", log_file)

    banner(log, "APhylo r2 — branch model 聚合")

    run_root = get_run_root(cfg)
    summary_dir = ensure_dir(run_root / SUMMARY_DIRNAME)

    inputs = cfg.get("inputs") or {}
    pep2cds_map = inputs.get("pep2cds_map", "")
    if not pep2cds_map:
        raise ValueError("[ERR] config.yaml 缺少 inputs.pep2cds_map")

    pep2cds_info = load_pep2cds_info(Path(str(pep2cds_map)).expanduser())

    codeml_dir = PROJECT_ROOT / str(paths.get("codeml_dir", "results/04_codeml"))
    bt_dir = PROJECT_ROOT / str(paths.get("bt_dir", "results/02_bt"))
    sets_dir = codeml_dir / "sets"
    order_dir = bt_dir / "order"

    fg_tips_cache: Dict[str, set] = {}
    og_order_cache: Dict[str, List[Tuple[str, str]]] = {}

    units = collect_all_units(run_root)
    if not units:
        log.error(f"[ERR] 未发现任何 r1 结果目录：{run_root / 'runs'}")
        return 2

    log.info(f"[INIT] RUN_ROOT={run_root}")
    log.info(f"[INIT] units={len(units)}")
    log.info(f"[INIT] sets_dir={sets_dir}")
    log.info(f"[INIT] order_dir={order_dir}")

    rows_all: List[Dict[str, Any]] = []
    rows_failed: List[Dict[str, Any]] = []
    rows_fgmap: List[Dict[str, Any]] = []

    for og, fg, one_mlc, two_mlc in units:
        bundle = build_foreground_gene_bundle(
            og=og,
            fg=fg,
            sets_dir=sets_dir,
            order_dir=order_dir,
            pep2cds_info=pep2cds_info,
            fg_tips_cache=fg_tips_cache,
            og_order_cache=og_order_cache,
            log=log
        )

        fg_map_row = {
            "OG": og,
            "foreground": fg,
            "foreground_species": bundle["foreground_species"],
            "foreground_gene_ids": bundle["foreground_gene_ids"],
            "foreground_gene_names": bundle["foreground_gene_names"],
        }
        rows_fgmap.append(fg_map_row)

        tree_path = run_root / "runs" / og / fg / "two" / "species_tree.nwk"
        n_tips = ""
        n_fg_tips = ""
        if tree_path.is_file():
            tree_txt = tree_path.read_text(encoding="utf-8", errors="ignore")
            tips = parse_tips(tree_txt)
            n_tips = str(len(tips))
            n_fg_tips = str(sum(1 for x in tips if "#1" in x))

        row = {
            "OG": og,
            "foreground": fg,
            "status": "ok",
            "n_tips": n_tips,
            "n_fg_tips": n_fg_tips,
            "foreground_species": bundle["foreground_species"],
            "foreground_gene_ids": bundle["foreground_gene_ids"],
            "foreground_gene_names": bundle["foreground_gene_names"],
            "lnL_one": "",
            "lnL_two": "",
            "LRT": "",
            "df": "1",
            "P": "",
            "Q": "",
            "omega_one": "",
            "omega_foreground": "",
            "omega_background": "",
            "delta_omega": "",
            "omega_fg_gt_bg": "",
            "omega_fg_gt_1": "",
            "is_testable": "FALSE",
            "omega_parse_method": "",
            "note": "",
        }

        fail_reason = ""

        if not one_mlc.is_file():
            fail_reason = f"missing_one_mlc:{one_mlc}"
        elif not two_mlc.is_file():
            fail_reason = f"missing_two_mlc:{two_mlc}"
        else:
            one_txt = one_mlc.read_text(encoding="utf-8", errors="ignore")
            two_txt = two_mlc.read_text(encoding="utf-8", errors="ignore")

            lnL_one = last_lnl(one_txt)
            lnL_two = last_lnl(two_txt)
            omega_one = parse_omega_one(one_txt)
            omega_bg, omega_fg, parse_method = parse_omega_bg_fg(two_txt)

            row["omega_parse_method"] = parse_method

            if lnL_one is None:
                fail_reason = "parse_fail:lnL_one"
            elif lnL_two is None:
                fail_reason = "parse_fail:lnL_two"
            elif omega_one is None:
                fail_reason = "parse_fail:omega_one"
            elif omega_bg is None or omega_fg is None:
                fail_reason = "parse_fail:omega_bg_or_fg"
            else:
                lrt = 2.0 * (lnL_two - lnL_one)
                if lrt < 0:
                    lrt = 0.0
                pval = chi2_df1_sf(lrt)
                delta = omega_fg - omega_bg

                row["lnL_one"] = f"{lnL_one:.6f}"
                row["lnL_two"] = f"{lnL_two:.6f}"
                row["LRT"] = f"{lrt:.6f}"
                row["P"] = f"{pval:.6g}"
                row["omega_one"] = f"{omega_one:.6g}"
                row["omega_foreground"] = f"{omega_fg:.6g}"
                row["omega_background"] = f"{omega_bg:.6g}"
                row["delta_omega"] = f"{delta:.6g}"
                row["omega_fg_gt_bg"] = "TRUE" if omega_fg > omega_bg else "FALSE"
                row["omega_fg_gt_1"] = "TRUE" if omega_fg > 1.0 else "FALSE"
                row["is_testable"] = "TRUE"

        if fail_reason:
            row["status"] = "failed"
            row["note"] = fail_reason
            rows_failed.append({
                "OG": og,
                "foreground": fg,
                "stage": "aggregate",
                "status": "failed",
                "reason": fail_reason,
                "one_mlc": str(one_mlc),
                "two_mlc": str(two_mlc),
            })

        rows_all.append(row)

    # 只对可检验行做 FDR
    testable_idx = [i for i, r in enumerate(rows_all) if r["is_testable"] == "TRUE" and r["P"] != ""]
    testable_p = [float(rows_all[i]["P"]) for i in testable_idx]
    testable_q = bh_fdr(testable_p)

    for i, q in zip(testable_idx, testable_q):
        rows_all[i]["Q"] = f"{q:.6g}"

    # 不可检验行 Q 填空
    for r in rows_all:
        if r["Q"] == "":
            r["Q"] = ""

    # 写主表
    out_all = summary_dir / ALL_TSV
    head_all = [
        "OG",
        "foreground",
        "status",
        "n_tips",
        "n_fg_tips",
        "foreground_species",
        "foreground_gene_ids",
        "foreground_gene_names",
        "lnL_one",
        "lnL_two",
        "LRT",
        "df",
        "P",
        "Q",
        "omega_one",
        "omega_foreground",
        "omega_background",
        "delta_omega",
        "omega_fg_gt_bg",
        "omega_fg_gt_1",
        "is_testable",
        "omega_parse_method",
        "note",
    ]
    with open(out_all, "w", encoding="utf-8") as w:
        w.write("\t".join(head_all) + "\n")
        for r in rows_all:
            w.write("\t".join(str(r.get(k, "")) for k in head_all) + "\n")

    # 写失败表
    out_fail = summary_dir / FAILED_TSV
    head_fail = ["OG", "foreground", "stage", "status", "reason", "one_mlc", "two_mlc"]
    with open(out_fail, "w", encoding="utf-8") as w:
        w.write("\t".join(head_fail) + "\n")
        for r in rows_failed:
            w.write("\t".join(str(r.get(k, "")) for k in head_fail) + "\n")

    # 写 foreground_gene_map
    out_fgmap = summary_dir / FG_MAP_TSV
    head_fgmap = [
        "OG",
        "foreground",
        "foreground_species",
        "foreground_gene_ids",
        "foreground_gene_names",
    ]
    with open(out_fgmap, "w", encoding="utf-8") as w:
        w.write("\t".join(head_fgmap) + "\n")
        for r in rows_fgmap:
            w.write("\t".join(str(r.get(k, "")) for k in head_fgmap) + "\n")

    write_done(summary_dir / DONE_FILENAME)

    n_ok = sum(1 for r in rows_all if r["is_testable"] == "TRUE")
    n_fail = len(rows_failed)

    log.info(f"[DONE] 写出：{out_all}")
    log.info(f"[DONE] 写出：{out_fail}")
    log.info(f"[DONE] 写出：{out_fgmap}")
    log.info(f"[STAT] total={len(rows_all)} testable={n_ok} failed={n_fail}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断")
        sys.exit(130)
