#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
11_cafe5_prepare_inputs.py —— 为 CAFE5 准备 family.tsv & 超时钟树

核心职责：
  1) 读取原始基因家族矩阵（OrthoFinder 等产生的 TSV）；
  2) 按超时钟树中的物种叶节点顺序，清洗并重排 family：
       - 删除 Total 等汇总列；
       - 物种名按 alias_map 归一；
       - 只保留树上出现的物种；
  3) 输出符合 CAFE5 要求的 family.tsv：
       - 第 1 列：Desc（描述列，占位，用 "null"）
       - 第 2 列：Orthogroup/FamilyID（家族 ID）
       - 后续列：各物种拷贝数（顺序与树一致）；
  4) 写出“清洁树”：叶名按 alias_map & TimeTree 尾缀清理到标准写法；
  5) 写入 .prep.done 哨兵。
"""

from __future__ import annotations
import sys, io, re, logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable, Optional
import yaml

DEFAULT_CONFIG = "config.yaml"

# ====================== 通用工具：配置、日志、保障 ======================

def _expand_publish_placeholders(obj, publish_dir: str):
    """把对象中出现的 <publish_dir> 占位符替换为真路径。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    """加载 config.yaml，并展开 publish_dir 占位符。"""
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    """确保目录存在。"""
    p.mkdir(parents=True, exist_ok=True)
    return p

def need_file(p: Path, what: str) -> Path:
    """确保文件存在，否则报错。"""
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """日志写文件 + 同步屏幕；并保持 stdout/stderr 实时刷新。"""
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name); lg.setLevel(level); lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt); fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout);     sh.setFormatter(fmt); sh.setLevel(level)
    lg.addHandler(fh); lg.addHandler(sh)

    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x): 
            self.s.write(x); self.s.flush(); 
            return len(x)
    sys.stdout = _Flush(sys.stdout); sys.stderr = _Flush(sys.stderr)
    return lg

def banner(logger: logging.Logger, text: str):
    """在日志中打印 banner。"""
    bar = "=" * max(10, len(text) + 2)
    logger.info(bar); logger.info(f" {text} "); logger.info(bar)

def write_done(path: Path):
    """写哨兵文件（空文件即可）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

# ====================== 专用工具：树和物种名处理 ======================

def extract_tree_leaves(nwk_text: str) -> List[str]:
    """
    从 Newick 字符串中提取叶节点标签（保持首次出现顺序，去重）。
    不依赖复杂解析；通过删除分支长度与括号等，再按逗号拆分。
    """
    # 去掉 :branch_length 部分
    s = re.sub(r":[0-9eE+\-\.]+", "", nwk_text)
    # 去掉括号、分号和空白
    s = re.sub(r"[\(\)\s;]", "", s)
    toks = [t for t in s.split(",") if t]
    seen, leaves = set(), []
    for t in toks:
        if t not in seen:
            seen.add(t)
            leaves.append(t)
    return leaves

def normalize_timetree_suffix(name: str) -> str:
    """
    清理 TimeTree 叶名的引号编号尾缀：
      例：Crassostrea_gigas'13' -> Crassostrea_gigas
          Pomacea_canaliculata'35''34''43' -> Pomacea_canaliculata
    """
    name2 = re.sub(r"(?:'[0-9]+')+$", "", name)
    name2 = name2.replace("'", "")
    return name2

def build_alias_resolver(alias_cfg: Dict[str, Any]):
    """
    构建别名解析器，兼容两种写法：
      a) {canonical: [alias1, alias2, ...]}
      b) {alias: canonical}
    返回：f(name) -> canonical_name
    """
    if not alias_cfg:
        return lambda x: x
    a2c: Dict[str, str] = {}
    for k, v in alias_cfg.items():
        if isinstance(v, list):
            canonical = str(k)
            a2c[canonical] = canonical
            for a in v:
                a2c[str(a)] = canonical
        elif isinstance(v, str):
            a2c[str(k)] = str(v)
            a2c[str(v)] = str(v)
    return lambda x: a2c.get(x, x)

# ====================== 专用工具：TSV 读取 ======================

def parse_tsv_header(path: Path) -> List[str]:
    """读取 TSV 第一行并拆分为列名列表。"""
    first = path.read_text(encoding="utf-8").splitlines()[0]
    return first.rstrip("\n\r").split("\t")

def iter_tsv_rows(path: Path):
    """
    逐行读取 TSV（跳过首行），返回 (row_id, [counts...])。
    row_id 一般是 Orthogroup / FamilyID。
    """
    with path.open("r", encoding="utf-8") as f:
        _ = f.readline()  # 跳过表头
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            yield parts[0], parts[1:]

def safe_int(x: str) -> int:
    """计数字段安全转 int，异常按 0 处理（不改变总体尺度）。"""
    try:
        return int(x)
    except Exception:
        return 0

# ====================== 主流程 ======================

def main():
    cfg   = load_config()
    paths = cfg["paths"]
    inputs = cfg["inputs"]
    cafe  = cfg.get("cafe5", {})

    logs_dir = Path(paths["logs_dir"])
    # 为避免动皇上的路径，这里沿用旧命名（08），只是文件名不影响后续流程
    LOG_FILE = logs_dir / "08_cafe5_prepare_inputs.log"
    log = get_logger("aphylo.08", LOG_FILE)
    banner(log, "APhylo 08 — CAFE5 输入准备")

    # 若关闭 CAFE5，直接写 prep.done 并退出
    if not cafe.get("enable", True):
        log.info("cafe5.enable=false —— 跳过本步")
        write_done(Path(paths["cafe_run_dir"]) / ".prep.done")
        return

    # ------- 1. 定位输入 -------
    fam_key = "family" if "family" in inputs else ("family_tsv" if "family_tsv" in inputs else None)
    if not fam_key:
        raise KeyError("[ERR] 缺少 inputs.family（或 family_tsv）")
    family = need_file(Path(inputs[fam_key]), "家族矩阵")
    utree  = need_file(Path(inputs["ultrametric_tree"]), "超时钟树")

    outdir = ensure_dir(Path(paths["cafe_run_dir"]) / "input")

    # ------- 2. 构建物种别名解析器 -------
    alias_map = (cfg.get("species") or {}).get("alias_map", {})
    resolve = build_alias_resolver(alias_map)

    # ------- 3. 解析并标准化树叶集合 -------
    tree_text = Path(utree).read_text(encoding="utf-8")
    leaf_order_raw = extract_tree_leaves(tree_text)
    leaf_order: List[str] = []
    for x in leaf_order_raw:
        canon = resolve(normalize_timetree_suffix(x.strip()))
        leaf_order.append(canon)
    leaf_set = set(leaf_order)

    if not leaf_order:
        raise ValueError("[ERR] 超时钟树中未解析到任何叶节点，请检查 ultrametric_tree")

    log.info(f"[TREE] 叶节点数：{len(leaf_order)}；示例：{', '.join(leaf_order[:5])}")

    # ------- 4. 读取 family 表头并清洗列名 -------
    header = parse_tsv_header(family)  # 例如 ["Orthogroup", sp1, sp2, ..., "Total"]
    if not header or len(header) < 2:
        raise ValueError("[ERR] family.tsv 表头列数异常（<2）")

    firstcol = header[0]
    if firstcol not in ("Orthogroup", "Orthogroups", "OG", "Family", "FamilyID"):
        log.warning(f"[WARN] 非典型家族ID列名：{firstcol}，将仍按“家族ID列”处理。")

    # 剩余列视为物种列 + 可能的 Total 列
    data_cols_raw = header[1:]

    # 去掉常见的 Total 列名
    total_like = {"Total", "TotalGenes", "Total_Genes", "Total_genes"}
    data_cols_no_total: List[str] = [c for c in data_cols_raw if c not in total_like]

    if not data_cols_no_total:
        raise ValueError("[ERR] 除家族ID外未检测到任何物种列，请检查 family.tsv 表头。")

    # 对物种列做别名归一
    col_canon: List[str] = [resolve(c.strip()) for c in data_cols_no_total]

    # 映射：canonical_species -> 这些列的索引（可能多个列同属一个物种，后续会相加）
    idx_by_sp: Dict[str, List[int]] = {}
    for idx, sp in enumerate(col_canon):
        idx_by_sp.setdefault(sp, []).append(idx)

    # ------- 5. 检查：树上的物种是否都在 family 中出现 -------
    missing_in_family = sorted(leaf_set - set(col_canon))
    if missing_in_family:
        msg = ", ".join(missing_in_family[:10])
        raise RuntimeError(
            "[ERR] family.tsv 中缺少以下树上物种的列，请检查物种命名与 alias_map："
            f"{msg}（共 {len(missing_in_family)} 个）"
        )

    # 最终使用的物种列顺序：严格按树的叶顺序来（如果树上有但 family 中没有，前面已报错）
    target_species: List[str] = [sp for sp in leaf_order if sp in idx_by_sp]

    log.info(f"[FAMILY] 原始物种列数：{len(data_cols_no_total)}；"
             f"树上物种数：{len(leaf_order)}；"
             f"最终使用物种列数：{len(target_species)}")

    # ------- 6. 写出新的 family.tsv（带 Desc 列） -------
    out_family = outdir / Path(family).name
    n_rows = 0
    kept_rows = 0

    with out_family.open("w", encoding="utf-8") as w:
        # 表头：Desc + 家族ID列名统一用 Orthogroup（仅列名，内部不改 ID）
        w.write("Desc\tOrthogroup\t" + "\t".join(target_species) + "\n")

        for og, vals_all in iter_tsv_rows(family):
            n_rows += 1
            # 与 data_cols_no_total 对齐：只取前 len(data_cols_no_total) 个值
            if len(vals_all) < len(data_cols_no_total):
                # 不足时右侧补零（非常保守的做法）
                vals = vals_all + ["0"] * (len(data_cols_no_total) - len(vals_all))
            else:
                vals = vals_all[:len(data_cols_no_total)]

            # 每个物种的拷贝数 = 属于该物种的所有列之和
            counts_per_sp: List[str] = []
            for sp in target_species:
                idx_list = idx_by_sp[sp]
                s = 0
                for idx in idx_list:
                    s += safe_int(vals[idx])
                counts_per_sp.append(str(s))

            # 不做额外过滤，所有 OG 原样写出
            w.write("null\t" + str(og) + "\t" + "\t".join(counts_per_sp) + "\n")
            kept_rows += 1

    log.info(f"[OUT] 写出 family：{out_family} （共 {kept_rows} 行；原始 {n_rows} 行）")

    # ------- 7. 写出“清洁树”到 input 目录 -------
    clean_tree = tree_text
    for raw in leaf_order_raw:
        canon = resolve(normalize_timetree_suffix(raw.strip()))
        if raw != canon:
            # 简单替换：在超时钟树中用标准名替换原始叶标（一般不会产生歧义）
            clean_tree = clean_tree.replace(raw, canon)

    out_tree = outdir / "utree_for_cafe.nwk"
    out_tree.write_text(clean_tree, encoding="utf-8")
    log.info(f"[OUT] 写出清洁树：{out_tree}")

    # ------- 8. 写 sentinel -------
    write_done(Path(paths["cafe_run_dir"]) / ".prep.done")
    log.info("[DONE] CAFE5 输入准备完成")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)

