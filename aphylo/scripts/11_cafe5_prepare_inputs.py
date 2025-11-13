#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_cafe5_prepare_inputs.py — 准备 CAFE5 输入（保持 I/O 契约不变）
更新要点（与既有流水线风格一致，零破坏）：
1) 移除 family.tsv 的 Total 列；日志分别报告“原始列数(含 Total)”与“有效列数(剔除 Total)”，避免误判为脚本出错。
2) 统一物种别名（config.species.alias_map），并清理 TimeTree 叶名中的引号编号尾缀（如 '13'、'35''34'…），两侧一致化。
3) 仅保留“树∩矩阵”的物种列，按树叶顺序重排；别名合并产生重复时逐行求和。
4) 若树中存在但矩阵缺失的物种：明确报错，列出缺失清单；不补零以避免伪信号。
5) 输出：将清洗后的 family.tsv 与“清洁树”副本写入 paths.cafe_run_dir/input/；原始发布文件不改；完成后落 .prep.done。
"""
from __future__ import annotations
import sys, io, logging, re
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable
import yaml

DEFAULT_CONFIG = "config.yaml"

# ================= 基础工具 =================
def _expand_publish_placeholders(obj, publish_dir: str):
    """展开 config 中的 <publish_dir> 占位符（保持既有风格）。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    """读取主配置；若存在 publish_dir，则展开 inputs.* 中的占位符。"""
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def need_file(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """文件 + 屏幕双通道日志；与 aphylo 其他脚本风格一致。"""
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name); lg.setLevel(level); lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt); fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout);     sh.setFormatter(fmt); sh.setLevel(level)
    lg.addHandler(fh); lg.addHandler(sh)
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout = _Flush(sys.stdout); sys.stderr = _Flush(sys.stderr)
    return lg

def banner(logger: logging.Logger, text: str):
    bar = "=" * max(10, len(text) + 2)
    logger.info(bar); logger.info(f" {text} "); logger.info(bar)

def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

# ================= 业务工具 =================
def extract_tree_leaves(nwk_text: str) -> List[str]:
    """
    解析 Newick 叶节点名称（轻量实现）：
    1) 去掉分支长度(:0.123)；2) 去括号/空白/分号；3) 按逗号切分；4) 过滤空串；5) 去重保持顺序。
    """
    s = re.sub(r":[0-9eE+\-\.]+", "", nwk_text)
    s = re.sub(r"[\(\)\s;]", "", s)
    toks = [t for t in s.split(",") if t]
    seen, leaves = set(), []
    for t in toks:
        if t not in seen:
            seen.add(t); leaves.append(t)
    return leaves

def normalize_timetree_suffix(name: str) -> str:
    """
    清理 TimeTree 叶名的引号编号尾缀：
    例：Crassostrea_gigas'13' -> Crassostrea_gigas
        Pomacea_canaliculata'35''34''43' -> Pomacea_canaliculata
    """
    # 去除结尾处重复出现的 '数字'
    name2 = re.sub(r"(?:'[0-9]+')+$", "", name)
    # 去除可能残留的单引号
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

def parse_tsv_header(path: Path) -> List[str]:
    """读取 TSV 第一行并拆分为列名列表（不去除首列）。"""
    first = path.read_text(encoding="utf-8").splitlines()[0]
    return first.rstrip("\n\r").split("\t")

def iter_tsv_rows(path: Path):
    """逐行读取 TSV（跳过首行），返回 (row_id, [counts...])。"""
    with path.open("r", encoding="utf-8") as f:
        _ = f.readline()
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            yield parts[0], parts[1:]

def safe_int(x: str) -> int:
    """计数字段安全转 int，异常按 0 处理（并不改变数据规模）。"""
    try:
        return int(x)
    except Exception:
        return 0

# ================= 主流程 =================
def main():
    cfg = load_config()
    paths = cfg["paths"]; inputs = cfg["inputs"]; cafe = cfg.get("cafe5", {})
    logs_dir = Path(paths["logs_dir"])
    LOG_FILE = logs_dir / "08_cafe5_prepare_inputs.log"
    log = get_logger("aphylo.08", LOG_FILE)
    banner(log, "APhylo 08 — CAFE5 输入准备")

    # 允许关闭本步
    if not cafe.get("enable", True):
        log.info("cafe5.enable=false —— 跳过本步")
        write_done(Path(paths["cafe_run_dir"]) / ".prep.done")
        return

    # 定位输入
    fam_key = "family" if "family" in inputs else ("family_tsv" if "family_tsv" in inputs else None)
    if not fam_key:
        raise KeyError("[ERR] 缺少 inputs.family（或 family_tsv）")
    family = need_file(Path(inputs[fam_key]), "家族矩阵")
    utree  = need_file(Path(inputs["ultrametric_tree"]), "超时钟树")

    outdir = ensure_dir(Path(paths["cafe_run_dir"]) / "input")

    # 别名与标准化器
    alias_map = (cfg.get("species") or {}).get("alias_map", {})
    resolve = build_alias_resolver(alias_map)

    # ---------- 解析并标准化树叶集合 ----------
    tree_text = Path(utree).read_text(encoding="utf-8")
    leaf_order_raw = extract_tree_leaves(tree_text)
    # 清理 TimeTree 尾缀 + 去引号 + 别名标准化
    leaf_order = [resolve(normalize_timetree_suffix(x.strip())) for x in leaf_order_raw]
    leaf_set = set(leaf_order)

    # ---------- 读取家族矩阵表头并清洗列名 ----------
    header = parse_tsv_header(family)  # 例如 ["Orthogroup", sp1, sp2, ..., "Total"]
    firstcol = header[0]
    # 兼容 OrthoFinder 的少量变体：Orthogroups -> Orthogroup
    if firstcol not in ("Orthogroup", "OG", "Orthogroups"):
        raise ValueError(f"[ERR] family.tsv 首列应为 Orthogroup（或 Orthogroups/OG），当前为：{firstcol}")
    if firstcol == "Orthogroups":
        header[0] = "Orthogroup"

    raw_species_all = header[1:]  # 含 Total
    # 剔除 Total（大小写无关；去前后空白）
    dropped_total = [c for c in raw_species_all if c.strip().lower() == "total"]
    species_wo_total = [c for c in raw_species_all if c.strip().lower() != "total"]

    # 应用别名映射（矩阵侧）
    renamed_pairs = []  # (old, new)
    species_norm = []
    for c in species_wo_total:
        c0 = c.strip()
        c1 = resolve(c0)
        if c1 != c0:
            renamed_pairs.append((c0, c1))
        species_norm.append(c1)

    # 统计信息（更准确的口径）
    log.info(f"[Info] 原始物种列数（含 Total，不含 Orthogroup）: {len(raw_species_all)}")
    log.info(f"[Info] 有效物种列数（剔除 Total 后）: {len(species_norm)}")
    if dropped_total:
        log.info(f"[Clean] 移除汇总列: {dropped_total}")
    if renamed_pairs:
        log.info(f"[Clean] 别名标准化（矩阵列）{len(renamed_pairs)} 项，示例前 5 条: {renamed_pairs[:5]}")

    # 记录矩阵中“树外”的列（将被丢弃，不进入 CAFE）
    extra_in_matrix = [c for c in species_norm if c not in leaf_set]
    if extra_in_matrix:
        log.info(f"[Clean] 丢弃树外物种列（不进入 CAFE）: {extra_in_matrix}")

    kept_in_matrix = [c for c in species_norm if c in leaf_set]

    # 缺失检查：树中存在但矩阵缺失（应用别名 + TimeTree 清理后）
    missing_in_matrix = [s for s in leaf_order if s not in set(kept_in_matrix)]
    if missing_in_matrix:
        raise ValueError(f"[ERR] 超时钟树中的物种在 family.tsv 中缺失（请在 config.species.alias_map 补别名或检查发布包）：{missing_in_matrix}")

    # ---------- 构建目标列顺序：严格按树叶顺序 ----------
    target_cols = list(leaf_order)

    # 构建“标准名 -> 原始索引列表”（用于别名合并求和）
    idx_map: Dict[str, List[int]] = {}
    for idx, cname in enumerate(species_norm):
        if cname in leaf_set:  # 仅收集树内物种
            idx_map.setdefault(cname, []).append(idx)

    # ---------- 逐行清洗并输出新 family.tsv ----------
    out_family = Path(outdir) / Path(family).name
    with out_family.open("w", encoding="utf-8") as w:
        # 写新的表头：Orthogroup + 按树顺序排列的物种列
        w.write("Orthogroup\t" + "\t".join(target_cols) + "\n")
        for og, row_vals in iter_tsv_rows(family):
            # row_vals 与 species_wo_total 对齐；聚合到 target_cols
            # 为所有目标物种准备 0
            agg: Dict[str, int] = {sp: 0 for sp in target_cols}
            # 对于每个标准物种，把映射到它的所有原始列求和
            for sp, id_list in idx_map.items():
                s = 0
                for j in id_list:
                    if j < len(row_vals):
                        s += safe_int(row_vals[j])
                agg[sp] = s
            # 按目标顺序写出
            w.write(og + "\t" + "\t".join(str(agg[sp]) for sp in target_cols) + "\n")

    # ---------- 写入“清洁树”副本到 input/（不改原文件） ----------
    # 仅清理 TimeTree 尾缀与引号，其余保持原状；别名标准化仅用于匹配，不改树文本
    # 出于可追溯性，这里只做“去引号编号”的清洁版副本，文件名与原名保持一致
    leaves_cleaned_text = tree_text
    # 清除所有叶名中的 '数字' 尾缀（多次重复也清除）
    leaves_cleaned_text = re.sub(r"(?:'[0-9]+')+", "", leaves_cleaned_text)
    out_tree = Path(outdir) / Path(utree).name
    out_tree.write_text(leaves_cleaned_text, encoding="utf-8")

    # ---------- 总结与收尾 ----------
    log.info(f"[OK] 输出家族矩阵：{out_family}")
    log.info(f"[OK] 输出超时钟树：{out_tree}（已移除 TimeTree 引号编号尾缀）")
    log.info(f"[Summary] 目标物种数量：{len(target_cols)}；被丢弃的矩阵列数量（含 Total 与树外列）：{len(dropped_total) + len(extra_in_matrix)}")
    write_done(Path(paths["cafe_run_dir"]) / ".prep.done")
    log.info("[DONE] CAFE5 输入准备完成")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)

