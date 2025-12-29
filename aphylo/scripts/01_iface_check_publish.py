#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_iface_check_publish.py —— 发布包接口体检（修正版）
修复要点：
1) 正确解析 Newick，只提取“叶节点标签（tip labels）”，忽略分支长度与内部节点支持值；
2) 与严格 SCO 的 AA-MSA 表头逐一比对：要求“仅物种名（不含‘|’）且覆盖集合=物种树叶名集合”；
3) 屏幕与日志双写；所有目录/后缀/路径从 config.yaml 读取（零猜测）。

依赖：
- Python 3.8+
- PyYAML（yaml）
"""

from __future__ import annotations
import sys, io, logging, re
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable
import yaml

DEFAULT_CONFIG = "config.yaml"

# ------------------ 基础工具 ------------------
def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True); return p

def need_dir(p: Path, what: str):
    p = Path(p)
    if not p.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p

def need_file(p: Path, what: str):
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
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

# ------------------ FASTA 与文件列举 ------------------
def read_fasta_headers(path: Path) -> List[str]:
    return [line[1:].strip() for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]

def list_ogs(msa_dir: Path, suffix: str) -> List[Tuple[str, Path]]:
    files = sorted(msa_dir.glob(f"OG*{suffix}"))
    pat = re.compile(r"^(OG\d+)")
    out = []
    for p in files:
        m = pat.match(p.name)
        if not m:
            raise ValueError(f"[ERR] 非法文件名（应以 OG 开头）：{p.name}")
        out.append((m.group(1), p))
    return out

# ------------------ 关键修复：严格解析 Newick 叶节点 ------------------
def parse_newick_tips(nwk_text: str) -> List[str]:
    """
    仅提取“叶节点标签（tip labels）”：
    - 忽略分支长度（出现在 ':' 之后的数字/科学计数法）；
    - 忽略内部节点标签/支持值（出现在 ')' 之后的 token，如 100、100/100、98.1/97 等）；
    - 兼容带引号的标签（'Label with spaces'）。

    实现思路（有限状态解析）：
      * 维护上一个分隔符 prev_delim ∈ {'(', ')', ',', ':', ';', None}
      * 收到一个 token 时：
          - 若 prev_delim == ')' ：这是内部节点标签 ⇒ 丢弃
          - 若 prev_delim == ':' ：这是分支长度 ⇒ 丢弃
          - 否则视为叶标签，追加到列表
      * token 的“纯数字 / 数字/数字”也会落在上述丢弃规则中（通常在 ')' 或 ':' 之后），
        但为保险再做一次数值/比值的兜底过滤。
    """
    tips: List[str] = []
    token = []
    in_quote = False
    prev_delim = None  # 最近一次遇到的结构字符：'(', ')', ',', ':', ';'

    def commit_token():
        nonlocal token, prev_delim, tips
        if not token:
            return
        s = "".join(token).strip()
        token = []
        if not s:
            return
        # 过滤典型的支持值/分支长度：纯数字或 数字/数字（含小数）
        if re.fullmatch(r"[\d.]+(?:/[\d.]+)?", s):
            return
        # 依据上一个分隔符判断类型
        if prev_delim == ")":
            # 内部节点标签（支持值/命名）——忽略
            return
        if prev_delim == ":":
            # 分支长度 —— 忽略
            return
        # 其余情况是叶标签
        tips.append(s)

    i = 0
    n = len(nwk_text)
    while i < n:
        ch = nwk_text[i]
        if in_quote:
            if ch == "'":
                # 结束引号 token
                in_quote = False
                commit_token()
            else:
                token.append(ch)
        else:
            if ch == "'":
                # 开始引号 token
                in_quote = True
            elif ch in ("(", ")", ",", ":", ";"):
                # 碰到结构分隔符，先提交当前 token，再更新 prev_delim
                commit_token()
                prev_delim = ch
            elif ch.isspace():
                # 空白也可作为 token 结束
                commit_token()
            else:
                token.append(ch)
        i += 1
    # 收尾
    commit_token()
    return sorted(set(tips))

# ------------------ 主流程 ------------------
def main():
    cfg = load_config()
    paths = cfg["paths"]; inputs = cfg["inputs"]
    logs_dir = Path(paths["logs_dir"])
    LOG_FILE = logs_dir / "01_iface_check_publish.log"
    log = get_logger("aphylo.01_iface", LOG_FILE)

    banner(log, "APhylo 01 —— 发布包接口体检（修正版）")

    iface_dir = ensure_dir(Path(paths["iface_dir"]))

    # —— 固定输入契约（零猜测）——
    sco_msa_dir    = need_dir(Path(inputs["sco_msa_dir"]), "严格SCO对齐目录")
    sco_msa_suffix = inputs["sco_msa_suffix"]
    pep2cds_map    = need_file(Path(inputs["pep2cds_map"]), "pep→cds映射表")
    species_tree   = need_file(Path(inputs["species_tree"]), "物种树(NWK)")

    fam_key = "family" if "family" in inputs else ("family_tsv" if "family_tsv" in inputs else None)
    if not fam_key:
        raise KeyError("[ERR] 缺少 inputs.family（或 family_tsv）")
    family_tsv = need_file(Path(inputs[fam_key]), "基因家族矩阵TSV")

    cafe = cfg.get("cafe5", {})
    utree_path = inputs.get("ultrametric_tree")
    cafe = cfg.get("cafe5", {})
    utree_path = inputs.get("ultrametric_tree")
    if utree_path:
        p = Path(utree_path)
        if p.is_file():
            ultrametric_tree = p
        else:
            log.warning(f"[WARN] 检测到需要{p}，但文件不存在；将由 MCMCTree 产出，本步骤忽略。")
            ultrametric_tree = None
    else:
        log.warning("[WARN] 未提供 inputs.ultrametric_tree；按方案A由 MCMCTree 产出，本步骤忽略。")
        ultrametric_tree = None

    # —— 解析物种树的“叶集合”（修复点）——
    nwk_text = Path(species_tree).read_text(encoding="utf-8")
    leaves = set(parse_newick_tips(nwk_text))
    if not leaves:
        raise ValueError("[ERR] 未能从物种树解析到叶节点标签，请检查 species_tree.nwk")

    # —— 遍历 OG 并体检 ——
    ogs = list_ogs(sco_msa_dir, sco_msa_suffix)
    if not ogs:
        raise RuntimeError("[ERR] 未发现 OG 文件，请检查发布包与后缀设定")

    bad_headers = []     # 含有“|”的非法表头
    cov_mismatch = []    # 覆盖集合与树不一致
    for og, p in ogs:
        heads = read_fasta_headers(p)
        # 表头必须仅为“物种名”（不能带 SeqID 等）
        for h in heads:
            if "|" in h:
                bad_headers.append((og, h))
        if set(heads) != leaves:
            cov_mismatch.append((og,
                                 sorted(set(heads) - leaves),
                                 sorted(leaves - set(heads))))
    if bad_headers:
        example = bad_headers[:5]
        raise ValueError(f"[ERR] 发现非物种名表头（含 '|'），示例：{example}")
    if cov_mismatch:
        # 严格 SCO 必须完全覆盖
        example = cov_mismatch[:3]
        raise ValueError(f"[ERR] OG 覆盖与树不一致（严格 SCO 应为全覆盖），示例：{example}")

    # —— 写简要报告（也同步打印）——
    rep = iface_dir / "iface_report.tsv"
    with rep.open("w", encoding="utf-8") as w:
        w.write("key\tvalue\n")
        w.write(f"sco_msa_dir\t{sco_msa_dir}\n")
        w.write(f"sco_msa_suffix\t{sco_msa_suffix}\n")
        w.write(f"pep2cds_map\t{pep2cds_map}\n")
        w.write(f"species_tree\t{species_tree}\n")
        w.write(f"ultrametric_tree\t{ultrametric_tree if ultrametric_tree else ''}\n")
        w.write(f"family_tsv\t{family_tsv}\n")
        w.write(f"n_og\t{len(ogs)}\n")
        w.write(f"n_leaves\t{len(leaves)}\n")

    log.info("=== 接口报告（iface_report.tsv）===")
    for line in rep.read_text(encoding="utf-8").splitlines():
        log.info(line)

    write_done(iface_dir / ".iface.done")
    log.info("[DONE] 接口验收完成（Newick 解析修复版）")

# ------------------ 入口 ------------------
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)

