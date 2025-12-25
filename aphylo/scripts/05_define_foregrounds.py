#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
05_define_foregrounds.py — 生成前景集合（terminals-only）并输出打标树 FG.nwk（#1）
产物（写死，符合 aphylo 流水线）：
  - 前景清单：results/04_codeml/sets/{FG}.list
  - 打标树：  results/04_codeml/sets/{FG}.nwk   （把 {FG}.list 里的每个叶子分支标为 #1）
  - 哨兵：    results/04_codeml/.fgsets.done
约束：
  - 仅使用 terminals 模式（每个叶子的“外连分支”打 #1，允许多个 #1、允许非单系）。
  - 不做自动剪枝；树叶名需与 {FG}.list 中物种名逐字一致（大小写、下划线）。
新增（可选）：
  - codeml.add_phylip_header: 是否为 {FG}.nwk 首行写入 “N 1”。默认 False，开启后按叶子数写 N。
"""

from __future__ import annotations
import sys, io, logging, re
from pathlib import Path
from typing import Dict, Any, List
import yaml
import subprocess  # 调用 gotree

DEFAULT_CONFIG = "config.yaml"

# ================= 基础工具（与项目其余脚本一致） =================
def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(config_path)
    if not p.exists(): raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub: cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True); return p

def need_file(p: Path, what: str):
    p = Path(p)
    if not p.is_file(): raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
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
    bar = "=" * max(10, len(text)+2); logger.info(bar); logger.info(f" {text} "); logger.info(bar)

def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

# ============== 利用 gotree 去除分支长度（供 CODEML 使用） ==============
def clear_branch_lengths_with_gotree(nwk: str, logger: logging.Logger) -> str:
    """
    使用 `gotree brlen clear` 去除分支长度，只保留树拓扑。
    输入和输出都是 Newick 字符串，不在磁盘上额外产生文件。
    """
    try:
        proc = subprocess.run(
            ["gotree", "brlen", "clear"],
            input=nwk.encode("utf-8"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
    except FileNotFoundError:
        raise RuntimeError("[ERR] 未找到 gotree，可执行文件不在 PATH 中，请先安装/加载 gotree")

    if proc.returncode != 0:
        err_msg = proc.stderr.decode("utf-8", errors="ignore")
        logger.error(err_msg.strip() or "[gotree] 无更多错误信息")
        raise RuntimeError(f"[ERR] gotree brlen clear 失败，返回码 {proc.returncode}")

    out = proc.stdout.decode("utf-8", errors="ignore").strip()
    if not out:
        raise RuntimeError("[ERR] gotree brlen clear 返回空输出，请检查输入物种树是否为合法 Newick 格式")

    # 确保以分号结尾，给后续解析/软件一个稳定格式
    out = out.strip()
    if not out.endswith(";"):
        out = out.rstrip(";") + ";"
    return out

# ================= Newick 处理（terminals-only 打标） =================
def extract_leaf_names(nwk: str) -> List[str]:
    """
    提取树的叶子名：在 '(' 或 ',' 之后出现的 token，且后面紧跟 ':', ')', ',', ';' 之一。
    """
    nwk = nwk.replace("#1", "")  # 清理既有标签，避免干扰
    pat = re.compile(r'(?<=\(|,)\s*([^\s():,;]+)\s*(?=[:),;])')
    names = pat.findall(nwk)
    seen=set(); out=[]
    for n in names:
        if n not in seen:
            seen.add(n); out.append(n)
    return out

def label_terminals(nwk: str, tips: List[str]) -> str:
    """
    把指定叶子的外连分支标成 #1（terminals）。
    通过在叶子名后追加 '#1' 实现，例如 'Homo_sapiens' -> 'Homo_sapiens#1'。
    """
    base = nwk.replace("#1", "")  # 幂等化
    for tip in tips:
        pat = re.compile(rf'(?<=\(|,)\s*({re.escape(tip)})\s*(?=[:),;])')
        base, n = pat.subn(r'\1#1', base, count=1)
        if n == 0:
            raise KeyError(f"[ERR] 无法在树上定位叶子：{tip}（命名可能与树不一致）")
    return base

# ============== 写入工具（可选 PHYLIP 头） ==============
def write_nwk_with_optional_phylip_header(
    out_path: Path,
    nwk_line: str,
    add_header: bool,
    n_leaves: int,
    logger: logging.Logger
):
    """
    将 Newick 写入 out_path。
    - add_header=True 时：首行写入 "N 1"，第二行写 Newick。
    - 自动幂等：若文件首行已是数字数字（^\d+\s+\d+$），不重复加头。
    """
    ensure_dir(out_path.parent)

    # 确保只有一个分号结尾并换行
    nwk_line = (nwk_line.rstrip() if nwk_line else "")
    if not nwk_line.endswith(";"):
        nwk_line = nwk_line.rstrip(";") + ";"

    content = nwk_line + "\n"

    if add_header:
        header = f"{int(n_leaves)} 1"
        logger.info(f"[Info] PHYLIP header: {header}")
        content = header + "\n" + content

    out_path.write_text(content, encoding="utf-8")

# ================= 主流程 =================
def main():
    cfg = load_config()
    paths  = cfg["paths"]; inputs = cfg["inputs"]; codeml = cfg.get("codeml", {})

    logs_dir = Path(paths["logs_dir"]); LOG_FILE = logs_dir / "05_define_foregrounds.log"
    log = get_logger("aphylo.05", LOG_FILE)
    banner(log, "APhylo 05 — 定义前景集合（terminals-only 打标）")

    # 产物目录：自动创建（根与 sets）
    codeml_root = ensure_dir(Path(paths["codeml_dir"]))
    sets_dir    = ensure_dir(codeml_root / "sets")

    # 基准物种树（未打标）
    tree_path = need_file(Path(inputs["species_tree"]), "物种树 species_tree.nwk")
    raw_nwk = tree_path.read_text(encoding="utf-8").strip()
    raw_nwk = raw_nwk.rstrip(";") + ";"  # 保证只有一个分号结尾

    # 使用 gotree 去除分支长度，仅保留拓扑（更贴合 CODEML 教程做法）
    log.info("[Info] 使用 gotree brlen clear 去除物种树分支长度，仅保留拓扑供 CODEML 使用")
    nwk = clear_branch_lengths_with_gotree(raw_nwk, log)

    # 叶子集合 & 别名映射（可选）
    leaves = extract_leaf_names(nwk)
    alias  = cfg.get("species", {}).get("alias_map", {})

    # 前景集合（来自 config.codeml.foreground_sets）
    fgs = codeml.get("foreground_sets", [])
    if not fgs:
        raise ValueError("[ERR] codeml.foreground_sets 为空（至少提供一个前景集合）")

    # 是否写入 PHYLIP header（默认 False，不影响既有行为）
    add_header: bool = bool(codeml.get("add_phylip_header", False))
    if add_header:
        log.info("[Info] 将为 {FG}.nwk 写入 PHYLIP 头（N 1）")

    # 逐前景生成 {FG}.list 与 {FG}.nwk（terminals-only 打标）
    for fg in fgs:
        name: str = fg["name"]
        tips_in: List[str] = fg.get("tips", [])
        if not tips_in:
            raise ValueError(f"[ERR] 前景 {name} 的 tips 为空")

        tips_mapped: List[str] = []
        missing: List[str] = []
        for t in tips_in:
            t2 = alias.get(t, t)
            if t2 in leaves:
                tips_mapped.append(t2)
            else:
                missing.append(t)

        if missing:
            miss_str = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
            raise KeyError(f"[ERR] 前景 {name} 中以下物种不在树叶名中：{miss_str}")

        # 写 {FG}.list（标准化后的叶名，一行一个）
        (sets_dir / f"{name}.list").write_text("\n".join(tips_mapped) + "\n", encoding="utf-8")

        # 生成 terminals-only 打标树
        nwk_fg = label_terminals(nwk, tips_mapped)

        # 写 {FG}.nwk —— 按需加 PHYLIP header
        out_nwk = sets_dir / f"{name}.nwk"
        write_nwk_with_optional_phylip_header(
            out_path=out_nwk,
            nwk_line=nwk_fg,
            add_header=add_header,
            n_leaves=len(leaves),
            logger=log
        )

        log.info(f"前景集合：{name} → {len(tips_mapped)} tips；已生成 {name}.list 与 {name}.nwk（#1 打标完成，分支长度已清除）")

    # 哨兵
    write_done(codeml_root / ".fgsets.done")
    log.info("[DONE] 前景集合与打标树已全部生成 —— 06 将使用 sets/{FG}.nwk 作为 treefile")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)

