#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_clean_tree_and_noblize.py —— 定根+净化（保留分支长度）并生成“无分支长度”的两行树
职责非常单一：只做 3 件事
  1) gotree 外群定根 + 清支持值/注释 + 规范 newick（得到 rooted_tree，保留分支长度）
  2) 写两行的 species_calib.trees（第一行 N 1，第二行为保留分支长度的 newick）
  3) 生成无分支长度的 species_calib.nobl.trees（第一行 N 1，第二行无任何 :length）

严格口径：
  - 不读取/不输出 calibration TSV
  - 不回填 B(L,U)；您后续需要时，直接在 .nobl.trees 第二行手动插入 'B(L,U)'
  - 外群对不上 → 硬失败（不做 midpoint 兜底）
  - 所有路径从 config.yaml:mcmctree 读取；不猜测、不扫描目录

config.yaml 期望键（位于 mcmctree 段）：
  mcmctree:
    unrooted_tree: "../phylo/results/publish/aphylo_ready/species_tree.nwk"
    rooted_tree:   "../phylo/results/publish/aphylo_ready/species_tree_rooted.nwk"
    source_tree:   "results/06_cafe/mcmctree/species_calib.trees"
    tree:          "results/06_cafe/mcmctree/species_calib.nobl.trees"
    outgroups:     ["Nematostella_vectensis"]
    outdir:        "results/06_cafe/mcmctree"

用法：
  直接运行（本脚本不接受命令行参数）
    python 08_clean_tree_and_noblize.py
"""

from __future__ import annotations
from pathlib import Path
from datetime import datetime
import subprocess, sys, re, csv

# ===================== 顶部参数（可被 config.yaml:mcmctree 覆盖） =====================
CONFIG_PATH       = Path("config.yaml")                          # 主配置
OUTDIR_DEFAULT    = Path("results/06_cafe/mcmctree")             # 输出目录（用于 .done 与默认日志）
LOG_PATH          = Path("logs/08_clean_tree_and_noblize.log")   # 日志路径（固定）
# =====================================================================================

# ----------------------------- YAML 读取（ruamel 优先） -----------------------------
def load_yaml_dict(p: Path):
    """读取 YAML 为 dict；优先 ruamel.yaml，失败回退到 PyYAML；读取失败返回 {}。"""
    if not p.exists():
        return {}
    try:
        from ruamel.yaml import YAML
        y = YAML(typ='rt'); y.preserve_quotes = True
        return y.load(p.read_text(encoding="utf-8")) or {}
    except Exception:
        try:
            import yaml
            return yaml.safe_load(p.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}

# -------------------------------- 简易日志器 --------------------------------
class Logger:
    """同时写屏幕与文件；文件换行统一为 LF"""
    def __init__(self, p: Path):
        p.parent.mkdir(parents=True, exist_ok=True)
        self.f = p.open("w", encoding="utf-8", newline="\n")
        self.write(f"[{datetime.now()}] 08_clean_tree_and_noblize 启动")

    def write(self, s: str):
        print(s)
        self.f.write(s + ("\n" if not s.endswith("\n") else ""))
        self.f.flush()

    def close(self):
        self.write(f"[{datetime.now()}] 完成"); self.f.close()

# -------------------------------- 子进程封装 --------------------------------
def run_capture(cmd, logger: Logger, input_text: str | None = None) -> str:
    """
    执行命令并捕获 stdout/stderr；非零返回码即硬失败。
    为避免日志过长，stdout 只保留前 2000 字节做摘要。
    """
    logger.write(f"[CMD] {' '.join(cmd)}")
    p = subprocess.run(
        cmd,
        input=input_text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if p.returncode != 0:
        if p.stderr:
            logger.write("[STDERR]\n" + p.stderr.strip())
        raise SystemExit(f"[ERR] 命令失败：{' '.join(cmd)}（详见日志）")
    if p.stderr.strip():
        logger.write("[STDERR]\n" + p.stderr.strip())
    if p.stdout:
        s = p.stdout
        if len(s) > 2000:
            s = s[:2000] + "...\n[截断]"
        logger.write("[STDOUT]\n" + s)
    return p.stdout

# --------------------------- Newick 解析（仅为去长度服务） ---------------------------
class Node:
    __slots__=("name","children","parent","length")
    def __init__(self, name=None, length=None):
        self.name=name; self.children=[]; self.parent=None; self.length=length
    def is_leaf(self): return not self.children

def _skip_ws(s, i):
    n=len(s)
    while i<n and s[i].isspace(): i+=1
    return i

def _parse_name(s, i):
    start=i; n=len(s)
    while i<n and s[i] not in "(),:;": i+=1
    return s[start:i].strip(), i

def _parse_number_after_colon(s, i):
    i=_skip_ws(s, i)
    m=re.match(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', s[i:])
    if not m: return None, i
    num=float(m.group(1)); i += len(m.group(1))
    return num, i

def parse_newick_recursive(s: str, i: int=0):
    i=_skip_ws(s, i)
    if i>=len(s): raise ValueError("Newick 解析到达末尾")
    if s[i] == '(':
        i+=1
        children=[]
        while True:
            child, i = parse_newick_recursive(s, i)
            children.append(child); child.parent=None
            i=_skip_ws(s, i)
            if i>=len(s): raise ValueError("缺少 ')'")
            if s[i]==',':
                i+=1; continue
            elif s[i]==')':
                i+=1; break
            else:
                raise ValueError(f"非法字符：{s[i]}")
        i=_skip_ws(s, i)
        name=None
        if i<len(s) and s[i] not in ":,);":
            name, i = _parse_name(s, i)
        i=_skip_ws(s, i)
        length=None
        if i<len(s) and s[i]==':':
            i+=1; length, i = _parse_number_after_colon(s, i)
        node=Node(name=name, length=length)
        for ch in children: ch.parent=node; node.children.append(ch)
        return node, i
    else:
        name, i = _parse_name(s, i)
        i=_skip_ws(s, i)
        length=None
        if i<len(s) and s[i]==':':
            i+=1; length, i = _parse_number_after_colon(s, i)
        return Node(name=name, length=length), i

def parse_newick(nwk: str) -> Node:
    """将 newick 文本解析为树；先剔除已有的 'B(...)' 注释以免干扰解析。"""
    s = re.sub(r"'B\([^']*\)'", "", nwk.strip())
    root, i = parse_newick_recursive(s, 0)
    i=_skip_ws(s, i)
    if i < len(s) and s[i]==';': i+=1
    while root and getattr(root, "parent", None): root=root.parent
    return root

def iter_leaves(n: Node):
    if n.is_leaf(): 
        yield n
    else:
        for ch in n.children:
            yield from iter_leaves(ch)

def leaf_names(n: Node): 
    return [x.name for x in iter_leaves(n)]

def to_newick_no_lengths(node: Node) -> str:
    """把树重新渲染为不带分支长度的新 newick（仅拓扑+叶名，不写内部名/长度）。"""
    if node.is_leaf():
        return node.name
    inner=",".join(to_newick_no_lengths(ch) for ch in node.children)
    return f"({inner})"

# ----------------------------------- 主流程 -----------------------------------
def main():
    log = Logger(LOG_PATH)
    try:
        # 1) 读取配置
        cfg = load_yaml_dict(CONFIG_PATH)
        m = (cfg.get("mcmctree") or {})
        # 必填键（不猜路径）
        UNROOTED = Path(m.get("unrooted_tree", ""))     # 未定根树
        ROOTED   = Path(m.get("rooted_tree", ""))       # 定根后 NWK（保留长度）
        SRC2L    = Path(m.get("source_tree", ""))       # 两行树（有长度）
        NOBL2L   = Path(m.get("tree", ""))              # 两行树（无长度）
        OUTDIR   = Path(m.get("outdir", OUTDIR_DEFAULT))
        OUTGROUP = list(m.get("outgroups", []))         # 外群列表（至少 1 个）

        # 基本断言
        if not UNROOTED or not ROOTED or not SRC2L or not NOBL2L:
            raise SystemExit("[ERR] config.yaml:mcmctree 中必须提供 unrooted_tree/rooted_tree/source_tree/tree/outdir/outgroups")
        OUTDIR.mkdir(parents=True, exist_ok=True)

        # gotree 检查
        if not shutil_which("gotree"):
            raise SystemExit("[ERR] 未找到 gotree，请先安装：mamba install -y -c bioconda gotree")

        if not UNROOTED.exists():
            raise SystemExit(f"[ERR] 未找到未定根树：{UNROOTED}")

        if not OUTGROUP:
            raise SystemExit("[ERR] 未提供外群（mcmctree.outgroups），本脚本不做 midpoint 兜底")

        log.write("[PATH] unrooted_tree = " + str(UNROOTED))
        log.write("[PATH] rooted_tree   = " + str(ROOTED))
        log.write("[PATH] source_tree   = " + str(SRC2L))
        log.write("[PATH] nobl_tree     = " + str(NOBL2L))
        log.write("[INFO] outgroups     = " + ", ".join(OUTGROUP))

        # 2) 定根+净化（保留长度）
        # 2.1 外群定根（外群匹配失败应硬退出）
        log.write("[1/5] 外群定根")
        rerooted = run_capture(["gotree","reroot","outgroup","-i",str(UNROOTED),"-o","-"]+OUTGROUP, log)

        # 2.2 清支持值
        log.write("[2/5] 清除内部支持值标签")
        cleared = run_capture(["gotree","support","clear","-i","-","-o","-"], log, input_text=rerooted)

        # 2.3 清注释（若版本不支持则跳过，保 WARN）
        log.write("[3/5] 清除注释（若不支持则跳过）")
        try:
            cleaned = run_capture(["gotree","comment","clear","-i","-","-o","-"], log, input_text=cleared)
        except SystemExit:
            log.write("[WARN] gotree 不支持 comment clear；已跳过")
            cleaned = cleared

        # 2.4 规范 newick，写 rooted_tree（保留分支长度）
        log.write("[4/5] 规范 newick 并写 rooted_tree（保留分支长度）")
        final_nwk = run_capture(["gotree","reformat","newick","-i","-","-o","-"], log, input_text=cleaned).strip()
        if not final_nwk.endswith(";"):
            final_nwk += ";"
        ROOTED.parent.mkdir(parents=True, exist_ok=True)
        ROOTED.write_text(final_nwk + "\n", encoding="utf-8")
        log.write(f"[OK] 已写：{ROOTED}")

        # 2.5 统计叶数 N（用 gotree stats tips）
        tips_text = run_capture(["gotree","stats","tips","-i",str(ROOTED)], log)
        tip_names = parse_tips_from_stats(tips_text)
        if not tip_names:
            raise SystemExit("[ERR] 无法解析叶名（gotree stats tips 输出为空）")
        N = len(tip_names)
        log.write(f"[OK] 叶子数 N = {N}")
        # 写两行有长度的 .trees
        SRC2L.parent.mkdir(parents=True, exist_ok=True)
        SRC2L.write_text(f"{N} 1\n{final_nwk}\n", encoding="utf-8")
        log.write(f"[OK] 已写两行（有长度）：{SRC2L}")
        # 哨兵
        (OUTDIR/".clean.done").write_text(str(datetime.now())+"\n", encoding="utf-8")

        # 3) 去分支长度（渲染为无长度 newick）
        log.write("[5/5] 生成无分支长度的两行树")
        root = parse_newick(final_nwk)           # 解析有长度的文本
        nobl_nwk = to_newick_no_lengths(root)    # 重建无长度的 newick
        nobl_nwk = re.sub(r"\s+", " ", nobl_nwk).strip()
        if not nobl_nwk.endswith(";"):
            nobl_nwk += ";"

        # 自检：确保没有 ':数字'
        if re.search(r":[0-9]", nobl_nwk):
            raise SystemExit("[ERR] 生成的无分支长度树仍含分支长度，请反馈日志排查。")

        NOBL2L.parent.mkdir(parents=True, exist_ok=True)
        NOBL2L.write_text(f"{N} 1\n{nobl_nwk}\n", encoding="utf-8")
        log.write(f"[OK] 已写两行（无长度）：{NOBL2L}")
        (OUTDIR/".nobl.done").write_text(str(datetime.now())+"\n", encoding="utf-8")

        log.write("[DONE] 全流程完成：定根+净化（有长度）→ 两行有长度 → 两行无长度")

    finally:
        log.close()

# -------------------------- 工具函数：which 与 tips 解析 --------------------------
def shutil_which(prog: str) -> bool:
    """跨平台检测可执行是否存在（不直接导入 shutil.which 以避免旧 Python 差异）"""
    try:
        import shutil as _sh
        return _sh.which(prog) is not None
    except Exception:
        return False

def parse_tips_from_stats(tips_text: str) -> list[str]:
    """
    兼容 gotree stats tips 的两种输出：
      a) 表格模式（首行含列名，至少含 'name'）
      b) 纯名单（一行一个名字）
    返回：按出现顺序的叶名列表（去重）
    """
    lines = [ln.strip() for ln in tips_text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split()
    names = []
    if any(tok.lower()=="name" for tok in header):
        idx = [i for i,t in enumerate(header) if t.lower()=="name"][0]
        for ln in lines[1:]:
            parts = ln.split()
            if len(parts) > idx:
                names.append(parts[idx])
    else:
        for ln in lines:
            parts = ln.split()
            names.append(parts[-1])
    # 去重保序
    seen=set(); ordered=[]
    for n in names:
        if n not in seen:
            seen.add(n); ordered.append(n)
    return ordered

# ----------------------------------- 入口 -----------------------------------
if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        # 让错误码/消息直观体现到终端，同时保持日志已写入
        raise
    except Exception as e:
        print(f"[FATAL] 未捕获异常：{e}", file=sys.stderr)
        sys.exit(99)

