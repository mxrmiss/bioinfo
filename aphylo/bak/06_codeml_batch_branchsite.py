#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
06_codeml_batch_branchsite.py —— 批量 codeml（branch-site，ALT/NULL；多线程由 config.yaml 控制）

输入（按 PDF 命名）：
  results/03_qc/keep_og.list
  results/03_codon/codon_msa/OG*.codon.fna
  results/04_codeml/sets/<FG>.nwk              # 05 产物（终端仅打 #1）
  templates/branch_site_alt.ctl
  templates/branch_site_null.ctl

输出（逐 OG×FG）：
  results/04_codeml/raw/OG/FG/alt/{seq.codon.fna,species_tree.nwk,alt.ctl,mlc.txt,...}
  results/04_codeml/raw/OG/FG/null/{seq.codon.fna,species_tree.nwk,null.ctl,mlc.txt,...}
  results/04_codeml/.codeml.done

稳态与容错：
  - 清理前景树内部标签/注释，仅保留末端物种名与 #1 打标
  - 物种集合比较时忽略 #1
  - 将 PAML 软失败（退出码!=0 但 mlc.txt 含 lnL）计为 soft_ok，不中断
"""

from __future__ import annotations
import sys, io, os, re, logging, subprocess, shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any, List, Tuple
import yaml

# ===== 固定路径（与 PDF 命名一致）=====
CONFIG_PATH = "config.yaml"
SEQ_DIR     = Path("results/03_codon/codon_msa")
QC_DIR      = Path("results/03_qc")
RAW_ROOT    = Path("results/04_codeml/raw")
SETS_DIR    = Path("results/04_codeml/sets")
LOG_FILE    = Path("logs/06_codeml_batch_branchsite.log")

# 断点续跑策略
FORCE_RERUN = False  # True=强制覆盖已完成；False=能跳过就跳过

# ========== 基础工具 ==========
def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(path: str = CONFIG_PATH) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists(): raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub: cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def need_file(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_file(): raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(logfile: Path, level=logging.INFO) -> logging.Logger:
    ensure_dir(logfile.parent)
    lg = logging.getLogger("aphylo.06"); lg.handlers.clear(); lg.setLevel(level)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt); fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout);     sh.setFormatter(fmt); sh.setLevel(level)
    lg.addHandler(fh); lg.addHandler(sh)
    # 立即刷新 stdout/stderr，便于 tail
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s=s
        def write(self, x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout=_Flush(sys.stdout); sys.stderr=_Flush(sys.stderr)
    return lg

def render_tpl(path: Path, rep: Dict[str, str]) -> str:
    txt = path.read_text(encoding="utf-8")
    for k, v in rep.items():
        txt = txt.replace(f"{{{{{k}}}}}", str(v))
    return txt

def has_lnl(mlc: Path) -> bool:
    """mlc.txt 内是否已有 lnL 行（判定软失败时使用）"""
    if not mlc.is_file(): return False
    try:
        with mlc.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("lnL("):  # PAML 的 lnL 行
                    return True
    except Exception:
        return False
    return False

def run_cmd_capture(cmd: List[str], cwd: Path, log: logging.Logger) -> int:
    """运行 codeml，流式打印输出，返回退出码；限制底层库线程=1"""
    log.info(f"CMD[{cwd}]: " + " ".join(map(str, cmd)))
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    p = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    for line in p.stdout: log.info(line.rstrip())
    return p.wait()

# ========== 树处理与集合校验 ==========
_TOK_TIP = re.compile(r'(?<=\(|,)\s*([^()\s:,;]+)\s*(?=[:),;])')

def sanitize_foreground_tree(raw_text: str) -> Tuple[str, Dict[str,int]]:
    """
    清理 Newick：
      - 去掉方括号注释
      - 去内部/分支标签（仅保留末端名与 #1）
      - 保持分号与括号完整
    返回：清理后的文本 + 统计信息
    """
    s0 = raw_text.strip()
    stats = {"removed_brackets":0, "removed_internal_labels":0}

    # 去掉方括号注释
    before = len(s0)
    s = re.sub(r'\[.*?\]', '', s0)
    stats["removed_brackets"] = before - len(s)

    # 去掉内部节点/分支上的名字（形如 )Name: 或 )Name, 或 ): 长度）
    def _rm_internal(m):
        stats["removed_internal_labels"] += 1
        return m.group(1)
    s = re.sub(r'(\))\s*[^():,;\s]+\s*(?=[:),;])', _rm_internal, s)

    # 压缩多余空白
    s = re.sub(r'\s+', ' ', s).strip()

    # 基本合法性
    if not s.endswith(';'): s += ';'
    if s.count('(') != s.count(')'):
        # 尽力修；若仍不匹配就抛错
        raise ValueError("[ERR] 物种树括号不匹配，无法用于 codeml")

    return s, stats

def parse_tips(nwk_text: str) -> List[str]:
    return _TOK_TIP.findall(nwk_text)

def strip_mark(name: str) -> str:
    """忽略 #1 打标：用于集合比较"""
    return name.replace("#1","")

def tips_equal_fasta_tree(fa_species: List[str], nwk_text: str) -> Tuple[bool, List[str], List[str]]:
    tree_tips = [strip_mark(x) for x in parse_tips(nwk_text)]
    fa_set, tr_set = set(fa_species), set(tree_tips)
    only_fa  = sorted(fa_set - tr_set)
    only_tr  = sorted(tr_set - fa_set)
    return (len(only_fa)==0 and len(only_tr)==0), only_fa, only_tr

# ========== 单个 OG×FG 任务 ==========
def prepare_and_run_one(og: str, fgname: str, codeml_bin: str,
                        alt_tpl: Path, null_tpl: Path,
                        log: logging.Logger) -> str:
    # 源文件
    seq_src  = need_file(SEQ_DIR / f"{og}.codon.fna", f"{og} codon MSA")
    tree_src = need_file(SETS_DIR / f"{fgname}.nwk", f"{fgname} 打标树")

    # 读取并“消毒”树
    raw_tree = tree_src.read_text(encoding="utf-8")
    clean_tree, stats = sanitize_foreground_tree(raw_tree)
    log.info(f"TREE_SANITIZE[{og}/{fgname}] removed_brackets={stats['removed_brackets']} "
             f"removed_internal_labels={stats['removed_internal_labels']} length={len(clean_tree)}")

    # FASTA 物种名（>Species）
    fa_species = []
    with seq_src.open("r", encoding="utf-8", errors="ignore") as f:
        name=None
        for line in f:
            if line.startswith(">"):
                name=line[1:].strip()
                fa_species.append(name)

    # 校验集合（忽略 #1）
    ok_set, only_fa, only_tr = tips_equal_fasta_tree(fa_species, clean_tree)
    if not ok_set:
        raise RuntimeError(f"[ERR] 物种集合不一致：FASTA={len(set(fa_species))} 与 TREE={len(set([strip_mark(x) for x in parse_tips(clean_tree)]))}\n"
                           f"FA_ONLY={only_fa}\nTREE_ONLY={only_tr}")

    # 目标目录
    work_root = ensure_dir(RAW_ROOT / og / fgname)
    alt_d = ensure_dir(work_root / "alt")
    nul_d = ensure_dir(work_root / "null")

    # 写入输入文件（拷贝 & 消毒树）
    for d in (alt_d, nul_d):
        shutil.copyfile(seq_src,  d / "seq.codon.fna")
        (d / "species_tree.nwk").write_text(clean_tree, encoding="utf-8")

    # 写 ctl（模板变量只用相对路径，增强可移植性）
    rep = {"seqfile": "seq.codon.fna", "treefile": "species_tree.nwk", "outfile": "mlc.txt"}
    (alt_d / "alt.ctl").write_text(render_tpl(alt_tpl, rep),  encoding="utf-8")
    (nul_d / "null.ctl").write_text(render_tpl(null_tpl, rep), encoding="utf-8")

    # 断点续跑逻辑
    alt_mlc, nul_mlc = alt_d / "mlc.txt", nul_d / "mlc.txt"
    if alt_mlc.is_file() and nul_mlc.is_file() and not FORCE_RERUN:
        if has_lnl(alt_mlc) and has_lnl(nul_mlc):
            return f"SKIP {og}/{fgname}"

    # 若 ALT 已有 lnL、NULL 无 → 只跑 NULL
    # 否则顺序跑 ALT → NULL
    soft_flags = []

    if not (alt_mlc.is_file() and has_lnl(alt_mlc)) or FORCE_RERUN:
        rc_alt = run_cmd_capture([codeml_bin, "alt.ctl"], alt_d, log)
        if rc_alt != 0 and has_lnl(alt_mlc):
            log.info(f"SOFT_OK[{og}/{fgname}/alt] rc={rc_alt} but mlc.txt has lnL()")
            soft_flags.append("alt")
        elif rc_alt != 0 and not has_lnl(alt_mlc):
            raise RuntimeError(f"[ERR] codeml 退出码 {rc_alt} @ {alt_d}")

    # NULL
    rc_null = run_cmd_capture([codeml_bin, "null.ctl"], nul_d, log)
    if rc_null != 0 and has_lnl(nul_mlc):
        log.info(f"SOFT_OK[{og}/{fgname}/null] rc={rc_null} but mlc.txt has lnL()")
        soft_flags.append("null")
    elif rc_null != 0 and not has_lnl(nul_mlc):
        raise RuntimeError(f"[ERR] codeml 退出码 {rc_null} @ {nul_d}")

    if soft_flags:
        return f"SOFT {og}/{fgname} {'+'.join(soft_flags)}"
    return f"DONE {og}/{fgname}"

# ========== 主流程 ==========
def main():
    cfg = load_config()
    log = get_logger(LOG_FILE)

    log.info("=======================")
    log.info(" APhylo 06 — 批量 codeml")
    log.info("=======================")

    # 模板路径
    alt_tpl  = need_file(Path(cfg["codeml"]["alt_template"]),  "ALT 模板")
    null_tpl = need_file(Path(cfg["codeml"]["null_template"]), "NULL 模板")
    keep     = need_file(QC_DIR / "keep_og.list", "保留 OG 列表")

    ogs = [x.strip() for x in keep.read_text(encoding="utf-8").splitlines() if x.strip()]
    fgs = sorted([p.stem for p in SETS_DIR.glob("*.nwk")])
    if not fgs:
        raise FileNotFoundError("[ERR] 未发现前景树：results/04_codeml/sets/*.nwk，请先运行 05 步。")

    codeml_bin = cfg.get("binaries", {}).get("codeml", "codeml")
    threads_cfg = int(cfg.get("codeml", {}).get("threads", 1))
    cpu = max(1, (os.cpu_count() or 1))
    threads = max(1, min(threads_cfg, cpu))

    log.info(f"计划运行：OG={len(ogs)} × FG={len(fgs)} × 2 (ALT/NULL)")
    log.info(f"并行任务数（config.yaml: codeml.threads）= {threads}")

    ensure_dir(RAW_ROOT)
    tasks: List[Tuple[str, str]] = [(og, fg) for og in ogs for fg in fgs]

    ok = soft_ok = skip = fail = 0
    if threads == 1:
        for og, fg in tasks:
            try:
                msg = prepare_and_run_one(og, fg, codeml_bin, Path(alt_tpl), Path(null_tpl), log)
                if   msg.startswith("SKIP"):  skip   += 1
                elif msg.startswith("SOFT"):  soft_ok+= 1
                else:                         ok     += 1
            except Exception as e:
                log.error(f"[FAIL] {og}/{fg} -> {e}"); fail += 1
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            fut2key = {ex.submit(prepare_and_run_one, og, fg, codeml_bin, Path(alt_tpl), Path(null_tpl), log):(og, fg)
                       for og, fg in tasks}
            for fut in as_completed(fut2key):
                og, fg = fut2key[fut]
                try:
                    msg = fut.result()
                    if   msg.startswith("SKIP"):  skip   += 1
                    elif msg.startswith("SOFT"):  soft_ok+= 1
                    else:                         ok     += 1
                except Exception as e:
                    log.error(f"[FAIL] {og}/{fg} -> {e}"); fail += 1

    log.info(f"[SUMMARY] ok={ok}  soft_ok={soft_ok}  skip={skip}  fail={fail}")
    if fail == 0:
        Path("results/04_codeml/.codeml.done").touch()
        log.info("[DONE] codeml 批量运行完成；写入哨兵 results/04_codeml/.codeml.done")
    else:
        sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n"); sys.exit(2)

