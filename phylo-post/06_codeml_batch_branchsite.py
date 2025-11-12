#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
06_codeml_batch_branchsite.py —— 两阶段批跑：先全量 ALT，再全量 NULL
设计要点（与皇上要求一致）：
  1) 不做 OG 级剪枝，不做有效性判定，不做软通过；不中断不早停。
  2) 每个 {OG, FG, alt|null} 独立目录，调用 codeml 前只清理产物文件，避免 "but mlc.txt has lnL()".
  3) 模板仅复制，不改写；模板内固定：seqfile=seq.codon.fna, treefile=species_tree.nwk, outfile=mlc.txt。
  4) Phase-A: 跑完所有 OG 的 ALT；Phase-B: 再跑所有 OG 的 NULL。NULL 与 ALT 无依赖门槛。
  5) 并发度与 codeml 路径从 config.yaml 读取；其余路径与 aphylo 既有契约保持一致。
"""

from __future__ import annotations
import os, sys, io, shutil, subprocess, logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Tuple

# ================= 固定路径（对齐 aphylo 约定） =================
PROJECT_ROOT = Path(".").resolve()
CONFIG_PATH  = PROJECT_ROOT / "config.yaml"

SEQ_DIR      = PROJECT_ROOT / "results/03_codon/codon_msa"      # 密码子对齐根（固定）
SETS_DIR     = PROJECT_ROOT / "results/04_codeml/sets"          # 前景树 .nwk（固定）
KEEP_OG_LIST = PROJECT_ROOT / "results/03_qc/keep_og.list"      # OG 清单（固定）

RAW_ROOT     = PROJECT_ROOT / "results/04_codeml/raw"           # 产物根：raw/OG/FG/{alt,null}
LOG_ROOT     = PROJECT_ROOT / "logs/06_branchsite"              # 日志根：logs/06_branchsite/OG/FG/{alt,null}.log
TEMPLATE_ALT = PROJECT_ROOT / "templates/branch_site_alt.ctl"   # 模板（仅复制）
TEMPLATE_NULL= PROJECT_ROOT / "templates/branch_site_null.ctl"

# ================= 顶部参数（可被 config.yaml 覆盖） =================
CODEML_BIN   = "codeml"            # 可在 config.yaml.binaries.codeml 覆盖
THREADS      = 8                   # 可在 config.yaml.codeml.threads 覆盖；表示并发 OG 数
ECHO_STREAM  = True                # 是否同时将子进程输出打印到屏幕
CLEAN_PATTERNS = [                 # 仅清理产物，避免 "mlc 已存在" 导致 codeml 退出
    "mlc.txt", "rst", "rst1", "rst2", "rub", "lnf",
    "2ML.t", "2ML.dN", "2ML.dS", "2ML.out", "2ML.omega",
    "2ML.trees", "2ML.siteclass", "2ML.siteclassN",
    "*.tmp", "tmp*"
]

# ================= 辅助函数 =================
def load_config(path: Path) -> Dict:
    """读取 config.yaml；失败时返回 {}，仅用于覆盖 codeml 路径与并发度。"""
    try:
        import yaml
    except Exception:
        return {}
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}
    except Exception:
        return {}
    return cfg if isinstance(cfg, dict) else {}

def apply_overrides(cfg: Dict):
    """从 config.yaml 读取 codeml 可执行路径与并发度；其他路径不覆盖。"""
    global CODEML_BIN, THREADS
    b = (cfg.get("binaries") or {})
    if isinstance(b, dict) and b.get("codeml"):
        CODEML_BIN = str(b["codeml"])
    c = (cfg.get("codeml") or {})
    if isinstance(c, dict) and c.get("threads") is not None:
        try:
            THREADS = int(c["threads"])
        except Exception:
            pass
    # 兜底限制
    cpu = max(1, (os.cpu_count() or 1))
    THREADS = max(1, min(THREADS, cpu))

def need_file(p: Path, what: str) -> Path:
    """校验必须存在的文件。"""
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def get_logger() -> logging.Logger:
    """同时写日志文件与控制台，方便 tail -f；立即刷新。"""
    ensure_dir(LOG_ROOT)
    lg = logging.getLogger("aphylo.06")
    lg.handlers.clear()
    lg.setLevel(logging.INFO)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s",
                            "%Y-%m-%d %H:%M:%S")
    # 主日志（总览）
    main_log = LOG_ROOT / "06_codeml_batch_branchsite.log"
    fh = logging.FileHandler(main_log, encoding="utf-8")
    fh.setFormatter(fmt); fh.setLevel(logging.INFO)
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt); sh.setLevel(logging.INFO)
    lg.addHandler(fh); lg.addHandler(sh)
    # 立即刷新 stdout/stderr，便于实时观察
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout = _Flush(sys.stdout); sys.stderr = _Flush(sys.stderr)
    return lg

def clean_outputs(d: Path):
    """仅清理产物，不动 seq/tree/ctl。"""
    for pat in CLEAN_PATTERNS:
        for p in d.glob(pat):
            try:
                if p.is_file():
                    p.unlink(missing_ok=True)
                else:
                    shutil.rmtree(p, ignore_errors=True)
            except Exception:
                pass

def copy_inputs_and_templates(og: str, fg: str, dest: Path):
    """为单元目录就位输入与模板（仅复制，不改写）。"""
    ensure_dir(dest)
    # 1) 复制密码子对齐到固定文件名
    seq_src = need_file(SEQ_DIR / f"{og}.codon.fna", f"{og} codon MSA")
    shutil.copy2(seq_src, dest / "seq.codon.fna")
    # 2) 复制前景树到固定文件名（不做消毒/剪枝）
    tree_src = need_file(SETS_DIR / f"{fg}.nwk", f"{fg} 前景树")
    shutil.copy2(tree_src, dest / "species_tree.nwk")
    # 3) 复制 ctl 模板（仅复制，不改）
    if dest.name == "alt":
        tpl = need_file(TEMPLATE_ALT, "ALT 模板")
        shutil.copy2(tpl, dest / "branch_site_alt.ctl")
    else:
        tpl = need_file(TEMPLATE_NULL, "NULL 模板")
        shutil.copy2(tpl, dest / "branch_site_null.ctl")

def run_streaming(cmd: List[str], cwd: Path, logfile: Path, tag: str) -> int:
    """流式运行 codeml，既打印也写入独立日志文件，返回退出码。"""
    ensure_dir(logfile.parent)
    env = os.environ.copy()
    # 强制单线程，避免并发爆栈/数值不稳
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    with open(logfile, "a", encoding="utf-8") as lf:
        lf.write(f"[CMD] (cwd={cwd}) $ {' '.join(cmd)}\n"); lf.flush()
        p = subprocess.Popen(
            cmd, cwd=str(cwd), env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, universal_newlines=True, bufsize=1
        )
        assert p.stdout is not None
        for line in p.stdout:
            line = line.rstrip("\n")
            if ECHO_STREAM:
                print(f"[{tag}] {line}")
            lf.write(line + "\n"); lf.flush()
        p.wait()
        lf.write(f"[EXIT] code={p.returncode}\n"); lf.flush()
        return int(p.returncode)

def run_one(og: str, fg: str, model: str) -> Tuple[str, str, str, int]:
    """
    运行单个 {OG, FG, model} 单元。
    model: 'alt' | 'null'
    返回: (og, fg, model, rc)
    """
    unit_root = RAW_ROOT / og / fg
    dest = unit_root / model
    # 日志路径
    log_path = LOG_ROOT / og / fg / f"{model}.log"
    # 创建并清理产物
    ensure_dir(unit_root); ensure_dir(dest)
    clean_outputs(dest)
    # 就位输入与模板
    copy_inputs_and_templates(og, fg, dest)
    # 运行 codeml
    ctl_name = "branch_site_alt.ctl" if model == "alt" else "branch_site_null.ctl"
    rc = run_streaming([CODEML_BIN, ctl_name], dest, log_path, f"{og}|{fg}|{model.upper()}")
    return og, fg, model, rc

def list_fgs() -> List[str]:
    """列出所有前景集合（来自 sets 目录的 .nwk 文件名）。"""
    if not SETS_DIR.exists():
        return []
    fgs = []
    for p in SETS_DIR.glob("*.nwk"):
        if p.is_file():
            fgs.append(p.stem)
    return sorted(set(fgs))

def read_ogs() -> List[str]:
    """读取 OG 清单（严格使用 keep_og.list）。"""
    need_file(KEEP_OG_LIST, "保留 OG 列表")
    out: List[str] = []
    with open(KEEP_OG_LIST, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s)
    return sorted(set(out))

# ================= 主流程 =================
def phase_run(ogs: List[str], fg: str, model: str, threads: int, lg: logging.Logger) -> Dict[str, int]:
    """
    以线程池并发跑一个前景 FG 的单个模型（alt 或 null）。
    不中断：全部提交并等待；最终返回计数汇总。
    """
    ok = fail = 0
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        fut2key = {ex.submit(run_one, og, fg, model): (og, fg) for og in ogs}
        for fut in as_completed(fut2key):
            og, _fg = fut2key[fut]
            try:
                _og, _fg2, _md, rc = fut.result()
                if rc == 0:
                    ok += 1
                    lg.info(f"[OK]    {og}|{fg}|{model.upper()} rc=0")
                else:
                    fail += 1
                    lg.info(f"[FAIL]  {og}|{fg}|{model.upper()} rc={rc}")
            except Exception as e:
                fail += 1
                lg.info(f"[FAIL]  {og}|{fg}|{model.upper()} EXC={e.__class__.__name__}: {e}")
    return {"ok": ok, "fail": fail, "total": ok + fail}

def main():
    lg = get_logger()
    # 配置覆盖
    cfg = load_config(CONFIG_PATH); apply_overrides(cfg)
    # 基础检查
    need_file(TEMPLATE_ALT,  "ALT 模板")
    need_file(TEMPLATE_NULL, "NULL 模板")
    if shutil.which(CODEML_BIN) is None:
        lg.error(f"[ERR] 未找到可执行文件：{CODEML_BIN}"); sys.exit(2)
    if not SEQ_DIR.exists():
        lg.error(f"[ERR] 密码子对齐根不存在：{SEQ_DIR}"); sys.exit(2)
    fgs = list_fgs()
    if not fgs:
        lg.error(f"[ERR] 未发现前景树：{SETS_DIR}/*.nwk；请先运行 05 步。"); sys.exit(2)
    ogs = read_ogs()
    if not ogs:
        lg.error(f"[ERR] keep_og.list 为空：{KEEP_OG_LIST}"); sys.exit(2)

    lg.info("=======================")
    lg.info(" APhylo 06 — 两阶段批跑（先 ALT 再 NULL）")
    lg.info("=======================")
    lg.info(f"[INIT] OG={len(ogs)}  FG={len(fgs)}  threads={THREADS}")
    lg.info(f"[PATH] SEQ_DIR={SEQ_DIR}")
    lg.info(f"[PATH] SETS_DIR={SETS_DIR}")
    lg.info(f"[PATH] RAW_ROOT={RAW_ROOT}")
    lg.info(f"[BIN ] codeml={CODEML_BIN}")

    grand = {"alt": {"ok":0,"fail":0,"total":0}, "null": {"ok":0,"fail":0,"total":0}}

    # 逐前景执行：Phase-A ALT → Phase-B NULL
    for fg in fgs:
        lg.info(f"\n[FG] >>> {fg}  (Phase-A: ALT)")
        r_alt = phase_run(ogs, fg, "alt", THREADS, lg)
        grand["alt"]["ok"]   += r_alt["ok"]
        grand["alt"]["fail"] += r_alt["fail"]
        grand["alt"]["total"]+= r_alt["total"]

        lg.info(f"\n[FG] >>> {fg}  (Phase-B: NULL)")
        r_null = phase_run(ogs, fg, "null", THREADS, lg)
        grand["null"]["ok"]   += r_null["ok"]
        grand["null"]["fail"] += r_null["fail"]
        grand["null"]["total"]+= r_null["total"]

    lg.info("\n[SUMMARY]")
    lg.info(f"  ALT:  ok={grand['alt']['ok']:4d}  fail={grand['alt']['fail']:4d}  total={grand['alt']['total']:4d}")
    lg.info(f"  NULL: ok={grand['null']['ok']:4d}  fail={grand['null']['fail']:4d}  total={grand['null']['total']:4d}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断"); sys.exit(130)