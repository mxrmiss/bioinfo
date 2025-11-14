#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— 运行 CAFE5（按 config.cafe5.models）

本版要点（在不改 11 的前提下）：
  - 输出扁平化：加 `-o .`，产物直接落在 results/06_cafe/models/<model>/ 下；
  - 自动自救：遇到 “Failed to initialize any reasonable values”/建议语，解析
    “Families with largest size differentials:” 区块，剔除这些 OG，生成
    family.autofixN.tsv 并重跑，轮数由 config.cafe5.max_autofix_rounds 控制；
  - 成功判定：仅当返回码=0 且日志未出现上述失败提示，才视为成功；
  - 失败处理：自动修正用尽仍失败则报错退出（不写 .cafe.done、不打印 DONE）；
  - 新增：
      * 每模型写入 models/<model>/run.log（与全局日志并存）；
      * 非最终轮生成的结果文件移动到 models/<model>/_roundN_tmp/，避免 13 聚合旧轮；
      * 被剔除的 OG 记录在 models/<model>/autofix_removed_roundN.tsv。
"""

from __future__ import annotations
import sys, io, logging, subprocess, re, shutil
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable
import yaml

DEFAULT_CONFIG = "config.yaml"

# ====================== 通用工具 ======================

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

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def need_dir(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_dir(): raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p

def need_file(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_file(): raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """文件日志 + 屏幕实时输出；stdout/stderr 立刻刷新。"""
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
    path.parent.mkdir(parents=True, exist_ok=True); Path(path).touch()

# ====================== CAFE5 相关工具 ======================

FAILED_INIT_PAT = re.compile(r"Failed to initialize any reasonable values", re.I)
ADVICE_PAT = re.compile(r"You may want to try removing the top few families with the largest difference", re.I)
BLOCK_HEAD_PAT = re.compile(r"Families with largest size differentials\s*:", re.I)
OG_LINE_PAT = re.compile(r"^\s*(OG\d+)\s*:\s*(\d+)\s*$", re.I)

def detect_cafe_errors_fatal(stdout: str, model: str, log_file: Path) -> None:
    """
    文本级“致命错误”检测：一旦发现直接抛异常（不交给自动修正流程）。
    不拦截初始化失败（交由外层处理）。
    """
    text = stdout or ""
    if re.search(r"was not found in gene family", text, flags=re.I):
        raise RuntimeError(f"[ERR] 输出含 'was not found in gene family' —— 树/矩阵物种不一致；模型：{model}")
    if re.search(r"unrecognized option", text, flags=re.I) or re.search(r"Unrecognized parameter", text, flags=re.I):
        raise RuntimeError(f"[ERR] 输出含 'unrecognized option/parameter' —— 版本不支持该参数；模型：{model}。详见 {log_file}")
    if re.search(r"Failed to open", text, flags=re.I):
        raise RuntimeError(f"[ERR] 输出含 'Failed to open' —— family/tree 路径无效；模型：{model}。详见 {log_file}")

def need_autofix(stdout: str) -> bool:
    """出现初始化失败或官方建议语 → 需要自动修正。"""
    return bool(stdout and (FAILED_INIT_PAT.search(stdout) or ADVICE_PAT.search(stdout)))

def parse_extreme_ogs_from_stdout(stdout: str) -> List[str]:
    """
    回溯解析日志中的“极端家族”区块：
      - 找到最后一个 'Families with largest size differentials:'；
      - 向下收集连续的 'OGxxxxx: NNN' 行。
    """
    if not stdout: return []
    lines = stdout.splitlines()
    start_idx = -1
    for i in range(len(lines)-1, -1, -1):
        if BLOCK_HEAD_PAT.search(lines[i]):
            start_idx = i; break
    if start_idx < 0: return []
    ogs: List[str] = []
    for j in range(start_idx+1, len(lines)):
        m = OG_LINE_PAT.match(lines[j])
        if not m: break
        ogs.append(m.group(1))
    # 去重保序
    seen, uniq = set(), []
    for og in ogs:
        if og not in seen: seen.add(og); uniq.append(og)
    return uniq

def filter_family_drop_ogs(in_family: Path, ogs: List[str], out_family: Path, log: logging.Logger):
    """从 family.tsv 中剔除指定 OG，写出新矩阵。"""
    ogset = set(ogs)
    n_in = n_out = n_drop = 0
    with in_family.open("r", encoding="utf-8") as fin, out_family.open("w", encoding="utf-8") as fout:
        header = fin.readline(); 
        if not header: raise RuntimeError(f"[ERR] family.tsv 为空：{in_family}")
        fout.write(header)
        for line in fin:
            if not line.strip(): continue
            n_in += 1
            og = line.split("\t", 1)[0].strip()
            if og in ogset:
                n_drop += 1
                continue
            fout.write(line); n_out += 1
    log.info(f"[FILTER] {in_family.name}：数据行 {n_in}，剔除 {n_drop}，保留 {n_out} → {out_family.name}")

def list_dirnames(mdir: Path) -> set:
    """列出目录内当前文件名集合（用于比较前后差异）。"""
    return set(p.name for p in mdir.iterdir())

def move_round_outputs(mdir: Path, new_names: set, round_idx: int, keep_names: set, log: logging.Logger):
    """
    将非最终轮新产生的结果文件移入临时归档子目录，避免被 13 聚合。
    - new_names: 本轮新出现的文件名集合
    - keep_names: 需要保留在根目录的文件（如 run.log / family.autofixN.tsv）
    """
    dst = ensure_dir(mdir / f"_round{round_idx}_tmp")
    moved = 0
    for name in sorted(new_names):
        if name in keep_names: continue
        src = mdir / name
        try:
            src.rename(dst / name); moved += 1
        except Exception:
            pass
    log.info(f"[CLEAN] 第 {round_idx} 轮产生的 {moved} 个结果文件已移入：{dst}")

def append_runlog(runlog: Path, header: str, lines: List[str]):
    """把本轮 stdout 追加写入 models/<model>/run.log（带轮次分隔）。"""
    with runlog.open("a", encoding="utf-8") as w:
        w.write(f"\n===== {header} =====\n")
        for ln in lines:
            w.write(ln + ("\n" if not ln.endswith("\n") else ""))

def run_cafe_streaming(cmd: List[str], mdir: Path,
                       log: logging.Logger, model: str) -> List[str]:
    """
    以流式方式运行 CAFE5，返回 stdout 行列表（不含换行符）。
    同时把每行写入全局日志。
    """
    log.info(f"[MODEL {model}] CMD: " + " ".join(map(str, cmd)))
    proc = subprocess.Popen(
        cmd, cwd=str(mdir),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1, universal_newlines=True,
    )
    lines: List[str] = []
    log.info(f"[MODEL {model}] ---- CAFE5 输出开始 ----")
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip("\n")
        lines.append(line)
        log.info(f"[MODEL {model}] {line}")
    proc.wait()
    log.info(f"[MODEL {model}] ---- CAFE5 输出结束 ----")
    if proc.returncode != 0:
        raise RuntimeError(f"[ERR] CAFE5 退出码 {proc.returncode}（模型：{model}）")
    return lines

# ====================== 主流程 ======================

def main():
    cfg = load_config()
    paths = cfg["paths"]; cafe = cfg.get("cafe5", {}) or {}

    logs_dir = Path(paths["logs_dir"])
    LOG_FILE = logs_dir / "12_cafe5_run_models.log"
    log = get_logger("aphylo.12", LOG_FILE)
    banner(log, "APhylo 12 — CAFE5 批量模型（自救增强版 v2）")

    if not cafe.get("enable", True):
        log.info("cafe5.enable=false —— 跳过本步")
        write_done(Path(paths["cafe_run_dir"]) / ".cafe.done"); return

    threads = int(cafe.get("threads", 4))
    k = cafe.get("gamma_k", None)
    pval_thr = cafe.get("pvalue", None)
    models = cafe.get("models", ["global"])
    max_autofix_rounds = int(cafe.get("max_autofix_rounds", 0))  # 0=关闭自动修正

    cafe_bin = cfg.get("binaries", {}).get("cafe5", "cafe5")

    # 输入
    inp = need_dir(Path(paths["cafe_run_dir"]) / "input", "CAFE 输入目录")
    ts = sorted(inp.glob("*.tsv"))
    if not ts: raise FileNotFoundError("[ERR] 未找到 family TSV（*.tsv）")
    family0 = ts[0]
    nwks = sorted(inp.glob("*.nwk"))
    if not nwks: raise FileNotFoundError("[ERR] 未找到超时钟树（*.nwk）")
    nwk = nwks[0]

    family0_abs = family0.resolve(); nwk_abs = nwk.resolve()
    log.info(f"[IN] family: {family0_abs}")
    log.info(f"[IN] tree  : {nwk_abs}")
    if max_autofix_rounds > 0:
        log.info(f"[CFG] max_autofix_rounds = {max_autofix_rounds}")

    # 输出
    models_dir = ensure_dir(Path(paths["cafe_run_dir"]) / "models")

    all_ok = True
    for model in models:
        mname = str(model).replace(":", "_")
        mdir = ensure_dir(models_dir / mname)
        runlog = mdir / "run.log"  # 每模型独立 run.log
        log.info(f"[MODEL {model}] 输出目录：{mdir}")

        current_family = family0_abs
        success = False

        rounds = max(0, max_autofix_rounds)
        for r in range(0, rounds + 1):
            before = list_dirnames(mdir)

            cmd = [cafe_bin, "-i", str(current_family), "-t", str(nwk_abs), "-c", str(threads), "-o", "."]
            if k: cmd += ["-k", str(k)]
            if pval_thr is not None: cmd += ["-P", str(pval_thr)]

            # 跑一轮
            out_lines = run_cafe_streaming(cmd, mdir, log, model)
            stdout_text = "\n".join(out_lines)

            # 模型级 run.log 也记录（附上轮次抬头）
            append_runlog(runlog, f"MODEL={model} ROUND={r+1}", out_lines)

            # 致命错误检查（不含初始化失败）
            detect_cafe_errors_fatal(stdout_text, model, LOG_FILE)

            # 成功判定 / 自动修正
            if need_autofix(stdout_text):
                # 解析极端 OG
                ogs = parse_extreme_ogs_from_stdout(stdout_text)
                if ogs and r < rounds:
                    # 记录被删 OG
                    removed_tsv = mdir / f"autofix_removed_round{r+1}.tsv"
                    with removed_tsv.open("w", encoding="utf-8") as w:
                        w.write("OG\n"); w.write("\n".join(ogs) + "\n")

                    # 生成下一轮 family
                    next_family = mdir / f"family.autofix{r+1}.tsv"
                    filter_family_drop_ogs(current_family, ogs, next_family, log)
                    current_family = next_family.resolve()

                    # 把本轮新产物移走，避免 13 聚合旧轮
                    after = list_dirnames(mdir)
                    new_names = after - before
                    keep = {runlog.name, removed_tsv.name, next_family.name}
                    move_round_outputs(mdir, new_names, r+1, keep, log)

                    log.info(f"[MODEL {model}] 第 {r+1} 轮自动修正完成，进入下一轮。")
                    continue
                else:
                    log.error(f"[MODEL {model}] 检测到初始化失败，但已无剩余轮次或未解析出 OG。")
                    success = False
                    break
            else:
                success = True
                log.info(f"[MODEL {model}] CAFE5 运行完成（第 {r+1} 轮，无失败提示）。")
                break

        if not success:
            all_ok = False
            break  # 按你的要求：任一模型失败即整体退出

    if all_ok:
        write_done(Path(paths["cafe_run_dir"]) / ".cafe.done")
        log.info("[DONE] CAFE5 模型运行完成")
    else:
        sys.stderr.write("[ERR] CAFE5 自动修正用尽仍失败；请检查 family/树/注释质量。\n")
        sys.exit(2)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)

