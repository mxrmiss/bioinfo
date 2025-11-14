#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— 运行 CAFE5（按 config.cafe5.models）

本版要点（在不改 11 的前提下）：
  - 输出扁平化：加 `-o .`，产物直接落在 results/06_cafe/models/<model>/ 下；
  - 自动自救：遇到初始化失败/官方建议语，解析 “Families with largest size differentials:” 区块，
    剔除这些 OG，生成 family.autofixN.tsv 并重跑；轮数由 config.cafe5.max_autofix_rounds 控制；
  - 修复：按“任意空白分隔”解析 family 第一列，并做交叉校验；若需删除但未命中任何 OG，立即报错退出；
  - 开跑前“卫生”：
      * 归档旧产物到 _prev_run_YYYYmmdd-HHMMSS/；
      * run.log 轮转到 run.log.YYYYmmdd-HHMMSS，新 run.log 仅记录本批次；
  - 成功判定：仅当返回码=0 且日志未出现初始化失败提示，才视为成功；
  - 失败处理：自动修正用尽仍失败则报错退出（不写 .cafe.done、不打印 DONE）。

保持既有特性：
  - 全局日志：logs/12_cafe5_run_models.log（实时刷屏）；
  - 每模型 run.log（追加写本批次；老日志已轮转）；
  - 旧轮产物移入 _roundN_tmp/，根目录只留最终轮。
"""

from __future__ import annotations
import sys, io, logging, subprocess, re, shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Tuple
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

def banner(log: logging.Logger, text: str):
    bar = "=" * max(10, len(text)+2)
    log.info(bar); log.info(f" {text} "); log.info(bar)

def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True); Path(path).touch()

# ====================== CAFE5 相关工具 ======================

FAILED_INIT_PAT = re.compile(r"Failed to initialize any reasonable values", re.I)
ADVICE_PAT = re.compile(r"You may want to try removing the top few families with the largest difference", re.I)
BLOCK_HEAD_PAT = re.compile(r"Families with largest size differentials\s*:", re.I)
OG_LINE_PAT = re.compile(r"^\s*(OG\d+)\s*:\s*(\d+)\s*$", re.I)

def detect_cafe_errors_fatal(stdout: str, model: str, log_file: Path) -> None:
    """致命类“伪成功”错误：一旦发现直接抛异常（初始化失败交由外层处理）。"""
    text = stdout or ""
    if re.search(r"was not found in gene family", text, flags=re.I):
        raise RuntimeError(f"[ERR] 输出含 'was not found in gene family' —— 树/矩阵物种不一致；模型：{model}")
    if re.search(r"unrecognized option", text, flags=re.I) or re.search(r"Unrecognized parameter", text, flags=re.I):
        raise RuntimeError(f"[ERR] 输出含 'unrecognized option/parameter' —— 版本不支持该参数；模型：{model}（详见 {log_file}）")
    if re.search(r"Failed to open", text, flags=re.I):
        raise RuntimeError(f"[ERR] 输出含 'Failed to open' —— family/tree 路径无效；模型：{model}（详见 {log_file}）")

def need_autofix(stdout: str) -> bool:
    """出现初始化失败或官方建议语 → 需要自动修正。"""
    return bool(stdout and (FAILED_INIT_PAT.search(stdout) or ADVICE_PAT.search(stdout)))

def parse_extreme_ogs_from_stdout(stdout: str) -> List[str]:
    """回溯解析‘largest differentials’区块，收集连续的 OGxxxx: NNN 行。"""
    if not stdout: return []
    lines = stdout.splitlines()
    start_idx = -1
    for i in range(len(lines)-1, -1, -1):
        if BLOCK_HEAD_PAT.search(lines[i]): start_idx = i; break
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

def first_token(s: str) -> str:
    """取首个非空白 token（适配制表符或空格混用）。"""
    return re.split(r"\s+", s.strip(), maxsplit=1)[0] if s.strip() else ""

def count_hits_in_family(family_path: Path, ogs: List[str]) -> int:
    """统计 OG 列表在 family 中的命中数（按首 token 匹配）。"""
    if not ogs: return 0
    targets = set(ogs)
    hits = 0
    with family_path.open("r", encoding="utf-8", errors="ignore") as fin:
        _ = fin.readline()  # header
        for line in fin:
            if not line.strip(): continue
            og = first_token(line)
            if og in targets: hits += 1
    return hits

def filter_family_drop_ogs(in_family: Path, ogs: List[str], out_family: Path, log: logging.Logger):
    """
    从 family 文件剔除指定 OG 行，写出新矩阵。
    - 识别列：按“任意空白分隔”取首 token；
    - 写回：保留原始行文本，不改变文件分隔/格式。
    """
    ogset = set(ogs)
    n_in = n_out = n_drop = 0
    with in_family.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_family.open("w", encoding="utf-8") as fout:
        header = fin.readline()
        if not header: raise RuntimeError(f"[ERR] family.tsv 为空：{in_family}")
        fout.write(header)
        for line in fin:
            if not line.strip(): continue
            n_in += 1
            og = first_token(line)
            if og in ogset:
                n_drop += 1
                continue
            fout.write(line)
            n_out += 1
    log.info(f"[FILTER] {in_family.name}：数据行 {n_in}，剔除 {n_drop}，保留 {n_out} → {out_family.name}")

def list_dirnames(mdir: Path) -> set:
    """列出现存文件/目录名集合（用于比较新产物）。"""
    return set(p.name for p in mdir.iterdir())

def move_round_outputs(mdir: Path, new_names: set, round_idx: int, keep_names: set, log: logging.Logger):
    """
    将非最终轮新产生的文件移入临时归档子目录（避免被 13 聚合）。
    keep_names：本轮需要保留在根目录的文件名（如 run.log / family.autofixN.tsv / removed.tsv）。
    """
    dst = ensure_dir(mdir / f"_round{round_idx}_tmp")
    moved = 0
    for name in sorted(new_names):
        if name in keep_names: continue
        src = mdir / name
        try:
            shutil.move(str(src), str(dst / name)); moved += 1
        except Exception:
            pass
    log.info(f"[CLEAN] 第 {round_idx} 轮产生的 {moved} 个结果文件已移入：{dst}")

def append_runlog(runlog: Path, header: str, lines: List[str]):
    """把本轮 stdout 追加到 models/<model>/run.log（带轮次分隔头）。"""
    with runlog.open("a", encoding="utf-8") as w:
        w.write(f"\n===== {header} =====\n")
        for ln in lines:
            w.write(ln + ("\n" if not ln.endswith("\n") else ""))

def run_cafe_streaming(cmd: List[str], mdir: Path,
                       log: logging.Logger, model: str) -> List[str]:
    """流式运行 CAFE5，返回 stdout 行列表（同时写全局日志）。"""
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

def archive_previous_run(mdir: Path, log: logging.Logger):
    """
    模型开跑前：把模型目录现有内容整体归档到 _prev_run_YYYYmmdd-HHMMSS/，避免历史残留干扰。
    run.log 的轮转另行处理。
    """
    items = list(mdir.iterdir())
    if not items: 
        return
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dst = ensure_dir(mdir / f"_prev_run_{stamp}")
    moved = 0
    for p in items:
        # run.log 类在轮转函数里处理，这里不移动
        if p.name.startswith("run.log"):
            continue
        try:
            shutil.move(str(p), str(dst / p.name)); moved += 1
        except Exception:
            pass
    if moved:
        log.info(f"[ARCHIVE] 发现历史产物，已归档 {moved} 项到：{dst}")

def rotate_runlog(runlog: Path, log: logging.Logger):
    """若存在旧 run.log，则轮转为 run.log.YYYYmmdd-HHMMSS，并创建空的新 run.log。"""
    if runlog.exists():
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        bak = runlog.with_name(f"run.log.{stamp}")
        try:
            runlog.rename(bak)
            log.info(f"[LOG] 轮转 run.log → {bak.name}")
        except Exception:
            # 若 rename 失败则复制一份
            shutil.copy2(str(runlog), str(bak))
            runlog.unlink(missing_ok=True)
            log.info(f"[LOG] 复制轮转 run.log → {bak.name}")
    # 新建空 run.log
    runlog.touch()

# ====================== 主流程 ======================

def main():
    cfg = load_config()
    paths = cfg["paths"]; cafe = cfg.get("cafe5", {}) or {}

    logs_dir = Path(paths["logs_dir"])
    LOG_FILE = logs_dir / "12_cafe5_run_models.log"
    log = get_logger("aphylo.12", LOG_FILE)
    banner(log, "APhylo 12 — CAFE5 批量模型（自救增强版 v3）")

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
        runlog = mdir / "run.log"
        log.info(f"[MODEL {model}] 输出目录：{mdir}")

        # —— 开跑前：归档旧产物 + run.log 轮转
        archive_previous_run(mdir, log)
        rotate_runlog(runlog, log)

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
            # 致命错误检查（不含初始化失败）
            detect_cafe_errors_fatal(stdout_text, model, LOG_FILE)

            # 记录模型 run.log（按轮次分块）
            append_runlog(runlog, f"MODEL={model} ROUND={r+1}", out_lines)

            # 成功判定 / 自动修正
            if need_autofix(stdout_text):
                ogs = parse_extreme_ogs_from_stdout(stdout_text)

                if ogs and r < rounds:
                    # 交叉校验：当前 family 是否能命中这些 OG
                    hit_n = count_hits_in_family(current_family, ogs)
                    if hit_n == 0:
                        # 说明分隔符/列解析有问题，直接报错，避免“删除 0”继续循环
                        log.error(f"[MODEL {model}] 解析到 {len(ogs)} 个极端 OG，但在 {current_family.name} 中命中数为 0。"
                                  f"请检查 family 的分隔符/首列是否为 OG ID。")
                        success = False
                        break

                    # 记录被删 OG
                    removed_tsv = mdir / f"autofix_removed_round{r+1}.tsv"
                    with removed_tsv.open("w", encoding="utf-8") as w:
                        w.write("OG\n"); w.write("\n".join(ogs) + "\n")

                    # 生成下一轮 family
                    next_family = mdir / f"family.autofix{r+1}.tsv"
                    filter_family_drop_ogs(current_family, ogs, next_family, log)
                    current_family = next_family.resolve()

                    # 移走本轮新产物，避免被 13 聚合旧轮
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
            break  # 任一模型失败即整体退出

    if all_ok:
        write_done(Path(paths["cafe_run_dir"]) / ".cafe.done")
        log.info("[DONE] CAFE5 模型运行完成")
    else:
        sys.stderr.write("[ERR] CAFE5 自动修正用尽或删除失败；请检查 family/树/注释质量与分隔符问题。\n")
        sys.exit(2)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)