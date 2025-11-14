#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— 运行 CAFE5（教程增强版）
特性汇总（全部由 config.yaml:cafe5 控制；不开则与旧版一致）：
  1) 自动自救：捕获 “Failed to initialize any reasonable values” 与 “largest differentials” 名单，
     按第2列(OG)剔除极端家族，生成 family.autofixN.tsv，最多 N 轮（max_autofix_rounds）。
  2) 目录整洁：开跑前归档旧产物(_prev_run_时间)；run.log 轮转；CAFE 用 -o . 扁平输出。
  3) 两阶段大拷贝家族（对齐官方 2.2.4 & 3.1.2）：
       - Stage0：用阈值 copy_threshold 将 input/family.tsv 分成 primary 与 large；
       - Stage1：在 primary 估计 λ（含自救）；
       - Stage3：在 large 用 Stage1 的全局 λ 固定 -l 再跑（可选叠加误差模型）。
  4) 误差模型（对齐官方 3.4）：
       - error_model.mode=estimate 先估计误差文件，再对指定子集 apply_to 重跑带 -e；
       - error_model.mode=use 直接引用外部误差文件（-e<file>），用于指定子集。
  5) 多-λ（对齐官方 3.1.3）：
       - multi_lambda.enable=true 且提供 y_tree，则在 primary 同时跑 global 与 multi(-y)；
       - 无论是否比较，large 阶段一律使用“global-λ”固定 -l（与官方建议一致）。
  6) 产物结构（不加深无谓层级）：
       results/06_cafe/models/<model>/
         primary_global/   …… 主流程（global-λ）
         primary_multi/    …… 主流程（multi-λ，若启用）
         primary_global_e/ …… 带误差（若启用 estimate/use）
         primary_multi_e/  …… 带误差（若启用并适用）
         large/            …… 大拷贝补集（固定 -l）        [可有]
         large_e/          …… 大拷贝补集 + 误差模型       [可有]
         error_model/      …… 误差模型文件（estimate）     [可有]
         flags/high_fail_ogs.list …… 质量提示提取的高失败率家族
         sentinels/.done_* …… 分阶段完成哨兵
         family.autofix*.tsv、autofix_removed_round*.tsv、_round*_tmp/、_prev_run_* …… 自救轨迹与归档
"""

from __future__ import annotations
import sys, io, logging, subprocess, re, shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
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

# ====================== CAFE5 常量/解析 ======================

FAILED_INIT_PAT = re.compile(r"Failed to initialize any reasonable values", re.I)
ADVICE_PAT      = re.compile(r"You may want to try removing the top few families", re.I)
BLOCK_HEAD_PAT  = re.compile(r"Families with largest size differentials\s*:", re.I)
OG_LINE_PAT     = re.compile(r"^\s*(OG\d+)\s*:\s*(\d+)\s*$", re.I)
LAMBDA_PAT      = re.compile(r"\bLambda:\s*([0-9Ee+\-\.]+(?:\s*,\s*[0-9Ee+\-\.]+)*)")

def detect_cafe_errors_fatal(stdout: str, model: str):
    """致命型伪成功；初始化失败类交由 need_autofix 决定。"""
    text = stdout or ""
    if re.search(r"was not found in gene family", text, re.I):
        raise RuntimeError(f"[ERR] CAFE 输出含 'was not found in gene family' —— 树/矩阵物种不一致；模型：{model}")
    if re.search(r"unrecognized option", text, re.I) or re.search(r"Unrecognized parameter", text, re.I):
        raise RuntimeError(f"[ERR] CAFE 输出含 'unrecognized option/parameter' —— 版本不支持该参数；模型：{model}")
    if re.search(r"Failed to open", text, re.I):
        raise RuntimeError(f"[ERR] CAFE 输出含 'Failed to open' —— family/tree 路径无效；模型：{model}")

def need_autofix(stdout: str) -> bool:
    return bool(stdout and (FAILED_INIT_PAT.search(stdout) or ADVICE_PAT.search(stdout)))

def parse_extreme_ogs_from_stdout(stdout: str) -> List[str]:
    if not stdout: return []
    lines = stdout.splitlines()
    start = -1
    for i in range(len(lines)-1, -1, -1):
        if BLOCK_HEAD_PAT.search(lines[i]): start = i; break
    if start < 0: return []
    ogs: List[str] = []
    for j in range(start+1, len(lines)):
        m = OG_LINE_PAT.match(lines[j])
        if not m: break
        ogs.append(m.group(1))
    # 去重保序
    out, seen = [], set()
    for og in ogs:
        if og not in seen: seen.add(og); out.append(og)
    return out

def parse_lambda_from_results_txt(p: Path) -> Optional[float]:
    """从 *results.txt 抓 'Lambda: x' 的首个数（global-λ）。"""
    if not p.is_file(): return None
    m = LAMBDA_PAT.search(p.read_text(encoding="utf-8", errors="ignore"))
    if not m: return None
    s = m.group(1).split(",")[0].strip()
    try: return float(s)
    except Exception: return None

# ====================== family（列固定：第1=Desc，第2=OG） ======================

def extract_og_col2(line: str) -> str:
    fields = re.split(r"\s+", line.strip())
    return fields[1] if len(fields) >= 2 else ""

def count_hits_in_family_by_col2(family_path: Path, ogs: List[str]) -> int:
    targets = set(ogs); hits = 0
    with family_path.open("r", encoding="utf-8", errors="ignore") as fin:
        _ = fin.readline()
        for ln in fin:
            if not ln.strip(): continue
            if extract_og_col2(ln) in targets: hits += 1
    return hits

def filter_family_drop_ogs_by_col2(in_family: Path, ogs: List[str], out_family: Path, log: logging.Logger):
    ogset = set(ogs); n_in = n_drop = n_out = 0
    with in_family.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_family.open("w", encoding="utf-8") as fout:
        header = fin.readline()
        if not header: raise RuntimeError(f"[ERR] family.tsv 为空：{in_family}")
        fout.write(header)
        for ln in fin:
            if not ln.strip(): continue
            n_in += 1
            if extract_og_col2(ln) in ogset:
                n_drop += 1
                continue
            fout.write(ln); n_out += 1
    log.info(f"[FILTER] {in_family.name}：数据行 {n_in}，剔除 {n_drop}，保留 {n_out} → {out_family.name}")

def split_family_primary_large(in_family: Path, out_primary: Path, out_large: Path, threshold: int, log: logging.Logger) -> Tuple[int,int]:
    """将 family.tsv 以“任一物种计数≥threshold”为 large，其余为 primary。"""
    n_p = n_l = 0
    with in_family.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_primary.open("w", encoding="utf-8") as fp, \
         out_large.open("w", encoding="utf-8") as fl:
        header = fin.readline()
        if not header: raise RuntimeError(f"[ERR] family.tsv 为空：{in_family}")
        fp.write(header); fl.write(header)
        for ln in fin:
            if not ln.strip(): continue
            fields = re.split(r"\s+", ln.strip())
            # 第3列开始是各物种计数
            maxcnt = 0
            for x in fields[2:]:
                try:
                    v = int(x)
                except Exception:
                    try: v = int(float(x))
                    except Exception: v = 0
                if v > maxcnt: maxcnt = v
            if maxcnt >= threshold:
                fl.write(ln); n_l += 1
            else:
                fp.write(ln); n_p += 1
    log.info(f"[SPLIT] 阈值≥{threshold}：primary={n_p} 行，large={n_l} 行")
    return n_p, n_l

# ====================== 目录与日志辅助 ======================

def list_names(d: Path) -> set: return set(p.name for p in d.iterdir())

def move_round_outputs(mdir: Path, new_names: set, round_idx: int, keep: set, log: logging.Logger):
    dst = ensure_dir(mdir / f"_round{round_idx}_tmp")
    moved = 0
    for name in sorted(new_names):
        if name in keep: continue
        src = mdir / name
        try: shutil.move(str(src), str(dst / name)); moved += 1
        except Exception: pass
    log.info(f"[CLEAN] 第 {round_idx} 轮新产物 {moved} 个已移入：{dst}")

def append_runlog(runlog: Path, header: str, lines: List[str]):
    with runlog.open("a", encoding="utf-8") as w:
        w.write(f"\n===== {header} =====\n")
        for ln in lines:
            w.write(ln + ("\n" if not ln.endswith("\n") else ""))

def archive_previous_run(mdir: Path, log: logging.Logger):
    items = list(mdir.iterdir())
    if not items: return
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dst = ensure_dir(mdir / f"_prev_run_{stamp}")
    moved = 0
    for p in items:
        if p.name.startswith("run.log"): continue
        try: shutil.move(str(p), str(dst / p.name)); moved += 1
        except Exception: pass
    if moved: log.info(f"[ARCHIVE] 发现历史产物，已归档 {moved} 项 → {dst}")

def rotate_runlog(runlog: Path, log: logging.Logger):
    if runlog.exists():
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        bak = runlog.with_name(f"run.log.{stamp}")
        try: runlog.rename(bak); log.info(f"[LOG] 轮转 run.log → {bak.name}")
        except Exception:
            shutil.copy2(str(runlog), str(bak)); runlog.unlink(missing_ok=True)
            log.info(f"[LOG] 复制轮转 run.log → {bak.name}")
    runlog.touch()

def take_last_block(runlog: Path) -> str:
    if not runlog.is_file(): return ""
    text = runlog.read_text(encoding="utf-8", errors="ignore")
    parts = re.split(r"^=+\s*MODEL=.*?ROUND=.*?=+\s*$", text, flags=re.M)
    return (parts[-1].strip() if parts else text)

def extract_high_fail_ogs_from_runlog(runlog: Path) -> List[str]:
    """从 run.log 最后块抓 'failure rates >20%' 区段的 OG 列表。"""
    blk = take_last_block(runlog)
    if not blk: return []
    out = []
    grab = False
    for ln in blk.splitlines():
        if re.search(r"failure rates >\s*20% of the time", ln, re.I):
            grab = True; continue
        if grab:
            m = re.search(r"(OG\d+)", ln)
            if m: out.append(m.group(1))
            elif ln.strip() == "": break
    # 去重保序
    seen=set(); uniq=[]
    for og in out:
        if og not in seen: seen.add(og); uniq.append(og)
    return uniq

# ====================== 运行 CAFE5 的公共函数 ======================

def run_cafe_streaming(cmd: List[str], cwd: Path, log: logging.Logger, tag: str) -> List[str]:
    log.info(f"[{tag}] CMD: " + " ".join(map(str, cmd)))
    proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1, universal_newlines=True)
    lines: List[str] = []
    log.info(f"[{tag}] ---- CAFE5 输出开始 ----")
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip("\n"); lines.append(line)
        log.info(f"[{tag}] {line}")
    proc.wait()
    log.info(f"[{tag}] ---- CAFE5 输出结束 ----")
    if proc.returncode != 0:
        raise RuntimeError(f"[ERR] CAFE5 退出码 {proc.returncode}（{tag}）")
    return lines

def run_with_autofix(
    cafe_bin: str, workdir: Path, family: Path, nwk: Path,
    threads: int, k: Optional[int], pval: Optional[float],
    add_y: Optional[Path], add_e: Optional[str], efile: Optional[Path],
    max_rounds: int, log: logging.Logger, tag: str
) -> Tuple[Path, Path]:
    """
    在 workdir 下执行带自救的 CAFE5：
      - add_y: 分簇树路径（-y），或 None；
      - add_e: "estimate" | "use" | None；
      - efile: 误差文件路径（当 add_e=="use" 时有效）；
    返回： (最终使用的 family 路径, results.txt 路径)
    """
    runlog = workdir / "run.log"
    rotate_runlog(runlog, log)

    current_family = family.resolve()
    rounds = max(0, max_rounds)

    for r in range(0, rounds + 1):
        before = list_names(workdir)
        cmd = [cafe_bin, "-i", str(current_family), "-t", str(nwk), "-c", str(threads), "-o", "."]
        if k: cmd += ["-k", str(k)]
        if pval is not None: cmd += ["-P", str(pval)]
        if add_y: cmd += ["-y", str(add_y)]
        if add_e == "estimate":
            cmd += ["-e"]  # 按官方：估计误差模型
        elif add_e == "use" and efile:
            cmd += [f"-e{efile}"]

        out_lines = run_cafe_streaming(cmd, workdir, log, tag + f" ROUND={r+1}")
        stdout_text = "\n".join(out_lines)

        detect_cafe_errors_fatal(stdout_text, tag)
        append_runlog(runlog, f"{tag} ROUND={r+1}", out_lines)

        if need_autofix(stdout_text):
            if r >= rounds:
                log.error(f"[{tag}] 检测到初始化失败，但已无剩余轮次。")
                raise RuntimeError("[ERR] 自动修正用尽")
            ogs = parse_extreme_ogs_from_stdout(stdout_text)
            if not ogs:
                log.error(f"[{tag}] 初始化失败但未解析到极端 OG。")
                raise RuntimeError("[ERR] 未解析到极端 OG")

            hit_n = count_hits_in_family_by_col2(current_family, ogs)
            if hit_n == 0:
                log.error(f"[{tag}] 解析到 {len(ogs)} 个极端 OG，但在 {current_family.name} 第2列命中数为 0。")
                raise RuntimeError("[ERR] 解析/分隔异常，终止。")

            removed_tsv = workdir / f"autofix_removed_round{r+1}.tsv"
            with removed_tsv.open("w", encoding="utf-8") as w:
                w.write("OG\n"); w.write("\n".join(ogs) + "\n")

            next_family = workdir / f"family.autofix{r+1}.tsv"
            filter_family_drop_ogs_by_col2(current_family, ogs, next_family, log)
            current_family = next_family.resolve()

            after = list_names(workdir)
            new_names = after - before
            keep = {runlog.name, removed_tsv.name, next_family.name}
            move_round_outputs(workdir, new_names, r+1, keep, log)

            log.info(f"[{tag}] 第 {r+1} 轮自动修正完成，进入下一轮。")
            continue
        else:
            log.info(f"[{tag}] 收敛，无失败提示。")
            break

    # 返回 results.txt（Gamma_results.txt / Base_results.txt 之一）
    # 优先找 *results.txt（含“Result:”与“Lambda:”）
    cand = sorted(workdir.glob("*results.txt"))
    if not cand:
        raise RuntimeError(f"[ERR] 未找到 *results.txt 于 {workdir}")
    return current_family, cand[0].resolve()

# ====================== 主流程 ======================

def main():
    cfg   = load_config()
    paths = cfg["paths"]; cafe = cfg.get("cafe5", {}) or {}

    logs_dir = Path(paths["logs_dir"])
    LOG_FILE = logs_dir / "12_cafe5_run_models.log"
    log = get_logger("aphylo.12", LOG_FILE)
    banner(log, "APhylo 12 — CAFE5（教程增强版 · 两阶段/误差/多-λ）")

    if not cafe.get("enable", True):
        log.info("cafe5.enable=false —— 跳过本步")
        write_done(Path(paths["cafe_run_dir"]) / ".cafe.done"); return

    threads = int(cafe.get("threads", 4))
    k       = cafe.get("gamma_k", None)
    pval    = cafe.get("pvalue", None)
    models  = cafe.get("models", ["global"])
    max_autofix_rounds = int(cafe.get("max_autofix_rounds", 0))

    two_stage = cafe.get("two_stage_large", {}) or {}
    large_enable   = bool(two_stage.get("enable", False))
    copy_threshold = int(two_stage.get("copy_threshold", 100))
    primary_tag    = two_stage.get("primary_tag", "primary")
    large_tag      = two_stage.get("large_tag",   "large")

    em_cfg   = cafe.get("error_model", {}) or {}
    em_mode  = em_cfg.get("mode", "off")            # off | estimate | use
    em_file  = em_cfg.get("file", None)
    em_apply = set(em_cfg.get("apply_to", ["primary","large"]))

    ml_cfg   = cafe.get("multi_lambda", {}) or {}
    ml_enable = bool(ml_cfg.get("enable", False))
    ml_y     = ml_cfg.get("y_tree", None)
    ml_compare = bool(ml_cfg.get("compare_with_global", True))

    cafe_bin = cfg.get("binaries", {}).get("cafe5", "cafe5")

    # 输入
    inp   = need_dir(Path(paths["cafe_run_dir"]) / "input", "CAFE 输入目录")
    ts    = sorted(inp.glob("*.tsv"))
    nwks  = sorted(inp.glob("*.nwk"))
    if not ts:  raise FileNotFoundError("[ERR] 未找到 family TSV（*.tsv）")
    if not nwks:raise FileNotFoundError("[ERR] 未找到超时钟树（*.nwk）")
    family0 = ts[0].resolve()
    nwk     = nwks[0].resolve()
    log.info(f"[IN] family: {family0}")
    log.info(f"[IN] tree  : {nwk}")
    if max_autofix_rounds>0: log.info(f"[CFG] max_autofix_rounds = {max_autofix_rounds}")

    # 输出根
    models_dir = ensure_dir(Path(paths["cafe_run_dir"]) / "models")

    # 逐模型运行（当前多为 'global'，但保留接口）
    for model in models:
        mname = str(model).replace(":", "_")
        mdir  = ensure_dir(models_dir / mname)
        log.info(f"[MODEL {model}] 输出目录：{mdir}")
        archive_previous_run(mdir, log)

        # 子目录规划
        primary_dir         = ensure_dir(mdir / f"{primary_tag}_global")   # global-λ 主流程目录
        primary_dir_multi   = ensure_dir(mdir / f"{primary_tag}_multi")    if ml_enable else None
        primary_dir_global_e= ensure_dir(mdir / f"{primary_tag}_global_e") if em_mode!="off" and ("primary" in em_apply) else None
        primary_dir_multi_e = ensure_dir(mdir / f"{primary_tag}_multi_e")  if ml_enable and em_mode!="off" and ("primary" in em_apply) else None
        large_dir           = ensure_dir(mdir / f"{large_tag}")            if large_enable else None
        large_dir_e         = ensure_dir(mdir / f"{large_tag}_e")          if large_enable and em_mode!="off" and ("large" in em_apply) else None
        err_dir             = ensure_dir(mdir / "error_model")             if em_mode=="estimate" else None
        flags_dir           = ensure_dir(mdir / "flags")
        sent_dir            = ensure_dir(mdir / "sentinels")

        # Stage0：按阈值拆分 primary / large
        fam_primary = family0
        fam_large   = None
        if large_enable:
            fam_primary = mdir / "family.primary.tsv"
            fam_large   = mdir / "family.large.tsv"
            split_family_primary_large(family0, fam_primary, fam_large, copy_threshold, log)
            if fam_primary.stat().st_size == 0:
                raise RuntimeError("[ERR] primary 集为空，阈值过严？")

        # Stage1：在 primary 估计 global-λ（可含自救）
        fam_used_primary, res_primary = run_with_autofix(
            cafe_bin=cafe_bin, workdir=primary_dir, family=fam_primary, nwk=nwk,
            threads=threads, k=k, pval=pval,
            add_y=None, add_e=None, efile=None,
            max_rounds=max_autofix_rounds, log=log, tag=f"MODEL={model} PRIMARY-GLOBAL"
        )
        glambda = parse_lambda_from_results_txt(res_primary)
        if glambda is None:
            raise RuntimeError("[ERR] 未能从 primary_global 的 *results.txt 解析出 Lambda")

        # 提取 high-fail 名单
        hf = extract_high_fail_ogs_from_runlog(primary_dir / "run.log")
        if hf:
            (flags_dir/"high_fail_ogs.list").write_text("\n".join(hf)+"\n", encoding="utf-8")
            log.info(f"[FLAG] high-fail 家族数：{len(hf)}")

        # Stage1b：multi-λ（若启用）
        y_tree = Path(ml_y).resolve() if (ml_enable and ml_y) else None
        if ml_enable:
            need_compare_global = ml_compare
            # 真正的 multi 跑
            run_with_autofix(
                cafe_bin=cafe_bin, workdir=primary_dir_multi, family=fam_used_primary, nwk=nwk,
                threads=threads, k=k, pval=pval,
                add_y=y_tree, add_e=None, efile=None,
                max_rounds=max_autofix_rounds, log=log, tag=f"MODEL={model} PRIMARY-MULTI"
            )
            # 若比较，global 已有；若未比较但 large 需要 λ，我们也已有 glambda。

        # Stage2：误差模型
        efile_path: Optional[Path] = None
        if em_mode == "estimate":
            # 在 primary_global 估计误差模型（-e）
            _, _ = run_with_autofix(
                cafe_bin=cafe_bin, workdir=ensure_dir(mdir / "error_model"), family=fam_used_primary, nwk=nwk,
                threads=threads, k=k, pval=pval,
                add_y=None, add_e="estimate", efile=None,
                max_rounds=0, log=log, tag=f"MODEL={model} ERROR-MODEL"
            )
            # 捕获 *_error_model.txt
            cand = sorted((mdir/"error_model").glob("*_error_model.txt"))
            if not cand:
                raise RuntimeError("[ERR] 未找到 *_error_model.txt（误差模型估计失败？）")
            # 将第一个作为模型文件
            efile_path = cand[0].resolve()
            log.info(f"[EM] 误差模型文件：{efile_path}")

        elif em_mode == "use":
            if not em_file: raise RuntimeError("[ERR] error_model.mode=use 但未指定 file")
            efile_path = Path(em_file).resolve()
            need_file(efile_path, "误差模型文件")
            log.info(f"[EM] 使用外部误差模型：{efile_path}")

        # 若需要对 primary 叠加误差模型重跑
        if efile_path and ("primary" in em_apply):
            # global 带 -e
            run_with_autofix(
                cafe_bin=cafe_bin, workdir=primary_dir_global_e, family=fam_used_primary, nwk=nwk,
                threads=threads, k=k, pval=pval,
                add_y=None, add_e="use", efile=efile_path,
                max_rounds=max_autofix_rounds, log=log, tag=f"MODEL={model} PRIMARY-GLOBAL(e)"
            )
            # multi 带 -e（若启用）
            if ml_enable and primary_dir_multi_e:
                run_with_autofix(
                    cafe_bin=cafe_bin, workdir=primary_dir_multi_e, family=fam_used_primary, nwk=nwk,
                    threads=threads, k=k, pval=pval,
                    add_y=y_tree, add_e="use", efile=efile_path,
                    max_rounds=max_autofix_rounds, log=log, tag=f"MODEL={model} PRIMARY-MULTI(e)"
                )

        # Stage3：large 集 —— 固定 global-λ 跑（可叠加误差）
        if large_enable and fam_large and fam_large.stat().st_size>0:
            # large（无误差）
            run_with_autofix(
                cafe_bin=cafe_bin, workdir=large_dir, family=fam_large, nwk=nwk,
                threads=threads, k=k, pval=pval,
                add_y=None, add_e=None, efile=None,
                max_rounds=0, log=log, tag=f"MODEL={model} LARGE(-l)"
            )  # 先跑一遍拿到正常文件（有些版本 -l 也可直接跑）
            # 用 -l 固定 λ（用 glambda）
            # 注意：为稳妥，这里直接再跑一遍，真正固定 -l；（简化起见，复用 run_cafe_streaming）
            cmd = [cafe_bin, "-i", str(fam_large), "-t", str(nwk), "-c", str(threads), "-o", ".", "-l", f"{glambda}"]
            if k: cmd += ["-k", str(k)]
            if pval is not None: cmd += ["-P", str(pval)]
            out_lines = run_cafe_streaming(cmd, large_dir, log, f"MODEL={model} LARGE -l {glambda}")
            append_runlog(large_dir/"run.log", f"MODEL={model} LARGE(-l) FINAL", out_lines)

            # large_e（若误差模型对 large 生效）
            if efile_path and ("large" in em_apply) and large_dir_e:
                cmd_e = [cafe_bin, "-i", str(fam_large), "-t", str(nwk), "-c", str(threads), "-o", ".", f"-e{efile_path}", "-l", f"{glambda}"]
                if k: cmd_e += ["-k", str(k)]
                if pval is not None: cmd_e += ["-P", str(pval)]
                out_lines = run_cafe_streaming(cmd_e, large_dir_e, log, f"MODEL={model} LARGE_e -l {glambda}")
                append_runlog(large_dir_e/"run.log", f"MODEL={model} LARGE_e(-l) FINAL", out_lines)

        # 标记完成
        write_done(sent_dir/".done_primary")
        if large_enable and fam_large and fam_large.stat().st_size>0:
            write_done(sent_dir/".done_large")
        if em_mode=="estimate": write_done(sent_dir/".done_error_model")

    # 总体完成
    write_done(Path(paths["cafe_run_dir"]) / ".cafe.done")
    log.info("[DONE] CAFE5 教程增强版全部完成")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)