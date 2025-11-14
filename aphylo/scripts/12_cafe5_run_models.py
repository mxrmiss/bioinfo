#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— CAFE5 两阶段（教程增强版/误差/流式输出/自动清理）

使用约定（与 aphylo 配置对齐）：
- 所有路径与参数从 config.yaml 读取（不接受命令行参数）。
- 关键键：
  paths:
    cafe_run_dir: "results/06_cafe"          # 本模块根目录
    logs_dir    : "logs"                      # 项目日志目录（脚本主日志放这里）
  inputs:
    family_tsv     : "results/06_cafe/input/family.tsv"     # family 输入
    ultrametric_tree: "results/06_cafe/input/utree_for_cafe.nwk"  # 超时树
  cafes:
    threads: 30
    k: 3
    P: 0.05
    models: ["global"]
    max_autofix_rounds: 3
    two_stage_large:
      enable: true
      copy_threshold: 100                      # 任一物种计数>=阈值 → large
    error_model:
      mode: "off"                              # off | estimate | use
      file: null                               # mode=use 时指定；estimate 时忽略
      apply_to: "both"                         # primary | both
    clean_on_start: true                       # 开跑前自动清空 models/ 与脚本主日志
    clean_strict_safety: true                  # 清理路径安全校验
"""

from __future__ import annotations
import os, sys, re, shutil, subprocess, datetime, io
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
import yaml

DEFAULT_CONFIG = "config.yaml"
SCRIPT_LOG_NAME = "12_cafe5_run_models.log"

# ------------------------- 工具：配置与输出 -------------------------

def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj

def load_config(cfg_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(cfg_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        for key in ("inputs", "paths", "cafes"):
            if key in cfg:
                cfg[key] = _expand_publish_placeholders(cfg[key], str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def banner(text: str) -> str:
    bar = "=" * max(16, len(text) + 4)
    return f"{bar}\n  {text}\n{bar}\n"

def now_ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# ------------------------- 工具：日志/流式 tee -------------------------

class TeeLogger:
    """
    将文本同时写入“主日志文件”和“屏幕”，并确保 flush。
    """
    def __init__(self, logfile: Path):
        ensure_dir(logfile.parent)
        # 主日志文件：每次开跑会清空
        self.logfile = logfile
        self.fh = open(self.logfile, "w", encoding="utf-8")
        # 包个 stdout/stderr 包装器，确保 flush
        class _Flush(io.TextIOBase):
            def __init__(self, inner): self.inner = inner
            def write(self, x):
                if x:
                    self.inner.write(x); self.inner.flush()
                return len(x)
        self.stdout = _Flush(sys.stdout)
        self.stderr = _Flush(sys.stderr)

    def write(self, text: str):
        self.stdout.write(text)
        self.fh.write(text)
        self.fh.flush()

    def close(self):
        try: self.fh.close()
        except: pass

# ------------------------- 工具：安全清理 -------------------------

def safe_clean_models(models_dir: Path, script_log: Path, enable: bool, strict: bool, tee: TeeLogger):
    if not enable:
        tee.write(f"[{now_ts()}] [CLEAN] 跳过自动清理（clean_on_start=false）。\n")
        return
    if not models_dir.exists():
        tee.write(f"[{now_ts()}] [CLEAN] models/ 不存在，无需清理。\n")
    else:
        if strict:
            # 基础保护：路径名必须包含 06_cafe 与 models
            safe_ok = ("06_cafe" in str(models_dir)) and (models_dir.name == "models")
            # 进一步保护：目录下如有“非模型特征”则拒绝
            suspicious = any(p.name in (".git", "scripts", "config.yaml") for p in models_dir.rglob("*"))
            if not safe_ok or suspicious:
                tee.write(f"[{now_ts()}] [CLEAN] 保护触发：拒绝清理 {models_dir}\n")
                raise SystemExit(2)
        shutil.rmtree(models_dir, ignore_errors=True)
        tee.write(f"[{now_ts()}] [CLEAN] 已清空：{models_dir}\n")
    # 清空脚本主日志（即将被当前进程重建）
    if script_log.exists():
        try: script_log.unlink()
        except: pass
    tee.write(f"[{now_ts()}] [CLEAN] 已重置脚本主日志：{script_log.name}\n")

# ------------------------- 工具：family 读取/拆分/修正 -------------------------

def read_family(family_tsv: Path) -> Tuple[List[str], List[List[str]]]:
    lines = family_tsv.read_text(encoding="utf-8").splitlines()
    if not lines:
        return [], []
    header = re.split(r"\t", lines[0].rstrip("\n"))
    body = [re.split(r"\t", ln.rstrip("\n")) for ln in lines[1:] if ln.strip()]
    return header, body

def write_family(path: Path, header: List[str], rows: List[List[str]]):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

def split_primary_large(header: List[str], rows: List[List[str]], threshold: int) -> Tuple[List[List[str]], List[List[str]]]:
    # 约定：首列 Desc，第二列 OG ID，从第三列起为物种计数
    pri, lg = [], []
    for r in rows:
        # 容错：不足列数直接归 primary
        if len(r) < 3:
            pri.append(r); continue
        counts = []
        for x in r[2:]:
            try:
                counts.append(int(x))
            except:
                # NA/空白等，按 0 处理
                counts.append(0)
        if counts and max(counts) >= threshold:
            lg.append(r)
        else:
            pri.append(r)
    return pri, lg

def remove_extreme_ogs_from_family(header: List[str], rows: List[List[str]], og_list: List[str]) -> Tuple[List[List[str]], List[List[str]]]:
    """
    从 family 中删除第二列（OG ID）命中 og_list 的行；返回 (保留行, 被删行)
    """
    og_set = set(og_list)
    kept, removed = [], []
    for r in rows:
        og = r[1] if len(r) > 1 else ""
        (removed if og in og_set else kept).append(r)
    return kept, removed

# ------------------------- 工具：解析 run.log（稳健，逐行） -------------------------

def parse_extreme_ogs_from_runlog(runlog_path: Path) -> List[str]:
    """
    解析“最后一个 blocks”极端家族列表：
    ... 'Families with largest size differentials:' 之后，连续的 'OGxxxxxxx: N' 行
    """
    if not runlog_path.exists():
        return []
    lines = runlog_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    # 找到所有块的起点行号
    anchors = [i for i, ln in enumerate(lines) if "Families with largest size differentials" in ln]
    if not anchors:
        return []
    i = anchors[-1] + 1
    ogs: List[str] = []
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            break
        # 逐行提取 'OG\d+:\s*\d+'
        m = re.search(r"(OG\d+)\s*:\s*\d+", ln)
        if m:
            ogs.append(m.group(1))
            i += 1
            continue
        # 非 OG 格式即视为结束
        break
    return ogs

def parse_high_fail_from_runlog(runlog_path: Path) -> List[str]:
    """
    解析“失败率>20%”清单：定位‘had failure rates >20%’行后，每行 'OGxxxxx had XX failures'
    """
    if not runlog_path.exists():
        return []
    lines = runlog_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    anchors = [i for i, ln in enumerate(lines) if "had failure rates >20%" in ln]
    if not anchors:
        return []
    i = anchors[-1] + 1
    out: List[str] = []
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            break
        m = re.search(r"(OG\d+)\b.*\bhad\b", ln)
        if m:
            out.append(m.group(1))
            i += 1
            continue
        # 非 OG 行结束
        break
    return out

def parse_lambda_from_results(out_dir: Path) -> Optional[float]:
    """
    从 Gamma_results.txt 或 Base_results.txt 解析 'Lambda: <float>'
    """
    for name in ("Gamma_results.txt", "Base_results.txt"):
        p = out_dir / name
        if p.exists():
            txt = p.read_text(encoding="utf-8", errors="ignore")
            m = re.search(r"Lambda\s*:\s*([0-9Ee\+\-\.]+)", txt)
            if m:
                try:
                    return float(m.group(1))
                except:
                    pass
    return None

# ------------------------- 工具：执行 CAFE5，流式输出 + 日志轮转 -------------------------

def rotate_log(logfile: Path):
    if logfile.exists():
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        bak = logfile.with_name(f"{logfile.stem}_{ts}{logfile.suffix}")
        try:
            logfile.rename(bak)
        except:
            pass

def run_cafe5_stream(
    cafe_bin: str,
    cwd: Path,
    family_file: Path,
    tree_file: Path,
    threads: int,
    k: int,
    P: float,
    tee: TeeLogger,
    stage_title: str,
    error_mode: str = "off",
    error_file: Optional[Path] = None,
    apply_error: bool = False,
    lambda_fix: Optional[float] = None,
) -> None:
    """
    在 cwd 下执行 CAFE5（不使用 -o 子目录），所有产物直接写入 cwd。
    将 stdout/stderr 同步写入 cwd/run.log 与屏幕。
    """
    ensure_dir(cwd)
    runlog = cwd / "run.log"
    rotate_log(runlog)

    # 构建命令
    cmd = [cafe_bin, "-i", str(family_file), "-t", str(tree_file), "-c", str(threads), "-k", str(k), "-P", str(P)]
    em_note = ""
    if apply_error:
        if error_mode == "estimate":
            cmd.append("-e")
            em_note = "【误差模型】估计并应用（epsilon 校正计数偏差）"
        elif error_mode == "use" and error_file:
            cmd.append("-e" + str(error_file))
            em_note = f"【误差模型】应用外部文件：{error_file.name}"
        else:
            em_note = "【误差模型】未应用（配置/文件缺失）"
    else:
        em_note = "【误差模型】未应用（apply_to=primary，仅 primary 阶段启用）"

    if lambda_fix is not None:
        cmd.extend(["-l", str(lambda_fix)])
        em_note += f"；固定 λ = {lambda_fix}"

    # 标题与命令回显
    tee.write(banner(stage_title))
    tee.write(f"[{now_ts()}] [CMD] {' '.join(cmd)}\n")
    tee.write(f"[{now_ts()}] [INFO] {em_note}\n")

    # 以 cwd 作为工作目录，实时读取输出
    with open(runlog, "w", encoding="utf-8") as rfh:
        proc = subprocess.Popen(
            cmd, cwd=str(cwd),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, universal_newlines=True
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            tee.write(line)
            rfh.write(line)
        proc.wait()

    tee.write(f"[{now_ts()}] [INFO] ---- CAFE5 输出结束 ----\n")

# ------------------------- 主流程 -------------------------

def main():
    cfg = load_config()

    # 读取路径
    paths = cfg.get("paths", {})
    cafe_root = Path(paths.get("cafe_run_dir", "results/06_cafe")).resolve()
    logs_dir  = Path(paths.get("logs_dir", "logs")).resolve()
    models_dir = cafe_root / "models"
    flags_dir  = models_dir / "global" / "flags"  # 先占位，实际创建在模型循环里
    sent_dir   = models_dir / "global" / "sentinels"

    # 主日志 tee
    script_log = Path(logs_dir) / SCRIPT_LOG_NAME
    tee = TeeLogger(script_log)

    # 配置：输入、参数
    inputs = cfg.get("inputs", {})
    family_tsv = Path(inputs.get("family_tsv", cafe_root / "input" / "family.tsv")).resolve()
    tree_file  = Path(inputs.get("ultrametric_tree", cafe_root / "input" / "utree_for_cafe.nwk")).resolve()

    cafes = cfg.get("cafes", {})
    threads = int(cafes.get("threads", 30))
    k = int(cafes.get("k", 3))
    P = float(cafes.get("P", 0.05))
    models = cafes.get("models", ["global"])
    max_rounds = int(cafes.get("max_autofix_rounds", 3))

    two_stage = cafes.get("two_stage_large", {}) or {}
    do_large = bool(two_stage.get("enable", True))
    copy_threshold = int(two_stage.get("copy_threshold", 100))

    em = cafes.get("error_model", {}) or {}
    em_mode = str(em.get("mode", "off")).lower()          # off|estimate|use
    em_file = Path(em["file"]).resolve() if (em.get("file") and em.get("file") != "null") else None
    em_apply_to = str(em.get("apply_to", "both")).lower() # primary|both

    cafe_bin = shutil.which("cafe5") or "cafe5"

    # 欢迎与配置概览
    tee.write(banner("APhylo 12 — CAFE5（教程增强版/两阶段/误差/流式输出）"))
    tee.write(f"[{now_ts()}] [IN ] family: {family_tsv}\n")
    tee.write(f"[{now_ts()}] [IN ] tree  : {tree_file}\n")
    tee.write(f"[{now_ts()}] [CFG] threads={threads}, k={k}, P={P}, models={models}\n")
    tee.write(f"[{now_ts()}] [CFG] large: {'enable' if do_large else 'disable'}, copy_threshold={copy_threshold}\n")
    if em_mode == "off":
        tee.write(f"[{now_ts()}] [CFG] error model: OFF → 解释：关闭误差校正，λ 可能偏大；与官方教程 3.4 相反。\n")
    elif em_mode == "estimate":
        tee.write(f"[{now_ts()}] [CFG] error model: ESTIMATE, apply_to={em_apply_to} → 解释：估计并应用 epsilon（计数偏差概率），通常使 λ 更稳/更小。\n")
    elif em_mode == "use":
        tee.write(f"[{now_ts()}] [CFG] error model: USE({em_file}), apply_to={em_apply_to} → 解释：直接使用外部误差模型文件进行校正。\n")
    else:
        tee.write(f"[{now_ts()}] [CFG] error model: 未知模式（{em_mode}），视为 OFF。\n")
        em_mode = "off"

    # 自动清理
    safe_clean_models(models_dir, script_log, bool(cafes.get("clean_on_start", True)), bool(cafes.get("clean_strict_safety", True)), tee)

    # 基础检查
    if not family_tsv.exists():
        tee.write(f"[{now_ts()}] [ERR] family.tsv 不存在：{family_tsv}\n"); sys.exit(2)
    if not tree_file.exists():
        tee.write(f"[{now_ts()}] [ERR] ultrametric tree 不存在：{tree_file}\n"); sys.exit(2)

    # 仅支持 global（与现行使用一致）
    for model in models:
        if model != "global":
            tee.write(f"[{now_ts()}] [WARN] 暂仅实现 global；忽略模型 {model}\n")

    # ---------- 准备目录 ----------
    model_root = models_dir / "global"
    primary_dir = model_root / "primary_global"
    large_dir   = model_root / "large"
    flags_dir   = model_root / "flags"
    sent_dir    = model_root / "sentinels"
    ensure_dir(primary_dir); ensure_dir(large_dir); ensure_dir(flags_dir); ensure_dir(sent_dir)

    # ---------- 拆分 family：primary / large ----------
    tee.write(banner("SPLIT — 拆分 primary/large"))
    hdr, rows = read_family(family_tsv)
    if not hdr:
        tee.write(f"[{now_ts()}] [ERR] 空的 family 文件：{family_tsv}\n"); sys.exit(2)
    if not (len(hdr) >= 2 and hdr[0].lower() == "desc"):
        tee.write(f"[{now_ts()}] [WARN] family 首列非 Desc（当前：{hdr[0]}），仍继续按“首列描述/第二列OG”解析。\n")

    pri_rows, lg_rows = split_primary_large(hdr, rows, copy_threshold) if do_large else (rows, [])
    pri_path = model_root / "family.primary.tsv"
    lg_path  = model_root / "family.large.tsv"
    write_family(pri_path, hdr, pri_rows)
    write_family(lg_path, hdr, lg_rows)
    tee.write(f"[{now_ts()}] [INFO] 拆分完成：primary={len(pri_rows)} 行，large={len(lg_rows)} 行\n")

    # ---------- PRIMARY 阶段：支持多轮自动修正 ----------
    primary_family_path = pri_path
    lambda_est: Optional[float] = None
    high_fail_all: List[str] = []
    for r in range(1, max_rounds + 1):
        stage = f"MODEL=global PRIMARY-GLOBAL ROUND={r}"
        apply_error = (em_apply_to in ("primary", "both"))
        run_cafe5_stream(
            cafe_bin=cafe_bin,
            cwd=primary_dir,
            family_file=primary_family_path,
            tree_file=tree_file,
            threads=threads, k=k, P=P,
            tee=tee,
            stage_title=stage,
            error_mode=em_mode, error_file=em_file, apply_error=apply_error,
            lambda_fix=None
        )
        # 每轮后解析极端家族
        ogs = parse_extreme_ogs_from_runlog(primary_dir / "run.log")
        # 解析高失败率 OG（累计）
        high_fail = parse_high_fail_from_runlog(primary_dir / "run.log")
        if high_fail:
            high_fail_all.extend(high_fail)

        # 尝试读 λ
        if lambda_est is None:
            ltmp = parse_lambda_from_results(primary_dir)
            if ltmp is not None:
                lambda_est = ltmp

        if ogs and (r < max_rounds):
            # 从当前 family.primary.tsv 删除极端 OG 并继续下一轮
            _, cur_rows = read_family(primary_family_path)
            kept, removed = remove_extreme_ogs_from_family(hdr, cur_rows, ogs)
            # 记录删除清单
            removed_path = model_root / f"autofix_removed_round{r}.tsv"
            write_family(removed_path, hdr, removed)
            # 生成新 family 文件
            primary_family_path = model_root / f"family.autofix{r}.tsv"
            write_family(primary_family_path, hdr, kept)
            tee.write(f"[{now_ts()}] [CLEAN] 第 {r} 轮剔除 {len(removed)} 个极端家族 → {primary_family_path.name}，进入下一轮。\n")
            continue
        else:
            if ogs:
                tee.write(f"[{now_ts()}] [INFO] 已到最大轮数（{max_rounds}），停止自动修正。\n")
            else:
                tee.write(f"[{now_ts()}] [INFO] 未再发现极端家族，自动修正结束。\n")
            break

    # 汇总高失败率 OG（去重）
    if high_fail_all:
        uniq = sorted(set(high_fail_all))
        (flags_dir / "high_fail_ogs.list").write_text("\n".join(uniq) + "\n", encoding="utf-8")
        tee.write(f"[{now_ts()}] [INFO] 高失败率家族：{len(uniq)} 个 → {flags_dir/'high_fail_ogs.list'}\n")
    else:
        tee.write(f"[{now_ts()}] [INFO] 未发现“失败率>20%”家族。\n")

    # ---------- LARGE 阶段：只跑一次（不再剔除极端 OG） ----------
    if do_large and lg_rows:
        # large 阶段若拿到 λ，则固定 λ 复算（与官方教程一致）
        stage = "MODEL=global LARGE"
        apply_error = (em_apply_to == "both")
        run_cafe5_stream(
            cafe_bin=cafe_bin,
            cwd=large_dir,
            family_file=lg_path,
            tree_file=tree_file,
            threads=threads, k=k, P=P,
            tee=tee,
            stage_title=stage,
            error_mode=em_mode, error_file=em_file, apply_error=apply_error,
            lambda_fix=lambda_est
        )
        # large 阶段也记录高失败率 OG（单独文件）
        hf_large = parse_high_fail_from_runlog(large_dir / "run.log")
        if hf_large:
            (flags_dir / "high_fail_ogs_large.list").write_text("\n".join(sorted(set(hf_large))) + "\n", encoding="utf-8")
            tee.write(f"[{now_ts()}] [INFO] [LARGE] 高失败率家族：{len(set(hf_large))} 个 → {flags_dir/'high_fail_ogs_large.list'}\n")
    else:
        tee.write(f"[{now_ts()}] [INFO] LARGE 阶段已跳过（enable={do_large}, large 行数={len(lg_rows)}）。\n")

    # ---------- 收尾 ----------
    # 写 sentinel
    (sent_dir / ".done_primary_global").write_text(now_ts() + "\n", encoding="utf-8")
    (cafe_root / ".cafe.done").write_text(now_ts() + "\n", encoding="utf-8")

    tee.write(banner("CAFE5 模型运行完成"))
    tee.write(f"[{now_ts()}] [OK ] 输出目录：{model_root}\n")
    tee.write(f"[{now_ts()}] [OK ] 主要结果文件（示例）:\n")
    tee.write(f"      - {primary_dir/'Gamma_family_results.txt'}（家族显著性）\n")
    tee.write(f"      - {primary_dir/'Gamma_clade_results.txt'}（分支汇总扩张/收缩/快速）\n")
    tee.write(f"      - {primary_dir/'Gamma_report.cafe'}（完整报告）\n")
    if (large_dir/'Gamma_family_results.txt').exists():
        tee.write(f"      - {large_dir/'Gamma_family_results.txt'}（大家族复算）\n")
    if lambda_est is not None:
        tee.write(f"[{now_ts()}] [OK ] 估计到 λ = {lambda_est}\n")
    tee.write(f"[{now_ts()}] [OK ] 结束时间：{now_ts()}\n")
    tee.close()

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        raise
    except Exception as e:
        sys.stderr.write(f"[FATAL] {e}\n")
        sys.exit(2)