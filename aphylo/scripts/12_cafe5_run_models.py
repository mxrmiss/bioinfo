#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 12_cafe5_run_models.py
# 精简版：按 config 运行 CAFE5；清理旧产物；primary 支持极端 OG 自动剔除；large 只跑一次；全程流式输出。

from __future__ import annotations
import os, sys, re, shutil, subprocess, logging, io
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import yaml
from datetime import datetime

CFG = "config.yaml"

# ---------- 基础 ----------
def load_cfg(p: str) -> dict:
    cf = Path(p)
    if not cf.exists():
        sys.exit(f"[ERR] 未找到配置文件：{cf}")
    return yaml.safe_load(cf.read_text(encoding="utf-8")) or {}

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def rm_rf(p: Path):
    if p.is_dir():
        shutil.rmtree(p)
    elif p.is_file():
        p.unlink(missing_ok=True)

def setup_logger(logfile: Path) -> logging.Logger:
    ensure_dir(logfile.parent)
    lg = logging.getLogger("cafe12"); lg.setLevel(logging.INFO); lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s","%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout);          sh.setFormatter(fmt)
    lg.addHandler(fh); lg.addHandler(sh)
    # 让 print 也立刻刷
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s=s
        def write(self,x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout=_Flush(sys.stdout); sys.stderr=_Flush(sys.stderr)
    return lg

# ---------- family 拆分/筛除 ----------
def read_family_table(path: Path) -> Tuple[List[str], List[List[str]]]:
    rows = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not rows: return [], []
    header = [c.strip() for c in rows[0].split("\t")]
    body = [[c.strip() for c in ln.split("\t")] for ln in rows[1:] if ln.strip()]
    return header, body

def write_family_table(header: List[str], body: List[List[str]], path: Path):
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in body:
            f.write("\t".join(r) + "\n")

def split_primary_large(header: List[str], body: List[List[str]], threshold: int) -> Tuple[List[List[str]], List[List[str]]]:
    # 首列固定 Desc，第 2 列 Orthogroup，计数列从第 3 列起
    pri, lg = [], []
    for r in body:
        counts = []
        for c in r[2:]:
            try: counts.append(int(c))
            except: counts.append(0)
        if any(v >= threshold for v in counts):
            lg.append(r)
        else:
            pri.append(r)
    return pri, lg

def remove_ogs_from_family(body: List[List[str]], og_set: set) -> List[List[str]]:
    # 按 Orthogroup（第 2 列）过滤
    out=[]
    for r in body:
        if len(r)>=2 and r[1] in og_set: 
            continue
        out.append(r)
    return out

# ---------- 解析 CAFE5 日志中的极端 OG ----------
_PAT_BLOCK = re.compile(r"Families with largest size differentials:\s*(?:\r?\n)((?:.*OG\d+\s*:\s*\d+\s*\r?\n)+)", re.I)
_PAT_OG = re.compile(r"(OG\d+)\s*:\s*(\d+)")
def parse_extreme_ogs(stdout_text: str) -> List[Tuple[str,int]]:
    blks = _PAT_BLOCK.findall(stdout_text)
    if not blks: return []
    last = blks[-1]
    out=[]
    for m in _PAT_OG.finditer(last):
        og = m.group(1); val = int(m.group(2))
        out.append((og, val))
    return out

# ---------- 子进程流式执行 ----------
def run_stream(cmd: List[str], cwd: Path, tee_fp: Path, logger: logging.Logger) -> Tuple[int, str]:
    ensure_dir(tee_fp.parent)
    logger.info(f"[CMD] {' '.join(cmd)}  (cwd={cwd})")
    proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
    buf=[]
    with tee_fp.open("w", encoding="utf-8") as fw:
        for line in proc.stdout:
            sys.stdout.write(line)  # 流式回显
            fw.write(line)
            buf.append(line)
    code = proc.wait()
    out = "".join(buf)
    if code != 0:
        logger.error(f"[ERR] 进程退出码 {code}")
    return code, out

# ---------- 主流程 ----------
def main():
    cfg = load_cfg(CFG)

    # 取路径/参数（只按您当前配置结构，不做兼容分支）
    pths = cfg.get("paths", {}) or {}
    run_root = Path(pths["cafe_run_dir"])          # results/06_cafe
    logs_root = Path(pths["logs_dir"])             # e.g. results/02_logs 或 logs
    cafes = cfg.get("cafes", {}) or {}
    models = cafes.get("models", {}) or {}
    if not isinstance(models, dict) or not models:
        sys.exit("[ERR] config.cafes.models 为空或不是字典")

    cafe5_bin = cafes.get("cafe5_bin", "cafe5")
    max_rounds = int(cafes.get("max_autofix_rounds", 3))

    # large 阶段
    large_cfg = cafes.get("two_stage_large", {}) or {}
    large_enable = bool(large_cfg.get("enable", True))
    copy_threshold = int(large_cfg.get("copy_threshold", 100))
    tag_primary = large_cfg.get("primary_tag","primary")
    tag_large   = large_cfg.get("large_tag","large")

    # 清理旧产物：models + 12相关日志
    models_dir = run_root / "models"
    if models_dir.exists(): shutil.rmtree(models_dir)
    ensure_dir(models_dir)
    # 清理 12 的顶级日志
    for lf in ["12_cafe5_run_models.log"]:
        (logs_root/lf).unlink(missing_ok=True)

    logger = setup_logger(logs_root / "12_cafe5_run_models.log")
    logger.info("="*45); logger.info(" APhylo 12 — CAFE5 "); logger.info("="*45)

    # 遍历模型
    for model_name, mcfg in models.items():
        tree_path   = Path(mcfg["tree"])
        family_path = Path(mcfg["family"])
        threads     = int(mcfg.get("threads", 8))
        extra       = mcfg.get("cmd_extra", "").strip()

        model_dir = ensure_dir(models_dir / model_name)
        # 拆分输入
        header, body = read_family_table(family_path)
        if not header or header[0] != "Desc" or (len(header)<2 or header[1] not in ("Orthogroup","Orthogroups","OrthogroupID")):
            sys.exit("[ERR] family.tsv 首列必须为 Desc，第二列为 Orthogroup")
        pri_rows, lg_rows = split_primary_large(header, body, copy_threshold)

        # 写入拆分文件（位于模型根目录，便于查看）
        f_primary = model_dir / f"family.{tag_primary}.tsv"
        f_large   = model_dir / f"family.{tag_large}.tsv"
        write_family_table(header, pri_rows, f_primary)
        write_family_table(header, lg_rows, f_large)

        # ---------- primary_global：多轮剔除极端 OG ----------
        pg_dir   = ensure_dir(model_dir / "primary_global")
        pg_log   = pg_dir / "run.log"
        cur_in   = pg_dir / "family.tsv"
        shutil.copy2(f_primary, cur_in)

        round_idx = 1
        while True:
            logger.info(f"[{model_name} PRIMARY-GLOBAL ROUND-{round_idx}] CAFE5 输出开始")
            cmd = [cafe5_bin, "-i", str(cur_in), "-t", str(tree_path), "-c", str(threads)]
            if extra: cmd += extra.split()
            code, out = run_stream(cmd, cwd=pg_dir, tee_fp=pg_log, logger=logger)
            logger.info(f"[{model_name} PRIMARY-GLOBAL ROUND-{round_idx}] ---- CAFE5 输出结束 ----")

            # 解析是否需要下一轮
            ext = parse_extreme_ogs(out)
            need_fix = ("Failed to initialize any reasonable values" in out) and bool(ext)
            if need_fix and round_idx < max_rounds:
                ogs = {og for og,_ in ext}
                logger.info(f"[{model_name}] 解析到极端 OG {len(ogs)} 个，进入自动修正下一轮")
                # 从当前输入剔除这些 OG
                h, b = read_family_table(cur_in)
                b2 = remove_ogs_from_family(b, ogs)
                write_family_table(h, b2, cur_in)
                round_idx += 1
                continue
            else:
                if need_fix:
                    logger.error(f"[{model_name}] 已达最大修正轮数（{max_rounds}），停止修正。")
                break

        # ---------- large：只跑一次，不剔除 ----------
        if large_enable and len(lg_rows) > 0:
            lg_dir = ensure_dir(model_dir / "large")
            lg_log = lg_dir / "run.log"
            shutil.copy2(f_large, lg_dir / "family.tsv")
            logger.info(f"[{model_name} LARGE] CAFE5 输出开始")
            cmd = [cafe5_bin, "-i", "family.tsv", "-t", str(tree_path), "-c", str(threads)]
            if extra: cmd += extra.split()
            code, out = run_stream(cmd, cwd=lg_dir, tee_fp=lg_log, logger=logger)
            logger.info(f"[{model_name} LARGE] ---- CAFE5 输出结束 ----")
        else:
            logger.info(f"[{model_name} LARGE] 跳过（开关关闭或 large 家族数=0）")

    # 完成标记
    (run_root / ".cafe.done").write_text(datetime.now().strftime("%F %T")+"\n", encoding="utf-8")
    logger.info("收工，全部模型已完成。")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e)+"\n")
        sys.exit(2)