#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_mcmctree_run_and_publish.py —— 两阶段调度（verbatim），配置优先，缺配也能独立跑

特性：
  1) 优先读取 config.yaml:mcmctree（仅用于定位路径与日志回显，不改模板内容）。
  2) 无配置文件或无该小节时，使用脚本顶部默认参数也能独立运行。
  3) 进入 work_dir 后，按模板原样直跑：Stage1→复制 out.BV→in.BV→Stage2。
  4) 实时把子进程输出打印到屏幕并写入日志；失败时打包最小排查包。
  5) 不做 BDparas M|C 预检；不做 ESS/QC，这些交给流水线其他脚本。

依赖：
  - Python 标准库 + PyYAML（pip/conda 安装 pyyaml 即可）
  - 外部程序：mcmctree（PAML 4.10.9）
使用：
  - 直接运行，无命令行参数；需要时在脚本顶部修改默认参数。
"""

from __future__ import annotations
import os, sys, shutil, tarfile, time, datetime, subprocess
from pathlib import Path
from typing import Optional

# ===================== 顶部默认参数（无配置文件时生效） =====================
# 默认工作目录（包含 mcmctree.ctl、mcmctree2.ctl、concat.clean.phy、species_calib.nobl.trees）
DEFAULT_WORK_DIR: Path = Path("results/06_cafe/mcmctree")

# 模板文件名（位于 work_dir 内；脚本不会改动）
STAGE1_CTL_NAME: str = "mcmctree.ctl"
STAGE2_CTL_NAME: str = "mcmctree2.ctl"

# 阶段1生成/阶段2需要的 BV 文件名（固定）
OUT_BV_NAME: str = "out.BV"
IN_BV_NAME: str = "in.BV"

# 当 in.BV 已存在时是否覆盖
OVERWRITE_INBV: bool = True

# 是否写 .done 哨兵
WRITE_SENTINELS: bool = True

# 配置文件路径（会优先使用此文件；若不存在则走默认参数）
CONFIG_PATH: Path = Path("config.yaml")

# ===================== 下面内容一般不需要改动 =====================

def ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ts_compact() -> str:
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def which_exec(name: str) -> Optional[str]:
    return shutil.which(name)

def run_streaming(cmd: list[str], cwd: Path, log_file: Path, prefix: str) -> int:
    """实时流式运行外部命令：同时写屏幕与日志。"""
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [CMD] (cwd={cwd}) $ {' '.join(cmd)}\n")
        lf.flush()

        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            print(f"[{prefix}] {line}")
            lf.write(line + "\n")
            lf.flush()
        proc.wait()
        lf.write(f"[{ts()}] [EXIT] code={proc.returncode}\n")
        lf.flush()
        return int(proc.returncode)

def summarize_file(p: Path) -> str:
    if not p.exists():
        return f"{p} [MISSING]"
    size = p.stat().st_size
    mtime = datetime.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    return f"{p} [size={size} bytes, mtime={mtime}]"

def copy_outBV_to_inBV(work_dir: Path, overwrite: bool, logf: Path) -> None:
    src = work_dir / OUT_BV_NAME
    dst = work_dir / IN_BV_NAME
    with open(logf, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INFO] Preparing BV files: {summarize_file(src)}\n")
    if not src.exists():
        raise FileNotFoundError(f"未发现 {src}，Stage1 可能失败，无法继续 Stage2。")
    if dst.exists() and not overwrite:
        with open(logf, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [INFO] {dst} 已存在，按策略不覆盖。\n")
        print(f"[middle] {dst} 已存在（保留原样）。")
        return
    shutil.copy2(src, dst)
    with open(logf, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INFO] Copied {src.name} -> {dst.name}\n")
    print(f"[middle] 已复制 {src.name} -> {dst.name}")

def pack_debug_bundle(work_dir: Path, debug_dir: Path, tb_prefix: str, extra_files: list[Path], log_path: Path) -> Path:
    ensure_dir(debug_dir)
    tar_path = debug_dir / f"{tb_prefix}_{ts_compact()}.tar.gz"
    candidates = [
        work_dir / STAGE1_CTL_NAME,
        work_dir / STAGE2_CTL_NAME,
        work_dir / "out.txt",
        work_dir / "mcmc.txt",
        work_dir / OUT_BV_NAME,
        work_dir / IN_BV_NAME,
        log_path,
    ] + extra_files
    with tarfile.open(tar_path, "w:gz") as tf:
        for f in candidates:
            if f.exists():
                try:
                    tf.add(f, arcname=f.relative_to(work_dir.parent) if work_dir in f.parents else f.name)
                except Exception as e:
                    print(f"[debug] 打包 {f} 出错：{e}")
    return tar_path

def write_sentinel(work_dir: Path, name: str) -> None:
    if not WRITE_SENTINELS:
        return
    sentinel = work_dir / f".{name}.done"
    with open(sentinel, "w", encoding="utf-8") as f:
        f.write(f"{name} done at {ts()}\n")

def load_config(path: Path) -> Optional[dict]:
    """仅使用 PyYAML 读取 config.yaml；不存在时返回 None。"""
    if not path.exists():
        return None
    try:
        import yaml  # PyYAML
    except Exception as e:
        print(f"[warn] 未安装 PyYAML（pyyaml），将使用默认参数。详细：{e}")
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}
        return cfg
    except Exception as e:
        print(f"[warn] 读取 {path} 失败，将使用默认参数。详细：{e}")
        return None

def resolve_workdir_from_config(cfg: dict) -> Optional[Path]:
    """从 config.yaml:mcmctree.work_dir 解析工作目录；不存在则返回 None。"""
    mc = cfg.get("mcmctree") if isinstance(cfg, dict) else None
    wd = mc.get("work_dir") if isinstance(mc, dict) else None
    if isinstance(wd, str) and wd.strip():
        return Path(wd.strip())
    return None

def main() -> int:
    # 0) 读取配置（若存在，以配置为主；否则使用默认）
    cfg = load_config(CONFIG_PATH)
    if cfg is None:
        work_dir = DEFAULT_WORK_DIR
        print(f"[init] 未发现或未能解析 config.yaml，使用默认工作目录：{work_dir}")
        cfg_seq = "concat.clean.phy"     # 仅用于日志回显
        cfg_tree = "species_calib.nobl.trees"
    else:
        wd = resolve_workdir_from_config(cfg)
        work_dir = wd if wd is not None else DEFAULT_WORK_DIR
        mc = cfg.get("mcmctree", {}) if isinstance(cfg, dict) else {}
        cfg_seq = mc.get("seqfile", "concat.clean.phy")
        cfg_tree = mc.get("treefile", "species_calib.nobl.trees")
        print(f"[init] 使用配置文件：{CONFIG_PATH}")
        print(f"[init] mcmctree.work_dir = {work_dir}")
        print(f"[init] (仅回显) seqfile = {cfg_seq}, treefile = {cfg_tree}")

    # 1) 组织日志/调试目录（基于最终 work_dir）
    LOG_DIR = work_dir / "logs" / "10_mcmctree"
    DEBUG_DIR = work_dir / "_debug"
    DEBUG_TARBALL_PREFIX = "mcmctree_debug"

    ensure_dir(work_dir)
    ensure_dir(LOG_DIR)
    ensure_dir(DEBUG_DIR)
    log_path = LOG_DIR / f"run_{ts_compact()}.log"

    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INIT] work_dir={work_dir}\n")

    # 2) 可执行文件存在性
    mctree = which_exec("mcmctree")
    if not mctree:
        msg = "未在 PATH 中找到 mcmctree，可执行文件缺失。"
        print(f"[ERR] {msg}")
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [ERR] {msg}\n")
        return 2

    # 3) 事实性存在性检查（不改模板；不做 M|C 校验）
    required = [
        work_dir / STAGE1_CTL_NAME,
        work_dir / STAGE2_CTL_NAME,
        work_dir / "concat.clean.phy",          # 以模板默认名检查
        work_dir / "species_calib.nobl.trees",
    ]
    # 若配置文件给出了自定义文件名，仅用于日志回显（不强改模板）
    if cfg is not None:
        if cfg_seq and cfg_seq != "concat.clean.phy":
            print(f"[note] 配置中的 seqfile = {cfg_seq}（模板仍以自身设置为准）")
        if cfg_tree and cfg_tree != "species_calib.nobl.trees":
            print(f"[note] 配置中的 treefile = {cfg_tree}（模板仍以自身设置为准）")

    missing = [p for p in required if not p.exists()]
    if missing:
        print("[ERR] 下列必需文件缺失：")
        for p in missing:
            print(f"      - {p}")
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [ERR] Missing required files:\n")
            for p in missing: lf.write(f"  - {p}\n")
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] 已生成排查包：{tb}")
        return 3

    # 4) Stage1
    print(f"[Stage1] 开始：{STAGE1_CTL_NAME}")
    t0 = time.time()
    code1 = run_streaming([mctree, STAGE1_CTL_NAME], work_dir, log_path, prefix="S1")
    t1 = time.time()
    print(f"[Stage1] 结束，退出码={code1}，耗时={t1 - t0:.1f}s")
    if code1 != 0:
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] Stage1 失败，已生成排查包：{tb}")
        return code1

    # 5) out.BV -> in.BV
    try:
        copy_outBV_to_inBV(work_dir, OVERWRITE_INBV, logf=log_path)
    except Exception as e:
        print(f"[ERR] 中间步骤失败：{e}")
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] 已生成排查包：{tb}")
        return 4

    # 6) Stage2
    print(f"[Stage2] 开始：{STAGE2_CTL_NAME}")
    t2 = time.time()
    code2 = run_streaming([mctree, STAGE2_CTL_NAME], work_dir, log_path, prefix="S2")
    t3 = time.time()
    print(f"[Stage2] 结束，退出码={code2}，耗时={t3 - t2:.1f}s")
    if code2 != 0:
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] Stage2 失败，已生成排查包：{tb}")
        return code2

    # 7) 成功收尾
    print("[OK] 两阶段运行完成。产物摘要：")
    for f in ["out.txt", "mcmc.txt", OUT_BV_NAME, IN_BV_NAME]:
        print("     - " + summarize_file(work_dir / f))
    if WRITE_SENTINELS:
        write_sentinel(work_dir, "mcmctree_stage1")
        write_sentinel(work_dir, "mcmctree_stage2")
    print(f"[OK] 日志文件：{log_path}")
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断。")
        sys.exit(130)

