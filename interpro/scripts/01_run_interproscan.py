#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_run_interproscan.py
================================================================================
功能：
  - 读取 config.yaml
  - 批量运行 InterProScan（TSV 输出）
  - 支持 -dp / -goterms / -appl / -cpu / -T tmp 等参数
  - 流式输出：InterProScan 运行过程实时滚动到屏幕
  - 日志系统：每个输入一个独立日志文件，同时写一个主日志与 run_manifest

结果目录（工程化）：
  - OUTDIR: results/01_interproscan/{ver}.core/
  - 期望 TSV: OUTDIR/{basename(input)}.tsv
================================================================================
"""

import os
import sys
import csv
import shlex
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

CONFIG_YAML = "config.yaml"


def die(msg: str, code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def load_yaml_config(path: str) -> Dict[str, Any]:
    p = Path(path)
    if not p.is_file():
        die(f"找不到配置文件：{path}")

    try:
        import yaml  # type: ignore
    except Exception as e:
        die(
            "缺少依赖 PyYAML，无法读取 config.yaml。\n"
            "解决方式：在你的 conda 环境中安装 PyYAML。\n"
            f"原始错误：{e}"
        )

    with p.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        die("config.yaml 解析失败：顶层必须是字典")
    return cfg


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def setup_logger(log_file: Path, level: str = "INFO"):
    import logging

    level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    lv = level_map.get(level.upper(), logging.INFO)

    logger = logging.getLogger(str(log_file))
    logger.setLevel(lv)
    logger.propagate = False

    if logger.handlers:
        return logger

    fmt = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(lv)
    sh.setFormatter(fmt)

    ensure_dir(log_file.parent)
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(lv)
    fh.setFormatter(fmt)

    logger.addHandler(sh)
    logger.addHandler(fh)
    return logger


def fmt_tpl(s: str, mapping: Dict[str, Any]) -> str:
    try:
        return s.format_map(mapping)
    except Exception:
        return s


def resolve_inputs(cfg: Dict[str, Any]) -> List[Path]:
    sec = cfg.get("inputs", {}) or {}
    inputs = sec.get("inputs", None)

    out: List[Path] = []
    if isinstance(inputs, list) and inputs:
        for x in inputs:
            out.append(Path(str(x)))
        return out

    inputs_dir = Path(str(sec.get("inputs_dir", "data/query")))
    ext = str(sec.get("inputs_ext", ".faa"))

    if not inputs_dir.is_dir():
        die(f"inputs_dir 不存在或不是目录：{inputs_dir}")

    for p in sorted(inputs_dir.iterdir()):
        if p.is_file() and p.name.endswith(ext):
            out.append(p)

    if not out:
        die(f"在 {inputs_dir} 下没有找到 *{ext} 输入文件")
    return out


def basename_keep_ext(p: Path) -> str:
    return p.name


def run_command_stream(
    cmd: List[str],
    logger,
    log_prefix: str,
    cwd: Optional[Path] = None,
) -> int:
    logger.info(f"{log_prefix} CMD = {shlex.join(cmd)}")

    p = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    assert p.stdout is not None
    for line in p.stdout:
        logger.info(f"{log_prefix} {line.rstrip()}")

    p.wait()
    logger.info(f"{log_prefix} Exit code = {p.returncode}")
    return int(p.returncode)


def write_manifest_row(manifest_tsv: Path, row: Dict[str, str]) -> None:
    ensure_dir(manifest_tsv.parent)
    file_exists = manifest_tsv.is_file()

    cols = [
        "timestamp_start",
        "timestamp_end",
        "ver",
        "input_faa",
        "outdir",
        "expected_tsv",
        "cpu",
        "dp",
        "goterms",
        "appl",
        "cmd",
        "exit_code",
    ]

    with manifest_tsv.open("a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        if not file_exists:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in cols})


def main() -> None:
    cfg = load_yaml_config(CONFIG_YAML)

    log_cfg = cfg.get("logging", {}) or {}
    log_dir = Path(str(log_cfg.get("log_dir", "logs")))
    log_level = str(log_cfg.get("level", "INFO"))
    print_command = bool(log_cfg.get("print_command", True))

    ips_cfg = cfg.get("interproscan", {}) or {}
    ver = str(ips_cfg.get("ver", "5.76-107.0")).strip()
    cpu = int(ips_cfg.get("cpu", 20))
    tmp_dir = str(ips_cfg.get("tmp_dir", "tmp"))
    outdir_tpl = str(ips_cfg.get("outdir_tpl", "results/01_interproscan/{ver}.core"))
    ips_tpl = str(ips_cfg.get("ips_sh_tpl", "data/interproscan-{ver}/interproscan.sh"))
    fmt = str(ips_cfg.get("format", "TSV")).strip()

    dp = bool(ips_cfg.get("dp", True))
    goterms = bool(ips_cfg.get("goterms", False))

    appl = ips_cfg.get("appl", [])
    if not isinstance(appl, list) or not appl:
        die("config.yaml: interproscan.appl 必须是非空列表")

    mapping = {"ver": ver}
    outdir = Path(fmt_tpl(outdir_tpl, mapping))
    ips_sh = Path(fmt_tpl(ips_tpl, mapping))

    ensure_dir(outdir)

    main_log = log_dir / "01_run_interproscan" / ver / "main.log"
    logger = setup_logger(main_log, level=log_level)

    logger.info(f"[01] Start at {now_iso()}")
    logger.info(f"[01] VER={ver}")
    logger.info(f"[01] IPS_SH={ips_sh}")
    logger.info(f"[01] OUTDIR={outdir}")
    logger.info(f"[01] CPU={cpu}")
    logger.info(f"[01] DP={dp}")
    logger.info(f"[01] GOTERMS={goterms}")
    logger.info(f"[01] APPL={','.join([str(x) for x in appl])}")

    if not ips_sh.is_file():
        die(f"InterProScan 脚本不存在：{ips_sh}")

    inputs = resolve_inputs(cfg)
    logger.info(f"[01] N_inputs={len(inputs)}")

    manifest = outdir / "run_manifest.tsv"

    for i, in_faa in enumerate(inputs, start=1):
        if not in_faa.is_file() or in_faa.stat().st_size == 0:
            die(f"输入文件不存在或为空：{in_faa}")

        base = basename_keep_ext(in_faa)
        expected_tsv = outdir / f"{base}.tsv"

        per_log = log_dir / "01_run_interproscan" / ver / f"{base}.log"
        per_logger = setup_logger(per_log, level=log_level)

        per_logger.info(f"[01] ({i}/{len(inputs)}) INPUT={in_faa}")
        per_logger.info(f"[01] EXPECTED_TSV={expected_tsv}")

        cmd = [
            str(ips_sh),
            "-i",
            str(in_faa),
            "-f",
            fmt,
            "-d",
            str(outdir),
            "-T",
            tmp_dir,
            "-cpu",
            str(cpu),
            "-appl",
            ",".join([str(x) for x in appl]),
        ]
        if dp:
            cmd.append("-dp")
        if goterms:
            cmd.append("-goterms")

        if print_command:
            per_logger.info(f"[01] CMD_PREVIEW={shlex.join(cmd)}")

        t_start = now_iso()
        rc = run_command_stream(cmd, per_logger, log_prefix="[IPS]")
        t_end = now_iso()

        write_manifest_row(
            manifest,
            {
                "timestamp_start": t_start,
                "timestamp_end": t_end,
                "ver": ver,
                "input_faa": str(in_faa),
                "outdir": str(outdir),
                "expected_tsv": str(expected_tsv),
                "cpu": str(cpu),
                "dp": str(dp),
                "goterms": str(goterms),
                "appl": ",".join([str(x) for x in appl]),
                "cmd": shlex.join(cmd),
                "exit_code": str(rc),
            },
        )

        if rc != 0:
            per_logger.error(f"[01] InterProScan 运行失败：exit_code={rc}")
            per_logger.error(f"[01] 请查看日志：{per_log}")
            continue

        if not expected_tsv.is_file() or expected_tsv.stat().st_size == 0:
            per_logger.warning("[01] 未检测到期望 TSV 或 TSV 为空。")
            per_logger.warning(f"[01] expected_tsv={expected_tsv}")
        else:
            per_logger.info("[01] OK: TSV 已生成")
            per_logger.info(f"[01] TSV={expected_tsv}")

    logger.info(f"[01] Done at {now_iso()}")
    logger.info(f"[01] Manifest={manifest}")


if __name__ == "__main__":
    main()

