#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
08_enrich_module.py — 富集模块包装器（流式日志 + 配置优先）
功能：
  1) 读取项目根目录下的 config.yaml，解析富集相关参数与路径。
  2) 调用 R 脚本 scripts/08_g_enrich.R 完成 GO/KEGG 富集与基因表输出。
  3) 屏幕与日志同时输出（简洁流式），返回码非 0 则报错退出。

约定：
  - 工作目录 = 项目根目录（config.yaml 所在目录）。
  - 所有路径为相对路径（相对于项目根目录）。
  - 二进制路径从 config.yaml -> 本脚本默认，config.yaml 覆盖脚本。
"""

import sys, os, json, subprocess, textwrap
from pathlib import Path
import yaml
import datetime

# ========== 本脚本默认（会被 config.yaml 覆盖） ==========
DEFAULTS = {
    "paths": {
        "anno_dir":   "ref/annotations",
        "deg_dir":    "results/deg",
        "matrix_dir": "results/matrix",
        "enrich_dir": "results/enrich",
        "logs_dir":   "logs",
    },
    "binaries": {
        "rscript": "Rscript",
    }
}

def load_yaml(fp: Path) -> dict:
    with open(fp, "r", encoding="utf-8") as r:
        return yaml.safe_load(r) or {}

def deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out

def main():
    proj = Path.cwd()
    cfg_fp = proj / "config.yaml"
    if not cfg_fp.exists():
        print("[ERR] 缺少配置文件：config.yaml（请在项目根目录运行）", file=sys.stderr)
        sys.exit(1)

    # 载入配置：config.yaml 覆盖脚本默认
    cfg = deep_merge(DEFAULTS, load_yaml(cfg_fp))

    # 解析关键路径
    paths = cfg.get("paths", {})
    enrich_dir = proj / paths.get("enrich_dir", "results/enrich")
    logs_dir   = proj / paths.get("logs_dir",   "logs")
    logs_dir.mkdir(parents=True, exist_ok=True)
    enrich_dir.mkdir(parents=True, exist_ok=True)

    # 组装传入 R 的运行时环境（用环境变量，不污染命令行）
    payload = {
        "paths": {
            "anno_dir":   paths.get("anno_dir",   "ref/annotations"),
            "deg_dir":    paths.get("deg_dir",    "results/deg"),
            "matrix_dir": paths.get("matrix_dir", "results/matrix"),
            "enrich_dir": paths.get("enrich_dir", "results/enrich"),
        },
        "deg": cfg.get("deg", {}),
        "enrich": cfg.get("enrich", {}),
    }
    env = os.environ.copy()
    env["RNASEQ_CONFIG_JSON"] = json.dumps(payload, ensure_ascii=False)

    rscript = cfg.get("binaries", {}).get("rscript", "Rscript")
    r_code  = str(proj / "scripts" / "08_g_enrich.R")

    cmd = [rscript, r_code]

    # 日志文件
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_fp = logs_dir / f"08_g_enrich.{stamp}.log"

    print("CMD:", " ".join(cmd))
    print(f"[INFO] 日志：{log_fp}")

    # 流式执行
    with open(log_fp, "w", encoding="utf-8") as logf:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            cwd=str(proj), env=env, text=True, bufsize=1
        )
        for line in proc.stdout:
            print(line, end="")
            logf.write(line)
        ret = proc.wait()

    if ret != 0:
        print(f"[ERR] 外部命令失败，退出码={ret}", file=sys.stderr)
        sys.exit(ret)
    else:
        print(f"富集分析完成 ✅ 结果目录： {enrich_dir}")

if __name__ == "__main__":
    main()

