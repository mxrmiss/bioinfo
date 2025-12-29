#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_prepare_trimmed_alignment.py  —— trimal 裁剪单基因蛋白对齐 + 哨兵

逻辑（保持原有主流程）：
  1) 读取 OrthoFinder 单基因 AA 对齐，使用 trimal 裁剪到 paths.msa_trim_dir
  2) 写哨兵：msa_trim_dir/.done
"""

import os
import sys
import subprocess
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"


def load_yaml(p):
    import yaml
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main():
    # 载入配置
    if not os.path.exists(DEFAULT_CONFIG):
        print(f"[ERR] 找不到配置文件 {DEFAULT_CONFIG}", file=sys.stderr)
        sys.exit(2)
    cfg = load_yaml(DEFAULT_CONFIG)

    paths = cfg.get("paths", {})
    bins = cfg.get("binaries", {})
    reports_dir = paths.get("reports_dir", "results/reports")
    of_root = paths.get("orthofinder_results_dir", "results/orthofinder")
    msa_trim_dir = paths.get("msa_trim_dir", "results/msa_trim")
    of_suffix = cfg.get("input", {}).get("of_msa_suffix", ".fa")
    keep_list = os.path.join(reports_dir, "ogs_selected.list")

    # 锚定 OrthoFinder 结果
    results_txt = os.path.join(of_root, "RESULTS_DIR.txt")
    if not os.path.exists(results_txt):
        print(f"[ERR] 缺少 RESULTS_DIR.txt：{results_txt}", file=sys.stderr)
        sys.exit(3)
    with open(results_txt, "r", encoding="utf-8") as f:
        results_dir = f.readline().strip()

    msa_src = os.path.join(results_dir, "MultipleSequenceAlignments")
    if not os.path.isdir(msa_src):
        print(f"[ERR] 找不到 OF 对齐目录：{msa_src}", file=sys.stderr)
        sys.exit(4)

    # 读严格 SCO 列表
    if not os.path.exists(keep_list):
        print(f"[ERR] 缺少严格SCO清单：{keep_list}", file=sys.stderr)
        sys.exit(5)
    with open(keep_list, "r", encoding="utf-8") as f:
        ogs = [x.strip() for x in f if x.strip()]

    os.makedirs(msa_trim_dir, exist_ok=True)

    trimal = bins.get("trimal", "trimal")
    trimmed_cnt = 0

    for og in ogs:
        src = os.path.join(msa_src, og + of_suffix)
        if not os.path.exists(src):
            print(f"[ERR] 缺少 MSA 源文件：{src}", file=sys.stderr)
            sys.exit(6)

        dst = os.path.join(msa_trim_dir, f"{og}.trim.faa")

        # 调 trimal
        r = subprocess.run(
            [trimal, "-in", src, "-out", dst, "-automated1"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if r.returncode != 0 or (not os.path.exists(dst)) or os.path.getsize(dst) == 0:
            sys.stderr.write(r.stderr)
            print(f"[ERR] trimal 失败：{src}", file=sys.stderr)
            sys.exit(7)

        trimmed_cnt += 1

    # === 哨兵文件（在所有循环完成后写入）===
    Path(os.path.join(msa_trim_dir, ".done")).touch()

    print(f"[DONE] 对齐裁剪完成：{trimmed_cnt} 个 → {msa_trim_dir}")


if __name__ == "__main__":
    main()


