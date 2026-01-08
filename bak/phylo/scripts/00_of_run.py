#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_of_run.py —— 封装 OrthoFinder 运行（零猜测 + 可复现锚点）

修复要点（本次更新）：
1) 不再提前创建 -o 目标目录（results/orthofinder），OrthoFinder 3.x 要求该目录“不能预先存在”；
   我们只在运行前删除旧目录（若存在），但不再 os.makedirs()。
2) 运行成功后写入 RESULTS_DIR.txt（唯一锚点），并 touch .OF_DONE（便于脚本独立运行）。

使用方式：
- 直接在仓库根目录运行：python scripts/00_of_run.py
- 或交由 Snakefile 的 rule: orthofinder_run 调用（该 rule 也会 touch 哨兵，不冲突）
"""

import os
import sys
import glob
import yaml
import shutil
import subprocess

def which(x):
    # 兼容低版本 Python，不用“|”类型注解
    from shutil import which as _w
    return _w(x)

def load_config(path="config.yaml"):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def main():
    cfg = load_config()
    P   = cfg["paths"]
    IN  = cfg["input"]
    OF  = cfg["orthofinder"]

    proteome_dir = IN["proteome_dir"]                # 蛋白目录（FASTA，后缀 .faa/.fa/.fasta）
    out_root     = P["orthofinder_results_dir"]      # 固定输出根目录（例如 results/orthofinder）
    sentinel     = os.path.join(out_root, ".OF_DONE")

    # 1) 若存在旧目录，先删除；注意：不要提前创建 out_root，让 OrthoFinder 自己建
    if os.path.exists(out_root):
        shutil.rmtree(out_root)

    # 2) 组装 OrthoFinder 命令（线程/搜索引擎/树构建引擎等全部来自 config.yaml）
    threads  = int(OF.get("threads", 16))
    athreads = int(OF.get("analysis_threads", max(1, threads // 2)))
    search   = str(OF.get("search_engine", "") or "").strip()
    msa_and  = bool(OF.get("msa_and_trees", True))
    gengine  = OF.get("gene_tree_engine", None)

    cmd = [
        "orthofinder",
        "-f", proteome_dir,
        "-t", str(threads),
        "-a", str(athreads),
        "-o", out_root
    ]

    # 搜索引擎：diamond/mmmseqs/blast；若指定 diamond 但未安装，则交给 OrthoFinder 自选
    if search:
        if search.lower() == "diamond":
            if which("diamond"):
                cmd += ["-S", "diamond"]
            else:
                print("[WARN] config 指定 diamond 但 PATH 未找到，将交由 OrthoFinder 默认搜索引擎", file=sys.stderr)
        else:
            cmd += ["-S", search]

    # -M msa（让 OrthoFinder 对单基因做 MSA，用于后续修剪/拼接）
    if msa_and:
        cmd += ["-M", "msa"]

    # -T <engine>（可选 gene tree 引擎）
    if gengine:
        cmd += ["-T", str(gengine)]

    # 3) 运行
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"

    print("[CMD]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)

    # 4) 解析唯一的 Results_* 子目录，写锚点
    rs = sorted(glob.glob(os.path.join(out_root, "Results_*")))
    if not rs:
        print("[ERR] 运行完成但未发现 Results_* 目录，请检查 OrthoFinder 日志", file=sys.stderr)
        sys.exit(50)
    results_dir = rs[0]
    with open(os.path.join(out_root, "RESULTS_DIR.txt"), "w", encoding="utf-8") as f:
        f.write(results_dir + "\n")

    # 5) 写哨兵，便于脚本单独运行时也能被下游识别
    try:
        with open(sentinel, "w") as f:
            f.write("done\n")
    except FileNotFoundError:
        # 防御性：正常情况下 out_root 会被 OrthoFinder 创建
        os.makedirs(out_root, exist_ok=True)
        with open(sentinel, "w") as f:
            f.write("done\n")

if __name__ == "__main__":
    main()

