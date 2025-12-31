#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_iqtree_concat.py —— IQ-TREE 构树（始终注入 -p）
改动要点（本版相对上一版）：
1) 只要 supermatrix_dir 下存在有效的 partitions.txt，就转换为 partitions.raxml，
   并且 **无条件** 在 IQ-TREE 命令中注入 `-p partitions.raxml`。
2) 不再检查 config.yaml 的 iqtree_params.extra 中是否已有 -p/-spp；即使有也照样注入。
3) 若未检测到 partitions.txt 或没有有效分区，自动回落为单分区构树（不报错）。

说明：
- 我们内部的 partitions.txt 为三列：OG  start  end（1-based, 闭区间），IQ-TREE 不直接识别，
  必须转换为 RAxML 风格：`AA, OG = start-end`，故生成 partitions.raxml 并用 -p 传入。
- 所有参数仍集中在 config.yaml（线程、模型、bootstrap、alrt、extra 等），本脚本不读命令行参数。
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"

def load_yaml(p):
    import yaml
    with open(p, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def get_config():
    cfg = os.environ.get('PHYLO_CONFIG', DEFAULT_CONFIG)
    if not os.path.exists(cfg):
        print(f"[ERR] 未找到配置文件：{cfg}", file=sys.stderr); sys.exit(2)
    return load_yaml(cfg)

def raxml_from_txt(txt_path: Path, raxml_path: Path, datatype: str = "AA") -> int:
    """
    将三列表 partitions.txt 转为 RAxML 风格：
      输入行：OG <tab/space> start <tab/space> end
      输出行：AA, OG = start-end
    返回写出的有效分区数。
    """
    n = 0
    with open(txt_path, "r", encoding="utf-8") as f, \
         open(raxml_path, "w", encoding="utf-8") as g:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # 兼容任意空白/逗号分隔
            parts = line.replace(",", " ").split()
            if len(parts) < 3:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
            og = parts[0]
            try:
                a = int(parts[1]); b = int(parts[2])
            except ValueError:
                continue
            if a <= 0 or b < a:
                continue
            g.write(f"{datatype}, {og} = {a}-{b}\n")
            n += 1
    return n

def main():
    cfg = get_config()
    P = cfg["paths"]
    C = cfg["species_tree"]["concat"]
    bins = cfg.get("binaries", {})

    sup_dir = Path(P["supermatrix_dir"])
    tree_dir = Path(P["trees_dir"])
    tree_dir.mkdir(parents=True, exist_ok=True)

    sup = sup_dir / "supermatrix.faa"
    if not sup.is_file():
        print(f"[ERR] 缺少 supermatrix：{sup}", file=sys.stderr); sys.exit(3)

    # IQ-TREE 可执行
    exe = bins.get("iqtree")
    if not exe:
        exe = shutil.which("iqtree2") or shutil.which("iqtree")
    if not exe:
        print("[ERR] PATH 未找到 iqtree/iqtree2，且 binaries.iqtree 为空", file=sys.stderr); sys.exit(4)

    # 读取参数（均来自 config.yaml）
    thr = int(C.get("iqtree_threads", 8))
    prm = C.get("iqtree_params", {})
    model = str(prm.get("model", "MFP+MERGE"))
    boot = int(prm.get("bootstrap", 1000))
    alrt = int(prm.get("alrt", 1000))
    extra_raw = str(prm.get("extra", "")).strip()
    extra = extra_raw.split() if extra_raw else []

    # --- 分区文件自动处理：只要有 partitions.txt 且有效，就无条件注入 -p ---
    txt = sup_dir / "partitions.txt"
    raxml = sup_dir / "partitions.raxml"
    injected_q = False

    if txt.is_file():
        try:
            n = raxml_from_txt(txt, raxml, datatype="AA")  # 蛋白数据 → AA
            if n > 0:
                # 无条件注入 -p，不检查 extra 中是否已有 -q/-spp
                extra = ["-p", str(raxml)] + extra
                injected_q = True
                print(f"[INFO] partitions: 从 {txt.name} 转换 {n} 个分区 → {raxml.name}；已注入 -p")
            else:
                print(f"[WARN] {txt.name} 中无有效分区，按单分区处理。")
        except Exception as e:
            print(f"[WARN] 分区文件转换失败（继续单分区）：{e}")
    else:
        print("[INFO] 未检测到 partitions.txt，将作为单分区处理。")

    # --- 组装并运行 IQ-TREE ---
    pre = tree_dir / "species_tree"
    cmd = [
        exe, "-s", str(sup), "-pre", str(pre),
        "-m", model, "-B", str(boot), "-alrt", str(alrt), "-T", str(thr),
        *extra
    ]

    print("[CMD]", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print("[ERR] IQ-TREE 运行失败", file=sys.stderr); sys.exit(r.returncode)

    # 收尾：复制并打点
    treefile = pre.with_suffix(".treefile")
    nwk = tree_dir / "species_tree.nwk"
    if not treefile.is_file():
        print(f"[ERR] 未生成 treefile：{treefile}", file=sys.stderr); sys.exit(5)

    shutil.copyfile(treefile, nwk)
    Path(tree_dir / ".iqtree.done").touch()
    print(f"[DONE] 物种树完成：{nwk}（{'分区已启用' if injected_q else '单分区'}）")

if __name__ == "__main__":
    main()

