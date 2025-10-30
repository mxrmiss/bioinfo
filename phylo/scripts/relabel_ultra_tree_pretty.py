#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
等时树标签美化（仅改 tip 文本，不改树结构与分支长度）：
- 修复：不把 ':' 当作标签的一部分；确保分支长度原样保留；
- 统一物种名为 Genus_species（把 Genus.species → Genus_species）；
- 若出现重复名，自动追加 .2/.3/... 保证唯一；
- 最后做一次 sanity check：每个 ':' 后应是数字。
"""

import re
import argparse

def clean_tip(raw: str) -> str:
    s = raw.strip().strip("'\"")
    # 去常见后缀/字段
    s = re.sub(r'(_|\.)protein.*$', '', s)
    s = re.sub(r'(_|\.)evm_model.*$', '', s)
    s = re.sub(r'(_|\.)g\d+_t\d+.*$', '', s)
    s = re.sub(r'(_|\.)XP?_[0-9\.]+.*$', '', s)
    s = re.sub(r'(_|\.)GN=\S+.*$', '', s)
    s = re.sub(r'\|.*$', '', s)  # Species|SeqID -> Species
    # Genus_species 或 "Genus species" 或 Genus.species → 统一成 Genus_species
    if "_" in s:
        parts = s.split("_")
        if len(parts) >= 2:
            s2 = f"{parts[0]}_{parts[1]}"
        else:
            s2 = s.replace(".", "_").replace(" ", "_")
    elif " " in s:
        sp = s.split()
        s2 = f"{sp[0]}_{sp[1]}" if len(sp) >= 2 else s.replace(".", "_").replace(" ", "_")
    else:
        s2 = s.replace(".", "_")
    # 仅保留字母/数字/下划线/连字符/点（点仅用于去重后缀 .2 等）
    s2 = re.sub(r'[^A-Za-z0-9_\.-]', '', s2)
    return s2

def relabel_newick(text: str) -> str:
    # 匹配 '(' 或 ',' 后的 tip 名称；后面必须紧跟 ':' 和数字（保持分支长度原样）
    pat = re.compile(r'([\(,])\s*([A-Za-z0-9._\-|]+)\s*(?=:\s*[0-9eE\+\-\.])')

    used = {}
    def _repl(m):
        prefix, name = m.group(1), m.group(2)
        new = clean_tip(name)
        # 确保唯一性（重复时追加 .2/.3/...）
        base, i = new, 2
        while new in used and used[new] != name:
            new = f"{base}.{i}"; i += 1
        used[new] = name
        return f"{prefix}{new}"

    new_text = pat.sub(_repl, text)

    # sanity check：每个 ':' 后都应是数字/科学计数法字符
    if re.search(r':\s*[^0-9eE\+\.-]', new_text):
        raise SystemExit("ERROR: 冒号后检测到非数字字符，疑似破坏了分支长度。")
    return new_text

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", required=True, help="输入 newick")
    ap.add_argument("--out", required=True, help="输出 newick（仅标签变更）")
    args = ap.parse_args()

    with open(args.__dict__["in"], "r") as f:
        text = f.read()
    new_text = relabel_newick(text)
    with open(args.out, "w") as f:
        f.write(new_text)

if __name__ == "__main__":
    main()

