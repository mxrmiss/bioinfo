#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
表头规范化脚本（增强稳健版）：
- 自动从 --proteome-dir 推断三套物种名集合（原始/点转下划线/.protein转_protein）
- 输入：--indir（裁剪后的 MSA 目录，*.trim.faa）
- 输出：--outdir（规范化后目录，>Species|SeqID）
- 匹配规则：按“最长前缀优先”做字面量前缀匹配，避免动态正则的环境差异
- 兜底与日志：未匹配 → _protein 优先，否则取首个下划线前；记录 WARN，不崩溃
"""

import argparse
import sys
import os
from pathlib import Path

def parse_args():
    ap = argparse.ArgumentParser(description="规范化 MSA 表头到 >Species|SeqID（自动推断物种集合）")
    ap.add_argument("--indir", required=True, help="输入目录（如 results/trim，包含 *.trim.faa）")
    ap.add_argument("--outdir", required=True, help="输出目录（如 results/trim_norm）")
    ap.add_argument("--proteome-dir", required=True, help="蛋白质目录（如 data/proteomes 或过滤后目录）")
    ap.add_argument("--workers", type=int, default=1, help="并行度（可选，默认1；建议用 Snakemake 控制并发）")
    return ap.parse_args()

def log(msg):
    sys.stderr.write(msg.rstrip() + "\n")

def stem_without_ext(p: Path) -> str:
    """去掉常见扩展名（.gz/.fasta/.faa/.fa），返回文件基名。"""
    name = p.name
    for suf in (".gz", ".fasta", ".faa", ".fa"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name

def build_species_sets_from_proteomes(prot_dir: Path):
    """
    从蛋白目录收集三套物种集合：
      1) 原始：去扩展名
      2) 变体1：将 . 替换为 _
      3) 变体2：将 .protein 替换为 _protein
    返回：keys（按长度降序去重合并后的列表），empty_flag
    """
    if not prot_dir.is_dir():
        return [], True
    raw = set()
    alt = set()
    pro = set()
    for p in sorted(prot_dir.iterdir()):
        if not p.is_file():
            continue
        s = stem_without_ext(p)
        raw.add(s)
        alt.add(s.replace(".", "_"))
        pro.add(s.replace(".protein", "_protein"))
    merged = set()
    merged.update(raw)
    merged.update(alt)
    merged.update(pro)
    keys = sorted(merged, key=lambda x: (-len(x), x))  # 最长前缀优先
    empty_flag = (len(keys) == 0)
    return keys, empty_flag

def parse_seqid(raw_header: str, sid: str) -> str:
    """
    提取 SeqID：
    - 若含 'protein_'：取其后连续 token（去掉空格）
    - 否则：取第一个 token；若以 sid_ 开头，剥掉该前缀
    """
    line = raw_header.rstrip("\r\n")
    if "protein_" in line:
        tmp = line.split("protein_", 1)[1]
        seqid = tmp.split()[0]
    else:
        t = line[1:] if line.startswith(">") else line
        tok = t.split()[0]
        seqid = tok
        prefix = sid + "_"
        if seqid.startswith(prefix):
            seqid = seqid[len(prefix):]
    return seqid

def normalize_header(raw_header: str, keys, empty_flag: bool):
    """
    返回 (sid, seqid, warn_msg)
    - sid：物种名
    - seqid：序列 ID
    - warn_msg：若发生兜底，返回告警文本；否则为 None
    """
    line = raw_header.rstrip("\r\n")
    t = line[1:] if line.startswith(">") else line
    sid = ""

    if not empty_flag:
        for k in keys:
            if t.startswith(k):
                sid = k
                break

    warn_msg = None
    if not sid:
        if "_protein_" in t:
            sid = t.split("_protein", 1)[0] + "_protein"
        else:
            sid = t.split("_", 1)[0]
        warn_msg = f"[WARN] 未匹配到已知物种，用兜底物种: {sid} | {line}"

    seqid = parse_seqid(line, sid)
    return sid, seqid, warn_msg

def atomic_write_text(path: Path, content_iter):
    """原子写入：先 .tmp 再 rename。"""
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as w:
        for s in content_iter:
            w.write(s)
    os.replace(tmp, path)

def process_file(in_path: Path, out_path: Path, keys, empty_flag: bool):
    """处理单个 *.trim.faa 文件，返回 (n_seq, n_warn)。"""
    n_seq = 0
    n_warn = 0

    def gen():
        nonlocal n_seq, n_warn
        with in_path.open("r", encoding="utf-8", errors="replace") as r:
            for raw in r:
                if raw.startswith(">"):
                    sid, seqid, warn = normalize_header(raw, keys, empty_flag)
                    if warn:
                        log(warn)
                        n_warn += 1
                    yield f">{sid}|{seqid}\n"
                    n_seq += 1
                else:
                    # 删除空白字符（空格、制表符、回车等）
                    yield "".join(ch for ch in raw if ch not in " \t\r")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(out_path, gen())
    return n_seq, n_warn

def main():
    args = parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir)
    pdir = Path(args.proteome_dir)

    if not indir.is_dir():
        log(f"[ERR] 输入目录不存在：{indir}")
        return 2
    outdir.mkdir(parents=True, exist_ok=True)

    keys, empty_flag = build_species_sets_from_proteomes(pdir)
    if empty_flag:
        log(f"[WARN] 未从 {pdir} 推断到任何物种名，将完全依赖兜底策略。")

    files = sorted(indir.glob("*.trim.faa"))
    if not files:
        log(f"[WARN] 未找到 *.trim.faa：{indir}")

    total_seq = total_warn = nfiles = 0
    for f in files:
        out = outdir / f.name
        try:
            n_seq, n_warn = process_file(f, out, keys, empty_flag)
        except Exception as e:
            log(f"[ERR] 处理失败：{f} | {e}")
            return 41
        total_seq += n_seq
        total_warn += n_warn
        nfiles += 1

    log(f"[INFO] 规范化完成：文件 {nfiles}，序列表头 {total_seq}，兜底警告 {total_warn} → {outdir}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

