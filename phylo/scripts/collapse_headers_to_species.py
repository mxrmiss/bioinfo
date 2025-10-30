#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把规范化后的 MSA 表头从 >Species|SeqID 收敛为 >Species；非表头行去空白。
幂等、原子写、失败即报错。
"""
from pathlib import Path
import argparse, sys, os

def log(msg): sys.stderr.write(msg.rstrip()+"\n")

def atomic_write_text(path: Path, it):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as w:
        for s in it: w.write(s)
    os.replace(tmp, path)

def process_one(src: Path, dst: Path):
    n=0
    def gen():
        nonlocal n
        with src.open("r", encoding="utf-8", errors="replace") as r:
            for line in r:
                if line.startswith(">"):
                    s=line.rstrip("\r\n")
                    i=s.find("|")
                    if i!=-1: s=s[:i]
                    yield s+"\n"; n+=1
                else:
                    yield "".join(ch for ch in line if ch not in " \t\r")
    dst.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(dst, gen())
    return n

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    args=ap.parse_args()
    indir, outdir = Path(args.indir), Path(args.outdir)
    if not indir.is_dir(): log(f"[ERR] 输入不存在: {indir}"); return 2
    outdir.mkdir(parents=True, exist_ok=True)
    files=sorted(indir.glob("*.faa"))
    total=nfiles=0
    for f in files:
        o=outdir/f.name
        try:
            cnt=process_one(f,o)
        except Exception as e:
            log(f"[ERR] 处理失败: {f} | {e}"); return 41
        total+=cnt; nfiles+=1
    log(f"[INFO] 收敛完成：文件 {nfiles}，表头 {total} → {outdir}")
    return 0

if __name__=="__main__":
    sys.exit(main())

