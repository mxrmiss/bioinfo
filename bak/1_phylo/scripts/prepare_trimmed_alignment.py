#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
遍历 sco.list 按顺序裁剪 MSA：
- 从 RESULTS_PATH.txt → MultipleSequenceAlignments 目录
- 对每个 OG 选择 og.fa / og.fasta
- 运行 trimal -automated1 输出到 outdir/OG.trim.faa
- 缺失或裁剪为空 → 生成空占位文件（保持 DAG 连续）
- 原子写、失败即报错；统计有效/总计
"""
import argparse, sys, subprocess, os
from pathlib import Path

def log(m): sys.stderr.write(m.rstrip()+"\n")

def read_results_path(path_file: Path) -> Path:
    return Path(path_file.read_text(encoding="utf-8").strip())

def load_ogs(list_file: Path):
    return [x.strip() for x in list_file.read_text(encoding="utf-8").splitlines() if x.strip()]

def run_trimal(trimal, src, dst_tmp):
    cmd=[trimal, "-in", str(src), "-out", str(dst_tmp), "-automated1"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--results-file", required=True)
    ap.add_argument("--og-list", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--trimal", default="trimal")
    args=ap.parse_args()

    resdir = read_results_path(Path(args.results_file))
    msadir = resdir/"MultipleSequenceAlignments"
    if not msadir.is_dir(): log(f"[ERR] 未找到 MSA 目录: {msadir}"); return 2

    ogs = load_ogs(Path(args.og_list))
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    ok=0; total=0
    for og in ogs:
        total+=1
        src = msadir/f"{og}.fa"
        if not src.exists(): src = msadir/f"{og}.fasta"
        dst = outdir/f"{og}.trim.faa"
        if not src.exists():
            log(f"[WARN] 缺少 MSA：{og} —— 生成空占位")
            dst.write_text("", encoding="utf-8"); continue
        try:
            tmp = dst.with_suffix(dst.suffix+".tmp")
            run_trimal(args.trimal, src, tmp)
            # 若 trimal 产物为空，保留空占位
            data = tmp.read_text(encoding="utf-8") if tmp.exists() else ""
            if data.strip():
                os.replace(tmp, dst); ok+=1
            else:
                if tmp.exists(): os.remove(tmp)
                dst.write_text("", encoding="utf-8")
                log(f"[WARN] 裁剪后为空：{og} —— 保留空占位")
        except subprocess.CalledProcessError as e:
            log(f"[WARN] trimal 失败：{og} —— 保留空占位 | {e.stderr.decode(errors='ignore')[:200]}")
            if tmp.exists(): os.remove(tmp)
            dst.write_text("", encoding="utf-8")
            continue

    if ok==0:
        log("[ERR] 无任何有效 MSA，流程中止"); return 4
    log(f"[INFO] 裁剪完成：有效 {ok} / 总计 {total}")
    return 0

if __name__=="__main__":
    sys.exit(main())

