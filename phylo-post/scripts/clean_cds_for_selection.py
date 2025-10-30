#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
clean_cds_for_selection.py
严格用于“正选择/密码子分析”之前对 CDS 的清洗。

运行位置：phylo-post 根目录（由 phylo-post 的 Snakefile 调用）
输入目录（只读）：../phylo/data/cds/   （来自上游 phylo）
可选蛋白：       ../phylo/data/proteomes/（做 ID 交集，保证一致）
输出目录：       results/cds_codon/{stem}.fna
临时目录：       results/.tmp_cds_clean/{stem}/  （每物种完成即清理，结束后删根目录）

依赖：seqkit, awk, gzip
参数：都在顶部，可手动修改，不走命令行参数。
"""

import os, sys, re, gzip, shutil, subprocess
from pathlib import Path
from typing import List, Optional

# ==================== 顶部参数（可手动改） ====================

# 上游输入（相对 phylo-post 根目录）
IN_CDS_DIR   = "../phylo/data/cds"
IN_PEP_DIR   = "../phylo/data/proteomes"   # 与蛋白 ID 取交集的数据源
USE_PROT_IDS = True                        # 推荐 True，确保与蛋白一一对应

# 本模块输出与临时
OUT_DIR      = "results/cds_codon"         # 清洗后的 CDS 放这里
TMP_ROOT     = "results/.tmp_cds_clean"    # 临时目录（会被清理）

# 严格阈值（正选择友好，若需要可再收紧）
CDS_MIN_LEN_NT            = 150       # 最短 CDS（nt）
CDS_MAX_N_FRAC            = 0.05      # N 比例上限
CDS_REQUIRE_MULTIPLE_OF_3 = True      # 是否必须为 3 的倍数
CDS_ALLOW_TERMINAL_STOP   = True      # 翻译后允许末尾单个 '*'
CDS_FORBID_INTERNAL_STOP  = True      # 禁止内部 '*'

# 支持的后缀（按顺序匹配，先 .gz）
FA_EXTS  = [".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"]
PEP_EXTS = [".faa.gz", ".faa", ".fa.gz", ".fa", ".fasta.gz", ".fasta"]

# 依赖工具
TOOLS = ["seqkit", "awk", "gzip"]

# ==================== 工具函数 ====================

def need(cmd: str):
    if subprocess.call(f"command -v {cmd} >/dev/null 2>&1", shell=True) != 0:
        print(f"[ERR] 未找到依赖命令：{cmd}", file=sys.stderr)
        sys.exit(127)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def is_gz(p: Path) -> bool:
    return str(p).endswith(".gz")

def gunzip_to(src: Path, dst: Path) -> Path:
    if not is_gz(src):
        return src
    ensure_dir(dst.parent)
    with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    return dst

def find_by_stem(stem: str, directory: Path, exts: List[str]) -> Optional[Path]:
    for e in exts:
        f = directory / f"{stem}{e}"
        if f.is_file():
            return f
    return None

def list_stems_from_dir(d: Path, exts: List[str]) -> List[str]:
    stems = set()
    for f in d.iterdir():
        if not f.is_file(): continue
        for e in exts:
            if f.name.endswith(e):
                stems.add(f.name[:-len(e)])
                break
    return sorted(stems)

def bash_out(cmd: str) -> str:
    return subprocess.getoutput(cmd)

def run(cmd: str):
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        print(f"[ERR] 命令失败：{cmd}", file=sys.stderr)
        sys.exit(res.returncode)

def cleanup_species_tmp(tmpdir: Path):
    if tmpdir.exists():
        shutil.rmtree(tmpdir, ignore_errors=True)

# ==================== 核心处理 ====================

def clean_one(stem: str, cds_in: Path, pep_in: Optional[Path], tmp_root: Path, out_dir: Path):
    tmpdir = tmp_root / stem
    ensure_dir(tmpdir)

    # 解压到临时（若为 .gz）
    cds_plain = gunzip_to(cds_in, tmpdir / (cds_in.name.replace(".gz","")))
    out = out_dir / f"{stem}.fna"

    # 第1阶段：字符合法性 + N% + 3倍数 + 最短长度
    tmp1 = tmpdir / "cds.qc1.fna"
    awk1 = r"""{
      id=$1; len=$2; seq=toupper($3);
      gsub(/U/,"T",seq);
      # 仅允许 IUPAC DNA
      if (seq ~ /[^ACGTRYKMSWBDHVN]/) next;
      # N 比例
      t=seq; gsub(/[^N]/,"",t); ncnt=length(t);
      if (len>0 && ncnt/len > %MAXN%) next;
      # 3 的倍数
      if (%REQ3% && (len%3)!=0) next;
      # 最短长度
      if (len < %MINLEN%) next;
      print ">" id "\n" seq
    }"""
    awk1 = awk1.replace("%MAXN%", str(CDS_MAX_N_FRAC))
    awk1 = awk1.replace("%REQ3%", "1" if CDS_REQUIRE_MULTIPLE_OF_3 else "0")
    awk1 = awk1.replace("%MINLEN%", str(CDS_MIN_LEN_NT))
    run(f"seqkit fx2tab -nl -i -s '{cds_plain}' | awk '{awk1}' | seqkit seq -w 0 > '{tmp1}'")

    # 第2阶段：翻译并剔除内部 *
    tmp_ids = tmpdir / "keep.ids"
    if CDS_FORBID_INTERNAL_STOP:
        pep_tmp = tmpdir / "cds.trans.pep"
        run(f"seqkit translate -w 0 '{tmp1}' > '{pep_tmp}'")
        awk2 = r"""{
          id=$1; pep=$3;
          t=pep; gsub(/\*/,"",t); sc=length(pep)-length(t);
          if (%ALLOW_TAIL% && sc==1) { if (substr(pep,length(pep),1)=="*") print id; }
          else if (sc==0) print id;
        }"""
        awk2 = awk2.replace("%ALLOW_TAIL%", "1" if CDS_ALLOW_TERMINAL_STOP else "0")
        run(f"seqkit fx2tab -nl -i -s '{pep_tmp}' | awk '{awk2}' > '{tmp_ids}'")
    else:
        # 不检查内部 * 时，用所有 ID
        run(f"seqkit fx2tab -n '{tmp1}' | awk '{{print $1}}' > '{tmp_ids}'")

    # 可选：与蛋白 IDs 求交集，保证与蛋白集一一对应
    if USE_PROT_IDS and pep_in:
        pep_plain = gunzip_to(pep_in, tmpdir / (pep_in.name.replace(".gz","")))
        prot_ids = tmpdir / "prot.ids"
        run(f"seqkit fx2tab -n '{pep_plain}' | awk '{{print $1}}' | sort -u > '{prot_ids}'")
        run(f"sort -u '{tmp_ids}' | comm -12 - '{prot_ids}' > '{tmp_ids}.x'")
        shutil.move(f"{tmp_ids}.x", tmp_ids)

    # 生成最终清洗后的 CDS（到 results/cds_codon）
    ensure_dir(out_dir)
    run(f"seqkit grep -f '{tmp_ids}' '{tmp1}' > '{out}'")

    # 数量统计（单行计数）
    n_out = int(bash_out(f"awk '/^>/'" + r"{c++} END{print c+0}" + f" '{out}'"))
    if n_out == 0:
        # 清理本物种临时文件再报错退出
        cleanup_species_tmp(tmpdir)
        print(f"[ERR] {stem}: 清洗后为 0 条。建议放宽阈值（MAX_N/3倍数/内部*）或检查源 CDS 与蛋白是否匹配。", file=sys.stderr)
        sys.exit(3)
    else:
        print(f"[DONE] {stem} cds_clean={n_out}")
        # 本物种完成后清理临时目录
        cleanup_species_tmp(tmpdir)

def main():
    # 依赖检测
    for t in TOOLS: need(t)

    in_cds = Path(IN_CDS_DIR)
    in_pep = Path(IN_PEP_DIR)
    out_dir = Path(OUT_DIR)
    tmp_root = Path(TMP_ROOT)

    if not in_cds.is_dir():
        print(f"[ERR] 输入目录不存在：{in_cds}", file=sys.stderr)
        sys.exit(1)

    ensure_dir(out_dir)
    ensure_dir(tmp_root)

    # 物种清单来自 ../phylo/data/cds
    stems = list_stems_from_dir(in_cds, FA_EXTS)
    if not stems:
        print("[ERR] 未在上游目录找到任何 CDS 文件。", file=sys.stderr)
        sys.exit(1)

    for stem in stems:
        cds = find_by_stem(stem, in_cds, FA_EXTS)
        pep = find_by_stem(stem, in_pep, PEP_EXTS) if USE_PROT_IDS else None
        if not cds:
            print(f"[WARN] 跳过 {stem}: 未找到 CDS 文件")
            continue
        if USE_PROT_IDS and not pep:
            print(f"[WARN] {stem}: 未找到对应蛋白文件，按仅CDS清洗继续（可能与蛋白不同步）")
        clean_one(stem, cds, pep, tmp_root, out_dir)

    # 全部完成后，再尝试删除临时根目录
    if Path(TMP_ROOT).exists():
        shutil.rmtree(TMP_ROOT, ignore_errors=True)

if __name__ == "__main__":
    main()

