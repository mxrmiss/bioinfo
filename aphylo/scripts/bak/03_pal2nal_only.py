#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_pal2nal_only.py —— Codon MSA（方案B最终版：抽蛋白→MAFFT→trimAl(automated1)+backtrans→codon alignment）

核心策略（皇上最终确认）：
1) 不使用 OrthoFinder 的 MSA 作为输入；仅用其 strict_sco_msa 目录的 OG 文件来枚举 OG 列表；
2) 每个 OG：从 ordered_cds.fna 解析 key → 从 proteomes 抽对应蛋白 → MAFFT 做蛋白对齐；
3) 用 trimAl -automated1 在蛋白对齐上做修剪，并用 -backtrans 将修剪同步回译到 CDS，直接输出 codon MSA；
4) 某个 OG 失败：记录日志并 continue 下一个 OG（不中断全局）；
5) 每次运行删除旧的 03 产物与临时目录；总日志覆盖写；每 OG 单独日志。

输入（默认路径可在脚本顶部修改）：
- OG 列表：phylo/results/publish/aphylo_ready/strict_sco_msa/OG*.raw.fa（仅用于枚举 OG）
- CDS：aphylo/results/02_bt/tmp/OGxxxx.ordered_cds.fna
- Proteome：phylo/data/proteomes/*.fa|*.faa|*.fasta

输出：
- aphylo/results/03_codon/codon_msa/OGxxxx.codon.fna          # 给 codeml 的 codon alignment（已修剪）
- aphylo/results/03_codon/pep_trimal/OGxxxx.pep.trimal.fa      # 修剪后的蛋白对齐（用于位点→蛋白坐标映射）
- aphylo/results/03_codon/colnumbering/OGxxxx.colnumbering.txt # 修剪后列→修剪前列映射
- aphylo/results/03_codon/.pal2nal.done                        # 哨兵文件（名称沿用旧习惯）
- aphylo/logs/03_pal2nal.log                                   # 总日志（覆盖）
- aphylo/logs/03_pal2nal/pal2nal_OGxxxx.log                    # 每 OG 日志（逐次覆盖）

依赖：
- mafft
- trimal
"""

from __future__ import annotations

import re
import sys
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Set, Tuple


# ==========================
# 皇上在这里改参数（不走命令行）
# ==========================

MSA_DIR = Path("~/project/phylo/results/publish/aphylo_ready/strict_sco_msa").expanduser()
MSA_SUFFIX = ".raw.fa"

CDS_DIR = Path("~/project/aphylo/results/02_bt/tmp").expanduser()
CDS_SUFFIX = ".ordered_cds.fna"

PROTEOME_DIR = Path("~/project/phylo/data/proteomes").expanduser()
PROTEOME_SUFFIXES = (".fa", ".faa", ".fasta")

OUT_CODON_DIR = Path("~/project/aphylo/results/03_codon/codon_msa").expanduser()
OUT_PEP_TRIM_DIR = Path("~/project/aphylo/results/03_codon/pep_trimal").expanduser()
OUT_COLNUM_DIR = Path("~/project/aphylo/results/03_codon/colnumbering").expanduser()

TMP_DIR = Path("~/project/aphylo/results/03_codon/_tmp_pal2nal").expanduser()
SENTINEL = Path("~/project/aphylo/results/03_codon/.pal2nal.done").expanduser()

LOG_TOTAL = Path("~/project/aphylo/logs/03_pal2nal.log").expanduser()
LOG_OG_DIR = Path("~/project/aphylo/logs/03_pal2nal").expanduser()

MAFFT = "mafft"
MAFFT_THREADS = 20

TRIMAL = "trimal"
TRIMAL_MODE = "automated1"

TRIM_TERMINAL_STOP = True
TRIM_TERMINAL_STAR_IN_PEP = True

# 为了防止修剪过度导致信息量极低：可设为 0 关闭
TRIMAL_MIN_ALIGN_LEN_AA = 30

PROGRESS_EVERY = 50


# ==========================
# 工具函数
# ==========================

def _now_banner() -> str:
    return "========== 03 CODON MSA (MAFFT -> trimAl automated1 + backtrans) =========="

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _rm_any(p: Path) -> None:
    if p.is_symlink() or p.is_file():
        p.unlink(missing_ok=True)
    elif p.is_dir():
        shutil.rmtree(p, ignore_errors=True)

def _write_text(p: Path, s: str) -> None:
    _ensure_dir(p.parent)
    with p.open("w", encoding="utf-8") as w:
        w.write(s)
        w.flush()

def _touch(p: Path) -> None:
    _ensure_dir(p.parent)
    p.touch()

def _key_from_header(h: str) -> str:
    """
    统一的 key 解析规则：
    - 取第一个 token
    - 如果包含 |，取最后一个 | 后的部分
    """
    h = (h or "").strip().split()[0]
    if "|" in h:
        return h.split("|")[-1]
    return h

def _read_fasta_iter(path: Path):
    """轻量 FASTA 解析器：yield (header_first_token, seq)"""
    h = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq_chunks)
                h = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if h is not None:
            yield h, "".join(seq_chunks)

def _stop_trim_nt(nt: str) -> str:
    """仅去掉末端 stop（TAA/TAG/TGA），不动内部"""
    nt = re.sub(r"\s+", "", nt.upper().replace("U", "T"))
    if not TRIM_TERMINAL_STOP:
        return nt
    if len(nt) >= 3 and len(nt) % 3 == 0 and nt[-3:] in {"TAA", "TAG", "TGA"}:
        return nt[:-3]
    return nt

def _trim_terminal_star_aa(aa: str) -> str:
    """仅去掉末端 *，不动内部"""
    aa = re.sub(r"\s+", "", aa)
    if not TRIM_TERMINAL_STAR_IN_PEP:
        return aa
    return aa[:-1] if aa.endswith("*") else aa

def _count_fasta_headers(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    n = 0
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n

def _headers_set(path: Path) -> Set[str]:
    s: Set[str] = set()
    if not path.exists() or path.stat().st_size == 0:
        return s
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip().split()[0]
                s.add(_key_from_header(h))
    return s

def _first_aln_seq_len(path: Path) -> int:
    for _h, seq in _read_fasta_iter(path):
        return len(re.sub(r"\s+", "", seq))
    return 0


# ==========================
# 主流程
# ==========================

def main() -> int:
    # ---- 基础检查
    if not MSA_DIR.is_dir():
        print(f"[ERR] msa_dir 不存在：{MSA_DIR}", file=sys.stderr, flush=True)
        return 2
    if not CDS_DIR.is_dir():
        print(f"[ERR] cds_dir 不存在：{CDS_DIR}", file=sys.stderr, flush=True)
        return 2
    if not PROTEOME_DIR.is_dir():
        print(f"[ERR] proteome_dir 不存在：{PROTEOME_DIR}", file=sys.stderr, flush=True)
        return 2

    # ---- 每次运行：删除旧产物（03 产物 & 临时文件 & OG 日志），总日志覆盖
    _rm_any(OUT_CODON_DIR)
    _rm_any(OUT_PEP_TRIM_DIR)
    _rm_any(OUT_COLNUM_DIR)
    _rm_any(TMP_DIR)
    _rm_any(SENTINEL)

    _rm_any(LOG_OG_DIR)
    _ensure_dir(LOG_OG_DIR)

    _ensure_dir(LOG_TOTAL.parent)
    total_log_fp = LOG_TOTAL.open("w", encoding="utf-8")

    def tlog(msg: str, also_stdout: bool = False) -> None:
        total_log_fp.write(msg.rstrip("\n") + "\n")
        total_log_fp.flush()
        if also_stdout:
            print(msg, flush=True)

    # ---- Banner + 参数回显（屏幕克制，只打印少量）
    tlog(_now_banner(), also_stdout=True)
    tlog(f"msa_dir       = {MSA_DIR}", also_stdout=True)
    tlog(f"cds_dir       = {CDS_DIR}", also_stdout=True)
    tlog(f"proteome_dir  = {PROTEOME_DIR}")
    tlog(f"out_codon_dir = {OUT_CODON_DIR}", also_stdout=True)
    tlog(f"out_pep_trim  = {OUT_PEP_TRIM_DIR}", also_stdout=True)
    tlog(f"out_colnum    = {OUT_COLNUM_DIR}", also_stdout=True)
    tlog(f"tmp_dir       = {TMP_DIR}")
    tlog(f"total_log     = {LOG_TOTAL}", also_stdout=True)
    tlog(f"og_log_dir    = {LOG_OG_DIR}", also_stdout=True)
    tlog(f"mafft         = {MAFFT}")
    tlog(f"mafft_threads = {MAFFT_THREADS}")
    tlog(f"trimal        = {TRIMAL}")
    tlog(f"trimal_mode   = {TRIMAL_MODE}")
    tlog(f"trimal_min_aln_len_aa = {TRIMAL_MIN_ALIGN_LEN_AA}")
    tlog(f"trim_terminal_stop    = {TRIM_TERMINAL_STOP}")
    tlog(f"trim_terminal_star_pep= {TRIM_TERMINAL_STAR_IN_PEP}")
    tlog("======================================")

    # ---- 枚举 OG（仅用于 OG 列表）
    og_files = sorted(MSA_DIR.glob(f"OG*{MSA_SUFFIX}"))
    og_pat = re.compile(r"^(OG\d+)")
    ogs: List[str] = []
    for p in og_files:
        m = og_pat.match(p.name)
        if m:
            ogs.append(m.group(1))
    if not ogs:
        tlog(f"[ERR] 未发现 OG 文件：{MSA_DIR}/OG*{MSA_SUFFIX}", also_stdout=True)
        total_log_fp.close()
        return 3
    tlog(f"total_ogs     = {len(ogs)}", also_stdout=True)

    # ---- 收集每个 OG 的 key 列表，并汇总全局 key 集合
    og2keys: Dict[str, List[str]] = {}
    all_need_keys: Set[str] = set()
    og_missing_cds: List[str] = []
    for og in ogs:
        cds_path = CDS_DIR / f"{og}{CDS_SUFFIX}"
        if not cds_path.is_file():
            og_missing_cds.append(og)
            continue
        keys: List[str] = []
        seen: Set[str] = set()
        for h, _seq in _read_fasta_iter(cds_path):
            key = _key_from_header(h)
            if key and key not in seen:
                keys.append(key)
                seen.add(key)
        og2keys[og] = keys
        all_need_keys.update(keys)

    if og_missing_cds:
        tlog(f"[WARN] 缺少 CDS 文件的 OG 数：{len(og_missing_cds)} sample={og_missing_cds[:5]}")

    # ---- 扫描 proteomes 一次：只抽取 all_need_keys 对应的蛋白序列（按 key 匹配）
    proteome_files: List[Path] = []
    for suf in PROTEOME_SUFFIXES:
        proteome_files.extend(sorted(PROTEOME_DIR.glob(f"*{suf}")))
    if not proteome_files:
        tlog(f"[ERR] proteome_dir 未发现 fasta：{PROTEOME_DIR} suffix={PROTEOME_SUFFIXES}", also_stdout=True)
        total_log_fp.close()
        return 4

    key2pep: Dict[str, str] = {}
    for fp in proteome_files:
        for h, seq in _read_fasta_iter(fp):
            key = _key_from_header(h)
            if key in all_need_keys and key not in key2pep:
                key2pep[key] = _trim_terminal_star_aa(seq)
        if len(key2pep) == len(all_need_keys):
            break

    # ---- 正式逐 OG：写 pep.raw.fa + cds.key.fna → MAFFT → trimAl(backtrans) 输出 codon + pep.trimal + colnumbering
    _ensure_dir(OUT_CODON_DIR)
    _ensure_dir(OUT_PEP_TRIM_DIR)
    _ensure_dir(OUT_COLNUM_DIR)
    _ensure_dir(TMP_DIR)

    ok = 0
    fail = 0
    skipped = 0

    for i, og in enumerate(ogs, start=1):
        og_log = LOG_OG_DIR / f"pal2nal_{og}.log"
        with og_log.open("w", encoding="utf-8") as w:
            def olog(msg: str) -> None:
                w.write(msg.rstrip("\n") + "\n")
                w.flush()

            cds_in = CDS_DIR / f"{og}{CDS_SUFFIX}"
            if not cds_in.is_file():
                olog(f"[SKIP] 缺少 CDS：{cds_in}")
                skipped += 1
                continue

            keys = og2keys.get(og, [])
            if not keys:
                olog(f"[SKIP] CDS 中未解析到任何序列：{cds_in}")
                skipped += 1
                continue

            missing_keys = [k for k in keys if k not in key2pep]
            if missing_keys:
                olog(f"[ERROR] proteome 缺少蛋白 key：n={len(missing_keys)} sample={missing_keys[:5]}")
                fail += 1
                continue

            pep_raw = TMP_DIR / f"{og}.pep.raw.fa"
            pep_mafft = TMP_DIR / f"{og}.pep.mafft.fa"
            cds_key = TMP_DIR / f"{og}.cds.key.fna"

            out_codon = OUT_CODON_DIR / f"{og}.codon.fna"
            out_pep_trimal = OUT_PEP_TRIM_DIR / f"{og}.pep.trimal.fa"
            out_colnum = OUT_COLNUM_DIR / f"{og}.colnumbering.txt"

            expected_keys: Set[str] = set(keys)

            # 1) 写 pep.raw.fa（header=key）
            with pep_raw.open("w", encoding="utf-8") as pr:
                for k in keys:
                    pr.write(f">{k}\n")
                    s = key2pep[k]
                    for j in range(0, len(s), 60):
                        pr.write(s[j:j+60] + "\n")

            if _count_fasta_headers(pep_raw) != len(keys):
                olog(f"[ERROR] pep.raw.fa 序列数异常：expect={len(keys)} got={_count_fasta_headers(pep_raw)} path={pep_raw}")
                fail += 1
                continue

            # 2) 写 cds.key.fna（header=key；可选去末端 stop）
            with cds_key.open("w", encoding="utf-8") as ck:
                for h, seq in _read_fasta_iter(cds_in):
                    key = _key_from_header(h)
                    nt = _stop_trim_nt(seq)
                    ck.write(f">{key}\n")
                    for j in range(0, len(nt), 60):
                        ck.write(nt[j:j+60] + "\n")

            if _count_fasta_headers(cds_key) != len(keys):
                olog(f"[ERROR] cds.key.fna 序列数异常：expect={len(keys)} got={_count_fasta_headers(cds_key)} path={cds_key}")
                fail += 1
                continue

            # 3) MAFFT
            mafft_cmd = [MAFFT, "--auto"]
            if MAFFT_THREADS and int(MAFFT_THREADS) > 0:
                mafft_cmd += ["--thread", str(int(MAFFT_THREADS))]
            mafft_cmd += [str(pep_raw)]
            olog(f"[INFO] mafft_cmd={' '.join(mafft_cmd)}")

            with pep_mafft.open("w", encoding="utf-8") as out_f:
                r = subprocess.run(
                    mafft_cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    text=True
                )
            if r.returncode != 0 or (not pep_mafft.exists()) or pep_mafft.stat().st_size == 0:
                olog("[ERROR] MAFFT 失败")
                if r.stderr:
                    olog("[MAFFT STDERR]")
                    olog(r.stderr)
                fail += 1
                continue

            # 4) trimAl：输出 codon alignment（修剪+回译同步）
            tr_back_cmd = [
                TRIMAL,
                "-in", str(pep_mafft),
                "-out", str(out_codon),
                "-fasta",
                f"-{TRIMAL_MODE}",
                "-backtrans", str(cds_key),
                "-ignorestopcodon",
                "-keepheader",
            ]
            olog(f"[INFO] trimal_backtrans_cmd={' '.join(tr_back_cmd)}")

            rbt = subprocess.run(
                tr_back_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if rbt.returncode != 0 or (not out_codon.exists()) or out_codon.stat().st_size == 0:
                olog("[ERROR] trimAl(backtrans) 失败或输出为空")
                if rbt.stdout:
                    olog("[TRIMAL STDOUT]")
                    olog(rbt.stdout)
                if rbt.stderr:
                    olog("[TRIMAL STDERR]")
                    olog(rbt.stderr)
                fail += 1
                continue

            # 5) trimAl：输出修剪后的 AA alignment（用于位点映射）
            tr_pep_cmd = [
                TRIMAL,
                "-in", str(pep_mafft),
                "-out", str(out_pep_trimal),
                "-fasta",
                f"-{TRIMAL_MODE}",
                "-keepheader",
            ]
            olog(f"[INFO] trimal_pep_cmd={' '.join(tr_pep_cmd)}")

            rpep = subprocess.run(
                tr_pep_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if rpep.returncode != 0 or (not out_pep_trimal.exists()) or out_pep_trimal.stat().st_size == 0:
                olog("[ERROR] trimAl(pep) 失败或输出为空")
                if rpep.stdout:
                    olog("[TRIMAL STDOUT]")
                    olog(rpep.stdout)
                if rpep.stderr:
                    olog("[TRIMAL STDERR]")
                    olog(rpep.stderr)
                fail += 1
                continue

            # 6) trimAl：输出列映射（修剪后列 → 修剪前列）
            tr_col_cmd = [
                TRIMAL,
                "-in", str(pep_mafft),
                "-out", "/dev/null",
                "-fasta",
                f"-{TRIMAL_MODE}",
                "-colnumbering",
            ]
            olog(f"[INFO] trimal_colnumbering_cmd={' '.join(tr_col_cmd)}")

            rcol = subprocess.run(
                tr_col_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if rcol.returncode != 0:
                olog("[ERROR] trimAl(colnumbering) 失败")
                if rcol.stderr:
                    olog("[TRIMAL STDERR]")
                    olog(rcol.stderr)
                fail += 1
                continue
            _write_text(out_colnum, rcol.stdout or "")

            # ---- 质量检查（保守而不苛刻）
            n_out = _count_fasta_headers(out_codon)
            if n_out != len(keys):
                olog(f"[ERROR] codon 输出序列数异常：expect={len(keys)} got={n_out} path={out_codon}")
                fail += 1
                continue

            pep_keys = _headers_set(out_pep_trimal)
            if pep_keys != expected_keys:
                olog("[ERROR] pep.trimal.fa 的 header key 集合不一致")
                olog(f"        expected={len(expected_keys)} got={len(pep_keys)}")
                diff1 = sorted(list(expected_keys - pep_keys))[:10]
                diff2 = sorted(list(pep_keys - expected_keys))[:10]
                olog(f"        missing_in_pep sample={diff1}")
                olog(f"        extra_in_pep   sample={diff2}")
                fail += 1
                continue

            aln_len = _first_aln_seq_len(out_pep_trimal)
            if int(TRIMAL_MIN_ALIGN_LEN_AA) > 0 and aln_len < int(TRIMAL_MIN_ALIGN_LEN_AA):
                olog(f"[ERROR] 修剪后 AA 对齐长度过短：len={aln_len} < {TRIMAL_MIN_ALIGN_LEN_AA} path={out_pep_trimal}")
                fail += 1
                continue

            # ---- 成功
            ok += 1

        if i % PROGRESS_EVERY == 0 or i == len(ogs):
            tlog(f"[PROGRESS] {i}/{len(ogs)} ok={ok} fail={fail} skipped={skipped}", also_stdout=True)

    _touch(SENTINEL)
    tlog(f"[DONE] ok={ok} fail={fail} skipped={skipped}", also_stdout=True)
    tlog(f"[DONE] sentinel={SENTINEL}", also_stdout=True)
    tlog(f"[DONE] total_log={LOG_TOTAL}", also_stdout=True)
    tlog(f"[DONE] og_log_dir={LOG_OG_DIR}", also_stdout=True)

    total_log_fp.close()
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

