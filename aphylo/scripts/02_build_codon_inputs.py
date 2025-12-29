#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_pal2nal_only.py —— 仅使用 PAL2NAL 生成 codon MSA（geneID 直接匹配版）

核心约定（皇上当前版本）：
  1) 蛋白 MSA（OGxxxx.raw.fa）的 FASTA header 为 geneID（与 CDS header 一致）
  2) 每个 OG 的 CDS FASTA（OGxxxx.ordered_cds.fna）的 FASTA header 也为 geneID
  3) 因此 PAL2NAL 使用 ID_MATCH：pep 与 nuc 通过相同 ID 自动配对（顺序无关）

重要修复（皇上本次定位到的真实根因）：
  - 部分基因ID存在同义写法差异（例如 gnl|WGS_XXXX| vs gnl|WGS:XXXX|）
  - 本脚本会在喂给 PAL2NAL 前对 header 做一致化归一（最小且安全）：
      将 '^gnl|WGS[_:]' 统一为 'gnl|WGS:'
    这样 AA 与 CDS 即便原始写法不同，也能通过同一 ID 匹配。

运行行为（皇上要求）：
  - 每次运行都删除旧结果产物并重新生成，不做备份：
      * out_codon_dir（results/03_codon/codon_msa）
      * tmp_dir（results/03_codon/_tmp_pal2nal）
      * og_log_dir（logs/03_pal2nal）
      * sentinel（results/03_codon/.pal2nal.done）
    总日志 total_log 覆盖写（不追加、不备份）。

输出（按皇上“少屏幕输出”规格）：
  - 屏幕仅打印：banner + 进度（每 50 OG）+ DONE 汇总
  - 细节全部写入：total_log 与每 OG log

依赖：
  - perl
  - pal2nal.pl（放在 aphylo/scripts/ 下）
"""

from __future__ import annotations

import io
import re
import sys
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple, Set


# =============================================================
# 参数区（不从命令行读；皇上需要改就改这里）
# =============================================================
APHYLO_DIR = Path(__file__).resolve().parent.parent

# 输入：蛋白对齐（geneID header）
MSA_DIR = APHYLO_DIR.parent / "phylo" / "results" / "publish" / "aphylo_ready" / "strict_sco_msa"
MSA_SUFFIX = ".raw.fa"

# 输入：每个 OG 的 CDS FASTA（geneID header）
CDS_DIR = APHYLO_DIR / "results" / "02_bt" / "tmp"
CDS_SUFFIX = ".ordered_cds.fna"

# 输出：codon MSA
OUT_CODON_DIR = APHYLO_DIR / "results" / "03_codon" / "codon_msa"

# 临时目录：写入“归一化后的 AA / 归一化且去末端 stop 的 CDS”
TMP_DIR = APHYLO_DIR / "results" / "03_codon" / "_tmp_pal2nal"

# 日志
LOGS_DIR = APHYLO_DIR / "logs"
TOTAL_LOG = LOGS_DIR / "03_pal2nal.log"
OG_LOG_DIR = LOGS_DIR / "03_pal2nal"

# sentinel
SENTINEL = APHYLO_DIR / "results" / "03_codon" / ".pal2nal.done"

# pal2nal
PAL2NAL_PL = APHYLO_DIR / "scripts" / "pal2nal.pl"
PERL = "perl"

# PAL2NAL 参数
CODON_TABLE = 1
OUTPUT_FMT = "fasta"
USE_NOGAP = False
USE_NOMISMATCH = False

# 末端 stop 处理（只剪最后一个终止密码子，不做其它修剪）
TRIM_TERMINAL_STOP = True

# 每次运行删除旧产物（皇上要求）
CLEAN_RUN = True

# 屏幕输出频率
PROGRESS_EVERY = 50


# =============================================================
# 工具函数
# =============================================================
def _flushify_std():
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x):
            self.s.write(x)
            self.s.flush()
            return len(x)
    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def rm_any(p: Path):
    if not p.exists():
        return
    if p.is_dir():
        shutil.rmtree(p)
    else:
        p.unlink()


def list_ogs(msa_dir: Path, suffix: str) -> List[Tuple[str, Path]]:
    files = sorted(msa_dir.glob(f"OG*{suffix}"))
    pat = re.compile(r"^(OG\d+)")
    out = []
    for p in files:
        m = pat.match(p.name)
        if m:
            out.append((m.group(1), p))
    return out


def normalize_id(s: str) -> str:
    """
    仅做最小且安全的归一化：
      - 把 gnl|WGS_XXXX| 与 gnl|WGS:XXXX| 统一为 gnl|WGS:XXXX|
    其他 ID 不动。
    """
    return re.sub(r"^gnl\|WGS[_:]", "gnl|WGS:", s)


def read_fasta_records(path: Path) -> List[Tuple[str, str]]:
    recs: List[Tuple[str, str]] = []
    h = None
    buf: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if h is not None:
                recs.append((h, "".join(buf)))
            h = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line)
    if h is not None:
        recs.append((h, "".join(buf)))
    return recs


def write_fasta_records(path: Path, recs: List[Tuple[str, str]]):
    with path.open("w", encoding="utf-8") as w:
        for h, seq in recs:
            w.write(f">{h}\n")
            s = "".join(seq.split()).strip()
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + "\n")


def trim_terminal_stop_codon(nt: str) -> str:
    s = "".join(nt.split()).upper().replace("U", "T")
    if len(s) >= 3 and (len(s) % 3 == 0):
        if s[-3:] in ("TAA", "TAG", "TGA"):
            return s[:-3]
    return s


def tlog_write(fp, msg: str):
    fp.write(msg + "\n")
    fp.flush()


def oglog_write(fp, msg: str):
    fp.write(msg + "\n")
    fp.flush()


def screen(msg: str):
    print(msg)


def uniq_check(ids: List[str]) -> Tuple[bool, List[str]]:
    """
    检查是否存在重复 ID（归一化后可能出现合并冲突）
    返回 (ok, duplicates)
    """
    seen = set()
    dups = []
    for x in ids:
        if x in seen:
            dups.append(x)
        else:
            seen.add(x)
    return (len(dups) == 0, dups)


# =============================================================
# 主流程
# =============================================================
def main():
    _flushify_std()
    ensure_dir(LOGS_DIR)

    if CLEAN_RUN:
        rm_any(OUT_CODON_DIR)
        rm_any(TMP_DIR)
        rm_any(OG_LOG_DIR)
        rm_any(SENTINEL)

    ensure_dir(OUT_CODON_DIR)
    ensure_dir(TMP_DIR)
    ensure_dir(OG_LOG_DIR)
    ensure_dir(SENTINEL.parent)

    # 屏幕 banner（简短）
    screen("========== 03 PAL2NAL ONLY ==========")
    screen(f"msa_dir       = {MSA_DIR}")
    screen(f"cds_dir       = {CDS_DIR}")
    screen(f"out_codon_dir = {OUT_CODON_DIR}")
    screen(f"total_log     = {TOTAL_LOG}")
    screen(f"og_log_dir    = {OG_LOG_DIR}")

    with TOTAL_LOG.open("w", encoding="utf-8") as TLOG:
        # 总日志写更全，但不刷屏
        tlog_write(TLOG, "========== 03 PAL2NAL ONLY ==========")
        tlog_write(TLOG, f"script        = {Path(__file__).resolve()}")
        tlog_write(TLOG, f"aphylo_dir     = {APHYLO_DIR}")
        tlog_write(TLOG, f"msa_dir        = {MSA_DIR}")
        tlog_write(TLOG, f"msa_suffix     = {MSA_SUFFIX}")
        tlog_write(TLOG, f"cds_dir        = {CDS_DIR}")
        tlog_write(TLOG, f"cds_suffix     = {CDS_SUFFIX}")
        tlog_write(TLOG, f"out_codon_dir  = {OUT_CODON_DIR}")
        tlog_write(TLOG, f"tmp_dir        = {TMP_DIR}")
        tlog_write(TLOG, f"total_log      = {TOTAL_LOG}")
        tlog_write(TLOG, f"og_log_dir     = {OG_LOG_DIR}")
        tlog_write(TLOG, f"pal2nal.pl     = {PAL2NAL_PL}")
        tlog_write(TLOG, f"perl           = {PERL}")
        tlog_write(TLOG, f"codon_table    = {CODON_TABLE}")
        tlog_write(TLOG, f"output_fmt     = {OUTPUT_FMT}")
        tlog_write(TLOG, f"use_nogap      = {USE_NOGAP}")
        tlog_write(TLOG, f"use_nomismatch = {USE_NOMISMATCH}")
        tlog_write(TLOG, f"trim_terminal_stop = {TRIM_TERMINAL_STOP}")
        tlog_write(TLOG, "=====================================")

        if not MSA_DIR.is_dir():
            raise FileNotFoundError(f"[ERR] msa_dir 不存在：{MSA_DIR}")
        if not CDS_DIR.is_dir():
            raise FileNotFoundError(f"[ERR] cds_dir 不存在：{CDS_DIR}")
        if not PAL2NAL_PL.is_file():
            raise FileNotFoundError(f"[ERR] pal2nal.pl 不存在：{PAL2NAL_PL}")

        ogs = list_ogs(MSA_DIR, MSA_SUFFIX)
        tlog_write(TLOG, f"total_ogs      = {len(ogs)}")
        if not ogs:
            raise RuntimeError("[ERR] 未发现任何 OG 蛋白对齐文件，请检查发布包与 msa_suffix")

        ok = 0
        fail = 0
        skipped = 0

        for idx, (og, msa_path) in enumerate(ogs, 1):
            og_log = OG_LOG_DIR / f"pal2nal_{og}.log"
            out_path = OUT_CODON_DIR / f"{og}.codon.fna"
            cds_path = CDS_DIR / f"{og}{CDS_SUFFIX}"

            # 每 OG 日志覆盖写（不刷屏）
            with og_log.open("w", encoding="utf-8") as OLOG:
                oglog_write(OLOG, f"[INFO] og={og}")
                oglog_write(OLOG, f"[INFO] msa={msa_path}")
                oglog_write(OLOG, f"[INFO] cds={cds_path}")
                oglog_write(OLOG, f"[INFO] out={out_path}")

                if not cds_path.is_file():
                    oglog_write(OLOG, "[SKIP] 缺少 CDS 文件")
                    skipped += 1
                    continue

                # 读 AA / CDS，做 header 归一化；CDS 同时做末端 stop 去除
                try:
                    aa_recs_raw = read_fasta_records(msa_path)
                    cds_recs_raw = read_fasta_records(cds_path)
                    if not aa_recs_raw:
                        raise RuntimeError("AA MSA FASTA 为空或无 header")
                    if not cds_recs_raw:
                        raise RuntimeError("CDS FASTA 为空或无 header")

                    aa_recs = []
                    changed_aa = 0
                    aa_ids_norm = []
                    for h, seq in aa_recs_raw:
                        hn = normalize_id(h)
                        if hn != h:
                            changed_aa += 1
                        aa_recs.append((hn, seq))
                        aa_ids_norm.append(hn)

                    cds_recs = []
                    changed_cds = 0
                    cds_ids_norm = []
                    for h, seq in cds_recs_raw:
                        hn = normalize_id(h)
                        if hn != h:
                            changed_cds += 1
                        nt = seq
                        if TRIM_TERMINAL_STOP:
                            nt = trim_terminal_stop_codon(nt)
                        cds_recs.append((hn, nt))
                        cds_ids_norm.append(hn)

                    oglog_write(OLOG, f"[INFO] normalize_changed_aa={changed_aa} normalize_changed_cds={changed_cds}")

                    ok1, dups1 = uniq_check(aa_ids_norm)
                    ok2, dups2 = uniq_check(cds_ids_norm)
                    if (not ok1) or (not ok2):
                        oglog_write(OLOG, "[ERROR] 归一化后出现重复 ID（可能导致匹配歧义）")
                        oglog_write(OLOG, f"[ERROR] aa_dups n={len(dups1)} sample={dups1[:10]}")
                        oglog_write(OLOG, f"[ERROR] cds_dups n={len(dups2)} sample={dups2[:10]}")
                        fail += 1
                        continue

                    set_aa: Set[str] = set(aa_ids_norm)
                    set_cds: Set[str] = set(cds_ids_norm)
                    miss_in_cds = sorted(set_aa - set_cds)
                    extra_in_cds = sorted(set_cds - set_aa)
                    if miss_in_cds or extra_in_cds:
                        oglog_write(OLOG, "[ERROR] geneID 集合不一致（已做归一化）：PAL2NAL 将无法按 ID 匹配")
                        oglog_write(OLOG, f"[ERROR] in_AA_not_in_CDS n={len(miss_in_cds)} sample={miss_in_cds[:10]}")
                        oglog_write(OLOG, f"[ERROR] in_CDS_not_in_AA n={len(extra_in_cds)} sample={extra_in_cds[:10]}")
                        fail += 1
                        continue

                    # 写临时文件（归一化后的 AA / CDS）
                    aa_tmp = TMP_DIR / f"{og}.aa.norm.fa"
                    cds_tmp = TMP_DIR / f"{og}.cds.clean.norm.fna"
                    write_fasta_records(aa_tmp, aa_recs)
                    write_fasta_records(cds_tmp, cds_recs)

                except Exception as e:
                    oglog_write(OLOG, f"[ERROR] 预处理失败：{e}")
                    fail += 1
                    continue

                # PAL2NAL（ID_MATCH：不改 header，不重排）
                cmd = [
                    PERL, str(PAL2NAL_PL),
                    str(aa_tmp),
                    str(cds_tmp),
                    "-output", OUTPUT_FMT,
                    "-codontable", str(CODON_TABLE),
                ]
                if USE_NOGAP:
                    cmd.append("-nogap")
                if USE_NOMISMATCH:
                    cmd.append("-nomismatch")

                oglog_write(OLOG, f"[INFO] cmd={' '.join(cmd)}")

                try:
                    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                except Exception as e:
                    oglog_write(OLOG, f"[ERROR] 运行 pal2nal 失败：{e}")
                    fail += 1
                    continue

                if r.stderr.strip():
                    oglog_write(OLOG, "")
                    oglog_write(OLOG, "[PAL2NAL STDERR]")
                    oglog_write(OLOG, r.stderr.rstrip())

                out_text = (r.stdout or "").strip()
                if r.returncode != 0:
                    oglog_write(OLOG, f"[ERROR] pal2nal returncode={r.returncode}")
                    fail += 1
                    continue
                if not out_text.startswith(">"):
                    oglog_write(OLOG, "[ERROR] pal2nal 输出不含任何 FASTA header，疑似空输出或被判不一致")
                    fail += 1
                    continue

                out_path.write_text(r.stdout, encoding="utf-8")
                ok += 1

            if (idx % PROGRESS_EVERY == 0) or (idx == len(ogs)):
                screen(f"[PROGRESS] {idx}/{len(ogs)} ok={ok} fail={fail} skipped={skipped}")
                tlog_write(TLOG, f"[PROGRESS] {idx}/{len(ogs)} ok={ok} fail={fail} skipped={skipped}")

        SENTINEL.touch()

        screen(f"[DONE] ok={ok} fail={fail} skipped={skipped}")
        screen(f"[DONE] sentinel={SENTINEL}")
        screen(f"[DONE] total_log={TOTAL_LOG}")
        screen(f"[DONE] og_log_dir={OG_LOG_DIR}")

        tlog_write(TLOG, f"[DONE] ok={ok} fail={fail} skipped={skipped}")
        tlog_write(TLOG, f"[DONE] sentinel={SENTINEL}")
        tlog_write(TLOG, f"[DONE] total_log={TOTAL_LOG}")
        tlog_write(TLOG, f"[DONE] og_log_dir={OG_LOG_DIR}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)

