#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07c_beb_map.py —— BEB 位点映射与 3D 标注清单导出

功能：
1) 读取 results/05_cmlagg/D_beb_sites.tsv（BEB 位点长表）
2) 映射到：
   - results/03_codon/pep_trimal/OGxxxx.pep.trimal.fa（修剪后蛋白对齐）
   - results/03_codon/colnumbering/OGxxxx.colnumbering.txt（修剪后列 -> MAFFT 列映射）
   - results/03_codon/_tmp_pal2nal/OGxxxx.pep.mafft.fa（修剪前 MAFFT 蛋白对齐）
   - phylo/data/proteomes/*（原始蛋白，用于核对 raw_aa）
3) 输出两张表（覆盖写）：
   A) results/05_cmlagg/D_beb_sites_mapped.tsv
      —— 每个位点一行（可追溯、可调试）
   B) results/05_cmlagg/D_beb_sites_mapped_gene.tsv
      —— 精简版：一个 (OG, foreground, gene_id) 一行；忽略 gap；只保留合并字段
         sites_3d_label 形如：Y39(0.998);H40(0.998);T44(0.992)

用法：
- 不走命令行参数；只需在“皇上在这里改参数”区域按需改路径/开关
- 运行：python3 scripts/07c_beb_map.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ==========================
# 皇上在这里改参数（不走命令行）
# ==========================

# —— 项目根目录：默认脚本位于 aphylo/scripts/ 下，因此根目录是 scripts 的上一级
APHYLO_ROOT = Path(__file__).resolve().parent.parent

# —— 输入：06 聚合输出的位点长表
INPUT_TSV = APHYLO_ROOT / "results/05_cmlagg/D_beb_sites.tsv"

# —— 03 输出：修剪后蛋白对齐、列映射、修剪前 MAFFT 蛋白对齐
PEP_TRIMAL_DIR = APHYLO_ROOT / "results/03_codon/pep_trimal"
COLNUM_DIR = APHYLO_ROOT / "results/03_codon/colnumbering"
MAFFT_DIR = APHYLO_ROOT / "results/03_codon/_tmp_pal2nal"  # 里面应有 OGxxxx.pep.mafft.fa

# —— 输出：映射后的长表（覆盖）
OUTPUT_MAPPED_TSV = APHYLO_ROOT / "results/05_cmlagg/D_beb_sites_mapped.tsv"

# —— 输出：精简基因表（覆盖）
OUTPUT_GENE_TSV = APHYLO_ROOT / "results/05_cmlagg/D_beb_sites_mapped_gene.tsv"

# —— 是否用原始 proteome 核对 raw_aa（建议 True）
DO_PROTEOME_CHECK = True

# phylo 项目默认与 aphylo 同级：~/project/phylo
PHYLO_ROOT = APHYLO_ROOT.parent / "phylo"
PROTEOME_DIR = PHYLO_ROOT / "data/proteomes"
PROTEOME_SUFFIXES = (".fa", ".faa", ".fasta")

# —— 打印进度（行数很多时用）
PROGRESS_EVERY = 20000

# —— sites_3d_label 里 posterior 保留的小数位（建议 3 位：0.998）
POST_DECIMALS = 3


# ==========================
# 基础工具函数
# ==========================

def normalize_id(x: str) -> str:
    """
    统一 ID：
    - 取第一个 token
    - 若有 |，取最后一个 | 后面
    """
    x = (x or "").strip().split()[0]
    if "|" in x:
        return x.split("|")[-1]
    return x


def read_fasta_as_dict(fp: Path) -> Dict[str, str]:
    """
    读取 FASTA，返回 {header_first_token: sequence}（sequence 保留 '-'）
    """
    d: Dict[str, str] = {}
    if not fp.exists():
        return d
    key = None
    buf: List[str] = []
    with fp.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if key is not None:
                    d[key] = "".join(buf)
                key = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if key is not None:
            d[key] = "".join(buf)
    return d


def split_gene_ids(s: str) -> List[str]:
    """
    foreground_gene_ids 可能是：
    - 单个：Sco01g11060.1
    - 多个：A,B,C 或 A;B;C 或 A|B|C
    """
    s = (s or "").strip()
    if not s:
        return []
    parts = re.split(r"[;,| ]+", s)
    out: List[str] = []
    for p in parts:
        p = p.strip()
        if p:
            out.append(p)
    # 去重保持顺序
    seen = set()
    uniq: List[str] = []
    for x in out:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return uniq


def load_proteomes_index(proteome_dir: Path, suffixes: Tuple[str, ...]) -> Dict[str, str]:
    """
    扫描 PROTEOME_DIR，构建 {normalized_id: aa_sequence} 索引。
    """
    idx: Dict[str, str] = {}
    files: List[Path] = []
    for suf in suffixes:
        files.extend(sorted(proteome_dir.glob(f"*{suf}")))
    for fp in files:
        d = read_fasta_as_dict(fp)
        for hid, seq in d.items():
            k = normalize_id(hid)
            if k and k not in idx:
                idx[k] = seq.replace("\r", "").replace("\n", "")
    return idx


def get_raw_aa_at(seq: str, residue_idx_1based: int) -> Optional[str]:
    """
    在原始蛋白序列中取第 residue_idx_1based 位（1-based）
    """
    if residue_idx_1based <= 0 or residue_idx_1based > len(seq):
        return None
    return seq[residue_idx_1based - 1]


def aln_char_at(seq: str, col_1based: int) -> str:
    """
    alignment 序列在 col_1based 列的字符（1-based）
    """
    if not seq or col_1based <= 0 or col_1based > len(seq):
        return ""
    return seq[col_1based - 1]


def residue_index_from_alignment(seq: str, col_1based: int) -> Optional[int]:
    """
    在 alignment 序列中，计算 col_1based 对应的“去 gap residue index”(1-based)
    - 若该列为 gap 或越界，返回 None
    """
    if not seq or col_1based <= 0 or col_1based > len(seq):
        return None
    c = seq[col_1based - 1]
    if c in {"-", "."}:
        return None
    pre = seq[:col_1based].replace("-", "").replace(".", "")
    return len(pre)


def find_seq_by_id(aln: Dict[str, str], gene_id_norm: str) -> Tuple[Optional[str], str]:
    """
    在 alignment dict 中找到与 gene_id_norm 匹配的序列：
    - 先直接 key 命中
    - 再扫描 normalize_id(hid) 命中
    返回 (seq, hit_id)
    """
    if gene_id_norm in aln:
        return aln[gene_id_norm], gene_id_norm
    for hid, seq in aln.items():
        if normalize_id(hid) == gene_id_norm:
            return seq, hid
    return None, ""


def parse_colnumbering_strict(fp: Path, target_len: int) -> Optional[List[int]]:
    """
    解析 trimAl 的 -colnumbering 输出，得到“修剪后列 -> 修剪前（MAFFT）列”的映射。

    关键点（这也是你现在 fg_aa_trim vs fg_aa_mafft 对不上的根因）：
    - trimAl 的 colnumbering 文件里真正可靠的是 `ColumnsMap`。
    - `ColumnsMap` 会包含：
        * -1：该列被 trimAl 剪掉了
        * 0-based 的列号：该列在“修剪前（MAFFT 对齐）”中的位置
    - 你的 BEB site 是“修剪后对齐”的 1-based 列号，因此：
        mafft_col_1based = (ColumnsMap 过滤掉 -1 后取第 site 个值) + 1

    返回：
      colmap0[i] = 修剪后第 (i+1) 列对应修剪前（MAFFT）第几列（0-based）
    """
    if (not fp.exists()) or fp.stat().st_size == 0:
        return None
    if target_len <= 0:
        return None

    txt = fp.read_text(encoding="utf-8", errors="ignore")
    p = txt.find("ColumnsMap")
    if p < 0:
        return None

    # 从 ColumnsMap 开始向后，把所有整数（包含 -1）都捞出来
    nums = [int(x) for x in re.findall(r"-?\d+", txt[p:])]
    if not nums:
        return None

    kept = [n for n in nums if n != -1]

    # 正常情况下 kept 的长度应当等于修剪后对齐长度（target_len）
    if len(kept) != target_len:
        return None

    return kept


def relpath(p: Path) -> str:
    """
    输出表里统一用相对 aphylo 根目录的路径，方便跨机器复现
    """
    try:
        return str(p.relative_to(APHYLO_ROOT))
    except Exception:
        return str(p)


def fmt_post(s: str, decimals: int = 3) -> str:
    """
    posterior 格式化为固定小数位（默认 3 位）
    """
    try:
        x = float(s)
        return f"{x:.{decimals}f}"
    except Exception:
        return s


# ==========================
# 主程序
# ==========================

def main() -> int:
    if not INPUT_TSV.exists():
        print(f"[ERR] INPUT_TSV not found: {INPUT_TSV}", file=sys.stderr, flush=True)
        return 2
    if not PEP_TRIMAL_DIR.is_dir():
        print(f"[ERR] PEP_TRIMAL_DIR not found: {PEP_TRIMAL_DIR}", file=sys.stderr, flush=True)
        return 2
    if not COLNUM_DIR.is_dir():
        print(f"[ERR] COLNUM_DIR not found: {COLNUM_DIR}", file=sys.stderr, flush=True)
        return 2
    if not MAFFT_DIR.is_dir():
        print(f"[ERR] MAFFT_DIR not found: {MAFFT_DIR}", file=sys.stderr, flush=True)
        return 2

    proteome_index: Optional[Dict[str, str]] = None
    if DO_PROTEOME_CHECK:
        if not PROTEOME_DIR.is_dir():
            print(f"[WARN] PROTEOME_DIR not found, proteome check disabled: {PROTEOME_DIR}", flush=True)
        else:
            print(f"[INFO] Loading proteomes index: {PROTEOME_DIR}", flush=True)
            proteome_index = load_proteomes_index(PROTEOME_DIR, PROTEOME_SUFFIXES)
            print(f"[INFO] Proteomes indexed: n={len(proteome_index)}", flush=True)

    # 缓存（避免重复读同一个 OG 的 fasta）
    og2pep_trim: Dict[str, Dict[str, str]] = {}
    og2pep_mafft: Dict[str, Dict[str, str]] = {}
    og2colmap: Dict[str, List[int]] = {}
    og2trim_len: Dict[str, int] = {}

    # 精简表聚合容器：key=(OG, foreground, gene_id)
    # value=list of (residue_index_raw:int, raw_aa:str, post_float:float, post_str_fmt:str)
    gene_sites: Dict[Tuple[str, str, str], List[Tuple[int, str, float, str]]] = {}

    with INPUT_TSV.open("r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n")
        if not header:
            print(f"[ERR] Empty file: {INPUT_TSV}", file=sys.stderr, flush=True)
            return 3

        cols = header.split("\t")
        col2i = {c: i for i, c in enumerate(cols)}
        required = ["OG", "foreground", "site", "aa", "post", "foreground_gene_ids"]
        missing = [c for c in required if c not in col2i]
        if missing:
            print(f"[ERR] Missing columns: {missing}", file=sys.stderr, flush=True)
            print(f"[ERR] Found columns: {cols}", file=sys.stderr, flush=True)
            return 3

        # A) 长表输出表头（与你现有列一致/兼容）
        mapped_cols = [
            "OG", "foreground", "site", "post",
            "gene_id",
            "ref_aa",
            "input_aa",
            "fg_aa_trim",
            "mafft_col",
            "fg_aa_mafft",
            "residue_index_raw",
            "raw_aa",
            "match_input_ref",
            "match_raw_fg",
            "status",
            "pep_trimal_fa",
            "pep_mafft_fa",
            "colnumbering_txt",
        ]

        OUTPUT_MAPPED_TSV.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_GENE_TSV.parent.mkdir(parents=True, exist_ok=True)

        with OUTPUT_MAPPED_TSV.open("w", encoding="utf-8") as w:
            w.write("\t".join(mapped_cols) + "\n")

            n_in = 0
            n_out = 0

            stat_gap = 0
            stat_ok = 0
            stat_ok_match = 0

            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                n_in += 1
                parts = line.split("\t")
                if len(parts) < len(cols):
                    continue

                og = parts[col2i["OG"]].strip()
                fg = parts[col2i["foreground"]].strip()
                site_str = parts[col2i["site"]].strip()
                input_aa = parts[col2i["aa"]].strip()
                post_raw = parts[col2i["post"]].strip()
                gene_ids_raw = parts[col2i["foreground_gene_ids"]].strip()

                if not og:
                    continue
                try:
                    site = int(site_str)
                except Exception:
                    continue

                gene_ids = split_gene_ids(gene_ids_raw)
                if not gene_ids:
                    gene_ids = [""]

                pep_trim_fp = PEP_TRIMAL_DIR / f"{og}.pep.trimal.fa"
                pep_mafft_fp = MAFFT_DIR / f"{og}.pep.mafft.fa"
                col_fp = COLNUM_DIR / f"{og}.colnumbering.txt"

                # 读/缓存 pep_trimal
                if og not in og2pep_trim:
                    og2pep_trim[og] = read_fasta_as_dict(pep_trim_fp) if pep_trim_fp.exists() else {}
                    if og2pep_trim[og]:
                        ref_id0 = next(iter(og2pep_trim[og].keys()))
                        og2trim_len[og] = len(og2pep_trim[og][ref_id0])
                    else:
                        og2trim_len[og] = 0

                # 读/缓存 pep_mafft
                if og not in og2pep_mafft:
                    og2pep_mafft[og] = read_fasta_as_dict(pep_mafft_fp) if pep_mafft_fp.exists() else {}

                # 读/缓存 colmap（严格解析：ColumnsMap）
                if og not in og2colmap:
                    target_len = og2trim_len.get(og, 0)
                    cm = parse_colnumbering_strict(col_fp, target_len) if target_len > 0 else None
                    og2colmap[og] = cm if cm is not None else []

                pep_trim = og2pep_trim[og]
                pep_mafft = og2pep_mafft[og]
                colmap = og2colmap[og]

                # 参考序列：pep_trimal 第一条
                ref_seq = ""
                if pep_trim:
                    ref_id = next(iter(pep_trim.keys()))
                    ref_seq = pep_trim[ref_id]

                ref_aa = aln_char_at(ref_seq, site) if ref_seq else ""

                # colmap 返回的是 0-based 的“修剪前列号”
                mafft_col0: Optional[int] = None
                if colmap and 1 <= site <= len(colmap):
                    mafft_col0 = colmap[site - 1]
                # 转成 1-based，供 aln_char_at / residue_index_from_alignment 使用
                mafft_col1: Optional[int] = (mafft_col0 + 1) if mafft_col0 is not None else None

                # posterior 统一格式化（用于 label）
                post_fmt = fmt_post(post_raw, POST_DECIMALS)
                post_float = None
                try:
                    post_float = float(post_raw)
                except Exception:
                    post_float = float("nan")

                for gid in gene_ids:
                    gid_norm = normalize_id(gid)
                    status = "OK"

                    if (not pep_trim_fp.exists()) and status == "OK":
                        status = "NO_PEP_TRIMAL"
                    if (not pep_mafft_fp.exists()) and status == "OK":
                        status = "NO_PEP_MAFFT"
                    if ((not col_fp.exists()) or (not colmap)) and status == "OK":
                        status = "NO_COLNUMBERING"

                    fg_seq_trim, _hit_trim = find_seq_by_id(pep_trim, gid_norm)
                    fg_aa_trim = ""
                    if fg_seq_trim is None:
                        if status == "OK":
                            status = "NO_HIT_IN_TRIMAL"
                    else:
                        fg_aa_trim = aln_char_at(fg_seq_trim, site)

                    fg_seq_mafft, _hit_mafft = find_seq_by_id(pep_mafft, gid_norm)
                    fg_aa_mafft = ""
                    residue_index_raw = ""
                    raw_aa = ""
                    mafft_col = ""

                    if mafft_col1 is None:
                        if status == "OK":
                            status = "SITE_OUT_OF_COLMAP"
                    else:
                        mafft_col = str(mafft_col1)  # 输出给人看的列号用 1-based
                        if fg_seq_mafft is None:
                            if status == "OK":
                                status = "NO_HIT_IN_MAFFT"
                        else:
                            fg_aa_mafft = aln_char_at(fg_seq_mafft, mafft_col1)
                            ridx = residue_index_from_alignment(fg_seq_mafft, mafft_col1)
                            if ridx is None:
                                if status == "OK":
                                    status = "GAP_AT_SITE"
                            else:
                                residue_index_raw = str(ridx)
                                if proteome_index is not None and gid_norm in proteome_index:
                                    raa = get_raw_aa_at(proteome_index[gid_norm], ridx)
                                    raw_aa = "" if raa is None else raa

                    # match_input_ref：输入 aa 是否等于 ref_aa（仅供参考）
                    match_input_ref = "NA"
                    if input_aa and ref_aa and ref_aa not in {"-", "."}:
                        match_input_ref = "TRUE" if input_aa.upper() == ref_aa.upper() else "FALSE"

                    # match_raw_fg：raw_aa 是否等于 fg_aa_mafft（闭环核对）
                    match_raw_fg = "NA"
                    if raw_aa and fg_aa_mafft and fg_aa_mafft not in {"-", "."}:
                        match_raw_fg = "TRUE" if raw_aa.upper() == fg_aa_mafft.upper() else "FALSE"

                    if status == "GAP_AT_SITE":
                        stat_gap += 1
                    if status == "OK":
                        stat_ok += 1
                    if status == "OK" and match_raw_fg == "TRUE" and residue_index_raw and raw_aa:
                        stat_ok_match += 1

                    # 写长表一行
                    row = {
                        "OG": og,
                        "foreground": fg,
                        "site": str(site),
                        "post": post_fmt,  # 长表也输出格式化后的 post，便于一致阅读
                        "gene_id": gid_norm,
                        "ref_aa": ref_aa,
                        "input_aa": input_aa,
                        "fg_aa_trim": fg_aa_trim,
                        "mafft_col": mafft_col,
                        "fg_aa_mafft": fg_aa_mafft,
                        "residue_index_raw": residue_index_raw,
                        "raw_aa": raw_aa,
                        "match_input_ref": match_input_ref,
                        "match_raw_fg": match_raw_fg,
                        "status": status,
                        "pep_trimal_fa": relpath(pep_trim_fp),
                        "pep_mafft_fa": relpath(pep_mafft_fp),
                        "colnumbering_txt": relpath(col_fp),
                    }
                    w.write("\t".join(row.get(c, "") for c in mapped_cols) + "\n")
                    n_out += 1

                    # 聚合精简表（只收“可标 3D 的有效位点”）
                    if (
                        status == "OK"
                        and match_raw_fg == "TRUE"
                        and residue_index_raw
                        and raw_aa
                    ):
                        try:
                            ridx_int = int(residue_index_raw)
                        except Exception:
                            ridx_int = -1
                        if ridx_int > 0:
                            key = (og, fg, gid_norm)
                            gene_sites.setdefault(key, []).append(
                                (ridx_int, raw_aa, float(post_raw) if post_float == post_float else float("nan"), post_fmt)
                            )

                if PROGRESS_EVERY and n_in % int(PROGRESS_EVERY) == 0:
                    print(f"[PROGRESS] lines={n_in} mapped_rows={n_out} gene_groups={len(gene_sites)}", flush=True)

    # 写精简表
    gene_cols = ["OG", "foreground", "gene_id", "n_sites", "max_post", "sites_3d_label"]
    with OUTPUT_GENE_TSV.open("w", encoding="utf-8") as wg:
        wg.write("\t".join(gene_cols) + "\n")

        # 为了输出稳定：按 OG, foreground, gene_id 排序
        for (og, fg, gid) in sorted(gene_sites.keys()):
            items = gene_sites[(og, fg, gid)]
            # 去重：同一个 residue_index_raw 可能因重复行出现，保留最高 post
            best_by_ridx: Dict[int, Tuple[str, float, str]] = {}
            for ridx, aa, postf, postfmt in items:
                if ridx not in best_by_ridx:
                    best_by_ridx[ridx] = (aa, postf, postfmt)
                else:
                    _, old_postf, old_postfmt = best_by_ridx[ridx]
                    # 选择更大的 posterior（若遇到 NaN，优先保留非 NaN）
                    if (old_postf != old_postf) and (postf == postf):
                        best_by_ridx[ridx] = (aa, postf, postfmt)
                    elif (postf == postf) and (old_postf == old_postf) and (postf > old_postf):
                        best_by_ridx[ridx] = (aa, postf, postfmt)
                    else:
                        # 保留原有
                        best_by_ridx[ridx] = (best_by_ridx[ridx][0], old_postf, old_postfmt)

            # 按 residue_index_raw 升序
            ridx_sorted = sorted(best_by_ridx.keys())
            tokens: List[str] = []
            max_post = float("-inf")
            for ridx in ridx_sorted:
                aa, postf, postfmt = best_by_ridx[ridx]
                tokens.append(f"{aa}{ridx}({postfmt})")
                if postf == postf:  # 非 NaN
                    if postf > max_post:
                        max_post = postf

            n_sites = len(tokens)
            max_post_str = ""
            if max_post == float("-inf"):
                max_post_str = ""
            else:
                max_post_str = fmt_post(str(max_post), POST_DECIMALS)

            wg.write("\t".join([
                og, fg, gid,
                str(n_sites),
                max_post_str,
                ";".join(tokens),
            ]) + "\n")

    print("[DONE] Wrote:", relpath(OUTPUT_MAPPED_TSV), flush=True)
    print("[DONE] Wrote:", relpath(OUTPUT_GENE_TSV), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

