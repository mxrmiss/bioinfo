#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
11_build_psg_subset_package.py —— 用物种子集重建 aphylo_ready 发布包（覆盖 10 脚本产物版）

核心功能：
  1) 从 OrthoFinder 的 Orthogroups.tsv 中，根据 config.positive_selection.species_subset
     筛选“对子集物种严格单拷贝（strict SCO）”的 OG 列表；
  2) 为这些 OG 构建只包含子集物种的 AA 对齐（优先复用 species-only 对齐，缺失时从 OF MSA+trimal 重建）；
  3) 按照 10_publish_aphylo_ready.py 的发布契约，额外发布 strict_sco_msa/OGxxxx.raw.fa：
     - 来源是 OrthoFinder 原始 MSA；
     - 只保留子集物种；
     - header 统一转成 Gene/ProteinID（即 Species|ID 的右侧部分）；
  4) 基于 strict_sco_msa，生成列掩码（colmask/OGXXXX.colmask），长度=对齐列数，全为 '1'；
  5) 从全物种树 species_tree.nwk 剪枝得到子集物种树 species_tree.nwk（文件名保持不变）；
  6) 直接利用 Orthogroups.tsv 中的 protein_id + pep2cds.tsv，构建 pep2cds_resolved.tsv；
     未能解析的条目写入 pep2cds_unresolved.tsv，但不会终止脚本；
  7) 基于 Orthogroups.GeneCount.tsv 生成子集物种的 family.tsv；
  8) 输出 manifest.json + QC_report.txt + .done 哨兵。

重要行为：
  * 发布目录直接使用 config.paths.publish_dir（通常是 results/publish/aphylo_ready）。
  * 若 OVERWRITE_PSG_DIR=True，则会先删除原有 publish_dir，再重建一份“子集物种版”的 aphylo_ready。
  * 本脚本现在与 10 脚本保持相同的 raw.fa 发布逻辑，因此无需修改 aphylo/03。
"""

from __future__ import annotations

import sys
import re
import csv
import json
import yaml
import shutil
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Set


# ====================== 顶部可配置（不走命令行） ======================

DEFAULT_CONFIG = "config.yaml"
LOG_DIR = "logs"
LOG_FILE = "11_build_psg_subset_package.log"
OVERWRITE_PSG_DIR = True
SKIP_OG_ON_TRIMAL_FAIL = True
DEFAULT_RAW_MSA_SUFFIX = ".raw.fa"


# ====================== 通用工具：配置与日志 ======================


def load_cfg(cfg_path: str = DEFAULT_CONFIG) -> dict:
    """读取 YAML 配置文件。"""
    p = Path(cfg_path)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 找不到配置文件：{p}")
    with p.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}



def need_file(p: Path, msg: str):
    """确保文件存在，否则抛错。"""
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{p} —— {msg}")



def need_dir(d: Path, msg: str):
    """确保目录存在，否则抛错。"""
    if not d.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{d} —— {msg}")



def init_logger(name: str, log_dir: Path) -> logging.Logger:
    """
    初始化日志：
      - INFO 级别；
      - 同时输出到屏幕（stdout）和 logs/LOG_FILE。
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)

    fh = logging.FileHandler(log_dir / LOG_FILE, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)

    logger.addHandler(sh)
    logger.addHandler(fh)
    logger.propagate = False
    return logger


# ====================== 物种名归一化与 OrthoFinder 锚点 ======================


def canon_species(s: str, alias_map: Dict[str, str]) -> str:
    """
    物种名统一：
      - strip；
      - 处理形如 "Sinonovacula_constricta_Sinonovacula_constricta" 的重复；
      - 应用 alias_map。
    """
    s = s.strip()
    if not s:
        return s

    for _ in range(3):
        toks = s.split('_')
        if len(toks) % 2 != 0 or not toks:
            break
        half = len(toks) // 2
        left = '_'.join(toks[:half])
        right = '_'.join(toks[half:])
        if left == right and left:
            s = left
        else:
            break

    if alias_map:
        s = alias_map.get(s, s)
    return s



def read_results_anchor(of_root: Path) -> Path:
    """
    根据 RESULTS_DIR.txt 或唯一的 Results_* 子目录，确定 OrthoFinder 结果目录。
    优先读取 RESULTS_DIR.txt；否则在 of_root 下寻找唯一 Results_* 目录。
    """
    anchor = of_root / "RESULTS_DIR.txt"
    if anchor.is_file():
        txt = anchor.read_text(encoding="utf-8").strip()
        p = Path(txt)
        if p.is_dir():
            return p

    cand = sorted(of_root.glob("Results_*"))
    cand = [c for c in cand if c.is_dir()]
    if len(cand) == 1:
        return cand[0]

    raise FileNotFoundError(
        "[ERR] 无法确定唯一的 OrthoFinder 结果目录："
        "既没有有效的 RESULTS_DIR.txt，也没有唯一一个 Results_* 目录"
    )


# ====================== FASTA 读写 ======================


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """简单 FASTA 读取：返回[(header, seq)]，header 不含 '>'。"""
    records: List[Tuple[str, str]] = []
    name = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line:
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts)))
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq_parts)))
    return records



def write_fasta(records: List[Tuple[str, str]], path: Path, linewrap: int = 60):
    """FASTA 写出：records 为 [(header, seq)]；header 不含 '>'。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for h, s in records:
            f.write(f">{h}\n")
            if linewrap and linewrap > 0:
                for i in range(0, len(s), linewrap):
                    f.write(s[i:i + linewrap] + "\n")
            else:
                f.write(s + "\n")


# ====================== Orthogroups 子集严格 SCO 筛选 ======================


def parse_orthogroups_subset_strict_sco(
    tsv: Path,
    subset_species: List[str],
    alias_map: Dict[str, str],
    logger: logging.Logger
) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    从 Orthogroups.tsv 中筛选对子集物种严格单拷贝（strict SCO）的 OG。

    条件：
      - 对 subset_species 中每个物种：
          * 对应列非空；
          * 拆分 ',' / 空白 后恰好 1 条蛋白 ID。
      - 其它物种列不做约束。

    返回：
      - og_list: 满足条件的 OG ID 列表；
      - og_to_sp_pid: {OG: {Species_std: protein_id}}（Species 为归一化后物种名）。
    """
    need_file(tsv, "OrthoFinder Orthogroups.tsv 缺失")
    subset_canon = [canon_species(s, alias_map) for s in subset_species]
    subset_canon_set = set(subset_canon)

    logger.info("[S1] 从 Orthogroups.tsv 筛选对子集严格 SCO 的 OG ...")
    logger.info("     子集物种（归一化后）：%s", ", ".join(sorted(subset_canon_set)))

    og_to_sp_pid: Dict[str, Dict[str, str]] = {}
    ogs: List[str] = []

    with tsv.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        if not header or header[0] != "Orthogroup":
            raise RuntimeError(f"[ERR] Orthogroups.tsv 表头异常（首列应为 Orthogroup）：{tsv}")

        name_to_col: Dict[str, int] = {}
        for idx, raw_name in enumerate(header[1:], start=1):
            sp_std = canon_species(raw_name, alias_map)
            if not sp_std:
                continue
            name_to_col[sp_std] = idx

        missing = [s for s in subset_canon if s not in name_to_col]
        if missing:
            raise RuntimeError(
                "[ERR] 子集物种在 Orthogroups.tsv 列名中不存在：" + ", ".join(missing)
            )

        n_total = 0
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            parts = line.split("\t")
            og = parts[0].strip()
            if not og:
                continue

            ok = True
            sp_pid_map: Dict[str, str] = {}

            for sp in subset_canon:
                col_idx = name_to_col[sp]
                if col_idx >= len(parts):
                    ok = False
                    break
                cell = parts[col_idx].strip()
                if not cell:
                    ok = False
                    break
                tokens = [x for x in re.split(r"[,\s]+", cell) if x]
                if len(tokens) != 1:
                    ok = False
                    break
                sp_pid_map[sp] = tokens[0]

            if ok:
                ogs.append(og)
                og_to_sp_pid[og] = sp_pid_map

    logger.info("[S1] Orthogroups.tsv 总行数（不含表头）：%d", n_total)
    logger.info("[S1] 子集严格 SCO OG 数量：%d", len(ogs))
    return ogs, og_to_sp_pid


# ====================== 构建子集 MSA ======================


def build_subset_msas(
    ogs: List[str],
    subset_species: List[str],
    alias_map: Dict[str, str],
    species_only_dir: Path,
    of_msa_dir: Path,
    of_msa_suffix: str,
    trimal_bin: str,
    out_msa_dir: Path,
    logger: logging.Logger
) -> List[str]:
    """
    为给定 OG 列表构建只包含子集物种的 AA 对齐：

      1) 若存在 species-only 对齐（paths.species_collapse_dir/OG*.trim.faa）：
         - 读取后只保留 subset 物种对应的序列（header 为物种名），写入 out_msa_dir；
      2) 否则：
         - 从 OrthoFinder 多物种对齐（MultipleSequenceAlignments/OG*.fa）读取；
         - 调用 trimal 自动裁剪；
         - 根据 header 中的物种名（假设为 Species|SeqID）筛选 subset 物种；
         - 写入 out_msa_dir。

    返回：最终成功构建的 OG 列表（可能少于输入 ogs）。
    """
    subset_canon = [canon_species(s, alias_map) for s in subset_species]
    subset_set = set(subset_canon)

    out_msa_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_msa_dir.parent / "tmp_psg_msa"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    need_dir(of_msa_dir, "OrthoFinder MultipleSequenceAlignments 目录缺失")
    if species_only_dir.exists():
        logger.info("[S2] species-only 对齐目录存在，将优先复用：%s", species_only_dir)
    else:
        logger.info("[S2] species-only 对齐目录不存在：%s —— 将全部从 OF MSA 重建", species_only_dir)

    n_from_species_only = 0
    n_from_of_raw = 0
    kept_ogs: List[str] = []

    for og in ogs:
        out_path = out_msa_dir / f"{og}.trim.faa"
        if out_path.is_file():
            kept_ogs.append(og)
            continue

        src_species = species_only_dir / f"{og}.trim.faa"
        if src_species.is_file():
            recs = read_fasta(src_species)
            sp_to_seq: Dict[str, str] = {}
            for h, seq in recs:
                sp = canon_species(h.split("|", 1)[0], alias_map)
                if sp in subset_set:
                    sp_to_seq[sp] = seq
            missing = subset_set - set(sp_to_seq.keys())
            if missing:
                logger.warning(
                    "[S2] OG %s 在 species-only 对齐中缺少子集物种：%s —— 跳过该 OG",
                    og, ", ".join(sorted(missing))
                )
                continue
            out_recs = [(sp, sp_to_seq[sp]) for sp in subset_canon]
            write_fasta(out_recs, out_path)
            n_from_species_only += 1
            kept_ogs.append(og)
            continue

        src_msa = of_msa_dir / f"{og}{of_msa_suffix}"
        if not src_msa.is_file():
            logger.warning("[S2] 缺少 OG %s 的 OF MSA：%s —— 跳过该 OG", og, src_msa)
            continue

        tmp_trim = tmp_dir / f"{og}.trim.faa"
        cmd = [trimal_bin, "-in", str(src_msa), "-out", str(tmp_trim), "-automated1"]
        logger.info("[S2] 调 trimal 裁剪 OG %s：%s", og, " ".join(cmd))
        r = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if r.returncode != 0 or (not tmp_trim.is_file()) or tmp_trim.stat().st_size == 0:
            logger.error("[S2] trimal 失败，OG=%s；stderr:\n%s", og, r.stderr)
            if not SKIP_OG_ON_TRIMAL_FAIL:
                raise RuntimeError(f"[ERR] trimal 失败：OG={og}")
            logger.warning("[S2] 由于配置允许，将跳过 OG %s", og)
            continue

        recs = read_fasta(tmp_trim)
        sp_to_seq: Dict[str, str] = {}
        for h, seq in recs:
            token = h.split("|", 1)[0]
            sp = canon_species(token, alias_map)
            if sp in subset_set and sp not in sp_to_seq:
                sp_to_seq[sp] = seq

        missing = subset_set - set(sp_to_seq.keys())
        if missing:
            logger.warning(
                "[S2] trimal 后 OG %s 的对齐缺少子集物种：%s —— 跳过该 OG",
                og, ", ".join(sorted(missing))
            )
            continue

        out_recs = [(sp, sp_to_seq[sp]) for sp in subset_canon]
        write_fasta(out_recs, out_path)
        n_from_of_raw += 1
        kept_ogs.append(og)

    logger.info(
        "[S2] 子集 MSA 构建完成：共请求 %d 个 OG，成功 %d 个；其中复用 species-only=%d，OF+trimal=%d",
        len(ogs), len(kept_ogs), n_from_species_only, n_from_of_raw
    )
    return kept_ogs


# ====================== 发布 raw.fa（与 10 脚本契约一致） ======================


def _sanitize_id_contract(s: str) -> str:
    """
    与 10 脚本保持一致的最小化 ID 规整：
      - 去首尾空白；
      - 只取首 token；
      - 去掉外围多余引号。
    不改动内部字符，以免破坏真实 protein_id。
    """
    s = s.strip().split()[0] if s.strip() else ""
    s = s.strip("\"'")
    return s



def find_of_raw_msa_file(msa_dir: Path, ogid: str) -> Path:
    """在 OrthoFinder 的 MultipleSequenceAlignments 目录中定位 OG 的原始 AA MSA 文件。"""
    preferred = [f"{ogid}.fa", f"{ogid}.fasta", f"{ogid}.faa", f"{ogid}.aln.fa", f"{ogid}.aln.fasta"]
    for name in preferred:
        p = msa_dir / name
        if p.is_file():
            return p

    cands = sorted([p for p in msa_dir.glob(f"{ogid}*") if p.is_file()])
    if len(cands) == 1:
        return cands[0]
    if not cands:
        raise FileNotFoundError(f"[ERR] 未找到 OrthoFinder 原始 MSA：{msa_dir}/{ogid}*")

    common_ext = {".fa", ".faa", ".fasta"}
    cands2 = [p for p in cands if p.suffix.lower() in common_ext]
    if len(cands2) == 1:
        return cands2[0]

    raise RuntimeError(
        f"[ERR] OrthoFinder 原始 MSA 存在多个候选，无法唯一判定：{ogid}\n"
        f"      candidates={','.join([x.name for x in cands[:10]])}"
    )



def normalize_of_msa_headers_to_geneid_subset(
    of_msa: Path,
    subset_species: List[str],
    alias_map: Dict[str, str],
) -> List[Tuple[str, str]]:
    """
    读取 OrthoFinder 原始 MSA，仅保留子集物种，并把 header 从 Species|Gene/ProteinID
    转成 Gene/ProteinID。返回的记录顺序严格按 subset_species 排列。

    这里故意保持和 10 脚本相同的思路：
      - 来源必须是 OrthoFinder 原始 MSA，而不是 trim.faa 改名；
      - header 必须转成右侧的真实 Gene/ProteinID，供 aphylo 03 与 CDS 精确对应。
    """
    subset_canon = [canon_species(s, alias_map) for s in subset_species]
    subset_set = set(subset_canon)

    recs = read_fasta(of_msa)
    if not recs:
        raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 为空：{of_msa}")

    sp_to_gid_seq: Dict[str, Tuple[str, str]] = {}
    seen_species: Set[str] = set()

    for h, seq in recs:
        tok = h.split()[0]
        if "|" not in tok:
            raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA header 不含 '|'（无法解析 Species|ID）：{tok} —— {of_msa}")
        sp_raw, gid = tok.split("|", 1)
        sp = canon_species(sp_raw, alias_map)

        if sp not in subset_set:
            continue
        if sp in seen_species:
            raise RuntimeError(f"[ERR] strict-SCO 冲突：同一物种在原始 MSA 中出现多条序列：{sp} —— {of_msa}")
        seen_species.add(sp)

        gid = _sanitize_id_contract(gid)
        if not gid:
            raise RuntimeError(f"[ERR] Gene/ProteinID 为空：{tok} —— {of_msa}")
        sp_to_gid_seq[sp] = (gid, seq)

    missing = [sp for sp in subset_canon if sp not in sp_to_gid_seq]
    if missing:
        raise RuntimeError(
            f"[ERR] OrthoFinder 原始 MSA 缺少子集物种：{of_msa}\n"
            f"      missing={','.join(missing)}"
        )

    out_recs: List[Tuple[str, str]] = []
    lengths: Set[int] = set()
    seen_gid: Set[str] = set()
    for sp in subset_canon:
        gid, seq = sp_to_gid_seq[sp]
        if gid in seen_gid:
            raise RuntimeError(f"[ERR] Gene/ProteinID 重复：{gid} —— {of_msa}")
        seen_gid.add(gid)
        out_recs.append((gid, seq))
        lengths.add(len(seq))

    if len(lengths) != 1:
        raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 非等长对齐：{of_msa}")

    return out_recs



def publish_subset_raw_msa(
    ogs: List[str],
    subset_species: List[str],
    alias_map: Dict[str, str],
    of_msa_dir: Path,
    out_msa_dir: Path,
    raw_suffix: str,
    logger: logging.Logger,
) -> int:
    """
    为 kept_ogs 额外发布子集版 raw.fa。
    产物路径：strict_sco_msa/OGxxxx.raw.fa
    """
    need_dir(of_msa_dir, "OrthoFinder MultipleSequenceAlignments 缺失")

    n_published = 0
    for og in ogs:
        of_raw = find_of_raw_msa_file(of_msa_dir, og)
        raw_recs = normalize_of_msa_headers_to_geneid_subset(
            of_raw,
            subset_species,
            alias_map,
        )
        out_raw = out_msa_dir / f"{og}{raw_suffix}"
        write_fasta(raw_recs, out_raw, linewrap=60)
        n_published += 1

    logger.info("[S2b] 已发布 raw AA MSA：%d 个 OG（后缀=%s）", n_published, raw_suffix)
    return n_published


# ====================== 基于 MSA 构建 colmask ======================


def build_colmask_for_msas(
    msa_dir: Path,
    colmask_dir: Path,
    logger: logging.Logger
):
    """
    根据 strict_sco_msa 里的对齐，为每个 OG 生成一份列掩码：
      - 输入：msa_dir 下的 OG*.trim.faa；
      - 输出：colmask_dir/OGXXXX.colmask；
      - 掩码内容：长度 = 对齐列数（单条序列长度），全为 '1'。
    """
    colmask_dir.mkdir(parents=True, exist_ok=True)

    n_files = 0
    for faa in sorted(msa_dir.glob("OG*.trim.faa")):
        ogid = faa.stem.split(".")[0]
        recs = read_fasta(faa)
        if not recs:
            raise RuntimeError(f"[ERR] 空对齐文件：{faa}")

        first_seq = recs[0][1]
        L = len(first_seq)
        if L <= 0:
            raise RuntimeError(f"[ERR] 对齐长度为 0：{faa}")

        bits = "1" * L
        out_path = colmask_dir / f"{ogid}.colmask"
        out_path.write_text(bits + "\n", encoding="utf-8")
        n_files += 1

    (colmask_dir / ".done").touch()
    logger.info("[S3] colmask 已生成：%d 个 OG（目录：%s）", n_files, colmask_dir)


# ====================== 剪枝物种树 ======================


def prune_species_tree_to_subset(
    src_tree: Path,
    subset_species: List[str],
    alias_map: Dict[str, str],
    out_tree: Path,
    logger: logging.Logger
):
    """
    从全物种树中剪枝得到只包含子集物种的树。
    """
    from Bio import Phylo

    need_file(src_tree, "species_tree.nwk 缺失（全物种树）")

    subset_canon = [canon_species(s, alias_map) for s in subset_species]
    subset_set = set(subset_canon)

    logger.info("[S4] 从 %s 剪枝构建子集物种树 ...", src_tree)
    tree = Phylo.read(str(src_tree), "newick")

    present: Set[str] = set()
    for term in tree.get_terminals():
        sp_std = canon_species(term.name, alias_map)
        if sp_std in subset_set:
            term.name = sp_std
            present.add(sp_std)

    missing = subset_set - present
    if missing:
        logger.warning("[S4] 子集物种在物种树中未找到：%s", ", ".join(sorted(missing)))

    for term in list(tree.get_terminals()):
        sp_std = canon_species(term.name, alias_map)
        if sp_std not in subset_set:
            tree.prune(term)

    Phylo.write(tree, str(out_tree), "newick")
    logger.info("[S4] 子集物种树已写出：%s", out_tree)


# ====================== pep2cds 子集映射 ======================


def _strip_dbprefix(s: str) -> str:
    """去掉 NCBI 类前缀（gb|/ref| 等），保留后半段。"""
    _DB_PREFIXES_REGEX = r"(?:lcl|gb|emb|dbj|ref|sp|tr|pir|prf|pdb|pat|pgp|bbs|bbm|gim|gi|tpg|tpe|tpd|tpa|gnl)"
    return re.sub(rf"^{_DB_PREFIXES_REGEX}\|", "", s)



def _norm_first_token(s: str) -> str:
    """对首 token 做 ':' 与 '_' 的规范化。"""
    head, *rest = s.split("|", 1)
    head = re.sub(r"([A-Za-z]{2,}[A-Za-z0-9]*)[:_]+([A-Za-z0-9]{2,})", r"\1_\2", head)
    return head if not rest else head + "|" + rest[0]



def _norm_for_match(h: str) -> str:
    """结合去前缀与首 token 规范化，用于 pep2cds 映射中多键匹配。"""
    return _norm_first_token(_strip_dbprefix(h.strip()))



def read_pep2cds_global(tsv: Path, alias_map: Dict[str, str]) -> Dict[Tuple[str, str], str]:
    """
    读取 pep2cds 全局映射。
    """
    mp: Dict[Tuple[str, str], str] = {}

    def add_keys(sp: str, pid: str, cds: str):
        mp.setdefault((sp, pid), cds)
        pid_base = _norm_for_match(pid)
        mp.setdefault((sp, pid_base), cds)
        pid_u = pid.replace("|", "_")
        mp.setdefault((sp, pid_u), cds)
        mp.setdefault((sp, re.sub(r"_+", "_", pid_u)), cds)
        pid_base_u = pid_base.replace("|", "_")
        mp.setdefault((sp, pid_base_u), cds)
        mp.setdefault((sp, re.sub(r"_+", "_", pid_base_u)), cds)
        if "|" in pid:
            tail = pid.rsplit("|", 1)[-1]
            mp.setdefault((sp, tail), cds)
            tail_u = tail.replace("|", "_")
            mp.setdefault((sp, tail_u), cds)
            mp.setdefault((sp, re.sub(r"_+", "_", tail_u)), cds)

    need_file(tsv, "pep2cds 汇总表缺失")
    with tsv.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, [])
        hl = [x.lower() for x in header]
        try:
            i_sp = hl.index("species")
            i_pid = hl.index("protein_id")
            i_cds = hl.index("cds_id")
        except ValueError:
            if len(header) >= 3:
                i_sp, i_pid, i_cds = 0, 1, 2
            else:
                raise RuntimeError(f"[ERR] {tsv} 格式异常：既无表头也不足三列")

        for row in reader:
            if not row or len(row) < 3:
                continue
            sp = canon_species(row[i_sp].strip(), alias_map)
            pid = row[i_pid].strip()
            cds = row[i_cds].strip()
            if sp and pid and cds:
                add_keys(sp, pid, cds)

    if not mp:
        raise RuntimeError(f"[ERR] 映射为空：{tsv}")
    return mp



def build_pep2cds_subset(
    pub_dir: Path,
    pep2cds_tsv: Path,
    og_list: List[str],
    subset_species: List[str],
    alias_map: Dict[str, str],
    og_to_sp_pid: Dict[str, Dict[str, str]],
    logger: logging.Logger
) -> Tuple[int, int]:
    """
    构建 PSG 子集版的 pep2cds_resolved.tsv。
    """
    subset_canon = [canon_species(s, alias_map) for s in subset_species]
    sp_pid_to_cds = read_pep2cds_global(pep2cds_tsv, alias_map)
    logger.info("[S5] pep2cds 全局映射键数量：%d", len(sp_pid_to_cds))

    out_resolved = pub_dir / "pep2cds_resolved.tsv"
    out_unres = pub_dir / "pep2cds_unresolved.tsv"

    unresolved: List[Tuple[str, str, str, str]] = []
    n_out = 0

    with out_resolved.open("w", encoding="utf-8", newline="") as wf:
        writer = csv.writer(wf, delimiter="\t")
        writer.writerow(["OG", "Species", "protein_id", "cds_id"])

        for og in og_list:
            sp_pid_map = og_to_sp_pid.get(og, {})
            for sp in subset_canon:
                pid = sp_pid_map.get(sp)
                if not pid:
                    unresolved.append((og, sp, "", "missing protein_id in Orthogroups.tsv"))
                    continue

                cds = None
                for x in (pid, _norm_for_match(pid), pid.rsplit("|", 1)[-1] if "|" in pid else None):
                    if not x:
                        continue
                    cds = (
                        sp_pid_to_cds.get((sp, x))
                        or sp_pid_to_cds.get((sp, x.replace("|", "_")))
                        or sp_pid_to_cds.get((sp, re.sub(r"_+", "_", x.replace("|", "_"))))
                    )
                    if cds:
                        break

                if cds is None:
                    unresolved.append((og, sp, pid, "protein→cds not found in pep2cds.tsv"))
                    continue

                writer.writerow([og, sp, pid, cds])
                n_out += 1

    n_unres = len(unresolved)
    if n_unres > 0:
        with out_unres.open("w", encoding="utf-8", newline="") as wf:
            writer = csv.writer(wf, delimiter="\t")
            writer.writerow(["OG", "Species", "protein_id", "detail"])
            for og, sp, pid, detail in unresolved:
                writer.writerow([og, sp, pid, detail])
        logger.warning("[S5] pep2cds_unresolved.tsv 条目：%d —— 不终止脚本，但建议人工检查", n_unres)
    else:
        logger.info("[S5] 所有 (OG,Species) 条目均成功解析到 cds_id")

    return n_out, n_unres


# ====================== family.tsv 构建 ======================


def build_family_subset(
    gene_count_src: Path,
    og_list: List[str],
    subset_species: List[str],
    alias_map: Dict[str, str],
    out_tsv: Path,
    logger: logging.Logger
):
    """
    基于 Orthogroups.GeneCount.tsv 构建子集物种的 family.tsv。
    """
    need_file(gene_count_src, "Orthogroups.GeneCount.tsv 缺失")
    og_set = set(og_list)
    subset_canon = [canon_species(s, alias_map) for s in subset_species]
    subset_set = set(subset_canon)

    with gene_count_src.open("r", encoding="utf-8") as f_in, out_tsv.open("w", encoding="utf-8", newline="") as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        header = next(reader, [])
        if not header or header[0] != "Orthogroup":
            raise RuntimeError(f"[ERR] Orthogroups.GeneCount.tsv 表头异常：{gene_count_src}")

        col_idx: List[int] = [0]
        new_header: List[str] = ["Orthogroup"]

        for idx, raw_name in enumerate(header[1:], start=1):
            sp_std = canon_species(raw_name, alias_map)
            if sp_std in subset_set:
                col_idx.append(idx)
                new_header.append(sp_std)

        missing = subset_set - set(new_header[1:])
        if missing:
            logger.warning("[S6] 子集物种在 Orthogroups.GeneCount.tsv 中未找到：%s", ", ".join(sorted(missing)))

        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(new_header)

        n_rows = 0
        n_kept = 0
        for row in reader:
            n_rows += 1
            if not row:
                continue
            og = row[0].strip()
            if og not in og_set:
                continue
            n_kept += 1
            out_row = [row[i] if i < len(row) else "" for i in col_idx]
            writer.writerow(out_row)

    logger.info("[S6] family.tsv 构建完成：原始行数=%d，保留行数(OG in PSG集合)=%d，输出文件=%s", n_rows, n_kept, out_tsv)


# ====================== 主流程 ======================


def main():
    cfg = load_cfg()
    paths = cfg.get("paths", {}) or {}
    bins = cfg.get("binaries", {}) or {}
    species_cfg = cfg.get("species", {}) or {}
    alias_map = species_cfg.get("alias_map", {}) or {}
    ps_cfg = cfg.get("positive_selection", {}) or {}
    publish_cfg = cfg.get("publish", {}) or {}

    log_dir = Path(LOG_DIR)
    logger = init_logger("psg_subset", log_dir)

    logger.info("============================================================")
    logger.info("11_build_psg_subset_package —— 用子集物种重建 aphylo_ready 发布包")
    logger.info("============================================================")

    subset_species = ps_cfg.get("species_subset") or []
    if not subset_species:
        raise RuntimeError("[ERR] positive_selection.species_subset 为空，请在 config.yaml 中配置子集物种列表")

    pub_dir = Path(paths.get("publish_dir", "results/publish/aphylo_ready")).resolve()
    logger.info("[CFG] PSG 子集物种原始列表：%s", ", ".join(subset_species))
    logger.info("[CFG] 发布目录（将覆盖原有 aphylo_ready）：%s", pub_dir)

    of_root = Path(paths.get("orthofinder_results_dir", "results/orthofinder"))
    reports_dir = Path(paths.get("reports_dir", "results/reports"))
    species_only_dir = Path(paths.get("species_collapse_dir", "results/trim_norm_collapse_sco"))
    trees_dir = Path(paths.get("trees_dir", "results/trees"))
    species_tree = trees_dir / "species_tree.nwk"
    pep2cds_tsv = Path(paths.get("pep2cds_tsv", "data/maps/pep2cds.tsv"))

    of_results_dir = read_results_anchor(of_root)
    orthogroups_tsv = of_results_dir / "Orthogroups" / "Orthogroups.tsv"
    gene_count_src = of_results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    of_msa_dir = of_results_dir / "MultipleSequenceAlignments"

    of_msa_suffix = cfg.get("input", {}).get("of_msa_suffix", ".fa")
    trimal_bin = bins.get("trimal", "trimal")
    raw_suffix = str(publish_cfg.get("raw_msa_suffix", DEFAULT_RAW_MSA_SUFFIX)) or DEFAULT_RAW_MSA_SUFFIX

    if OVERWRITE_PSG_DIR and pub_dir.exists():
        logger.info("[INIT] 删除已有发布目录（原 10 脚本产物）：%s", pub_dir)
        shutil.rmtree(pub_dir)
    pub_dir.mkdir(parents=True, exist_ok=True)
    out_msa_dir = pub_dir / "strict_sco_msa"
    out_msa_dir.mkdir(parents=True, exist_ok=True)
    colmask_dir = pub_dir / "colmask"

    ogs_psg, og_to_sp_pid = parse_orthogroups_subset_strict_sco(
        orthogroups_tsv, subset_species, alias_map, logger
    )
    if not ogs_psg:
        raise RuntimeError("[FATAL] 未筛选到任何对子集严格 SCO 的 OG，无法继续构建发布包")

    kept_ogs = build_subset_msas(
        ogs_psg,
        subset_species,
        alias_map,
        species_only_dir,
        of_msa_dir,
        of_msa_suffix,
        trimal_bin,
        out_msa_dir,
        logger
    )
    if not kept_ogs:
        raise RuntimeError("[FATAL] 子集 MSA 构建后，成功的 OG 数量为 0，无法继续")

    with (pub_dir / "ogs_selected.list").open("w", encoding="utf-8") as f:
        for og in kept_ogs:
            f.write(og + "\n")
    with (pub_dir / "sco_filelist.txt").open("w", encoding="utf-8") as f:
        for og in kept_ogs:
            f.write(str((out_msa_dir / f"{og}.trim.faa").resolve()) + "\n")
    logger.info("[S2] ogs_selected.list & sco_filelist.txt 已写出")

    raw_published = publish_subset_raw_msa(
        ogs=kept_ogs,
        subset_species=subset_species,
        alias_map=alias_map,
        of_msa_dir=of_msa_dir,
        out_msa_dir=out_msa_dir,
        raw_suffix=raw_suffix,
        logger=logger,
    )

    build_colmask_for_msas(
        msa_dir=out_msa_dir,
        colmask_dir=colmask_dir,
        logger=logger
    )

    subset_tree = pub_dir / "species_tree.nwk"
    prune_species_tree_to_subset(
        species_tree,
        subset_species,
        alias_map,
        subset_tree,
        logger
    )

    n_resolved, n_unresolved = build_pep2cds_subset(
        pub_dir,
        pep2cds_tsv,
        kept_ogs,
        subset_species,
        alias_map,
        og_to_sp_pid,
        logger
    )

    family_psg = pub_dir / "family.tsv"
    build_family_subset(
        gene_count_src,
        kept_ogs,
        subset_species,
        alias_map,
        family_psg,
        logger
    )

    manifest = {
        "psg_subset_species_raw": subset_species,
        "psg_subset_species_canonical": [canon_species(s, alias_map) for s in subset_species],
        "n_ogs_strict_sco_subset": len(ogs_psg),
        "n_ogs_msa_kept": len(kept_ogs),
        "n_raw_msa_published": raw_published,
        "raw_msa_suffix": raw_suffix,
        "publish_dir": str(pub_dir),
        "orthofinder_results_dir": str(of_results_dir),
        "orthogroups_tsv": str(orthogroups_tsv),
        "gene_count_tsv": str(gene_count_src),
        "of_msa_dir": str(of_msa_dir),
        "of_msa_suffix": of_msa_suffix,
        "species_tree_source": str(species_tree),
        "species_tree_psg": str(subset_tree),
        "pep2cds_tsv": str(pep2cds_tsv),
        "n_pep2cds_resolved": n_resolved,
        "n_pep2cds_unresolved": n_unresolved,
        "reports_dir": str(reports_dir),
    }
    with (pub_dir / "manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    with (pub_dir / "QC_report.txt").open("w", encoding="utf-8") as f:
        f.write(f"[INFO] PSG subset species (canonical): {', '.join(manifest['psg_subset_species_canonical'])}\n")
        f.write(f"[INFO] strict SCO OG（子集定义）数量：{manifest['n_ogs_strict_sco_subset']}\n")
        f.write(f"[INFO] 严格 SCO MSA（成功构建）数量：{manifest['n_ogs_msa_kept']}\n")
        f.write(f"[INFO] raw AA MSA（发布）数量：{manifest['n_raw_msa_published']}\n")
        f.write(f"[INFO] pep2cds_resolved.tsv 行数：{n_resolved}\n")
        f.write(f"[INFO] pep2cds_unresolved.tsv 条目数：{n_unresolved}\n")
        f.write(f"[INFO] family.tsv 路径：{family_psg}\n")
        f.write(f"[INFO] species_tree.nwk 路径：{subset_tree}\n")

    (pub_dir / ".done").touch()
    (out_msa_dir / ".done").touch()

    logger.info("[DONE] PSG 子集版 aphylo_ready 发布包构建完成：%s", pub_dir)
    logger.info(
        "[STAT] OG: %d（子集严格 SCO）→ %d（成功构建 MSA）→ %d（发布 raw.fa）；pep2cds_resolved: %d；unresolved: %d",
        len(ogs_psg), len(kept_ogs), raw_published, n_resolved, n_unresolved
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        sys.exit(1)

