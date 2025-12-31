#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_publish_aphylo_ready.py —— 发布到 aphylo 可直接消费的“五件套”，稳健严格版
要点：
  * Newick 叶名稳健解析（忽略内部节点/支持值/分支长度，支持引号）。
  * SeqID 归一化：统一剥离 NCBI 数据库前缀（gb|/ref|/sp|/tr|/gnl|...），仅规范第一个 '|' 前 token 的 ':' 与 '_'。
  * 前缀白名单检测：pep2cds 与 CDS FASTA 出现未知数据库前缀时立即报错并给出示例。
  * pep2cds → 真实 FASTA 头：FASTA 索引同时使用“整行 header”与“首 token”两类键，并生成规范化/尾段候选，避免 410 类 miss。
产物契约保持不变：strict_sco_msa/、pep2cds_resolved.tsv、family.tsv、species_tree.nwk、清单与报告等。
"""

from __future__ import annotations
import os, re, sys, csv, json, yaml, shutil
from pathlib import Path
from typing import Dict, Tuple, List, Set, Optional

# ====================== 顶部可配置（不走命令行） ======================
OVERWRITE_PUBLISH_DIR = True                   # True=重建发布目录；False=原地覆盖
DEFAULT_PEP2CDS_TSV   = "data/maps/pep2cds.tsv"
STRICT_NORMALIZE_CDS  = True                   # True=未能规范化为真实 FASTA 头则硬失败
DEFAULT_CDS_DIR       = "data/cds"
DEFAULT_CDS_SUFFIX    = ".fna"
DEFAULT_RAW_MSA_SUFFIX = ".raw.fa"


# 数据库前缀白名单（基于 NCBI SeqID 规范，可在 config.publish.db_prefix_allowlist 扩展）
ALLOWED_DB_PREFIXES: Set[str] = {
    "lcl","gb","emb","dbj","ref","sp","tr","pir","prf","pdb","pat","pgp",
    "bbs","bbm","gim","gi","tpg","tpe","tpd","tpa","gnl"
}
_DB_PREFIXES_REGEX = r"(?:lcl|gb|emb|dbj|ref|sp|tr|pir|prf|pdb|pat|pgp|bbs|bbm|gim|gi|tpg|tpe|tpd|tpa|gnl)"

# =======================

def load_cfg(cfg_path: str = "config.yaml") -> dict:
    with open(cfg_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def need_file(p: Path, msg: str):
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{p} —— {msg}")

def need_dir(d: Path, msg: str):
    if not d.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{d} —— {msg}")

def copy_if_exists(src: Path, dst: Path):
    if src.is_file():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

def canon_species(s: str, alias_map: Dict[str,str]) -> str:
    """物种名统一：strip + alias 映射 + 防止“重复物种名拼接”"""
    s = s.strip()
    if not s:
        return s
    # 防御性：处理 "Sinonovacula_constricta_Sinonovacula_constricta" 这种偶然重复
    for _ in range(3):
        toks = s.split('_')
        if len(toks) % 2 != 0 or not toks:
            break
        half = len(toks)//2
        left = '_'.join(toks[:half]); right = '_'.join(toks[half:])
        if left == right and left:
            s = left
        else:
            break
    if alias_map:
        s = alias_map.get(s, s)
    return s

# ---------- 常用工具 ----------
def load_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def read_results_anchor(of_root: Path) -> Path:
    """
    根据 RESULTS_DIR.txt 或唯一的 Results_* 子目录，确定 Orthofinder 结果目录。
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

def og_id_from_name(name: str) -> str:
    m = re.match(r"^(OG\d+)", name)
    if not m:
        raise ValueError(f"[ERR] 文件名不含 OG 前缀：{name}")
    return m.group(1)

def read_fasta_headers(path: Path) -> List[str]:
    heads: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(">"):
            heads.append(line[1:].strip())
    return heads



def read_fasta_as_ordered_dict(path: Path) -> "OrderedDict[str,str]":
    """读取 FASTA 为 OrderedDict[header->seq]（去空白、保留顺序、全大写）"""
    from collections import OrderedDict
    d = OrderedDict()
    name = None
    buf: List[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    d[name] = re.sub(r"\s+", "", "".join(buf)).upper()
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            d[name] = re.sub(r"\s+", "", "".join(buf)).upper()
    return d

def write_fasta_ordered_dict(path: Path, d: "OrderedDict[str,str]", linewrap: int = 60) -> None:
    """写 FASTA（按顺序）"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as w:
        for h, s in d.items():
            w.write(f">{h}\n")
            if linewrap and linewrap > 0:
                for i in range(0, len(s), linewrap):
                    w.write(s[i:i+linewrap] + "\n")
            else:
                w.write(s + "\n")

def normalize_of_msa_headers_to_species(of_msa: Path,
                                       alias_map: Dict[str,str],
                                       leaves: Set[str]) -> "OrderedDict[str,str]":
    """把 OrthoFinder 的原始 MSA（header 常为 Species|Gene）转换为：header=规范物种名，序列保持对齐。"""
    from collections import OrderedDict
    recs = read_fasta_as_ordered_dict(of_msa)
    if not recs:
        raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 为空：{of_msa}")
    out = OrderedDict()
    for h, seq in recs.items():
        tok = h.split()[0]
        sp_raw = tok.split("|", 1)[0] if "|" in tok else tok
        sp = canon_species(sp_raw, alias_map)
        if sp in out:
            raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 出现重复物种名（header 归一后冲突）：{sp} —— {of_msa}")
        out[sp] = seq
    # 等长性检查
    Ls = {len(s) for s in out.values()}
    if len(Ls) != 1:
        raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 非等长对齐：{of_msa}")
    if leaves and set(out.keys()) != set(leaves):
        miss = sorted(set(leaves) - set(out.keys()))
        extra = sorted(set(out.keys()) - set(leaves))
        raise RuntimeError(
            f"[ERR] OrthoFinder 原始 MSA 物种覆盖与物种树不一致：{of_msa}\n"
            f"      missing={miss[:5]}... ({len(miss)})\n"
            f"      extra  ={extra[:5]}... ({len(extra)})"
        )
    return out


def normalize_of_msa_headers_to_geneid(of_msa: Path,
                                      alias_map: Dict[str,str],
                                      leaves: Set[str]) -> "OrderedDict[str,str]":
    """把 OrthoFinder 的原始 MSA（header 常为 Species|Gene/ProteinID）转换为：header=Gene/ProteinID。

    设计目标：
      - 让下游 PAL2NAL 可以用“真实蛋白 ID”去匹配 CDS（而不是物种名）；
      - 同时仍然严格检验这是 strict-SCO：每个物种只能出现一次，且覆盖集合=物种树叶集合。

    返回：
      OrderedDict[gene_or_protein_id] = aligned_aa_seq
    """
    from collections import OrderedDict
    recs = read_fasta_as_ordered_dict(of_msa)
    if not recs:
        raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 为空：{of_msa}")
    out = OrderedDict()
    seen_species: Set[str] = set()
    for h, seq in recs.items():
        tok = h.split()[0]
        if "|" not in tok:
            raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA header 不含 '|'（无法解析 Species|ID）：{tok} —— {of_msa}")
        sp_raw, gid = tok.split("|", 1)
        sp = canon_species(sp_raw, alias_map)
        if sp in seen_species:
            raise RuntimeError(f"[ERR] strict-SCO 冲突：同一物种出现多条序列：{sp} —— {of_msa}")
        seen_species.add(sp)

        gid = _sanitize_id_contract(gid.strip())
        if not gid:
            raise RuntimeError(f"[ERR] Gene/ProteinID 为空：{tok} —— {of_msa}")
        if gid in out:
            raise RuntimeError(f"[ERR] Gene/ProteinID 重复：{gid} —— {of_msa}")
        out[gid] = seq

    # 等长性检查
    Ls = {len(s) for s in out.values()}
    if len(Ls) != 1:
        raise RuntimeError(f"[ERR] OrthoFinder 原始 MSA 非等长对齐：{of_msa}")
    if leaves and set(seen_species) != set(leaves):
        miss = sorted(set(leaves) - set(seen_species))
        extra = sorted(set(seen_species) - set(leaves))
        raise RuntimeError(
            f"[ERR] OrthoFinder 原始 MSA 物种覆盖与物种树不一致：{of_msa}\n"
            f"      missing={miss[:5]}... ({len(miss)})\n"
            f"      extra  ={extra[:5]}... ({len(extra)})"
        )
    return out


def find_of_raw_msa_file(msa_dir: Path, ogid: str) -> Path:
    """在 OrthoFinder 的 MultipleSequenceAlignments 目录中定位 OG 的原始 AA MSA 文件（稳健查找）。"""
    # 优先常见命名
    preferred = [f"{ogid}.fa", f"{ogid}.fasta", f"{ogid}.faa", f"{ogid}.aln.fa", f"{ogid}.aln.fasta"]
    for name in preferred:
        p = msa_dir / name
        if p.is_file():
            return p
    # 退化：匹配所有以 ogid 开头的文件（避免不同版本后缀差异）
    cands = sorted([p for p in msa_dir.glob(f"{ogid}*") if p.is_file()])
    if len(cands) == 1:
        return cands[0]
    if not cands:
        raise FileNotFoundError(f"[ERR] 未找到 OrthoFinder 原始 MSA：{msa_dir}/{ogid}*")
    # 多个候选时：尽量挑选最短文件名且后缀在常见集合里
    common_ext = {".fa", ".faa", ".fasta"}
    cands2 = [p for p in cands if p.suffix.lower() in common_ext]
    if len(cands2) == 1:
        return cands2[0]
    # 仍然不唯一则报错，避免选错
    raise RuntimeError(
        f"[ERR] OrthoFinder 原始 MSA 存在多个候选，无法唯一判定：{ogid}\n"
        f"      candidates={','.join([x.name for x in cands[:10]])}"
    )
# ---------- Newick 叶标签解析（忽略内部节点/分支长度/支持值） ----------
def parse_newick_tips(nwk_text: str) -> List[str]:
    tips: List[str] = []
    token = []
    in_quote = False
    prev_delim = None

    def commit():
        nonlocal token, prev_delim, tips
        if not token:
            return
        s = "".join(token).strip()
        token.clear()
        if not s:
            return
        # 内部节点名或支持值通常紧跟在 ')' 或 ':' 之后；叶节点则出现在 '(' 或 ',' 之后
        if prev_delim in (")", ":"):
            # 可能是内部节点名或支持值（纯数字 / 数字/数字）
            if re.fullmatch(r"[\d.]+(?:/[\d.]+)?", s):
                return
        tips.append(s)

    for ch in nwk_text:
        if in_quote:
            if ch == "'":
                in_quote = False
                commit()
            else:
                token.append(ch)
        else:
            if ch == "'":
                in_quote = True
            elif ch in "(),:;":
                commit()
                prev_delim = ch
            elif ch.isspace():
                commit()
            else:
                token.append(ch)
    commit()
    return sorted(set(tips))

# ---------- SeqID 归一化 & 数据库前缀剥离 ----------
def _strip_dbprefix(s: str) -> str:
    return re.sub(rf"^{_DB_PREFIXES_REGEX}\|", "", s)

def _norm_first_token(s: str) -> str:
    head, *rest = s.split("|", 1)
    head = re.sub(r"([A-Za-z]{2,}[A-Za-z0-9]*)[:_]+([A-Za-z0-9]{2,})", r"\1_\2", head)
    return head if not rest else head + "|" + rest[0]

def _norm_for_match(h: str) -> str:
    return _norm_first_token(_strip_dbprefix(h.strip()))

# ---------- 前缀白名单检测 ----------
def _detect_unknown_db_prefixes(peps_tsv: Path, cds_dir: Path, cds_suffix: str, extra: Set[str]) -> None:
    allow = set(ALLOWED_DB_PREFIXES) | set(extra or set())
    seen: Dict[str, List[str]] = {}

    def collect(h: str, where: str):
        m = re.match(r"^([A-Za-z]{2,5})\|", h.strip())
        if m and m.group(1) not in allow:
            seen.setdefault(m.group(1), [])
            if len(seen[m.group(1)]) < 3:
                seen[m.group(1)].append(f"{where}: {h[:120]}")

    # pep2cds.tsv
    if peps_tsv.is_file():
        with open(peps_tsv, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, [])
            try:
                i_pid = [x.lower() for x in header].index("protein_id")
            except ValueError:
                i_pid = 1 if len(header) > 1 else -1
            if i_pid >= 0:
                for row in reader:
                    if not row:
                        continue
                    collect(row[i_pid], "pep2cds.tsv:protein_id")

    # CDS FASTA 头（每个物种抽一个）
    for fa in sorted(Path(cds_dir).glob(f"*{cds_suffix}")):
        try:
            with open(fa, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith(">"):
                        collect(line[1:].strip(), f"CDS:{fa.name}")
                        break
        except Exception:
            continue

    if seen:
        lines = [
            "[ERR] 发现未知的数据库前缀（形如 'xxx|' 开头），为避免静默错配，发布终止。",
            f"[INFO] 允许前缀白名单 = {sorted(allow)}；可在 config.publish.db_prefix_allowlist 扩展。"
        ]
        for pref, samples in seen.items():
            lines.append(f"[HINT] 前缀 '{pref}|' 的示例：")
            for s in samples:
                lines.append(f"       {s}")
        raise RuntimeError("\n".join(lines))

# ---------- pep2cds 全局映射 ----------
def read_pep2cds_global(tsv: Path, alias_map: Dict[str, str]) -> Dict[Tuple[str,str], str]:
    mp: Dict[Tuple[str,str], str] = {}

    def add_keys(sp: str, pid: str, cds: str):
        # 原样
        mp.setdefault((sp, pid), cds)
        # 规范化主形态
        pid_base = _norm_for_match(pid)
        mp.setdefault((sp, pid_base), cds)
        # 管道 -> 下划线（及压缩）
        pid_u = pid.replace("|", "_")
        mp.setdefault((sp, pid_u), cds)
        mp.setdefault((sp, re.sub(r"_+", "_", pid_u)), cds)
        pid_base_u = pid_base.replace("|", "_")
        mp.setdefault((sp, pid_base_u), cds)
        mp.setdefault((sp, re.sub(r"_+", "_", pid_base_u)), cds)
        # 尾段 tail
        if "|" in pid:
            tail = pid.rsplit("|", 1)[-1]
            mp.setdefault((sp, tail), cds)
            tail_u = tail.replace("|", "_")
            mp.setdefault((sp, tail_u), cds)
            mp.setdefault((sp, re.sub(r"_+", "_", tail_u)), cds)

    need_file(tsv, "pep2cds 汇总表缺失")
    with open(tsv, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, [])
        hl = [x.lower() for x in header]
        try:
            i_sp  = hl.index("species"); i_pid = hl.index("protein_id"); i_cds = hl.index("cds_id")
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

# ---------- Orthogroups 成员解析（用于 CAFE 全量家族映射） ----------
def _sanitize_id_contract(x: str) -> str:
    """
    与 extract_longest_from_gff.py 中保持一致的 ID “合同化清洗”规则：
    - 保留竖线 '|'（例如 gnl|WGS_AMQN|CAPTEDRAFT... 依然合法）
    - 仅允许 [A-Za-z0-9._:-|]，其它字符替换为 '_'
    - 去首尾空白
    """
    return re.sub(r"[^A-Za-z0-9._:\-\|]", "_", x.strip())


def parse_orthogroups_members(tsv: Path,
                              alias_map: Dict[str, str]) -> List[Tuple[str, str, str]]:
    """
    解析 OrthoFinder 的 Orthogroups.tsv：

    输入表典型结构：
        Orthogroup  SpeciesA  SpeciesB  ...
        OG0000001   id1,id2   id3      ...

    返回：
        List[(OG, Species, protein_id)]

    约束与约定：
      - Species 列名按 alias_map 规范化；
      - 每个单元格内允许以逗号/空格分隔多个蛋白；
      - protein_id 取去除“Species|”前缀后的部分，
        再应用 _sanitize_id_contract，与 pep2cds.tsv 的 protein_id 保持同一规范。
    """
    need_file(tsv, "Orthofinder Orthogroups.tsv 缺失")
    triplets: List[Tuple[str, str, str]] = []

    with open(tsv, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, [])
        if not header or header[0] != "Orthogroup":
            raise RuntimeError(f"[ERR] Orthogroups.tsv 首列应为 Orthogroup：{tsv}")

        species_cols = header[1:]
        if not species_cols:
            raise RuntimeError(f"[ERR] Orthogroups.tsv 中缺少物种列：{tsv}")

        for row in reader:
            if not row or not row[0].strip():
                continue
            og = row[0].strip()
            # 遍历每个物种列
            for idx_col, sp_raw in enumerate(species_cols, start=1):
                if idx_col >= len(row):
                    continue
                cell = row[idx_col].strip()
                if not cell:
                    continue
                sp = canon_species(sp_raw, alias_map)
                # 拆分多个蛋白 ID（以逗号和空白符为分隔）
                for token in re.split(r"[;,\s]+", cell):
                    token = token.strip()
                    if not token:
                        continue
                    # 去掉可能存在的 Species| 前缀，仅保留 SeqID 部分
                    if "|" in token:
                        _, seqid = token.split("|", 1)
                    else:
                        seqid = token
                    seqid = _sanitize_id_contract(seqid)
                    triplets.append((og, sp, seqid))

    if not triplets:
        raise RuntimeError(f"[ERR] 从 Orthogroups.tsv 未解析到任何 OG×物种×蛋白条目：{tsv}")
    return triplets

# ---------- FASTA 索引：同时按“整行 header”与“首 token”建键，并生成规范化别名 ----------
class FastaIndex:
    """
    根据 CDS FASTA 建立两类键 → 真实整行 header：
      - full header：'>' 后整行（去掉换行）
      - first token：整行头的首个非空白 token（遇到描述时只取 token）
    同时为 token 生成别名：去数据库前缀、统一第1 token 的 ':' 与 '_'、'|'→'_'、尾段 tail 等。
    """
    def __init__(self, path: Path):
        self.full_to_full: Dict[str, str] = {}
        self.key_to_full: Dict[str, str] = {}
        self._build(path)

    @staticmethod
    def _first_token(h: str) -> str:
        return h.split()[0] if h else h

    def _put_key(self, k: str, full: str):
        if not k: return
        self.key_to_full.setdefault(k, full)

    def _build(self, path: Path):
        name: Optional[str] = None
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    name = line[1:].rstrip("\n")
                    self.full_to_full[name] = name
                    tok = self._first_token(name)
                    if tok:
                        base = _norm_for_match(tok)
                        self._put_key(tok, name)
                        self._put_key(base, name)
                        # '|'→'_' 派生键
                        for k in (tok, base):
                            if not k:
                                continue
                            u = k.replace("|", "_")
                            self._put_key(u, name)
                            self._put_key(re.sub(r"_+", "_", u), name)
                        # 尾段 tail
                        if "|" in tok:
                            tail = tok.rsplit("|", 1)[-1]
                            self._put_key(tail, name)
                            tu = tail.replace("|", "_")
                            self._put_key(tu, name)
                            self._put_key(re.sub(r"_+", "_", tu), name)
                else:
                    continue

    def resolve(self, cds_id: str) -> Optional[str]:
        s = cds_id.strip()
        if not s:
            return None

        # 1) 完全等于某条整行 header
        if s in self.full_to_full:
            return s

        # 2) 尝试多种 token 别名
        cand: List[str] = []

        def add(x: str):
            if x and x not in cand:
                cand.append(x)

        add(self._first_token(s))
        base = _norm_for_match(s)
        add(self._first_token(base))
        for k in (s, base):
            if not k:
                continue
            u = k.replace("|", "_")
            add(u)
            add(re.sub(r"_+", "_", u))
        if "|" in s:
            tail = s.rsplit("|", 1)[-1]
            add(tail)
            tu = tail.replace("|", "_")
            add(tu)
            add(re.sub(r"_+", "_", tu))

        for k in cand:
            full = self.key_to_full.get(k)
            if full:
                return full
        return None

def main():
    cfg = load_cfg()
    P = cfg.get("paths", {}) or {}
    alias_map = (cfg.get("species", {}) or {}).get("alias_map", {}) or {}
    allow_extra = set((cfg.get("publish", {}) or {}).get("db_prefix_allowlist", []) or [])

    # 路径锚点
    pub = Path(P["publish_dir"]).resolve()
    reports_dir = Path(P["reports_dir"])
    species_only_dir = Path(P["species_collapse_dir"])
    trees_dir = Path(P.get("trees_dir", "results/trees"))
    of_root = Path(P["orthofinder_results_dir"])

    # 关键文件
    ogs_selected = reports_dir / "ogs_selected.list"
    sco_filelist = reports_dir / "sco_filelist.txt"
    og_species_protein = reports_dir / "og_species_protein.tsv"
    matrix_tsv = reports_dir / "matrix.tsv"
    species_tree = trees_dir / "species_tree.nwk"

    # CDS 仓库
    cds_dir = Path(cfg.get("input", {}).get("cds_dir", DEFAULT_CDS_DIR))
    cds_suffix = cfg.get("input", {}).get("cds_suffix", DEFAULT_CDS_SUFFIX)

    # OrthoFinder 与 family.tsv
    of_results_dir = read_results_anchor(of_root)
    gene_count_src = of_results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    need_file(gene_count_src, "OrthoFinder GeneCount 缺失")
    orthogroups_src = of_results_dir / "Orthogroups" / "Orthogroups.tsv"
    need_file(orthogroups_src, "OrthoFinder Orthogroups.tsv 缺失")

    # 物种叶集合
    leaves = set(parse_newick_tips(Path(species_tree).read_text(encoding="utf-8")))
    if not leaves:
        raise RuntimeError("[ERR] 无法从物种树解析叶节点标签（tip labels）")

    # 前缀白名单检测
    pep2cds_src = Path(P.get("pep2cds_tsv", DEFAULT_PEP2CDS_TSV))
    need_file(pep2cds_src, "缺少蛋白→CDS 全局映射（pep2cds.tsv）")
    _detect_unknown_db_prefixes(pep2cds_src, cds_dir, cds_suffix, allow_extra)

    # 重建发布目录
    if OVERWRITE_PUBLISH_DIR and pub.exists():
        shutil.rmtree(pub)
    pub.mkdir(parents=True, exist_ok=True)

    # ---------- A) 严格 SCO MSA 发布 ----------
    def load_kept_ogs_from_filelist(filelist: Path) -> List[Path]:
        files: List[Path] = []
        with open(filelist, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s:
                    files.append(Path(s))
        if not files:
            raise RuntimeError(f"[ERR] sco_filelist.txt 为空：{filelist}")
        return files

    cand_files = load_kept_ogs_from_filelist(sco_filelist)

    kept_files: List[Path] = []
    kept_ogs: List[str] = []
    excluded: List[Tuple[str,str]] = []

    for raw_path in cand_files:
        p = raw_path if raw_path.is_file() else (species_only_dir / raw_path.name)
        if not p.is_file():
            excluded.append((raw_path.name, "file_missing")); continue
        ogid = og_id_from_name(p.name)
        heads = [canon_species(h, alias_map) for h in read_fasta_headers(p)]
        if not heads:
            excluded.append((ogid, "empty_alignment")); continue
        if any('|' in h for h in heads):
            excluded.append((ogid, "non_species_headers")); continue
        if set(heads) != leaves:
            excluded.append((ogid, "coverage_mismatch")); continue
        kept_files.append(p); kept_ogs.append(ogid)

    # 发布 strict_sco_msa/
    sco_dir = pub / "strict_sco_msa"; sco_dir.mkdir(parents=True, exist_ok=True)
    for src in kept_files:
        shutil.copy2(src, sco_dir / src.name)

    # ---------- A2) 发布 OrthoFinder 原始 AA MSA（同目录；header 归一为 Gene/ProteinID） ----------
    # 说明：
    #   - downstream 的 aphylo 03b 会使用 “raw AA（Gene/ProteinID header）+ CDS” 来回译 codon；
    #   - 为了让 raw AA 的 header 与 CDS 的 header 能 1:1 对应，这里将 OrthoFinder 原始 header（常为 Species|Gene）
    #     在发布时统一转换为 Gene/ProteinID（即 '|' 右侧部分），并写入 strict_sco_msa/。
    #   - 产物命名固定为：OGxxxx{DEFAULT_RAW_MSA_SUFFIX}（默认 .raw.fa）
    of_msa_dir = of_results_dir / "MultipleSequenceAlignments"
    need_dir(of_msa_dir, "OrthoFinder MultipleSequenceAlignments 缺失")
    raw_suffix = str((cfg.get("publish", {}) or {}).get("raw_msa_suffix", DEFAULT_RAW_MSA_SUFFIX)) or DEFAULT_RAW_MSA_SUFFIX

    raw_published = 0
    for ogid in kept_ogs:
        of_raw = find_of_raw_msa_file(of_msa_dir, ogid)
        raw_norm = normalize_of_msa_headers_to_geneid(of_raw, alias_map, leaves)
        out_raw = sco_dir / f"{ogid}{raw_suffix}"
        write_fasta_ordered_dict(out_raw, raw_norm, linewrap=60)
        raw_published += 1

    # ---------- C) 复制其它固定件 ----------
    shutil.copy2(ogs_selected,  pub / "ogs_selected.list")
    shutil.copy2(sco_filelist,  pub / "sco_filelist.txt")
    copy_if_exists(og_species_protein, pub / "og_species_protein.tsv")
    copy_if_exists(matrix_tsv,        pub / "matrix.tsv")
    shutil.copy2(species_tree,        pub / "species_tree.nwk")
    shutil.copy2(gene_count_src,      pub / "family.tsv")

    include_ultra = bool(cfg.get("publish", {}).get("include_ultrametric_tree", False))
    if include_ultra:
        ultra_path_cfg = P.get("ultrametric_tree", "")
        ultra_path = Path(ultra_path_cfg) if isinstance(ultra_path_cfg, str) and ultra_path_cfg else None
        if ultra_path and ultra_path.is_file():
            shutil.copy2(ultra_path, pub / "species_tree_ultrametric.nwk")

    # ---------- D) 发布集合清单 ----------
    (pub / "ogs_selected.published.list").write_text("\n".join(kept_ogs) + "\n", encoding="utf-8")
    with open(pub / "sco_filelist.published.txt", "w", encoding="utf-8") as wf:
        for src in kept_files:
            wf.write(src.name + "\n")

    # ---------- E) pep2cds_resolved.tsv ----------
    def parse_og_species_protein(tsv: Path, only_ogs: Set[str], alias_map: Dict[str,str]) -> List[Tuple[str,str,str]]:
        triplets: List[Tuple[str,str,str]] = []
        with open(tsv, "r", encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            def idx(col): return header.index(col) if col in header else -1
            i_og = idx("OG"); i_sp = idx("Species"); i_id = idx("SeqID")
            if min(i_og, i_sp, i_id) < 0:
                raise RuntimeError(f"[ERR] {tsv} 缺少必要列（OG/Species/SeqID）")
            for line in f:
                if not line.strip(): continue
                parts = line.rstrip("\n").split("\t")
                og, sp_raw, pid = parts[i_og], parts[i_sp], parts[i_id]
                if og in only_ogs:
                    triplets.append((og, canon_species(sp_raw, alias_map), pid))
        if not triplets:
            raise RuntimeError("[ERR] og_species_protein.tsv 过滤后为空；请检查与发布集合是否匹配")
        return triplets

    kept_from_filelist = {og_id_from_name(p.name) for p in cand_files}
    og_sp_pid_all = parse_og_species_protein(og_species_protein, kept_from_filelist, alias_map)

    sp_pid_to_cds = read_pep2cds_global(pep2cds_src, alias_map)

    # 逐物种构建 FASTA 索引（含 token 别名）
    sp_to_index: Dict[str, FastaIndex] = {}
    def cds_index_for(sp: str) -> FastaIndex:
        sp_std = canon_species(sp, alias_map)
        if sp_std in sp_to_index:
            return sp_to_index[sp_std]
        fa = cds_dir / f"{sp_std}{DEFAULT_CDS_SUFFIX if not cds_suffix else cds_suffix}"
        need_file(fa, f"CDS FASTA 缺失（期望 {sp_std}{cds_suffix}；根目录 {cds_dir}）")
        idx = FastaIndex(fa)
        sp_to_index[sp_std] = idx
        return idx

    out_resolved = pub / "pep2cds_resolved.tsv"
    unresolved: List[Tuple[str,str,str,str]] = []
    n_out = 0
    keep_set = set(kept_ogs)

    with open(out_resolved, "w", encoding="utf-8", newline="") as w:
        writer = csv.writer(w, delimiter="\t")
        writer.writerow(["OG", "Species", "protein_id", "cds_id"])
        for og, sp, pid in og_sp_pid_all:
            if og not in keep_set:
                continue
            # protein → cds 原始映射（多候选）
            cds = None
            for x in (pid, _norm_for_match(pid), pid.rsplit("|",1)[-1] if "|" in pid else None):
                if not x: continue
                cds = (sp_pid_to_cds.get((sp, x)) or
                       sp_pid_to_cds.get((sp, x.replace("|","_"))) or
                       sp_pid_to_cds.get((sp, re.sub(r"_+", "_", x.replace("|","_")))))
                if cds: break
            if cds is None:
                unresolved.append((og, sp, pid, "NA(protein→cds not found)"))
                continue

            # 规范到真实 FASTA header（利用 token 索引）
            idx = cds_index_for(sp)
            real_hdr = idx.resolve(cds)
            if real_hdr is None:
                unresolved.append((og, sp, pid, cds))
                continue

            writer.writerow([og, sp, pid, real_hdr])
            n_out += 1

    n_unres = len(unresolved)
    if n_unres > 0:
        with open(pub / "pep2cds_unresolved.tsv", "w", encoding="utf-8") as wf:
            wf.write("OG\tSpecies\tprotein_id\tcds_id_raw\n")
            for og, sp, pid, cds in unresolved:
                wf.write(f"{og}\t{sp}\t{pid}\t{cds}\n")
        if STRICT_NORMALIZE_CDS:
            raise RuntimeError(
                f"[ERR] pep2cds 规范化失败：共有 {n_unres} 条无法在 CDS FASTA 中找到真实表头；"
                f"详见 {pub/'pep2cds_unresolved.tsv'}。\n"
                f"建议：检查该物种的 CDS 头写法与 pep2cds.tsv 的一致性，或在 config.publish.db_prefix_allowlist 扩展允许的数据库前缀。"
            )

    # ---------- E2) all_pep2cds_resolved.tsv（全 OG，用于 CAFE 扩张/收缩分析） ----------
    og_sp_pid_all_members = parse_orthogroups_members(orthogroups_src, alias_map)
    all_out = pub / "all_pep2cds_resolved.tsv"
    n_all_out = 0
    with open(all_out, "w", encoding="utf-8", newline="") as w_all:
        writer_all = csv.writer(w_all, delimiter="\t")
        writer_all.writerow(["OG", "Species", "protein_id", "cds_id"])
        for og, sp, pid in og_sp_pid_all_members:
            # protein → cds 原始映射（多候选；不对缺失条目强制报错，仅跳过）
            cds = None
            for x in (pid, _norm_for_match(pid), pid.rsplit("|",1)[-1] if "|" in pid else None):
                if not x:
                    continue
                cds = (sp_pid_to_cds.get((sp, x)) or
                       sp_pid_to_cds.get((sp, x.replace("|","_"))) or
                       sp_pid_to_cds.get((sp, re.sub(r"_+", "_", x.replace("|","_")))))
                if cds:
                    break
            if cds is None:
                continue

            idx = cds_index_for(sp)
            real_hdr = idx.resolve(cds)
            if real_hdr is None:
                continue

            writer_all.writerow([og, sp, pid, real_hdr])
            n_all_out += 1

    # ---------- F) 报告 ----------
    manifest = {
        "publish_dir": str(pub),
        "kept_msa_count": len(kept_files),
        "species_tree": str(pub / "species_tree.nwk"),
        "family_tsv": str(pub / "family.tsv"),
        "ogs_selected": str(pub / "ogs_selected.list"),
        "filelist": str(pub / "sco_filelist.txt"),
        "ogs_selected_published": str(pub / "ogs_selected.published.list"),
        "sco_filelist_published": str(pub / "sco_filelist.published.txt"),
        "has_matrix": (pub / "matrix.tsv").is_file(),
        "has_og_species_protein": (pub / "og_species_protein.tsv").is_file(),
        "pep2cds_resolved": str(out_resolved),
        "pep2cds_resolved_lines": n_out,
        "pep2cds_unresolved_count": n_unres,
        "all_pep2cds_resolved": str(all_out),
        "all_pep2cds_resolved_lines": n_all_out,
        "excluded_count": len(excluded),
        "cds_dir": str(cds_dir),
        "cds_suffix": cds_suffix,
        "strict_normalize_cds": STRICT_NORMALIZE_CDS,
        "db_prefix_allowlist": sorted(list(ALLOWED_DB_PREFIXES | allow_extra)),
    }
    with open(pub / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    with open(pub / "QC_report.txt", "w", encoding="utf-8") as f:
        f.write(f"[INFO] strict_sco_msa（发布）数量：{len(kept_files)}\n")
        f.write(f"[INFO] species_tree.nwk 是否存在：{(pub / 'species_tree.nwk').is_file()}\n")
        f.write(f"[INFO] family.tsv 是否存在：{(pub / 'family.tsv').is_file()}\n")
        f.write(f"[INFO] ogs_selected.list 是否存在：{(pub / 'ogs_selected.list').is_file()}\n")
        f.write(f"[INFO] sco_filelist.txt 是否存在：{(pub / 'sco_filelist.txt').is_file()}\n")
        f.write(f"[INFO] pep2cds_resolved.tsv 行数：{n_out}\n")
        f.write(f"[INFO] all_pep2cds_resolved.tsv 行数：{n_all_out}\n")
        if excluded:
            f.write(f"[WARN] 剔除 OG 数量：{len(excluded)}（详见 excluded_reason.tsv）\n")
        if n_unres > 0:
            f.write(f"[WARN] pep2cds_unresolved.tsv 条目：{n_unres} —— {'终止发布' if STRICT_NORMALIZE_CDS else '继续发布'}\n")

    if excluded:
        with open(pub / "excluded_reason.tsv", "w", encoding="utf-8") as w:
            w.write("OG\treason\n")
            for og, rs in excluded:
                w.write(f"{og}\t{rs}\n")

    # 哨兵
    (pub / '.done').touch()
    (pub / 'strict_sco_msa' / '.done').touch()

    print(f"[DONE] 发布完成：{pub}")
    print(f"[STAT] strict_sco_msa: {len(kept_files)}；pep2cds_resolved: {n_out}；unresolved: {n_unres}；all_pep2cds_resolved: {n_all_out}")
    if excluded:
        print(f"[STAT] 剔除 OG: {len(excluded)} —— {pub/'excluded_reason.tsv'}")

if __name__ == "__main__":
    main()

