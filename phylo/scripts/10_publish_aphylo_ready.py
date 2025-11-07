#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_publish_aphylo_ready.py —— 发布到 aphylo 可直接消费的“五件套”，稳健严格版
要点：
  * 只接受 AA 位串掩码（长度=AA 列数），若存在但不合规则报错；不做任何“猜测或兼容”。
  * Newick 叶名稳健解析（忽略内部节点/支持值/分支长度，支持引号）。
  * SeqID 归一化：统一剥离 NCBI 数据库前缀（gb|/ref|/sp|/tr|/gnl|...），仅规范第一个 '|' 前 token 的 ':' 与 '_'。
  * 前缀白名单检测：pep2cds 与 CDS FASTA 出现未知数据库前缀时立即报错并给出示例。
  * pep2cds → 真实 FASTA 头：FASTA 索引同时使用“整行 header”与“首 token”两类键，并生成规范化/尾段候选，避免 410 类 miss。
产物契约保持不变：strict_sco_msa/、colmask/、pep2cds_resolved.tsv、family.tsv、species_tree.nwk、清单与报告等。
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

# 数据库前缀白名单（基于 NCBI SeqID 规范，可在 config.publish.db_prefix_allowlist 扩展）
ALLOWED_DB_PREFIXES: Set[str] = {
    "lcl","gb","emb","dbj","ref","sp","tr","pir","prf","pdb","pat","pgp",
    "bbs","bbm","gim","gi","tpg","tpe","tpd","tpa","gnl"
}
_DB_PREFIXES_REGEX = r"(?:lcl|gb|emb|dbj|ref|sp|tr|pir|prf|pdb|pat|pgp|bbs|bbm|gim|gi|tpg|tpe|tpd|tpa|gnl)"
# =====================================================================

# ---------- 物种名规范化（与上游保持一致） ----------
def canon_species(sp_raw: str, alias_map: Optional[Dict[str, str]] = None) -> str:
    s = sp_raw.strip()
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

def read_results_anchor(of_root: Path) -> Path:
    anchor = of_root / "RESULTS_DIR.txt"
    if anchor.is_file():
        txt = anchor.read_text(encoding="utf-8").strip()
        p = Path(txt)
        if p.is_dir():
            return p
    cand = sorted(of_root.glob("Results_*"))
    if len(cand) == 1 and cand[0].is_dir():
        return cand[0]
    raise FileNotFoundError("[ERR] 无法确定唯一的 Orthofinder Results_* 目录；请检查 RESULTS_DIR.txt 或目录结构")

def load_kept_ogs_from_filelist(sco_filelist: Path) -> List[Path]:
    files: List[Path] = []
    with open(sco_filelist, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                files.append(Path(s))
    if not files:
        raise RuntimeError("[ERR] sco_filelist.txt 中没有任何条目")
    return files

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

# ---------- Newick 叶标签解析（忽略内部节点/分支长度/支持值） ----------
def parse_newick_tips(nwk_text: str) -> List[str]:
    tips: List[str] = []
    token = []
    in_quote = False
    prev_delim = None

    def commit():
        nonlocal token, prev_delim, tips
        if not token: return
        s = "".join(token).strip()
        token.clear()
        if not s: return
        if re.fullmatch(r"[\d.]+(?:/[\d.]+)?", s): return   # 纯数字 / 数字/数字
        if prev_delim in (")", ":"): return                # 内部节点/分支长度
        tips.append(s)

    for ch in nwk_text:
        if in_quote:
            if ch == "'":
                in_quote = False; commit()
            else:
                token.append(ch)
        else:
            if ch == "'":
                in_quote = True
            elif ch in ("(", ")", ",", ":", ";"):
                commit(); prev_delim = ch
            elif ch.isspace():
                commit()
            else:
                token.append(ch)
    commit()
    return sorted(set(tips))

# ---------- SeqID 归一化与前缀检测 ----------
def _strip_dbprefix(s: str) -> str:
    return re.sub(rf"^{_DB_PREFIXES_REGEX}\|", "", s)

def _norm_first_token(s: str) -> str:
    head, *rest = s.split("|", 1)
    head = re.sub(r"([A-Za-z]{2,}[A-Za-z0-9]*)[:_]+([A-Za-z0-9]{2,})", r"\1_\2", head)
    return head if not rest else head + "|" + rest[0]

def _norm_for_match(h: str) -> str:
    return _norm_first_token(_strip_dbprefix(h.strip()))

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
            lower = [x.lower() for x in header]
            i_pid = lower.index("protein_id") if "protein_id" in lower else (1 if len(header) > 1 else -1)
            i_cds = lower.index("cds_id")     if "cds_id"     in lower else (2 if len(header) > 2 else -1)
            for row in reader:
                if not row: continue
                if 0 <= i_pid < len(row): collect(row[i_pid], "pep2cds.protein_id")
                if 0 <= i_cds < len(row): collect(row[i_cds], "pep2cds.cds_id")

    # CDS FASTA（取每个物种一个示例）
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
        msg = ["[ERR] 发现未知的数据库前缀（形如 'xxx|' 开头）。为避免静默错配，发布终止。"]
        msg.append(f"[INFO] 白名单={sorted(allow)}；可在 config.publish.db_prefix_allowlist 扩展。")
        for k, ex in seen.items():
            msg.append(f"  - 前缀 '{k}|' 示例：")
            for s in ex: msg.append(f"      {s}")
        raise RuntimeError("\n".join(msg))

# ---------- pep2cds 读取：为 protein_id 生成稳健候选键 ----------
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
        # 最后一段 tail
        if "|" in pid:
            tail = pid.rsplit("|", 1)[-1]
            mp.setdefault((sp, tail), cds)
            tail_u = tail.replace("|", "_")
            mp.setdefault((sp, tail_u), cds)
            mp.setdefault((sp, re.sub(r"_+", "_", tail_u)), cds)

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
                if not line: continue
                if line.startswith(">"):
                    full = line[1:].strip()
                    self.full_to_full[full] = full

                    tok = self._first_token(full)
                    # 主键
                    self._put_key(tok, full)
                    # 规范化与派生键
                    base = _norm_for_match(tok)
                    self._put_key(base, full)
                    self._put_key(tok.replace("|","_"), full)
                    self._put_key(re.sub(r"_+", "_", tok.replace("|","_")), full)
                    self._put_key(base.replace("|","_"), full)
                    self._put_key(re.sub(r"_+", "_", base.replace("|","_")), full)
                    if "|" in tok:
                        tail = tok.rsplit("|", 1)[-1]
                        self._put_key(tail, full)
                        self._put_key(re.sub(r"_+", "_", tail.replace("|","_")), full)

    def resolve(self, cds_id: str) -> Optional[str]:
        # 1) 精确等于某条“整行 header”
        if cds_id in self.full_to_full:
            return cds_id
        # 2) 用候选键在 key_to_full 命中
        s = cds_id.strip()
        cand: List[str] = []
        def add(x: str):
            if x and x not in cand: cand.append(x)
        add(self._first_token(s))
        base = _norm_for_match(s)
        add(self._first_token(base))
        add(s.replace("|","_")); add(re.sub(r"_+", "_", s.replace("|","_")))
        add(base.replace("|","_")); add(re.sub(r"_+", "_", base.replace("|","_")))
        if "|" in s:
            tail = s.rsplit("|",1)[-1]
            add(tail); add(re.sub(r"_+", "_", tail.replace("|","_")))
        for k in cand:
            full = self.key_to_full.get(k)
            if full: return full
        return None

# ---------- 掩码：只接受 AA 位串 ----------
def _aa_length_from_faa(path: Path) -> int:
    L, began = 0, False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if began: break
                began = True; continue
            if began:
                L += len(line.strip().replace(" ", "").replace("\t", ""))
    return L

def _read_mask_bits_strict(path: Path) -> str:
    s = Path(path).read_text(encoding="utf-8", errors="ignore").strip()
    bits = "".join(ch for ch in s if ch in "01")
    if not bits:
        raise RuntimeError(f"[ERR] 掩码文件不包含 0/1 位串：{path}")
    if re.search(r"[^01\s]", s):
        # 出现除 0/1 以外的可见字符，视为非位串
        raise RuntimeError(f"[ERR] 掩码文件含非 0/1 字符：{path}")
    return bits

# ------------------ 主流程 ------------------
def main():
    cfg = load_cfg()
    P = cfg.get("paths", {}) or {}
    alias_map = (cfg.get("species", {}) or {}).get("alias_map", {}) or {}
    allow_extra = set((cfg.get("publish", {}) or {}).get("db_prefix_allowlist", []) or [])

    # 路径锚点
    pub = Path(P["publish_dir"]).resolve()
    reports_dir = Path(P["reports_dir"])
    species_only_dir = Path(P["species_collapse_dir"])
    colmask_dir = Path(P.get("colmask_dir", "results/colmask"))  # 可无；存在则严格校验
    trees_dir = Path(P.get("trees_dir", "results/trees"))
    of_root = Path(P["orthofinder_results_dir"])

    # 关键文件
    ogs_selected = reports_dir / "ogs_selected.list"
    sco_filelist = reports_dir / "sco_filelist.txt"
    og_species_protein = reports_dir / "og_species_protein.tsv"
    matrix_tsv = reports_dir / "matrix.tsv"
    species_tree = trees_dir / "species_tree.nwk"

    # CDS 仓库
    cds_dir = (
        Path(cfg.get("inputs", {}).get("cds_dir", "")) if "inputs" in cfg and "cds_dir" in cfg["inputs"]
        else Path(P.get("cds_dir", "")) if "cds_dir" in P
        else Path(DEFAULT_CDS_DIR)
    )
    cds_suffix = cfg.get("inputs", {}).get("cds_suffix", DEFAULT_CDS_SUFFIX)

    # 存在性检查
    need_file(ogs_selected, "严格 SCO 列表缺失")
    need_file(sco_filelist, "SCO 文件清单缺失")
    need_file(og_species_protein, "OG×Species×SeqID 对照缺失")
    need_file(species_tree, "物种树缺失")
    need_dir(species_only_dir, "species-only 目录缺失")
    need_dir(cds_dir, "CDS 仓库缺失")

    # OrthoFinder 与 family.tsv
    of_results_dir = read_results_anchor(of_root)
    gene_count_src = of_results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    need_file(gene_count_src, "OrthoFinder GeneCount 缺失")

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

    # ---------- A) 严格筛选 MSA ----------
    cand_files: List[Path] = load_kept_ogs_from_filelist(sco_filelist)
    kept_files: List[Path] = []
    kept_ogs:   List[str]  = []
    excluded:   List[Tuple[str,str]] = []

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

    # ---------- B) colmask（可无；若存在必须是 AA 位串且长度=AA） ----------
    pub_colmask = pub / "colmask"
    pub_colmask.mkdir(parents=True, exist_ok=True)
    mask_checked = 0
    if colmask_dir.is_dir():
        for ogid in kept_ogs:
            aa_path = sco_dir / f"{ogid}.trim.faa"
            mask_src = colmask_dir / f"{ogid}.colmask"
            if not mask_src.is_file():
                # 允许缺失（不强制）
                continue
            L = _aa_length_from_faa(aa_path)
            bits = _read_mask_bits_strict(mask_src)
            if len(bits) != L:
                raise RuntimeError(f"[ERR] 掩码长度 != AA 列数：{mask_src}（mask={len(bits)} vs AA={L}）")
            # 合规即复制
            shutil.copy2(mask_src, pub_colmask / mask_src.name)
            mask_checked += 1

    # ---------- C) 复制其它固定件 ----------
    shutil.copy2(ogs_selected,  pub / "ogs_selected.list")
    shutil.copy2(sco_filelist,  pub / "sco_filelist.txt")
    copy_if_exists(og_species_protein, pub / "og_species_protein.tsv")
    copy_if_exists(matrix_tsv,        pub / "matrix.tsv")
    shutil.copy2(species_tree,        pub / "species_tree.nwk")
    shutil.copy2(gene_count_src,      pub / "family.tsv")

    include_ultra = bool(cfg.get("publish", {}).get("include_ultrametric_tree", False))
    if include_ultra:
        ultra_path = Path(P.get("ultrametric_tree", "")) if isinstance(P.get("ultrametric_tree", ""), str) else None
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
    unresolved = []
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
        "has_colmask": (pub / "colmask").is_dir(),
        "has_matrix": (pub / "matrix.tsv").is_file(),
        "has_og_species_protein": (pub / "og_species_protein.tsv").is_file(),
        "pep2cds_resolved": str(out_resolved),
        "pep2cds_resolved_lines": n_out,
        "pep2cds_unresolved_count": n_unres,
        "excluded_count": len(excluded),
        "cds_dir": str(cds_dir),
        "cds_suffix": cds_suffix,
        "strict_normalize_cds": STRICT_NORMALIZE_CDS,
        "db_prefix_allowlist": sorted(list(ALLOWED_DB_PREFIXES | allow_extra)),
        "colmask_checked": mask_checked
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
        f.write(f"[INFO] colmask（AA 位串）严格校验通过：{mask_checked} 个\n")
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
    (pub / 'colmask' / '.done').touch()

    print(f"[DONE] 发布完成：{pub}")
    print(f"[STAT] strict_sco_msa: {len(kept_files)}；pep2cds_resolved: {n_out}；unresolved: {n_unres}")
    if excluded:
        print(f"[STAT] 剔除 OG: {len(excluded)} —— {pub/'excluded_reason.tsv'}")
    print(f"[STAT] colmask（AA 位串）严格校验通过：{mask_checked}")

if __name__ == "__main__":
    main()

