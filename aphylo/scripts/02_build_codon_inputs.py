#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_build_codon_inputs.py —— 构建回译（codon）所需的输入文件（最终版）

核心职责（不做对齐，仅做“准备材料”）：
  1) 读取发布包中的严格 SCO 蛋白质比对（AA-MSA），按“表头物种顺序”生成该 OG 的
     CDS 输入（ordered_cds.fna），用于后续 03 步 pal2nal-like 回译；
  2) 使用发布包中的 pep→cds 解析表（pep2cds_resolved.tsv），按 (OG, Species) 唯一映射到 cds_id；
  3) 从固定的 CDS 仓库（inputs.cds_dir/*.fna）提取序列，严格一对一写入；
  4) 生成对照表（order.tsv），便于后续审计与排错；
  5) 屏幕与日志双写，可复查每一步；写入 .inputs.done 哨兵。

关键修复与加固：
  - 【必修复】OG 解析：从文件名“^OG\\d+”提取（而非 Path.stem），彻底规避“.trim.faa”等双后缀误判；
  - 映射口径与 10 发布脚本完全一致：严格按 (OG, Species) 取唯一记录；
  - 全量清洗：统一 strip 空白、去 CR、物种按 alias_map 规范；
  - 映射健全性：同一 (OG, Species) 多条/缺失都会明确报错，并给出 5 条示例；
  - I/O 稳健：CDS FASTA 按物种仅索引一次（缓存），大幅减少磁盘开销。

仅依赖：PyYAML（yaml）
"""

from __future__ import annotations
import sys, io, re, logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Iterable, Optional
import yaml

DEFAULT_CONFIG = "config.yaml"

# ====================== 基础工具：配置、日志、文件保障 ======================

def _expand_publish_placeholders(obj, publish_dir: str):
    """把对象中出现的 <publish_dir> 占位符替换为真路径。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def need_dir(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p

def need_file(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """日志写文件 + 同步屏幕；并保持 stdout/stderr 实时刷新。"""
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name); lg.setLevel(level); lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8"); fh.setFormatter(fmt); fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout);     sh.setFormatter(fmt); sh.setLevel(level)
    lg.addHandler(fh); lg.addHandler(sh)

    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x): self.s.write(x); self.s.flush(); return len(x)
    sys.stdout = _Flush(sys.stdout); sys.stderr = _Flush(sys.stderr)
    return lg

def banner(logger: logging.Logger, text: str):
    bar = "=" * max(10, len(text) + 2)
    logger.info(bar); logger.info(f" {text} "); logger.info(bar)

def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

# ====================== 专用工具：FASTA/OG/物种名 ======================

def og_id_from_filename(name: str) -> str:
    """
    从文件名中提取 OG 基础名（严格 ^OG\\d+），例如：
      OG0006410.trim.faa -> OG0006410
      OG0006410.faa      -> OG0006410
    """
    m = re.match(r"^(OG\d+)", name)
    if not m:
        raise RuntimeError(f"[ERR] 无法从文件名解析 OG：{name}")
    return m.group(1)

def read_fasta_headers(path: Path) -> List[str]:
    """读取 FASTA 的表头（去掉 '>'），按出现顺序返回。"""
    heads: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                heads.append(line[1:].strip().replace("\r", ""))
    return heads

def fasta_index_by_id(path: Path) -> Dict[str, str]:
    """
    将给定 FASTA 文件构建为字典：{header_id -> seq}。
    要求 header_id = '>' 后的整行（不含空白和描述），与发布包 cds_id 一致。
    """
    idx: Dict[str, str] = {}
    name: Optional[str] = None
    seq_parts: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n").replace("\r", "")
            if not line:
                continue
            if line.startswith(">"):
                # 收尾上一个
                if name is not None:
                    idx[name] = "".join(seq_parts).upper()
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if name is not None:
            idx[name] = "".join(seq_parts).upper()
    if not idx:
        raise RuntimeError(f"[ERR] 空的 CDS FASTA：{path}")
    return idx

def canon_species(sp: str, alias_map: Dict[str, str]) -> str:
    """按 alias_map 统一物种名；未命中则原样返回。"""
    s = sp.strip().replace("\r", "")
    return alias_map.get(s, s)

# ====================== 读取发布包映射：pep2cds_resolved.tsv ======================

def load_pep2cds_resolved(path: Path, alias_map: Dict[str, str]) -> Dict[Tuple[str, str], Tuple[str, str]]:
    """
    读取发布包的 pep2cds_resolved.tsv：
      列：OG, Species, protein_id, cds_id
    返回：{(OG, Species) -> (protein_id, cds_id)}；若同键多条则报错。
    """
    mp: Dict[Tuple[str, str], Tuple[str, str]] = {}
    dups: List[Tuple[str, str, str, str]] = []

    with open(path, "r", encoding="utf-8") as f:
        header = f.readline()
        cols = [c.strip().lower() for c in header.strip().split("\t")]
        try:
            i_og = cols.index("og")
            i_sp = cols.index("species")
            i_pid = cols.index("protein_id")
            i_cds = cols.index("cds_id")
        except ValueError:
            raise RuntimeError(f"[ERR] {path} 缺少列（需含 OG/Species/protein_id/cds_id）")

        for line in f:
            if not line.strip():
                continue
            parts = [x.strip().replace("\r", "") for x in line.rstrip("\n").split("\t")]
            if len(parts) < max(i_og, i_sp, i_pid, i_cds) + 1:
                continue
            og = parts[i_og]
            sp = canon_species(parts[i_sp], alias_map)
            pid = parts[i_pid]
            cds = parts[i_cds]
            key = (og, sp)
            if key in mp and mp[key] != (pid, cds):
                dups.append((og, sp, pid, cds))
            else:
                mp[key] = (pid, cds)

    if dups:
        ex = "\n".join(f"  - {a}\t{b}\t{c}\t{d}" for a,b,c,d in dups[:5])
        raise RuntimeError(f"[ERR] pep2cds_resolved.tsv 存在同一 (OG,Species) 多条记录（应唯一）：\n{ex}\n... 共 {len(dups)} 条冲突")
    if not mp:
        raise RuntimeError(f"[ERR] pep2cds_resolved.tsv 读取为空：{path}")
    return mp

# ====================== 主流程 ======================

def main():
    cfg = load_config()
    I, P = cfg["inputs"], cfg["paths"]
    alias_map = (cfg.get("species", {}) or {}).get("alias_map", {}) or {}

    logs_dir = Path(P["logs_dir"])
    LOG_FILE = logs_dir / "02_build_codon_inputs.log"
    log = get_logger("aphylo.02_build_codon_inputs", LOG_FILE)

    banner(log, "APhylo 02 — 构建 codon 输入")

    # 固定输入契约（零猜测）
    msa_dir     = need_dir(Path(I["sco_msa_dir"]),       "严格 SCO 对齐目录")
    msa_suffix  = I["sco_msa_suffix"]
    colmask_ok  = bool(I.get("colmask_dir"))  # 本步不读取 colmask，仅提示
    pep2cds_tsv = need_file(Path(I["pep2cds_map"]),      "pep→cds 解析表（发布生成）")
    cds_root    = need_dir(Path(I["cds_dir"]),           "CDS 仓库根目录")
    cds_suffix  = I.get("cds_suffix", ".fna")

    bt_dir   = ensure_dir(Path(P["bt_dir"]))
    tmp_dir  = ensure_dir(bt_dir / "tmp")
    order_dir= ensure_dir(bt_dir / "order")

    # 列举 OG MSA
    og_msa_files: List[Path] = sorted(msa_dir.glob(f"OG*{msa_suffix}"))
    if not og_msa_files:
        raise RuntimeError(f"[ERR] 未在 {msa_dir} 找到任何 OG*{msa_suffix} 文件")

    # 读取 pep2cds 解析（与 10 脚本口径一致）
    map_og_sp_to_ids = load_pep2cds_resolved(pep2cds_tsv, alias_map)
    log.info(f"映射条目数：{len(map_og_sp_to_ids)}（(OG,Species) → (protein_id, cds_id)）")
    if colmask_ok:
        log.info("检测到 colmask：后续 03 将按列掩码 ×3 同步回译。")

    # 物种 → CDS fasta 路径缓存；物种 → (id→seq) 索引缓存
    sp_to_cds_fa: Dict[str, Path] = {}
    sp_to_index: Dict[str, Dict[str, str]] = {}

    def cds_index_for(sp: str) -> Dict[str, str]:
        """为给定物种建立/返回 CDS 索引。"""
        sp_std = canon_species(sp, alias_map)
        if sp_std in sp_to_index:
            return sp_to_index[sp_std]
        fa = sp_to_cds_fa.get(sp_std)
        if fa is None:
            fa = cds_root / f"{sp_std}{cds_suffix}"
            need_file(fa, f"CDS FASTA 缺失（期望 {sp_std}{cds_suffix}）")
            sp_to_cds_fa[sp_std] = fa
        idx = fasta_index_by_id(fa)
        sp_to_index[sp_std] = idx
        return idx

    # 处理每个 OG
    n_ok = 0
    bad_missing_map: List[str] = []
    bad_missing_seq: List[str] = []

    for msa_path in og_msa_files:
        og = og_id_from_filename(msa_path.name)
        heads = [canon_species(h, alias_map) for h in read_fasta_headers(msa_path)]
        if not heads:
            bad_missing_map.append(f"{og}\t(empty_alignment)")
            continue

        # 构建 (Species 顺序) → (protein_id, cds_id)
        pair_list: List[Tuple[str, str, str]] = []  # (Species, protein_id, cds_id)
        miss_this: List[str] = []
        for sp in heads:
            key = (og, sp)
            ids = map_og_sp_to_ids.get(key)
            if ids is None:
                miss_this.append(sp)
            else:
                pid, cds = ids
                pair_list.append((sp, pid, cds))

        if miss_this:
            # 明确报错并给出示例，便于快速修复
            ex = ", ".join(miss_this[:5])
            raise RuntimeError(
                f"[ERR] pep2cds 缺映射 (OG={og})：缺少物种 {len(miss_this)} 个（示例：{ex}）。\n"
                f"请检查发布包 pep2cds_resolved.tsv 是否包含对应 (OG,Species) 记录。"
            )

        # 拉取 CDS，并按 heads 顺序写入 ordered_cds.fna
        out_fna = tmp_dir / f"{og}.ordered_cds.fna"
        out_ord = order_dir / f"{og}.order.tsv"

        n_written = 0
        with open(out_fna, "w", encoding="utf-8") as w_fa, \
             open(out_ord, "w", encoding="utf-8") as w_tab:
            w_tab.write("OG\tSpecies\tprotein_id\tcds_id\n")
            for sp, pid, cds_id in pair_list:
                idx = cds_index_for(sp)
                seq = idx.get(cds_id)
                if seq is None:
                    bad_missing_seq.append(f"{og}\t{sp}\t{cds_id}")
                    # 不中断，先记账，最后统一报错（便于一次修齐）
                    continue
                w_fa.write(f">{cds_id}\n")
                # 每 60 列换行（方便肉眼查看；后续脚本不依赖这个宽度）
                for i in range(0, len(seq), 60):
                    w_fa.write(seq[i:i+60] + "\n")
                w_tab.write(f"{og}\t{sp}\t{pid}\t{cds_id}\n")
                n_written += 1

        if n_written != len(heads):
            # 若有缺失的 CDS 序列，汇总后报错
            pass
        else:
            n_ok += 1

    # 汇总与错误处理
    if bad_missing_seq:
        ex = "\n".join(bad_missing_seq[:5])
        raise RuntimeError(
            "[ERR] 在 CDS 仓库中找不到以下条目（示例最多 5 条）：\n"
            "OG\tSpecies\tcds_id\n" + ex +
            f"\n... 共 {len(bad_missing_seq)} 条缺失。"
        )

    write_done(bt_dir / ".inputs.done")
    log.info(f"[DONE] 构建完成：OG={n_ok}/{len(og_msa_files)} 全部生成 ordered_cds.fna 与 order.tsv")

# ====================== 入口 ======================

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # 同时写屏幕（已由 get_logger 重定向 stdout/stderr 为“实时刷新”）
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)

