#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量统一蛋白 FASTA 表头为 `Species|SeqID`（仅改表头，不改序列）

使用场景：
- 你的 data/proteomes/*.fa* 目前没有 `|` 分隔的物种名前缀（如 >evm.model..., >Sco01g...）
- 期望以“文件名”作为物种名来源（Cyclina_sinensis.faa -> Cyclina_sinensis）

特点：
- 支持 .fa/.faa/.fasta 以及 .gz 压缩
- 默认就地覆盖保存；可选创建带时间戳的备份目录
- 不依赖命令行参数，所有参数在脚本顶部填写

安全性：
- 仅修改 FASTA 的 `>` 头行；氨基酸序列逐行原样写回
- 序列名唯一性检查：同一文件中若出现重复 SeqID，会自动加后缀 `_1/_2/...`

作者：小茹（零猜测 · 可复现）
"""

from __future__ import annotations
import os, re, sys, gzip, shutil, time
from pathlib import Path
from typing import Iterable

# ===================== 用户可配置区（无需命令行） =====================
PROTEOME_DIR = "data/proteomes"   # 蛋白目录（会遍历 *.fa *.faa *.fasta 以及 *.gz）
INPLACE_SAVE = True               # True=就地覆盖；False=写到 OUT_DIR
OUT_DIR = "data/proteomes_renamed"  # 当 INPLACE_SAVE=False 时生效
MAKE_BACKUP = True                # 就地覆盖时，是否先把原文件归档到 backup 目录
BACKUP_ROOT = "backup/proteomes_headers"  # 备份根目录

# 若需要将文件名映射为更规范的物种名，可在此编写字典；不写则使用文件名 stem
# 例如：{"S_constricta": "Sinonovacula_constricta"}
SPECIES_CANON_MAP: dict[str, str] = {}

# 支持的扩展名（不区分大小写）
FA_EXTS = (".fa", ".faa", ".fasta")
# =====================================================================

def open_auto(path: Path, mode: str):
    """自动识别是否为 .gz 并返回文件对象；mode 可为 'rt' 或 'wt'（文本）"""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[arg-type]
    return open(path, mode, encoding="utf-8")

def sanitize_species(name: str) -> str:
    """物种名清洗：去掉空白并把空格/制表替换为下划线"""
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    return name

def sanitize_seqid(s: str) -> str:
    """
    SeqID 清洗（与 pep2cds.tsv 对齐）：
    - 仅取第一个空白前的 token；
    - **保留** 竖线 '|'（如 gnl|WGS_AMQN|...），避免与 pep2cds 不一致；
    - 仅允许 [A-Za-z0-9._:-|]，其余字符一律替换为 '_'。
    """
    s = s.strip()
    if not s:
        return "NA"
    token = s.split()[0]
    token = re.sub(r"[^A-Za-z0-9._:\-\|]", "_", token)
    return token or "NA"

def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    """简易 FASTA 迭代器：yield (header_without_gt, seq)；header 不含 '>'"""
    with open_auto(path, "rt") as fh:
        header = None
        chunks = []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield (header, "".join(chunks))
                header = line[1:].rstrip("\n\r")
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield (header, "".join(chunks))

def write_fasta(path: Path, recs: Iterable[tuple[str, str]]):
    """将 (header_without_gt, seq) 序列写出为 FASTA"""
    with open_auto(path, "wt") as out:
        for h, seq in recs:
            out.write(">")
            out.write(h)
            out.write("\n")
            out.write(seq)

def choose_output_path(in_path: Path, out_dir: Path, inplace: bool) -> Path:
    if inplace:
        return in_path
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / in_path.name

def do_backup(src: Path, backup_root: Path):
    """把原始文件按时间戳移动到 backup 目录，保持同名"""
    ts = time.strftime("%Y%m%d_%H%M%S")
    target_dir = backup_root / ts
    target_dir.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(target_dir / src.name))

def process_one_file(fa: Path, species_name: str, inplace: bool, out_dir: Path|None, do_bak: bool, bak_root: Path):
    """处理单个 FASTA：把每条 header 变为 Species|SeqID"""
    species = sanitize_species(species_name)
    # 输出路径（就地或 OUT_DIR）
    out_path = choose_output_path(fa, out_dir or fa.parent, inplace=False)  # 先写到临时路径
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    # 唯一性检查集合
    seen = set()
    new_records = []

    for raw_h, seq in iter_fasta(fa):
        seqid = sanitize_seqid(raw_h)
        new_h = f"{species}|{seqid}"
        # 避免同文件内重复
        if new_h in seen:
            i = 1
            while f"{new_h}_{i}" in seen:
                i += 1
            new_h = f"{new_h}_{i}"
        seen.add(new_h)
        new_records.append((new_h, seq))

    # 写 tmp，再原子替换
    write_fasta(tmp_path, new_records)

    if inplace:
        # 备份
        if do_bak:
            do_backup(fa, bak_root)
        # 覆盖到原路径
        tmp_path.replace(fa)
    else:
        # 输出到 OUT_DIR
        out_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path.replace(out_path)

def main():
    root = Path(PROTEOME_DIR)
    assert root.is_dir(), f"[ERR] 目录不存在：{root}"

    # 收集候选文件
    files = []
    for p in sorted(root.iterdir()):
        low = p.name.lower()
        if any(low.endswith(ext) for ext in FA_EXTS) or low.endswith(tuple(ext + ".gz" for ext in FA_EXTS)):
            files.append(p)
    assert files, f"[ERR] 在 {root} 未找到 *.fa/*.faa/*.fasta[.gz]"

    # 准备输出与备份
    out_dir = Path(OUT_DIR) if not INPLACE_SAVE else None
    bak_root = Path(BACKUP_ROOT)

    # 逐个处理
    n = 0
    for fa in files:
        stem = fa.name
        # 去掉多重扩展名（如 .faa.gz）
        s = stem
        for ext in (".gz",) + FA_EXTS:
            s = re.sub(re.escape(ext) + r"$", "", s, flags=re.IGNORECASE)
        file_stem = s  # 原始文件名去扩展
        species = SPECIES_CANON_MAP.get(file_stem, file_stem)

        print(f"[DO] {fa.name}  ->  Species={species}")
        process_one_file(
            fa=fa,
            species_name=species,
            inplace=INPLACE_SAVE,
            out_dir=out_dir,
            do_bak=MAKE_BACKUP and INPLACE_SAVE,
            bak_root=bak_root
        )
        n += 1

    print(f"[OK] 共处理 {n} 个 FASTA；模式：{'就地覆盖' if INPLACE_SAVE else '输出到 ' + str(out_dir)}。")
    if INPLACE_SAVE and MAKE_BACKUP:
        print(f"[INFO] 原始文件备份位置：{BACKUP_ROOT}/<时间戳>/")

if __name__ == "__main__":
    main()

