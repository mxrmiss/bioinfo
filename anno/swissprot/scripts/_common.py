#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import re
import sys
import hashlib
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def find_project_root(start: Optional[Path] = None) -> Path:
    """向上查找包含 config.yaml 的目录作为项目根目录"""
    if start is None:
        start = Path.cwd()
    start = start.resolve()
    cur = start
    while True:
        if (cur / "config.yaml").is_file():
            return cur
        if cur.parent == cur:
            raise FileNotFoundError("Cannot find project root containing config.yaml")
        cur = cur.parent


def read_yaml_config(root: Path) -> Dict[str, Any]:
    cfg_path = root / "config.yaml"
    try:
        import yaml  # type: ignore
    except Exception as ex:
        raise RuntimeError("PyYAML is required (import yaml failed). Please install pyyaml.") from ex
    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("config.yaml must be a YAML mapping (dict)")
    return cfg


def ensure_no_duplicates(items: List[str], name: str) -> None:
    if len(items) != len(set(items)):
        raise ValueError(f"{name} contains duplicates: {items}")


def clamp01(x: float) -> float:
    if x < 0.0:
        return 0.0
    if x > 1.0:
        return 1.0
    return x


def md5sum_file(p: Path) -> str:
    h = hashlib.md5()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def count_fasta_seqs(p: Path) -> int:
    n = 0
    with p.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def _read_current_rel_file(root: Path) -> Optional[str]:
    """
    读取 results/01_db/CURRENT_REL.txt 作为当前 REL 的权威来源。
    - 文件不存在/为空 -> 返回 None
    - 存在且内容合法 -> 返回 REL 字符串
    """
    cur = root / "results" / "01_db" / "CURRENT_REL.txt"
    if not cur.is_file():
        return None
    s = cur.read_text(encoding="utf-8", errors="replace").strip()
    if not s:
        return None
    # REL 通常类似 2025_04（保险起见只允许数字和下划线）
    if not re.fullmatch(r"[0-9_]+", s):
        raise ValueError(f"Invalid REL in {cur}: {s!r}")
    return s


def _validate_rel_points_to_db(root: Path, rel: str) -> None:
    """
    校验 CURRENT_REL.txt 指向的 REL 至少满足“项目内可用”的最低条件。
    这里不强制要求 dmnd 一定存在（因为有人可能先跑 db 产物清单，再构建库），
    但至少要求 results/01_db/{REL} 目录存在（01_db.py 会创建）。
    """
    rel_dir = root / "results" / "01_db" / rel
    if not rel_dir.exists():
        raise FileNotFoundError(
            f"CURRENT_REL points to missing directory: {rel_dir}. "
            f"Please run scripts/01_db.py or fix results/01_db/CURRENT_REL.txt"
        )


def parse_rel_from_01_db(root: Path) -> str:
    """
    解析当前 REL（UniProt/Swiss-Prot Release）。

    新规则（与你的 Snakefile 完全一致）：
    1) 若 results/01_db/CURRENT_REL.txt 存在且有效，则它是唯一权威来源，直接返回。
       即使 results/01_db/ 下存在多个历史 REL，也不再要求“唯一”。

    兼容兜底：
    2) 若 CURRENT_REL.txt 不存在/为空，则回退到旧逻辑：扫描 results/01_db/* 提取 REL，
       并要求解析出的 REL 唯一。
    """
    # 1) 权威来源：CURRENT_REL.txt
    cur_rel = _read_current_rel_file(root)
    if cur_rel:
        _validate_rel_points_to_db(root, cur_rel)
        return cur_rel

    # 2) 兜底：旧逻辑扫描（仅当没有 CURRENT_REL.txt 或为空时使用）
    db_root = root / "results" / "01_db"
    if not db_root.exists():
        raise FileNotFoundError("results/01_db not found. Please build DB first.")

    rels = set()

    for mf in db_root.glob("*/db_manifest.tsv"):
        try:
            with mf.open("r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if line.startswith("REL\t"):
                        rels.add(line.split("\t", 1)[1].strip())
                        break
        except Exception:
            pass

    if not rels:
        for rf in db_root.glob("*/reldate.txt"):
            try:
                with rf.open("r", encoding="utf-8", errors="replace") as f:
                    for line in f:
                        s = line.strip()

                        # 兼容新/旧格式：Release 不一定在行首
                        m = re.search(r"\bRelease\s+([0-9_]+)\b", s)
                        if m:
                            rels.add(m.group(1))
                            break

                        # 兜底：从 Swiss-Prot 行取
                        m2 = re.search(r"Swiss-Prot\s+Release\s+([0-9_]+)\b", s)
                        if m2:
                            rels.add(m2.group(1))
                            break
            except Exception:
                pass

    if len(rels) != 1:
        raise ValueError(
            f"REL must be unique when CURRENT_REL.txt is missing/empty. "
            f"Found: {sorted(rels)} in results/01_db/*"
        )
    rel = next(iter(rels))
    _validate_rel_points_to_db(root, rel)
    return rel


def get_species_list(cfg: Dict[str, Any]) -> List[str]:
    species = cfg.get("species_prefixes", [])
    if not isinstance(species, list) or not all(isinstance(x, str) and x for x in species):
        raise ValueError("config.yaml: species_prefixes must be a non-empty list of strings")
    ensure_no_duplicates(species, "config.yaml: species_prefixes")
    for s in species:
        if "/" in s:
            raise ValueError(f"Invalid species_id contains '/': {s}")
    return species


def get_only_species_from_env() -> Optional[str]:
    s = os.environ.get("SPECIES_ID", "").strip()
    return s if s else None


def resolve_db_dmnd(root: Path, rel: str) -> Tuple[Path, str]:
    """返回 (db_dmnd_path, db_prefix)"""
    expected = root / "data" / "db" / "diamond" / f"uniprot_sprot_{rel}.dmnd"
    if expected.is_file() and expected.stat().st_size > 0:
        return expected, str(expected.with_suffix(""))

    dmnd_dir = root / "data" / "db" / "diamond"
    dmnds = sorted(dmnd_dir.glob("*.dmnd"))
    dmnds = [p for p in dmnds if p.is_file() and p.stat().st_size > 0]
    if len(dmnds) == 1:
        p = dmnds[0]
        return p, str(p.with_suffix(""))
    if not dmnds:
        raise FileNotFoundError(f"No .dmnd found under {dmnd_dir}")
    raise ValueError(f"Multiple .dmnd found under {dmnd_dir}; expected only one or {expected.name}")


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

