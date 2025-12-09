#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ID 一致性体检脚本（仅检查，不修改任何原始文件）

功能概要：
1. 自动扫描 raw_data/protein 下的物种前缀；
2. 对每个物种：
   - 读取蛋白 FASTA 头信息，抽取蛋白 ID 集合；
   - 读取对应 GFF 文件（.gff / .gff3 / .gz），抽取属性列中的各种 ID（ID/Parent/protein_id 等）；
   - 对蛋白 ID 依次应用多种规则（原样 / 去掉 .数字 / 去掉 -RA 等），统计与 GFF ID 的匹配情况；
   - 输出每个规则的匹配数量、覆盖率，以及若干匹配示例，方便肉眼判断哪条规则是“真”。

产物目录：
  tmp/idcheck/
    idcheck_summary.tsv         # 总表：每个物种 × 每条规则的匹配统计
    <species>/pep_ids.txt       # 该物种蛋白 ID 原始列表（去重）
    <species>/gff_ids.txt       # 该物种 GFF 中抽取的 ID 集合（去重）
    <species>/rule_examples.tsv # 该物种每条规则的命中示例（便于人工确认）

使用方式：
  1. 将本脚本保存到 project/magic/scripts/idcheck_protein_gff.py
  2. cd ~/project/magic
  3. python scripts/idcheck_protein_gff.py

然后把 tmp/idcheck/idcheck_summary.tsv 以及某些物种的 rule_examples.tsv 内容贴给我看，
我们就可以精确写死每个物种的 ID 规则，后续共线性脚本就不再“瞎猜”。
"""

# ===================== 用户参数区（可按需修改） =====================

# 项目根目录（默认自动用脚本的上上级目录，如果你想指定就手动写死绝对路径）
PROJECT_ROOT = None  # 例如 "/home/mxrmiss/project/magic"

# 相对目录（相对于 PROJECT_ROOT）
PROTEIN_DIR = "raw_data/protein"
GFF_DIR = "raw_data/gff"
OUTPUT_DIR = "tmp/idcheck"

# 每个物种最多读取多少条蛋白 ID 用于统计（节省内存，通常几万条也没问题）
MAX_PROTEIN_IDS = 50000

# 为了展示规律，每条规则最多输出多少个“命中示例”
MAX_EXAMPLES_PER_RULE = 20

# 如果只想检查部分物种，可以在这里写前缀列表；留空列表则表示“自动检测全部”
SPECIES_WHITELIST = []  # 例如 ["Sinonovacula_constricta", "Sinonovacula_rivularis"]

# ===================== 下面是脚本主体，皇上不要改 =====================

import os
import re
import gzip
from pathlib import Path
from collections import defaultdict

def resolve_project_root() -> Path:
    """解析项目根目录：如果用户没填 PROJECT_ROOT，就用脚本所在位置的上上级目录。"""
    if PROJECT_ROOT:
        return Path(PROJECT_ROOT).resolve()
    # 本脚本一般放在 project/magic/scripts 下
    here = Path(__file__).resolve()
    return here.parent.parent


def open_maybe_gzip(path: Path, mode: str = "rt", encoding: str = "utf-8"):
    """
    根据文件后缀自动选择 open 或 gzip.open。
    mode 一般用 "rt"（文本读）。
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode=mode, encoding=encoding, errors="replace")
    else:
        return open(path, mode=mode, encoding=encoding, errors="replace")


def detect_species_from_protein(protein_dir: Path):
    """从 protein 目录推断物种前缀列表（去掉第一个点之后的所有部分）。"""
    species = set()
    for f in protein_dir.iterdir():
        if not f.is_file():
            continue
        name = f.name
        # 只考虑常见蛋白后缀
        if not any(name.endswith(ext) for ext in (".faa", ".fa", ".fasta", ".pep.fa", ".pep.fasta")):
            continue
        prefix = name.split(".")[0]
        species.add(prefix)
    return sorted(species)


def find_gff_for_species(gff_dir: Path, sp: str):
    """按优先级寻找与物种前缀匹配的 GFF 文件。找不到就返回 None。"""
    candidates = [
        gff_dir / f"{sp}.gff3",
        gff_dir / f"{sp}.gff",
        gff_dir / f"{sp}.gff3.gz",
        gff_dir / f"{sp}.gff.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def find_protein_for_species(protein_dir: Path, sp: str):
    """按优先级寻找蛋白文件。"""
    candidates = [
        protein_dir / f"{sp}.pep.fa",
        protein_dir / f"{sp}.pep.fasta",
        protein_dir / f"{sp}.faa",
        protein_dir / f"{sp}.fa",
        protein_dir / f"{sp}.fasta",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def load_protein_ids(protein_path: Path, max_ids: int):
    """读取蛋白 FASTA 头行，抽取 ID（取 '>' 后第一个空格前的 token），去重。"""
    ids = []
    seen = set()
    count = 0
    with open_maybe_gzip(protein_path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if not header:
                continue
            token = header.split()[0]
            if token not in seen:
                seen.add(token)
                ids.append(token)
                count += 1
                if count >= max_ids:
                    break
    return ids


ATTR_SPLIT_RE = re.compile(r"[;,]")  # 分割属性字段用


def parse_gff_ids(gff_path: Path):
    """
    解析 GFF，收集属性列中的各种 ID 值。
    我们主要关注的字段：ID, Parent, Name, gene_id, transcript_id, protein_id 等。
    返回：set(all_ids)
    """
    id_keys = {"ID", "Parent", "Name", "gene_id", "transcript_id", "protein_id", "proteinId", "mRNA", "transcript"}
    ids = set()
    with open_maybe_gzip(gff_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attr_field = parts[8]
            # 属性通常形如 "ID=xxx;Parent=yyy;Name=zzz"
            for item in ATTR_SPLIT_RE.split(attr_field):
                if "=" in item:
                    key, val = item.split("=", 1)
                    key = key.strip()
                    val = val.strip()
                    if key in id_keys and val:
                        # 有的字段可能有逗号分隔多个 ID
                        for v in val.split(","):
                            v = v.strip()
                            if v:
                                ids.add(v)
    return ids


def transform_id(raw_id: str, rule: str) -> str:
    """
    对单个 ID 应用转换规则。
    这里只是用于“尝试”，不会改动原始文件。
    """
    if rule == "raw":
        return raw_id
    # 去掉最后的 .数字 （如 gene.1）
    if rule == "strip_dot_number":
        return re.sub(r"\.[0-9]+$", "", raw_id)
    # 去掉 -R.. （如 Teg00005373-RA）
    if rule == "strip_dash_R":
        return re.sub(r"-R[A-Za-z0-9]+$", "", raw_id)
    # 去掉常见的 .t1 / .mRNA1 / .p1 这类后缀
    if rule == "strip_transcript_suffix":
        return re.sub(r"\.(?:t|T|mrna|mRNA|p|P)[0-9]+$", "", raw_id)
    # 组合：先去 dashR 再去 dot_number（防止顺序差异）
    if rule == "strip_dashR_then_dot":
        x = re.sub(r"-R[A-Za-z0-9]+$", "", raw_id)
        x = re.sub(r"\.[0-9]+$", "", x)
        return x
    # 默认：不变
    return raw_id


RULES = [
    "raw",
    "strip_dot_number",
    "strip_dash_R",
    "strip_transcript_suffix",
    "strip_dashR_then_dot",
]


def analyze_species(sp: str, protein_path: Path, gff_path: Path, out_dir: Path):
    """
    对单个物种进行 ID 一致性分析：
    - 读取蛋白 ID
    - 读取 GFF ID
    - 对每条规则计算覆盖率和示例
    返回：list[dict]，每个 dict 是一条规则的统计结果。
    """
    sp_dir = out_dir / sp
    sp_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] === 分析物种 {sp} ===")
    print(f"[INFO] 蛋白文件: {protein_path}")
    print(f"[INFO] GFF 文件: {gff_path}")

    # 1. 载入蛋白 ID
    pep_ids = load_protein_ids(protein_path, MAX_PROTEIN_IDS)
    n_pep = len(pep_ids)
    print(f"[INFO] 蛋白 ID 数量（去重，最多 {MAX_PROTEIN_IDS} 条）: {n_pep}")

    with open(sp_dir / "pep_ids.txt", "w", encoding="utf-8") as fh:
        for pid in pep_ids:
            fh.write(pid + "\n")

    # 2. 载入 GFF ID
    gff_ids = parse_gff_ids(gff_path)
    n_gff = len(gff_ids)
    print(f"[INFO] GFF 中抽取的 ID 数量: {n_gff}")

    with open(sp_dir / "gff_ids.txt", "w", encoding="utf-8") as fh:
        for gid in sorted(gff_ids):
            fh.write(gid + "\n")

    gff_id_set = gff_ids  # 已经是 set

    results = []
    examples_per_rule = defaultdict(list)

    # 3. 对每条规则进行匹配统计
    for rule in RULES:
        transformed = [transform_id(pid, rule) for pid in pep_ids]
        transformed_set = set(transformed)
        # 与 GFF ID 的交集
        matches = transformed_set & gff_id_set
        n_match = len(matches)
        coverage = n_match / n_pep if n_pep > 0 else 0.0

        # 选出若干匹配示例（记录原始 ID -> 转换 ID）
        if n_match > 0:
            match_set = set(matches)
            for raw_id, t_id in zip(pep_ids, transformed):
                if t_id in match_set:
                    examples_per_rule[rule].append((raw_id, t_id))
                    if len(examples_per_rule[rule]) >= MAX_EXAMPLES_PER_RULE:
                        break

        print(f"[INFO] 规则={rule:22s} 命中数={n_match:6d} 覆盖率={coverage:7.3%}")

        results.append(
            {
                "species": sp,
                "n_pep_ids": n_pep,
                "n_gff_ids": n_gff,
                "rule": rule,
                "n_match": n_match,
                "coverage": coverage,
            }
        )

    # 4. 写出该物种的规则示例表
    ex_path = sp_dir / "rule_examples.tsv"
    with open(ex_path, "w", encoding="utf-8") as fh:
        fh.write("rule\traw_protein_id\ttransformed_id\n")
        for rule in RULES:
            examples = examples_per_rule.get(rule, [])
            for raw_id, t_id in examples:
                fh.write(f"{rule}\t{raw_id}\t{t_id}\n")

    print(f"[INFO] 示例写出: {ex_path}")
    print()
    return results


def main():
    project_root = resolve_project_root()
    protein_dir = (project_root / PROTEIN_DIR).resolve()
    gff_dir = (project_root / GFF_DIR).resolve()
    out_dir = (project_root / OUTPUT_DIR).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] 项目根目录: {project_root}")
    print(f"[INFO] 蛋白目录:   {protein_dir}")
    print(f"[INFO] GFF 目录:    {gff_dir}")
    print(f"[INFO] 输出目录:    {out_dir}")
    print()

    if not protein_dir.exists():
        raise SystemExit(f"[ERROR] 蛋白目录不存在: {protein_dir}")
    if not gff_dir.exists():
        raise SystemExit(f"[ERROR] GFF 目录不存在: {gff_dir}")

    # 1. 物种列表
    if SPECIES_WHITELIST:
        species_list = SPECIES_WHITELIST
        print(f"[INFO] 使用用户指定物种列表: {species_list}")
    else:
        species_list = detect_species_from_protein(protein_dir)
        print(f"[INFO] 自动检测到蛋白物种数={len(species_list)}: {', '.join(species_list)}")
    print()

    all_results = []

    # 2. 逐物种体检
    for sp in species_list:
        protein_path = find_protein_for_species(protein_dir, sp)
        gff_path = find_gff_for_species(gff_dir, sp)

        if protein_path is None:
            print(f"[WARN] 物种 {sp} 找不到蛋白文件，跳过。")
            continue
        if gff_path is None:
            print(f"[WARN] 物种 {sp} 找不到 GFF 文件，跳过。")
            continue

        try:
            res = analyze_species(sp, protein_path, gff_path, out_dir)
            all_results.extend(res)
        except Exception as e:
            print(f"[ERROR] 分析物种 {sp} 时出错: {e}")
            continue

    # 3. 写出总汇总表
    summary_path = out_dir / "idcheck_summary.tsv"
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("species\tn_pep_ids\tn_gff_ids\trule\tn_match\tcoverage\n")
        for row in all_results:
            fh.write(
                f"{row['species']}\t{row['n_pep_ids']}\t{row['n_gff_ids']}\t"
                f"{row['rule']}\t{row['n_match']}\t{row['coverage']:.6f}\n"
            )

    print(f"[INFO] 全部物种规则统计写出: {summary_path}")
    print("[INFO] 完成。你可以先 cat 这个 summary 和某些物种的 rule_examples.tsv 发给我。")


if __name__ == "__main__":
    main()

