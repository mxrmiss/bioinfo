#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv
import re
from typing import Dict, List, Optional, Tuple
from openpyxl import Workbook


GENOMIC_DIR = Path("data/genomic")

QC_DIR = Path("results/qc")
QUAST_DIR = QC_DIR / "quast"
BUSCO_DIR = QC_DIR / "busco"
AGAT_DIR = QC_DIR / "agat"
SEQKIT_DIR = QC_DIR / "seqkit"

RESULTS_DIR = Path("results")
OUTPUT_DIR = RESULTS_DIR / "supp_table"

ANNOT_MATCH_TSV = RESULTS_DIR / "proteomes" / "annot_genome_match.tsv"

MATRIX_TSV_CANDIDATES = [
    RESULTS_DIR / "reports" / "matrix.tsv",
    RESULTS_DIR / "publish" / "aphylo_ready" / "matrix.tsv",
]

TREE_FILE_CANDIDATES = [
    RESULTS_DIR / "trees" / "species_tree.nwk",
    RESULTS_DIR / "publish" / "aphylo_ready" / "species_tree.nwk",
]

MAIN_TABLE_TSV = "Supplementary_Table_S1.tsv"
MAIN_TABLE_XLSX = "Supplementary_Table_S1.xlsx"
QC_TABLE_TSV = "supp_table_qc.tsv"
QC_TABLE_XLSX = "supp_table_qc.xlsx"

GENOME_SUFFIXES = (".fa", ".fa.gz", ".fna", ".fna.gz", ".fasta", ".fasta.gz")


def ensure_dir(path: Path) -> None:
    """确保目录存在。"""
    path.mkdir(parents=True, exist_ok=True)


def try_read_text(path: Path) -> str:
    """读取文本文件；失败时返回空字符串。"""
    if not path.exists():
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def read_tsv_dicts(path: Path) -> List[Dict[str, str]]:
    """读取 TSV 为字典列表。"""
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def write_tsv(rows: List[Dict[str, str]], out_path: Path, fieldnames: List[str]) -> None:
    """写出 TSV 文件。"""
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_xlsx(rows: List[Dict[str, str]], out_path: Path, fieldnames: List[str]) -> None:
    """写出 XLSX 文件。使用 write_only 模式提高速度。"""
    wb = Workbook(write_only=True)
    ws = wb.create_sheet(title="Sheet1")
    ws.append(fieldnames)

    for row in rows:
        ws.append([row.get(col, "") for col in fieldnames])

    wb.save(out_path)


def pick_first_existing(paths: List[Path]) -> Optional[Path]:
    """从候选路径中返回第一个存在的文件。"""
    for p in paths:
        if p.exists():
            return p
    return None


def parse_bool_to_yes_no(value: Optional[bool]) -> str:
    """布尔值转 Yes/No。"""
    if value is None:
        return ""
    return "Yes" if value else "No"


def format_percent(value: Optional[float], digits: int = 2) -> str:
    """格式化百分比。"""
    if value is None:
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return ""


def normalize_spaces(text: str) -> str:
    """压缩多余空白。"""
    return re.sub(r"\s+", " ", text.strip())


def strip_allowed_suffix(filename: str, allowed_suffixes: Tuple[str, ...]) -> Optional[str]:
    """去掉允许的后缀，返回物种名。"""
    for suffix in sorted(allowed_suffixes, key=len, reverse=True):
        if filename.endswith(suffix):
            return filename[:-len(suffix)]
    return None


def discover_species_from_genomic_dir() -> List[str]:
    """从 data/genomic 自动发现物种名。"""
    if not GENOMIC_DIR.exists():
        raise FileNotFoundError(f"[ERROR] 目录不存在：{GENOMIC_DIR}")

    species_list = []
    for p in sorted(GENOMIC_DIR.iterdir()):
        if not p.is_file():
            continue
        species = strip_allowed_suffix(p.name, GENOME_SUFFIXES)
        if species is not None:
            species_list.append(species)

    species_list = sorted(set(species_list))

    if not species_list:
        raise RuntimeError(f"[ERROR] 在 {GENOMIC_DIR} 中未发现任何可识别的 genome 文件")

    return species_list


def load_annot_match_map(path: Path) -> Dict[str, Dict[str, str]]:
    """读取 annot_genome_match.tsv。"""
    rows = read_tsv_dicts(path)
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        species = row.get("species", "").strip()
        if species:
            result[species] = row
    return result


def load_matrix_map(path: Path) -> Dict[str, Dict[str, str]]:
    """读取 matrix.tsv。"""
    rows = read_tsv_dicts(path)
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        species = row.get("species", "").strip()
        if species:
            result[species] = row
    return result


def parse_species_from_newick(path: Path) -> List[str]:
    """从 Newick 树中粗略提取叶节点物种名。"""
    text = try_read_text(path).strip()
    if not text:
        return []

    text = re.sub(r":[^,()]+", "", text)
    tokens = re.split(r"[(),;]", text)

    seen = set()
    species_list = []
    for token in tokens:
        name = token.strip()
        if name and name not in seen:
            seen.add(name)
            species_list.append(name)

    return species_list


def normalize_quast_key(key: str) -> str:
    """标准化 QUAST 指标名称。"""
    return normalize_spaces(key).lower()


def parse_quast_report_tsv(path: Path) -> Dict[str, str]:
    """读取 QUAST 的 report.tsv。"""
    result = {
        "assembly_size_bp": "",
        "gc_percent": "",
        "sequence_count": "",
        "n50_bp": "",
    }

    if not path.exists():
        return result

    rows: List[List[str]] = []
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row:
                rows.append(row)

    if len(rows) < 2:
        return result

    for row in rows:
        if len(row) < 2:
            continue

        key = normalize_quast_key(row[0])
        val = row[1].strip()

        if key in {"total length", "total length (>= 0 bp)"}:
            result["assembly_size_bp"] = val
        elif key in {"gc (%)", "gc"}:
            result["gc_percent"] = val
        elif key in {"# contigs", "# contigs (>= 0 bp)", "# scaffolds", "# sequences"}:
            result["sequence_count"] = val
        elif key == "n50":
            result["n50_bp"] = val

    return result


def load_quast_for_species(species: str) -> Dict[str, str]:
    """读取单个物种的 QUAST 结果。"""
    return parse_quast_report_tsv(QUAST_DIR / species / "report.tsv")


def find_busco_summary_file(species: str) -> Optional[Path]:
    """在 BUSCO 输出目录中查找 short_summary 文件。"""
    species_dir = BUSCO_DIR / species
    if not species_dir.exists():
        return None

    candidates = sorted(species_dir.glob("short_summary*.txt"))
    if candidates:
        return candidates[0]
    return None


def parse_busco_short_summary(path: Path) -> Dict[str, str]:
    """读取 BUSCO short_summary*.txt。"""
    result = {
        "complete_busco_percent": "",
        "single_copy_busco_percent": "",
        "duplicated_busco_percent": "",
        "fragmented_busco_percent": "",
        "missing_busco_percent": "",
    }

    text = try_read_text(path)
    if not text:
        return result

    m = re.search(
        r"C:(?P<C>[\d.]+)%\s*\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],\s*F:(?P<F>[\d.]+)%,\s*M:(?P<M>[\d.]+)%",
        text
    )
    if not m:
        return result

    result["complete_busco_percent"] = m.group("C")
    result["single_copy_busco_percent"] = m.group("S")
    result["duplicated_busco_percent"] = m.group("D")
    result["fragmented_busco_percent"] = m.group("F")
    result["missing_busco_percent"] = m.group("M")
    return result


def load_busco_for_species(species: str) -> Dict[str, str]:
    """读取单个物种的 BUSCO 结果。"""
    summary = find_busco_summary_file(species)
    if summary is None:
        return {
            "complete_busco_percent": "",
            "single_copy_busco_percent": "",
            "duplicated_busco_percent": "",
            "fragmented_busco_percent": "",
            "missing_busco_percent": "",
        }
    return parse_busco_short_summary(summary)


def parse_agat_stats_text(path: Path) -> Dict[str, str]:
    """精确解析 AGAT stats.txt 中 mrna 主统计区块的 gene/mrna/cds 数量。逐行扫描以提高速度。"""
    result = {
        "gene_count": "",
        "mrna_count": "",
        "cds_count": "",
    }

    if not path.exists():
        return result

    in_mrna_block = False

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            line_stripped = line.strip()
            line_lower = line_stripped.lower()

            # 进入主统计区块
            if re.match(r"^-+\s*mrna\s*-+$", line_lower):
                in_mrna_block = True
                continue

            # 进入 mrna 区块后，如果遇到新的大标题区块，则停止
            if in_mrna_block and re.match(r"^-+\s*[a-z0-9_]+\s*-+$", line_lower):
                break

            if not in_mrna_block:
                continue

            m_gene = re.match(r"^Number of gene\s+(\d+)\s*$", line_stripped, flags=re.IGNORECASE)
            if m_gene:
                result["gene_count"] = m_gene.group(1)
                if result["gene_count"] and result["mrna_count"] and result["cds_count"]:
                    break
                continue

            m_mrna = re.match(r"^Number of mrna\s+(\d+)\s*$", line_stripped, flags=re.IGNORECASE)
            if m_mrna:
                result["mrna_count"] = m_mrna.group(1)
                if result["gene_count"] and result["mrna_count"] and result["cds_count"]:
                    break
                continue

            m_cds = re.match(r"^Number of cds\s+(\d+)\s*$", line_stripped, flags=re.IGNORECASE)
            if m_cds:
                result["cds_count"] = m_cds.group(1)
                if result["gene_count"] and result["mrna_count"] and result["cds_count"]:
                    break
                continue

    return result


def load_agat_for_species(species: str) -> Dict[str, str]:
    """读取单个物种的 AGAT 统计结果。"""
    return parse_agat_stats_text(AGAT_DIR / species / "stats.txt")


def parse_seqkit_stats_tsv(path: Path) -> Dict[str, str]:
    """读取 SeqKit stats TSV。"""
    result = {
        "sequence_count": "",
        "assembly_size_bp": "",
        "gc_percent": "",
        "n50_bp": "",
        "record_count": "",
    }

    rows = read_tsv_dicts(path)
    if not rows:
        return result

    row = rows[0]

    def pick(*keys: str) -> str:
        for k in keys:
            if k in row and str(row[k]).strip() != "":
                return str(row[k]).strip()
        return ""

    result["sequence_count"] = pick("num_seqs", "num.seq", "num_seqs:")
    result["assembly_size_bp"] = pick("sum_len", "sum.len", "sum_len:")
    result["gc_percent"] = pick("GC(%)", "gc(%)", "GC%", "gc%")
    result["n50_bp"] = pick("N50", "n50")
    result["record_count"] = pick("num_seqs", "num.seq", "num_seqs:")

    return result


def load_seqkit_genome_for_species(species: str) -> Dict[str, str]:
    """读取单个物种 genome 的 SeqKit 统计。"""
    return parse_seqkit_stats_tsv(SEQKIT_DIR / f"{species}.genome.tsv")


def load_seqkit_protein_for_species(species: str) -> Dict[str, str]:
    """读取单个物种 protein 的 SeqKit 统计。"""
    return parse_seqkit_stats_tsv(SEQKIT_DIR / f"{species}.protein.tsv")


def merge_primary_metric(primary: str, secondary: str) -> str:
    """合并两个候选值，优先使用 primary。"""
    if str(primary).strip() != "":
        return str(primary).strip()
    if str(secondary).strip() != "":
        return str(secondary).strip()
    return ""


def derive_contig_scaffold_counts(sequence_count: str, assembly_level: str) -> Tuple[str, str]:
    """根据 assembly_level 将 sequence_count 映射到 contig/scaffold count。"""
    level = str(assembly_level).strip().lower()
    count = str(sequence_count).strip()

    if count == "":
        return "", ""
    if level == "contig":
        return count, ""
    if level in {"scaffold", "chromosome"}:
        return "", count
    return "", ""


def status_is_pass(match_status: str) -> bool:
    """判断 annot match status 是否属于通过状态。"""
    status_upper = str(match_status).strip().upper()
    if status_upper == "":
        return False
    if status_upper.startswith("PASS"):
        return True
    if status_upper in {"OK", "PASSED"}:
        return True
    return False


def collect_species_row(
    species: str,
    annot_match_map: Dict[str, Dict[str, str]],
    matrix_map: Dict[str, Dict[str, str]],
    species_in_tree: set
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """整合单个物种的信息。"""

    quast_info = load_quast_for_species(species)
    busco_info = load_busco_for_species(species)
    agat_info = load_agat_for_species(species)
    seqkit_genome = load_seqkit_genome_for_species(species)
    seqkit_protein = load_seqkit_protein_for_species(species)

    assembly_size_bp = merge_primary_metric(quast_info["assembly_size_bp"], seqkit_genome["assembly_size_bp"])
    gc_percent = merge_primary_metric(quast_info["gc_percent"], seqkit_genome["gc_percent"])
    sequence_count = merge_primary_metric(quast_info["sequence_count"], seqkit_genome["sequence_count"])
    n50_bp = merge_primary_metric(quast_info["n50_bp"], seqkit_genome["n50_bp"])
    protein_count = seqkit_protein["record_count"]

    annot_match = annot_match_map.get(species, {})
    matrix_info = matrix_map.get(species, {})

    match_cov = annot_match.get("cov_frac", "")
    match_status = annot_match.get("status", "")
    pep_kept = annot_match.get("pep_kept", "")
    annot_note = annot_match.get("note", "")

    nongap_sites = matrix_info.get("nongap_sites", "")
    occupancy_percent = matrix_info.get("occupancy_percent", "")

    included_in_matrix = species in matrix_map
    included_in_tree = species in species_in_tree if species_in_tree else None

    assembly_id = ""
    source = ""
    assembly_level = ""
    notes = ""

    contig_count, scaffold_count = derive_contig_scaffold_counts(sequence_count, assembly_level)
    longest_protein_count = pep_kept if str(pep_kept).strip() != "" else protein_count

    passed_phylo_qc = None
    if match_status or included_in_matrix:
        passed_phylo_qc = status_is_pass(match_status) and included_in_matrix

    if passed_phylo_qc is True and included_in_tree is False:
        passed_phylo_qc = False

    merged_notes = []
    if notes:
        merged_notes.append(notes)
    if annot_note:
        merged_notes.append(f"annot_match_note={annot_note}")
    final_notes = "; ".join(merged_notes)

    main_row = {
        "Species": species,
        "Assembly ID": assembly_id,
        "Data source": source,
        "Assembly level": assembly_level,
        "Assembly size (bp)": assembly_size_bp,
        "GC (%)": gc_percent,
        "Sequence count": sequence_count,
        "Contig count": contig_count,
        "Scaffold count": scaffold_count,
        "N50 (bp)": n50_bp,
        "Gene count": agat_info["gene_count"],
        "mRNA count": agat_info["mrna_count"],
        "CDS count": agat_info["cds_count"],
        "Protein count": protein_count,
        "Longest protein count": longest_protein_count,
        "Complete BUSCOs (%)": busco_info["complete_busco_percent"],
        "Single-copy BUSCOs (%)": busco_info["single_copy_busco_percent"],
        "Duplicated BUSCOs (%)": busco_info["duplicated_busco_percent"],
        "Fragmented BUSCOs (%)": busco_info["fragmented_busco_percent"],
        "Missing BUSCOs (%)": busco_info["missing_busco_percent"],
        "Annotation-genome match (%)": format_percent(float(match_cov) * 100.0 if str(match_cov).strip() != "" else None, 2),
        "Annot match status": match_status,
        "Phylogenomic nongap sites": nongap_sites,
        "Phylogenomic occupancy (%)": occupancy_percent,
        "Included in phylogenomic matrix": parse_bool_to_yes_no(included_in_matrix),
        "Included in final tree": parse_bool_to_yes_no(included_in_tree) if included_in_tree is not None else "",
        "Passed phylo QC": parse_bool_to_yes_no(passed_phylo_qc),
        "Notes": final_notes,
    }

    qc_issues = []

    if str(assembly_size_bp).strip() == "":
        qc_issues.append("missing_quast_or_seqkit_genome_size")
    if str(gc_percent).strip() == "":
        qc_issues.append("missing_quast_or_seqkit_gc")
    if str(n50_bp).strip() == "":
        qc_issues.append("missing_quast_n50")
    if str(agat_info["gene_count"]).strip() == "":
        qc_issues.append("missing_agat_gene_count")
    if str(agat_info["mrna_count"]).strip() == "":
        qc_issues.append("missing_agat_mrna_count")
    if str(agat_info["cds_count"]).strip() == "":
        qc_issues.append("missing_agat_cds_count")
    if str(busco_info["complete_busco_percent"]).strip() == "":
        qc_issues.append("missing_busco_summary")
    if str(protein_count).strip() == "":
        qc_issues.append("missing_seqkit_protein_count")
    if species not in annot_match_map:
        qc_issues.append("species_not_found_in_annot_genome_match")
    if species not in matrix_map:
        qc_issues.append("species_not_found_in_matrix")
    if species_in_tree and species not in species_in_tree:
        qc_issues.append("species_not_found_in_tree")

    qc_row = {
        "Species": species,
        "Has QUAST result": parse_bool_to_yes_no((QUAST_DIR / species / "report.tsv").exists()),
        "Has BUSCO result": parse_bool_to_yes_no(find_busco_summary_file(species) is not None),
        "Has AGAT result": parse_bool_to_yes_no((AGAT_DIR / species / "stats.txt").exists()),
        "Has SeqKit genome result": parse_bool_to_yes_no((SEQKIT_DIR / f"{species}.genome.tsv").exists()),
        "Has SeqKit protein result": parse_bool_to_yes_no((SEQKIT_DIR / f"{species}.protein.tsv").exists()),
        "Found in annot_genome_match": parse_bool_to_yes_no(species in annot_match_map),
        "Found in matrix": parse_bool_to_yes_no(species in matrix_map),
        "Found in species tree": parse_bool_to_yes_no(species in species_in_tree) if species_in_tree else "",
        "QC issues": "; ".join(qc_issues),
    }

    return main_row, qc_row


def main() -> None:
    """主程序。"""
    ensure_dir(OUTPUT_DIR)

    species_list = discover_species_from_genomic_dir()

    annot_match_map: Dict[str, Dict[str, str]] = {}
    matrix_map: Dict[str, Dict[str, str]] = {}
    species_in_tree: set = set()

    if ANNOT_MATCH_TSV.exists():
        annot_match_map = load_annot_match_map(ANNOT_MATCH_TSV)

    matrix_file = pick_first_existing(MATRIX_TSV_CANDIDATES)
    if matrix_file is not None:
        matrix_map = load_matrix_map(matrix_file)

    tree_file = pick_first_existing(TREE_FILE_CANDIDATES)
    if tree_file is not None:
        species_in_tree = set(parse_species_from_newick(tree_file))

    main_rows: List[Dict[str, str]] = []
    qc_rows: List[Dict[str, str]] = []

    for species in species_list:
        main_row, qc_row = collect_species_row(
            species=species,
            annot_match_map=annot_match_map,
            matrix_map=matrix_map,
            species_in_tree=species_in_tree,
        )
        main_rows.append(main_row)
        qc_rows.append(qc_row)

    main_fields = [
        "Species",
        "Assembly ID",
        "Data source",
        "Assembly level",
        "Assembly size (bp)",
        "GC (%)",
        "Sequence count",
        "Contig count",
        "Scaffold count",
        "N50 (bp)",
        "Gene count",
        "mRNA count",
        "CDS count",
        "Protein count",
        "Longest protein count",
        "Complete BUSCOs (%)",
        "Single-copy BUSCOs (%)",
        "Duplicated BUSCOs (%)",
        "Fragmented BUSCOs (%)",
        "Missing BUSCOs (%)",
        "Annotation-genome match (%)",
        "Annot match status",
        "Phylogenomic nongap sites",
        "Phylogenomic occupancy (%)",
        "Included in phylogenomic matrix",
        "Included in final tree",
        "Passed phylo QC",
        "Notes",
    ]

    qc_fields = [
        "Species",
        "Has QUAST result",
        "Has BUSCO result",
        "Has AGAT result",
        "Has SeqKit genome result",
        "Has SeqKit protein result",
        "Found in annot_genome_match",
        "Found in matrix",
        "Found in species tree",
        "QC issues",
    ]

    main_tsv = OUTPUT_DIR / MAIN_TABLE_TSV
    main_xlsx = OUTPUT_DIR / MAIN_TABLE_XLSX
    qc_tsv = OUTPUT_DIR / QC_TABLE_TSV
    qc_xlsx = OUTPUT_DIR / QC_TABLE_XLSX

    write_tsv(main_rows, main_tsv, main_fields)
    write_xlsx(main_rows, main_xlsx, main_fields)
    write_tsv(qc_rows, qc_tsv, qc_fields)
    write_xlsx(qc_rows, qc_xlsx, qc_fields)

    print("=" * 70)
    print("[DONE] 主表已生成：")
    print(main_tsv)
    print(main_xlsx)
    print()
    print("[DONE] QC 表已生成：")
    print(qc_tsv)
    print(qc_xlsx)
    print("=" * 70)


if __name__ == "__main__":
    main()
