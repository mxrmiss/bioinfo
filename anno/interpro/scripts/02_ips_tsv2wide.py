#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_ips_tsv2wide.py
================================================================================
功能：
  - 读取 config.yaml
  - 自动定位 01 生成的 TSV（results/01_interproscan/{ver}.core/{input}.tsv）
  - 输出宽表到 results/02_wide/{ver}.core/
  - 可选同时输出 TSV 和 XLSX
  - 可选把 GO 抽成一列（wide_table.include_go=true）

GO 抽取策略（稳，不依赖固定列号）：
  - 在每行的所有字段中扫描形如 GO:0000000 的 token
  - 去重、排序后聚合到 protein_id 级别
================================================================================
"""

import re
import sys
import csv
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple
from collections import defaultdict

CONFIG_YAML = "config.yaml"

PROGRESS_EVERY_LINES = 2_000_000
PROGRESS_EVERY_PROTEINS = 50_000

GO_RE = re.compile(r"GO:\d{7}")


def die(msg: str, code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def load_yaml_config(path: str) -> Dict[str, Any]:
    p = Path(path)
    if not p.is_file():
        die(f"找不到配置文件：{path}")

    try:
        import yaml  # type: ignore
    except Exception as e:
        die(
            "缺少依赖 PyYAML，无法读取 config.yaml。\n"
            "解决方式：在你的 conda 环境中安装 PyYAML。\n"
            f"原始错误：{e}"
        )

    with p.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        die("config.yaml 解析失败：顶层必须是字典")
    return cfg


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def setup_logger(log_file: Path, level: str = "INFO"):
    import logging

    level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    lv = level_map.get(level.upper(), logging.INFO)

    logger = logging.getLogger(str(log_file))
    logger.setLevel(lv)
    logger.propagate = False

    if logger.handlers:
        return logger

    fmt = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(lv)
    sh.setFormatter(fmt)

    ensure_dir(log_file.parent)
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(lv)
    fh.setFormatter(fmt)

    logger.addHandler(sh)
    logger.addHandler(fh)
    return logger


def fmt_tpl(s: str, mapping: Dict[str, Any]) -> str:
    try:
        return s.format_map(mapping)
    except Exception:
        return s


def resolve_inputs(cfg: Dict[str, Any]) -> List[Path]:
    sec = cfg.get("inputs", {}) or {}
    inputs = sec.get("inputs", None)

    out: List[Path] = []
    if isinstance(inputs, list) and inputs:
        for x in inputs:
            out.append(Path(str(x)))
        return out

    inputs_dir = Path(str(sec.get("inputs_dir", "data/query")))
    ext = str(sec.get("inputs_ext", ".faa"))

    if not inputs_dir.is_dir():
        die(f"inputs_dir 不存在或不是目录：{inputs_dir}")

    for p in sorted(inputs_dir.iterdir()):
        if p.is_file() and p.name.endswith(ext):
            out.append(p)

    if not out:
        die(f"在 {inputs_dir} 下没有找到 *{ext} 输入文件")
    return out


def basename_keep_ext(p: Path) -> str:
    return p.name


def norm_desc(s: str, sanitize: bool = True) -> str:
    s = (s or "").strip()
    if sanitize:
        s = s.replace("\t", " ").replace(";", ",")
    s = " ".join(s.split())
    return s


def locate_tsv(outdir: Path, input_faa: Path) -> Path:
    base = basename_keep_ext(input_faa)
    return outdir / f"{base}.tsv"


def extract_go_terms(row: List[str]) -> List[str]:
    out: List[str] = []
    for field in row:
        if "GO:" not in field:
            continue
        out.extend(GO_RE.findall(field))
    return out


def parse_ips_tsv(
    in_tsv: Path,
    target_analyses: List[str],
    sanitize_desc_flag: bool,
    include_go: bool,
    logger,
) -> Tuple[List[str], Dict[str, Dict[str, Dict[str, str]]], Dict[str, set]]:
    target_set = set(target_analyses)
    hits = defaultdict(lambda: defaultdict(dict))
    proteins_seen = set()

    go_by_protein: Dict[str, set] = defaultdict(set)

    bad_short = 0
    total = 0

    logger.info(f"[02] Reading TSV: {in_tsv}")

    with in_tsv.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            total += 1
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 6:
                bad_short += 1
                continue

            protein_id = (row[0] or "").strip()
            analysis = (row[3] or "").strip()
            accession = (row[4] or "").strip()
            desc = norm_desc(row[5], sanitize=sanitize_desc_flag)

            if not protein_id:
                continue

            proteins_seen.add(protein_id)

            if include_go:
                for go in extract_go_terms(row):
                    go_by_protein[protein_id].add(go)

            if analysis not in target_set:
                if total % PROGRESS_EVERY_LINES == 0:
                    logger.info(f"[02] TSV lines processed: {total:,} ; proteins_seen={len(proteins_seen):,}")
                continue
            if not accession:
                if total % PROGRESS_EVERY_LINES == 0:
                    logger.info(f"[02] TSV lines processed: {total:,} ; proteins_seen={len(proteins_seen):,}")
                continue

            old = hits[protein_id][analysis].get(accession, "")
            if len(desc) > len(old):
                hits[protein_id][analysis][accession] = desc

            if total % PROGRESS_EVERY_LINES == 0:
                logger.info(f"[02] TSV lines processed: {total:,} ; proteins_seen={len(proteins_seen):,}")

    logger.info(f"[02] TSV lines total: {total:,}")
    if bad_short:
        logger.warning(f"[02] Rows with <6 columns skipped: {bad_short:,}")

    proteins_sorted = sorted(proteins_seen)
    logger.info(f"[02] Proteins total: {len(proteins_sorted):,}")

    return proteins_sorted, hits, go_by_protein


def write_wide_tsv(
    out_tsv: Path,
    proteins: List[str],
    hits: Dict[str, Dict[str, Dict[str, str]]],
    go_by_protein: Dict[str, set],
    target_analyses: List[str],
    join_sep: str,
    include_go: bool,
    go_col_name: str,
    go_join_sep: str,
    logger,
) -> None:
    ensure_dir(out_tsv.parent)
    logger.info(f"[02] Writing wide TSV: {out_tsv}")

    header = ["protein_id"] + target_analyses
    if include_go:
        header.append(go_col_name)

    with out_tsv.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)

        for idx, pid in enumerate(proteins, start=1):
            row_out = [pid]
            for ana in target_analyses:
                acc2desc = hits.get(pid, {}).get(ana, {})
                if not acc2desc:
                    row_out.append("")
                    continue

                items = []
                for acc in sorted(acc2desc.keys()):
                    d = acc2desc[acc]
                    items.append(f"{acc}|{d}" if d else acc)
                row_out.append(join_sep.join(items))

            if include_go:
                gos = sorted(go_by_protein.get(pid, set()))
                row_out.append(go_join_sep.join(gos) if gos else "")

            w.writerow(row_out)

            if idx % PROGRESS_EVERY_PROTEINS == 0:
                logger.info(f"[02] TSV written proteins: {idx:,}/{len(proteins):,}")

    logger.info("[02] TSV done")


def write_wide_xlsx(
    out_xlsx: Path,
    proteins: List[str],
    hits: Dict[str, Dict[str, Dict[str, str]]],
    go_by_protein: Dict[str, set],
    target_analyses: List[str],
    join_sep: str,
    include_go: bool,
    go_col_name: str,
    go_join_sep: str,
    logger,
) -> None:
    ensure_dir(out_xlsx.parent)
    logger.info(f"[02] Writing wide XLSX: {out_xlsx}")

    try:
        from openpyxl import Workbook
    except Exception as e:
        logger.error(f"[02] openpyxl 导入失败，无法写 XLSX：{e}")
        return

    wb = Workbook()
    ws = wb.active
    ws.title = "wide"
    ws.freeze_panes = "B2"

    header = ["protein_id"] + target_analyses
    if include_go:
        header.append(go_col_name)
    ws.append(header)

    for idx, pid in enumerate(proteins, start=1):
        row_out = [pid]
        for ana in target_analyses:
            acc2desc = hits.get(pid, {}).get(ana, {})
            if not acc2desc:
                row_out.append("")
                continue

            items = []
            for acc in sorted(acc2desc.keys()):
                d = acc2desc[acc]
                items.append(f"{acc}|{d}" if d else acc)
            row_out.append(join_sep.join(items))

        if include_go:
            gos = sorted(go_by_protein.get(pid, set()))
            row_out.append(go_join_sep.join(gos) if gos else "")

        ws.append(row_out)

        if idx % PROGRESS_EVERY_PROTEINS == 0:
            logger.info(f"[02] XLSX written proteins: {idx:,}/{len(proteins):,}")

    wb.save(str(out_xlsx))
    logger.info("[02] XLSX done")


def main() -> None:
    cfg = load_yaml_config(CONFIG_YAML)

    log_cfg = cfg.get("logging", {}) or {}
    log_dir = Path(str(log_cfg.get("log_dir", "logs")))
    log_level = str(log_cfg.get("level", "INFO"))

    ips_cfg = cfg.get("interproscan", {}) or {}
    ver = str(ips_cfg.get("ver", "5.76-107.0")).strip()
    outdir_tpl = str(ips_cfg.get("outdir_tpl", "results/01_interproscan/{ver}.core"))
    outdir = Path(fmt_tpl(outdir_tpl, {"ver": ver}))

    wt_cfg = cfg.get("wide_table", {}) or {}
    wt_outdir_tpl = str(wt_cfg.get("outdir_tpl", "results/02_wide/{ver}.core"))
    wt_outdir = Path(fmt_tpl(wt_outdir_tpl, {"ver": ver}))

    write_tsv_flag = bool(wt_cfg.get("write_tsv", True))
    write_xlsx_flag = bool(wt_cfg.get("write_xlsx", False))

    include_go = bool(wt_cfg.get("include_go", False))
    go_col_name = str(wt_cfg.get("go_column_name", "GO"))
    go_join_sep = str(wt_cfg.get("go_join_sep", "; "))

    join_sep = str(wt_cfg.get("join_sep", "; "))
    sanitize_desc_flag = bool(wt_cfg.get("sanitize_desc", True))

    target_analyses = wt_cfg.get("target_analyses", None)
    if not isinstance(target_analyses, list) or not target_analyses:
        die("config.yaml: wide_table.target_analyses 必须是非空列表")

    xlsx_skip_if_too_many = bool(wt_cfg.get("xlsx_skip_if_too_many_rows", True))
    EXCEL_MAX_ROWS = 1_048_576

    main_log = log_dir / "02_ips_tsv2wide" / ver / "main.log"
    logger = setup_logger(main_log, level=log_level)

    logger.info(f"[02] Start at {now_iso()}")
    logger.info(f"[02] VER={ver}")
    logger.info(f"[02] OUTDIR(ips)={outdir}")
    logger.info(f"[02] OUTDIR(wide)={wt_outdir}")
    logger.info(f"[02] write_tsv={write_tsv_flag}")
    logger.info(f"[02] write_xlsx={write_xlsx_flag}")
    logger.info(f"[02] include_go={include_go}")
    logger.info(f"[02] target_analyses={','.join([str(x) for x in target_analyses])}")

    if not outdir.is_dir():
        die(f"InterProScan outdir 不存在：{outdir}（请先运行 01_run_interproscan.py）")

    inputs = resolve_inputs(cfg)
    logger.info(f"[02] N_inputs={len(inputs)}")
    ensure_dir(wt_outdir)

    for i, in_faa in enumerate(inputs, start=1):
        base = basename_keep_ext(in_faa)
        in_tsv = locate_tsv(outdir, in_faa)

        per_log = log_dir / "02_ips_tsv2wide" / ver / f"{base}.log"
        per_logger = setup_logger(per_log, level=log_level)

        per_logger.info(f"[02] ({i}/{len(inputs)}) INPUT={in_faa}")
        per_logger.info(f"[02] TSV={in_tsv}")

        if not in_tsv.is_file() or in_tsv.stat().st_size == 0:
            per_logger.error("[02] 找不到 TSV 或 TSV 为空，请先确认 01 是否成功生成输出。")
            continue

        proteins, hits, go_by_protein = parse_ips_tsv(
            in_tsv=in_tsv,
            target_analyses=target_analyses,
            sanitize_desc_flag=sanitize_desc_flag,
            include_go=include_go,
            logger=per_logger,
        )

        out_prefix = wt_outdir / f"{base}.wide"
        out_tsv = Path(str(out_prefix) + ".tsv")
        out_xlsx = Path(str(out_prefix) + ".xlsx")

        if write_tsv_flag:
            write_wide_tsv(
                out_tsv=out_tsv,
                proteins=proteins,
                hits=hits,
                go_by_protein=go_by_protein,
                target_analyses=target_analyses,
                join_sep=join_sep,
                include_go=include_go,
                go_col_name=go_col_name,
                go_join_sep=go_join_sep,
                logger=per_logger,
            )

        if write_xlsx_flag:
            if xlsx_skip_if_too_many and (len(proteins) + 1) > EXCEL_MAX_ROWS:
                per_logger.warning(
                    f"[02] proteins={len(proteins):,} 超过 Excel 行上限，按配置跳过写 XLSX。"
                )
            else:
                write_wide_xlsx(
                    out_xlsx=out_xlsx,
                    proteins=proteins,
                    hits=hits,
                    go_by_protein=go_by_protein,
                    target_analyses=target_analyses,
                    join_sep=join_sep,
                    include_go=include_go,
                    go_col_name=go_col_name,
                    go_join_sep=go_join_sep,
                    logger=per_logger,
                )

        per_logger.info("[02] Done one input")

    logger.info(f"[02] Done at {now_iso()}")


if __name__ == "__main__":
    main()

