#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny02_blocks_macro.py
—— 朴素模式 Step02：只生成相邻 pairs 的 anchors.simple（MINSPAN 唯一旋钮）

输入（硬位置）：
- raw_data/synteny_species_meta.tsv
- output/synteny_01_mcscan_catalog/<species>.geneorder.bed
- raw_data/proteins/<species>.(faa/fa/fasta/pep[.gz]) 自动探测

输出（硬位置）：
output/synteny_02_blocks_macro/
  raw_anchors/
    <A>__vs__<B>.anchors
    <A>__vs__<B>.anchors.new
    <A>__vs__<B>.anchors.simple
  work_pairs/<pair>/   jcvi 工作目录（便于排查）
  summaries/step02.pairs.summary.tsv

要点：
- 相邻链式：只跑 (s1,s2),(s2,s3)...
- 线程：只给 ortholog 传 --cpus
- screen：只传 --minspan + --simple
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

THREADS = 25
MINSPAN = 5
JCVI_PYTHON = "python"

META_TSV_REL = "raw_data/synteny_species_meta.tsv"
STEP01_DIR_REL = "output/synteny_01_mcscan_catalog"
PROTEIN_DIR_REL = "raw_data/proteins"
OUTPUT_DIR_REL = "output/synteny_02_blocks_macro"
LOG_DIR_REL = "logs"

PROTEIN_SUFFIX_CANDIDATES = [
    ".faa", ".fa", ".fasta", ".pep",
    ".faa.gz", ".fa.gz", ".fasta.gz", ".pep.gz",
]

# ============================================================
# 实现区（皇上勿改）
# ============================================================

import os
import csv
import gzip
import time
import shutil
import traceback
import subprocess
from pathlib import Path
from typing import List, Tuple, Iterable


def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


class Logger:
    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

    def _write(self, level: str, msg: str) -> None:
        line = f"{now_ts()} [{level}] {msg}"
        print(line, flush=True)
        with self.log_file.open("a", encoding="utf-8") as fw:
            fw.write(line + "\n")

    def info(self, msg: str) -> None:
        self._write("INFO", msg)

    def warn(self, msg: str) -> None:
        self._write("WARN", msg)

    def error(self, msg: str) -> None:
        self._write("ERROR", msg)


def die(logger: Logger, msg: str, code: int = 1) -> None:
    logger.error(msg)
    raise SystemExit(code)


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def clean_dir(path: Path, clean: bool, logger: Logger) -> None:
    if clean and path.exists():
        logger.info(f"Clean output dir: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fw:
        w = csv.writer(fw, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def read_meta(meta_tsv: Path) -> List[str]:
    species: List[str] = []
    with open_text_auto(meta_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != ["species_id", "group"]:
            raise ValueError(f"meta header must be: species_id<TAB>group, got: {header}")
        for ln, line in enumerate(fr, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"meta line {ln} must have 2 columns: {line}")
            sid = parts[0].strip()
            if not sid:
                raise ValueError(f"meta line {ln} species_id empty")
            species.append(sid)
    if len(species) < 2:
        raise ValueError("meta must contain at least 2 species")
    return species


def detect_protein_file(protein_dir: Path, species_id: str) -> Path:
    for suf in PROTEIN_SUFFIX_CANDIDATES:
        p = protein_dir / f"{species_id}{suf}"
        if p.exists():
            return p
    raise FileNotFoundError(f"Protein FASTA not found for {species_id} in {protein_dir}")


def run_cmd(logger: Logger, cmd: List[str], cwd: Path) -> None:
    logger.info("[CMD] " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = (r.stdout or "").rstrip("\n")
    if out:
        for ln in out.splitlines():
            logger.info(f"[JCVI] {ln}")
    if r.returncode != 0:
        raise RuntimeError(f"Command failed (exit={r.returncode})")


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        project = Path(__file__).resolve().parents[1]
    else:
        project = Path(str(pr)).expanduser().resolve()

    meta_tsv = project / META_TSV_REL
    step01_dir = project / STEP01_DIR_REL
    protein_dir = project / PROTEIN_DIR_REL
    out_dir = project / OUTPUT_DIR_REL
    log_dir = project / LOG_DIR_REL

    logger = Logger(log_dir / "synteny02_blocks_macro.log")
    t0 = time.time()

    logger.info("========== synteny02 — plain anchors.simple (adjacent, MINSPAN only) ==========")
    logger.info(f"PROJECT_ROOT={project}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"THREADS={THREADS}")
    logger.info(f"MINSPAN={MINSPAN}")
    logger.info(f"JCVI_PYTHON={JCVI_PYTHON}")

    if not meta_tsv.exists():
        die(logger, f"Missing meta: {meta_tsv}")
    if not step01_dir.exists():
        die(logger, f"Missing step01 dir: {step01_dir}")
    if not protein_dir.exists():
        die(logger, f"Missing protein dir: {protein_dir}")

    species = read_meta(meta_tsv)
    pairs = [(species[i], species[i + 1]) for i in range(len(species) - 1)]
    logger.info("META species order=" + " | ".join(species))
    logger.info(f"N_species={len(species)} N_pairs={len(pairs)}")

    clean_dir(out_dir, CLEAN_OUTPUT, logger)

    d_raw = out_dir / "raw_anchors"
    d_work = out_dir / "work_pairs"
    d_inputs = out_dir / "_jcvi_inputs"
    d_sum = out_dir / "summaries"
    for d in (d_raw, d_work, d_inputs, d_sum):
        d.mkdir(parents=True, exist_ok=True)

    for sid in species:
        bed_src = step01_dir / f"{sid}.geneorder.bed"
        if not bed_src.exists():
            die(logger, f"Missing geneorder bed: {bed_src}")
        pep_src = detect_protein_file(protein_dir, sid)

        bed_dst = d_inputs / f"{sid}.bed"
        pep_dst = d_inputs / f"{sid}.pep"
        for src, dst in ((bed_src, bed_dst), (pep_src, pep_dst)):
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            dst.symlink_to(src)

    sum_rows: List[List[str]] = []
    sum_header = ["pair_id", "a_species", "b_species", "anchors_simple_lines", "runtime_sec"]

    for a, b in pairs:
        pair_id = f"{a}__vs__{b}"
        wdir = d_work / pair_id
        wdir.mkdir(parents=True, exist_ok=True)

        for x in wdir.glob("*"):
            try:
                if x.is_dir():
                    shutil.rmtree(x)
                else:
                    x.unlink()
            except Exception:
                pass

        for sid in (a, b):
            for ext in (".bed", ".pep"):
                src = d_inputs / f"{sid}{ext}"
                dst = wdir / f"{sid}{ext}"
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                dst.symlink_to(src)

        t_pair = time.time()

        cmd1 = [
            JCVI_PYTHON, "-m", "jcvi.compara.catalog", "ortholog",
            "--dbtype", "prot",
            "--no_strip_names",
            a, b,
            f"--cpus={THREADS}",
        ]
        run_cmd(logger, cmd1, cwd=wdir)

        cmd2 = [
            JCVI_PYTHON, "-m", "jcvi.compara.synteny", "screen",
            f"--minspan={MINSPAN}",
            "--simple",
            f"{a}.{b}.anchors",
            f"{a}.{b}.anchors.new",
        ]
        run_cmd(logger, cmd2, cwd=wdir)

        anchors = wdir / f"{a}.{b}.anchors"
        anchors_new = wdir / f"{a}.{b}.anchors.new"
        anchors_simple = wdir / f"{a}.{b}.anchors.simple"
        if not anchors.exists() or not anchors_new.exists() or not anchors_simple.exists():
            die(logger, f"[{pair_id}] Missing outputs in work dir: {wdir}")

        shutil.copy2(anchors, d_raw / f"{pair_id}.anchors")
        shutil.copy2(anchors_new, d_raw / f"{pair_id}.anchors.new")
        shutil.copy2(anchors_simple, d_raw / f"{pair_id}.anchors.simple")

        n_simple = 0
        with open_text_auto(d_raw / f"{pair_id}.anchors.simple") as fr:
            for line in fr:
                if line and line.strip():
                    n_simple += 1

        rt = int(time.time() - t_pair)
        sum_rows.append([pair_id, a, b, str(n_simple), str(rt)])
        logger.info(f"[{pair_id}] anchors.simple lines={n_simple} runtime_sec={rt}")

    write_tsv(d_sum / "step02.pairs.summary.tsv", sum_header, sum_rows)
    logger.info(f"Done. runtime_sec={int(time.time()-t0)}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        pr = PROJECT_ROOT
        if pr is None or str(pr).strip() == "":
            project = Path(__file__).resolve().parents[1]
        else:
            project = Path(str(pr)).expanduser().resolve()
        lg = Logger(project / LOG_DIR_REL / "synteny02_blocks_macro.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

