#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny05_plot_karyotype.py
—— Step05：adjacent only + 不筛线 + dominant-ref 上色（以 S. constricta 为唯一参考）+ chr 级 flip + 只出 PDF

输入（硬位置）：
- raw_data/palette.tsv
- output/synteny_04_layout_tracks/layout_species_tracks.tsv
- output/synteny_03_order_seqids/seqids_species/<species>.seqids
- output/synteny_03_order_seqids/chr2ref/chr2ref.tsv
- output/synteny_03_order_seqids/chr2ref/chr_flip.tsv
- output/synteny_01_mcscan_catalog/<species>.geneorder.bed
- output/synteny_02_blocks_macro/raw_anchors/<A>__vs__<B>.anchors.simple
- output/synteny_00_chr_rename/chr_rename_<species>.tsv   （用于 flip 的 chr length）

输出（硬位置）：
output/synteny_05_plot/
  inputs/<species>.bed
  simple/<A>__vs__<B>.simple
  layout/seqids
  layout/karyotype.layout
  figures/karyotype.pdf
  summaries/step05.links.summary.tsv

要点：
- 连线：不做任何过滤、不截断
- 上色：完全按 ref_chr（来自 Step03 chr2ref.tsv）决定颜色
- flip：完全按 Step03 chr_flip.tsv 决定是否对 bed 坐标反转（chr 级）
- 防 KeyError：bed 内为每个 gene 复制多份 #HEX*gene（palette 全部颜色）
- 物种名大小/画布高度：放在用户参数区（尽量走 jcvi 原生参数；若不支持则不强行加，保证稳定）
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

JCVI_PYTHON = "python"
FONT_FAMILY = "Arial"
OUT_PDF = "karyotype.pdf"

FIGURE_WIDTH_INCH = 16.0
FIGURE_HEIGHT_INCH = None

# 物种名字号（若当前 jcvi 版本支持 --labelsize，则自动生效；不支持则忽略以保证稳定）
SPECIES_LABEL_SIZE_PT = 12

TRACK_XSTART = 0.10
TRACK_XEND = 0.92
LABEL_VA = "top"

NO_CIRCLES = True

PALETTE_TSV_REL = "raw_data/palette.tsv"
STEP00_DIR_REL = "output/synteny_00_chr_rename"
STEP01_DIR_REL = "output/synteny_01_mcscan_catalog"
STEP02_DIR_REL = "output/synteny_02_blocks_macro"
STEP03_DIR_REL = "output/synteny_03_order_seqids"
STEP04_DIR_REL = "output/synteny_04_layout_tracks"
OUTPUT_DIR_REL = "output/synteny_05_plot"
LOG_DIR_REL = "logs"

# ============================================================
# 实现区（皇上勿改）
# ============================================================

import csv
import gzip
import time
import shutil
import traceback
import subprocess
import re
from pathlib import Path
from typing import List, Dict, Tuple, Iterable, Optional


_ALLOWED_JCVI_FONTS = {"Helvetica", "Liberation Sans", "Palatino", "Schoolbook", "Arial"}
_CHR_NUM_RE = re.compile(r"^Chr0*([0-9]+)$")


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


def load_palette(palette_tsv: Path) -> List[str]:
    if not palette_tsv.exists():
        raise FileNotFoundError(f"Missing palette: {palette_tsv}")
    rows: List[Tuple[int, str]] = []
    with open_text_auto(palette_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != ["rank", "color_hex"]:
            raise ValueError(f"palette header must be: rank<TAB>color_hex, got: {header}")
        for ln, line in enumerate(fr, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            rk, hx = line.split("\t")
            rk = int(rk) if rk.isdigit() else 0
            hx = hx.strip()
            if rk > 0 and hx:
                rows.append((rk, hx))
    rows.sort(key=lambda x: x[0])
    pal = [hx for _, hx in rows]
    if len(pal) < 2:
        raise ValueError("palette must contain >=2 colors")
    return pal


def read_layout_tracks(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing layout tracks: {path}")
    need = ["order_index", "species_id", "display_label", "group", "y_center", "track_height", "chr_bar_color_hex", "label_color_hex"]
    rows: List[Dict[str, str]] = []
    with open_text_auto(path) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != need:
            raise ValueError(f"Bad header in {path.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(need):
                continue
            rows.append(dict(zip(need, parts)))
    rows.sort(key=lambda d: int(d["order_index"]))
    return rows


def validate_seqids_line(line: str) -> List[str]:
    raw = (line or "").strip()
    if not raw:
        raise ValueError("seqids line is empty")
    tokens = [t.strip() for t in raw.split(",") if t.strip()]
    if not tokens:
        raise ValueError("seqids has no tokens")
    return tokens


def build_seqids_file(track_species: List[str], step03_seqids_dir: Path, out_seqids: Path, logger: Logger) -> None:
    lines: List[str] = []
    for sid in track_species:
        p = step03_seqids_dir / f"{sid}.seqids"
        if not p.exists():
            raise FileNotFoundError(f"Missing step03 seqids for {sid}: {p}")
        line = p.read_text(encoding="utf-8", errors="replace").strip()
        tokens = validate_seqids_line(line)
        logger.info(f"[SEQIDS] {sid}: n_chr={len(tokens)}")
        lines.append(line)
    out_seqids.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_layout_file(tracks: List[Dict[str, str]], out_layout: Path, bed_dir: Path, simple_dir: Path) -> None:
    lines: List[str] = []
    rot = 0
    for d in tracks:
        sid = d["species_id"]
        y = float(d["y_center"])
        bar_hex = d["chr_bar_color_hex"]
        label_text = d["display_label"]
        bed = bed_dir / f"{sid}.bed"
        lines.append(f"{y:.6f}, {TRACK_XSTART:.2f}, {TRACK_XEND:.2f}, {rot}, {bar_hex}, {label_text}, {LABEL_VA}, {bed.as_posix()}, center")
    for i in range(len(tracks) - 1):
        a = tracks[i]["species_id"]
        b = tracks[i + 1]["species_id"]
        simple = simple_dir / f"{a}__vs__{b}.simple"
        lines.append(f"e, {i}, {i+1}, {simple.as_posix()}")
    out_layout.write_text("\n".join(lines) + "\n", encoding="utf-8")


def normalize_gene_id(raw: str) -> str:
    x = (raw or "").strip().strip('"').strip("'").strip()
    if not x:
        return ""
    if "*" in x and x.startswith("#") and len(x) >= 8:
        try:
            x = x.split("*", 1)[1]
        except Exception:
            pass
    if "|" in x:
        x = x.split("|")[-1]
    if x.startswith("rna-"):
        x = x[4:]
    return x.strip()


def load_gene_to_chr(step01_bed: Path) -> Dict[str, str]:
    mp: Dict[str, str] = {}
    with step01_bed.open("r", encoding="utf-8", errors="replace") as fr:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chr_, _, _, gid = parts[:4]
            gid = normalize_gene_id(gid)
            if gid:
                mp[gid] = chr_
    if not mp:
        raise ValueError(f"Empty gene->chr map from: {step01_bed}")
    return mp


def chr_to_num(chr_name: str) -> int:
    m = _CHR_NUM_RE.match(chr_name)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return 1
    return 1


def load_chr2ref(chr2ref_tsv: Path) -> Dict[str, Dict[str, str]]:
    if not chr2ref_tsv.exists():
        raise FileNotFoundError(f"Missing chr2ref: {chr2ref_tsv}")
    need = ["species_id", "chr", "ref_chr", "support"]
    out: Dict[str, Dict[str, str]] = {}
    with open_text_auto(chr2ref_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != need:
            raise ValueError(f"Bad header in chr2ref.tsv: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                continue
            sid, chr_, ref_chr, _sup = parts
            out.setdefault(sid, {})[chr_] = ref_chr
    return out


def load_chr_flip(chr_flip_tsv: Path) -> Dict[str, Dict[str, int]]:
    if not chr_flip_tsv.exists():
        raise FileNotFoundError(f"Missing chr_flip: {chr_flip_tsv}")
    need = ["species_id", "chr", "flip", "support"]
    out: Dict[str, Dict[str, int]] = {}
    with open_text_auto(chr_flip_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != need:
            raise ValueError(f"Bad header in chr_flip.tsv: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                continue
            sid, chr_, flip, _sup = parts
            try:
                out.setdefault(sid, {})[chr_] = 1 if int(flip) == 1 else 0
            except Exception:
                out.setdefault(sid, {})[chr_] = 0
    return out


def load_chr_lengths(step00_chr_rename_tsv: Path) -> Dict[str, int]:
    """
    从 Step00 的 chr_rename_<species>.tsv 读取 ChrNN -> length_bp（is_chromosome=yes）
    """
    if not step00_chr_rename_tsv.exists():
        raise FileNotFoundError(f"Missing step00 chr_rename: {step00_chr_rename_tsv}")
    need = ["species_id", "seqid_raw", "seqid_renamed", "rank", "length_bp", "is_chromosome"]
    out: Dict[str, int] = {}
    with open_text_auto(step00_chr_rename_tsv) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header != need:
            raise ValueError(f"Bad header in {step00_chr_rename_tsv.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                continue
            _sid, _raw, renamed, _rk, L, is_chr = parts
            if is_chr != "yes":
                continue
            if renamed == "NA":
                continue
            try:
                out[renamed] = int(L)
            except Exception:
                continue
    if not out:
        raise ValueError(f"Empty chr lengths from: {step00_chr_rename_tsv}")
    return out


def copy_bed_with_palette_dups_and_flip(
    step01_bed: Path,
    out_bed: Path,
    palette: List[str],
    flip_map_chr: Dict[str, int],
    chr_len_map: Dict[str, int],
) -> int:
    """
    写出 bed，并对 flip=1 的 chr 做坐标反转：
      new_start = L - old_end
      new_end   = L - old_start
    同时保留“为每个 gene 复制多份 #HEX*gene”以防 KeyError（保持原行为稳定）。
    最后按 (chr, start, gene) 排序，保证 jcvi 读入稳定。
    """
    lines: List[Tuple[str, int, int, str]] = []
    n_gene = 0

    with step01_bed.open("r", encoding="utf-8", errors="replace") as fr:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chr_, start_s, end_s, gid_raw = parts[:4]
            gid = normalize_gene_id(gid_raw)
            if not gid:
                continue
            try:
                s = int(float(start_s))
                e = int(float(end_s))
            except Exception:
                continue
            if e <= s:
                e = s + 1

            if flip_map_chr.get(chr_, 0) == 1:
                L = chr_len_map.get(chr_, None)
                if L is not None and L > 0:
                    ns = L - e
                    ne = L - s
                    if ne <= ns:
                        ne = ns + 1
                    if ns < 0:
                        ns = 0
                    if ne < 1:
                        ne = 1
                    s, e = ns, ne

            lines.append((chr_, s, e, gid))
            for hx in palette:
                lines.append((chr_, s, e, f"{hx}*{gid}"))
            n_gene += 1

    lines.sort(key=lambda x: (chr_to_num(x[0]), x[0], x[1], x[3]))

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with out_bed.open("w", encoding="utf-8", newline="") as fw:
        for chr_, s, e, gid in lines:
            fw.write("\t".join([chr_, str(s), str(e), gid]) + "\n")

    return n_gene


def ref_chr_to_color(ref_chr: str, palette: List[str]) -> str:
    n = chr_to_num(ref_chr)
    return palette[(n - 1) % len(palette)]


def gene_to_ref_color(
    species_id: str,
    gene_id: str,
    gene2chr: Dict[str, Dict[str, str]],
    chr2ref: Dict[str, Dict[str, str]],
    palette: List[str],
) -> str:
    gid = normalize_gene_id(gene_id)
    chr_ = gene2chr.get(species_id, {}).get(gid, "")
    ref_chr = chr2ref.get(species_id, {}).get(chr_, "NA")
    if not ref_chr or ref_chr == "NA":
        return palette[0]
    return ref_chr_to_color(ref_chr, palette)


def build_simple_dominant_ref_colored(
    a_species: str,
    src_simple: Path,
    out_simple: Path,
    gene2chr: Dict[str, Dict[str, str]],
    chr2ref: Dict[str, Dict[str, str]],
    palette: List[str],
) -> Tuple[int, int]:
    kept = 0
    dropped = 0
    out_simple.parent.mkdir(parents=True, exist_ok=True)

    with open_text_auto(src_simple) as fr, out_simple.open("w", encoding="utf-8", newline="") as fw:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                dropped += 1
                continue
            a1, a2, b1, b2, n, o = parts

            # 颜色只按 A 侧基因映射到 ref_chr 决定（dominant ref）
            hx = gene_to_ref_color(a_species, a1, gene2chr, chr2ref, palette)

            # 为保证 bed 中能找到对应 gid，这里统一用 normalize 后的 gene（与 Step01 一致）
            a1n = normalize_gene_id(a1)
            a2n = normalize_gene_id(a2)
            b1n = normalize_gene_id(b1)
            b2n = normalize_gene_id(b2)

            if not a1n:
                dropped += 1
                continue

            fw.write("\t".join([f"{hx}*{a1n}", f"{hx}*{a2n}", f"{hx}*{b1n}", f"{hx}*{b2n}", n, o]) + "\n")
            kept += 1

    return kept, dropped


def _estimate_height_inch(n_tracks: int) -> float:
    h = max(6.0, 0.66 * float(n_tracks))
    return round(h, 2)


def _karyotype_help_text(logger: Logger) -> str:
    cmd = [JCVI_PYTHON, "-m", "jcvi.graphics.karyotype", "-h"]
    try:
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return (r.stdout or "")
    except Exception as e:
        logger.warn(f"Cannot probe jcvi karyotype help: {e}")
        return ""


def run_jcvi_karyotype(seqids: Path, layout: Path, fig_dir: Path, n_tracks: int, logger: Logger) -> Path:
    fig_dir.mkdir(parents=True, exist_ok=True)
    if FONT_FAMILY not in _ALLOWED_JCVI_FONTS:
        die(logger, f"FONT_FAMILY illegal: {FONT_FAMILY}")

    w = float(FIGURE_WIDTH_INCH)
    h = float(FIGURE_HEIGHT_INCH) if FIGURE_HEIGHT_INCH is not None else _estimate_height_inch(n_tracks)

    cmd = [
        JCVI_PYTHON, "-m", "jcvi.graphics.karyotype",
        str(seqids), str(layout),
        "--format=pdf",
        "--notex",
        f"--font={FONT_FAMILY}",
        f"--figsize={int(round(w))}x{int(round(h))}",
    ]
    if NO_CIRCLES:
        cmd.append("--nocircles")

    # 尽量使用 jcvi 原生参数控制物种名字号（若当前版本支持）
    help_txt = _karyotype_help_text(logger)
    if help_txt:
        if "--labelsize" in help_txt and SPECIES_LABEL_SIZE_PT is not None:
            cmd.append(f"--labelsize={int(SPECIES_LABEL_SIZE_PT)}")

    logger.info("[CMD] " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(fig_dir), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = (r.stdout or "").rstrip("\n")
    if out:
        for ln in out.splitlines():
            if "pkg_resources is deprecated" in ln:
                continue
            logger.info(f"[JCVI] {ln}")
    if r.returncode != 0:
        raise RuntimeError(f"jcvi karyotype failed (exit={r.returncode})")

    pdf = None
    for p in sorted(fig_dir.glob("*.pdf")):
        pdf = p
        break
    if pdf is None or not pdf.exists():
        raise RuntimeError("jcvi did not produce any PDF")
    return pdf


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        project = Path(__file__).resolve().parents[1]
    else:
        project = Path(str(pr)).expanduser().resolve()

    palette_tsv = project / PALETTE_TSV_REL
    step00_dir = project / STEP00_DIR_REL
    step01_dir = project / STEP01_DIR_REL
    step02_dir = project / STEP02_DIR_REL
    step03_dir = project / STEP03_DIR_REL
    step04_dir = project / STEP04_DIR_REL
    out_dir = project / OUTPUT_DIR_REL
    log_dir = project / LOG_DIR_REL

    logger = Logger(log_dir / "synteny05_plot_karyotype.log")
    t0 = time.time()

    logger.info("========== synteny05 — dominant-ref color (S.constricta) + chr flip + no filtering ==========")
    logger.info(f"PROJECT_ROOT={project}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"JCVI_PYTHON={JCVI_PYTHON}")
    logger.info(f"FONT_FAMILY={FONT_FAMILY}")
    logger.info(f"SPECIES_LABEL_SIZE_PT={SPECIES_LABEL_SIZE_PT}")

    for p in (palette_tsv, step00_dir, step01_dir, step02_dir, step03_dir, step04_dir):
        if not p.exists():
            die(logger, f"Missing required input: {p}")

    palette = load_palette(palette_tsv)
    tracks = read_layout_tracks(step04_dir / "layout_species_tracks.tsv")
    track_species = [d["species_id"] for d in tracks]
    logger.info("Track order=" + " | ".join(track_species))

    # Step03 的 dominant-ref 映射与 flip
    chr2ref = load_chr2ref(step03_dir / "chr2ref" / "chr2ref.tsv")
    chr_flip = load_chr_flip(step03_dir / "chr2ref" / "chr_flip.tsv")

    clean_dir(out_dir, CLEAN_OUTPUT, logger)
    d_inputs = out_dir / "inputs"
    d_simple = out_dir / "simple"
    d_layout = out_dir / "layout"
    d_fig = out_dir / "figures"
    d_sum = out_dir / "summaries"
    for d in (d_inputs, d_simple, d_layout, d_fig, d_sum):
        d.mkdir(parents=True, exist_ok=True)

    # 预载 gene->chr（用于上色时把 gene 映射到 ref_chr）
    gene2chr: Dict[str, Dict[str, str]] = {}
    bed_sum = []
    for sid in track_species:
        src_bed = step01_dir / f"{sid}.geneorder.bed"
        if not src_bed.exists():
            die(logger, f"Missing Step01 bed: {src_bed}")
        gene2chr[sid] = load_gene_to_chr(src_bed)

    # 写 bed（按 chr_flip 做坐标反转；并维持 palette dups 行为）
    for sid in track_species:
        src_bed = step01_dir / f"{sid}.geneorder.bed"
        out_bed = d_inputs / f"{sid}.bed"

        step00_chr_rename = step00_dir / f"chr_rename_{sid}.tsv"
        chr_len_map = load_chr_lengths(step00_chr_rename)
        flip_map_chr = chr_flip.get(sid, {})

        n = copy_bed_with_palette_dups_and_flip(
            src_bed, out_bed, palette,
            flip_map_chr=flip_map_chr,
            chr_len_map=chr_len_map
        )
        bed_sum.append([sid, str(n), str(sum(1 for c, f in flip_map_chr.items() if f == 1))])
        logger.info(f"[BED] {sid}: genes={n} flip_chr={sum(1 for c,f in flip_map_chr.items() if f==1)}")

    write_tsv(d_sum / "step05.bed.summary.tsv", ["species_id", "n_genes", "n_flip_chr"], bed_sum)

    # 写 simple：dominant-ref 上色（按 A 侧 gene -> A chr -> ref_chr）
    link_rows: List[List[str]] = []
    for i in range(len(track_species) - 1):
        a = track_species[i]
        b = track_species[i + 1]
        pair_id = f"{a}__vs__{b}"

        src_simple = step02_dir / "raw_anchors" / f"{pair_id}.anchors.simple"
        if not src_simple.exists():
            die(logger, f"Missing Step02 raw anchors.simple: {src_simple}")

        out_simple = d_simple / f"{pair_id}.simple"
        kept, dropped = build_simple_dominant_ref_colored(
            a_species=a,
            src_simple=src_simple,
            out_simple=out_simple,
            gene2chr=gene2chr,
            chr2ref=chr2ref,
            palette=palette,
        )
        link_rows.append([pair_id, str(kept), str(dropped), src_simple.as_posix()])
        logger.info(f"[SIMPLE] {pair_id}: kept={kept} dropped={dropped}")

    write_tsv(d_sum / "step05.links.summary.tsv", ["pair_id", "kept_lines", "dropped_lines", "src_simple"], link_rows)

    # 写 seqids + layout
    seqids_path = d_layout / "seqids"
    layout_path = d_layout / "karyotype.layout"
    build_seqids_file(track_species, step03_dir / "seqids_species", seqids_path, logger)
    write_layout_file(tracks, layout_path, d_inputs, d_simple)

    pdf_tmp = run_jcvi_karyotype(seqids_path, layout_path, d_fig, n_tracks=len(track_species), logger=logger)
    pdf_out = d_fig / OUT_PDF
    if pdf_tmp != pdf_out:
        if pdf_out.exists():
            pdf_out.unlink()
        pdf_tmp.rename(pdf_out)

    logger.info(f"Output PDF: {pdf_out}")
    logger.info(f"Done. runtime_sec={int(time.time() - t0)}")


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
        lg = Logger(project / LOG_DIR_REL / "synteny05_plot_karyotype.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

