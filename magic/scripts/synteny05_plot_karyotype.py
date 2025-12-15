#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny05_plot_karyotype.py
—— Step05：最终出图（相邻连线 + 防爆阀 + 参考红灰横条 + 贯穿连线配色）

输入（硬位置）：
- output/synteny_04_layout_tracks/layout_species_tracks.tsv
- output/synteny_03_order_seqids/global_chr_style.tsv
- output/synteny_03_order_seqids/seqids_species/<species>.seqids
- output/synteny_01_mcscan_catalog/<species>.geneorder.bed
- output/synteny_01_mcscan_catalog/geneorder_index_<species>.tsv
- output/synteny_02_blocks_macro/anchors_filtered/<A>__vs__<B>.anchors.simple.filtered
- output/synteny_02_blocks_macro/pair_chr_weight/<A>__vs__<B>.chr_weight.tsv

输出（硬位置）：
output/synteny_05_plot/
  bed_noheader/
  simple/
  layout/
  figures/ karyotype.pdf
  summaries/ step05.links.summary.tsv

硬合同要点（本脚本强制执行）：
- edges 只写相邻（track i 到 track i+1）
- simple 只保留两端 color_id 相同且 >0 的块，并在 gene token 前加 #HEX* 前缀
- 防爆阀门：MAX_LINKS_PER_PAIR（每对相邻最多保留多少条块）
  截断规则：优先保留 chr-pair total_span_bp 更大的块；同权重保持输入顺序（稳定截断）
- 横条颜色 + 标签颜色：读取 Step04 的 chr_bar_color_hex / label_color_hex
- seqids 翻转编码：ChrNN-（严格）；严禁 ChrNNr（发现立即报错退出）
- 图形：只输出 PDF；并强制去除白色圆圈与染色体数字（强制 --nocircles）
- 字体：默认 Arial（与 jcvi karyotype 的 --font 选项兼容）
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

# 1) 项目根目录
# - 留空/None：自动用脚本所在位置的上两级作为项目根（.../magic）
# - 如果你把脚本拷到别处运行，可以手动填绝对路径
PROJECT_ROOT = None

# 2) 是否清空并重建 Step05 输出目录
# - True：每次运行都 rm -rf output/synteny_05_plot 再重建（推荐）
# - False：保留旧结果，覆盖同名文件
CLEAN_OUTPUT = True

# 3) 防爆阀：每对相邻物种最多绘制多少条块（越大越容易“杂线爆炸”）
# - 建议先 2000~10000 试探，再逐步加
MAX_LINKS_PER_PAIR = 25000

# 4) JCVI 运行入口（一般写 "python" 即可；也可写绝对路径）
# - 关键：要指向装了 jcvi 的那个 Python
JCVI_PYTHON = "python"

# 5) 输出 PDF 文件名（只会在 figures/ 下生成这一份）
OUT_PDF = "karyotype.pdf"

# 6) 字体（默认 Arial）
# - 说明：jcvi karyotype 的 --font 不是系统任意字体，而是固定枚举
# - 你若要改，只能从下面列表选一个：
#   Arial / Helvetica / Liberation Sans / Palatino / Schoolbook
FONT_FAMILY = "Arial"

# 7) 画布大小（英寸）
# - FIGURE_WIDTH_INCH：宽度建议 12~20
# - FIGURE_HEIGHT_INCH：
#     - None：根据物种数自动估算高度（更省心）
#     - 指定数值：强制高度（适合你要对齐论文版式时）
FIGURE_WIDTH_INCH = 14.0
FIGURE_HEIGHT_INCH = None

# 8) 轨道横向占比（0~1）
# - TRACK_XSTART：轨道起点（越小越靠左，留给物种名的空间越大）
# - TRACK_XEND  ：轨道终点（越大越靠右）
TRACK_XSTART = 0.10
TRACK_XEND = 0.92

# 9) 物种名的垂直对齐方式（传给 layout 的 va 字段）
# - 常用：top / center / bottom
LABEL_VA = "top"

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
from typing import Dict, List, Tuple, Optional, Iterable


_ALLOWED_JCVI_FONTS = {"Helvetica", "Liberation Sans", "Palatino", "Schoolbook", "Arial"}
_HEX_PREFIX_RE = re.compile(r"^#[0-9a-fA-F]{6}\*")
_HEX_COLOR_RE = re.compile(r"^#[0-9a-fA-F]{6}$")
_SEQ_TOKEN_RE = re.compile(r"^Chr\d+(?:-)?$")


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


def strip_color_prefix(token: str) -> str:
    token = token.strip()
    if _HEX_PREFIX_RE.match(token):
        return token.split("*", 1)[1]
    return token


def normalize_gene_id(raw: str) -> str:
    x = raw.strip().strip('"').strip("'").strip()
    if "|" in x:
        x = x.split("|")[-1]
    if x.startswith("rna-"):
        x = x[4:]
    return x


def validate_seqids_line(line: str) -> Tuple[List[str], List[str]]:
    """
    防线：允许 ChrNN / ChrNN-；禁止 ChrNNr；禁止空 token；禁止重复 base Chr
    返回：(tokens, base_tokens)
    """
    raw = (line or "").strip()
    if not raw:
        raise ValueError("seqids line is empty")

    tokens = [t.strip() for t in raw.split(",")]
    if any(t == "" for t in tokens):
        raise ValueError(f"seqids contains empty token (extra commas?): {raw}")

    base_seen = set()
    base_tokens: List[str] = []
    for t in tokens:
        if t.endswith("r") or t.endswith("R"):
            raise ValueError(f"seqids token uses 'r' (FORBIDDEN). Use '-' for flip: {t}")
        if not _SEQ_TOKEN_RE.match(t):
            raise ValueError(f"illegal seqids token: {t} (must be ChrNN or ChrNN-)")
        base = t[:-1] if t.endswith("-") else t
        if base in base_seen:
            raise ValueError(f"duplicated chromosome in seqids: {base}")
        base_seen.add(base)
        base_tokens.append(base)
    return tokens, base_tokens


def read_layout_tracks(path: Path) -> List[Dict[str, str]]:
    """
    Step04 layout_species_tracks.tsv
    header:
      order_index species_id display_label group y_center track_height chr_bar_color_hex label_color_hex
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing layout tracks: {path}")
    rows: List[Dict[str, str]] = []
    need = ["order_index", "species_id", "display_label", "group", "y_center", "track_height", "chr_bar_color_hex", "label_color_hex"]
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
            d = dict(zip(need, parts))
            rows.append(d)
    rows.sort(key=lambda d: int(d["order_index"]))
    return rows


def load_global_chr_style(path: Path) -> Dict[str, Dict[str, Tuple[int, str]]]:
    """
    global_chr_style.tsv
    header:
      species_id chr order_index orientation dominant_prev_chr dominant_weight_span_bp color_id color_hex
    返回：
      style[species][chr] = (color_id, color_hex)
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing global chr style: {path}")
    style: Dict[str, Dict[str, Tuple[int, str]]] = {}
    with open_text_auto(path) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        need = ["species_id", "chr", "order_index", "orientation", "dominant_prev_chr", "dominant_weight_span_bp", "color_id", "color_hex"]
        if header != need:
            raise ValueError(f"Bad header in {path.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 8:
                continue
            sid = parts[0]
            chr_ = parts[1]
            cid = int(parts[6]) if parts[6].isdigit() else 0
            chex = parts[7]
            style.setdefault(sid, {})[chr_] = (cid, chex)
    return style


def load_gene_index(step01_dir: Path, sid: str) -> Dict[str, str]:
    """
    geneorder_index_<sid>.tsv：gene_id -> chr
    """
    p = step01_dir / f"geneorder_index_{sid}.tsv"
    if not p.exists():
        raise FileNotFoundError(f"Missing gene index: {p}")
    mp: Dict[str, str] = {}
    with open_text_auto(p) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        if header[:2] != ["gene_id", "chr"]:
            raise ValueError(f"Bad header in {p.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gid = normalize_gene_id(strip_color_prefix(parts[0]))
            chr_ = parts[1]
            if gid:
                mp[gid] = chr_
    if not mp:
        raise ValueError(f"Empty gene index: {p}")
    return mp


def load_pair_chr_weight(step02_dir: Path, pair_id: str) -> Dict[Tuple[str, str], int]:
    """
    pair_chr_weight/<pair>.chr_weight.tsv
    header: pair_id a_chr b_chr total_span_bp n_blocks
    返回： (a_chr,b_chr) -> total_span_bp
    """
    p = step02_dir / "pair_chr_weight" / f"{pair_id}.chr_weight.tsv"
    if not p.exists():
        raise FileNotFoundError(f"Missing chr_weight: {p}")
    mp: Dict[Tuple[str, str], int] = {}
    with open_text_auto(p) as fr:
        header = fr.readline().rstrip("\n").split("\t")
        need = ["pair_id", "a_chr", "b_chr", "total_span_bp", "n_blocks"]
        if header != need:
            raise ValueError(f"Bad header in {p.name}: {header}")
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 5:
                continue
            a_chr, b_chr = parts[1], parts[2]
            w = int(parts[3]) if parts[3].isdigit() else 0
            mp[(a_chr, b_chr)] = w
    return mp


def prepare_bed_with_color_dups(
    sid: str,
    step01_dir: Path,
    chr_style: Dict[str, Tuple[int, str]],
    out_bed: Path,
    logger: Logger
) -> None:
    """
    bed_noheader/<sid>.bed
    - 先写原始 gene_id
    - 若该 gene 所在 chr 的 color_id>0 且 color_hex 合法(#RRGGBB)，则再写一条 #HEX*gene_id 的镜像记录
    """
    src = step01_dir / f"{sid}.geneorder.bed"
    if not src.exists():
        raise FileNotFoundError(f"Missing geneorder bed: {src}")

    n_in = 0
    n_out = 0
    n_pref = 0
    n_bad_hex = 0

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with src.open("r", encoding="utf-8", errors="replace") as fr, out_bed.open("w", encoding="utf-8", newline="") as fw:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            n_in += 1
            chr_, start, end, gid = parts[0], parts[1], parts[2], parts[3]
            gid0 = normalize_gene_id(strip_color_prefix(gid))

            fw.write("\t".join([chr_, start, end, gid0]) + "\n")
            n_out += 1

            if chr_ in chr_style:
                cid, chex = chr_style[chr_]
                if cid > 0 and chex != "NA":
                    if not _HEX_COLOR_RE.match(chex):
                        n_bad_hex += 1
                    else:
                        fw.write("\t".join([chr_, start, end, f"{chex}*{gid0}"]) + "\n")
                        n_out += 1
                        n_pref += 1

    if n_bad_hex > 0:
        logger.warn(f"[BED] {sid}: found {n_bad_hex} invalid color_hex (must be #RRGGBB); colored duplicates skipped for those")
    logger.info(f"[BED] {sid}: in_lines={n_in} out_lines={n_out} colored_dups={n_pref}")


def build_simple_for_pair(
    a: str,
    b: str,
    step02_dir: Path,
    gene_chr_a: Dict[str, str],
    gene_chr_b: Dict[str, str],
    style_a: Dict[str, Tuple[int, str]],
    style_b: Dict[str, Tuple[int, str]],
    pair_weight: Dict[Tuple[str, str], int],
    out_simple: Path,
    logger: Logger
) -> Tuple[int, int, int, int, int]:
    """
    anchors_filtered/<pair>.anchors.simple.filtered
    6列：a1 a2 b1 b2 n orient

    过滤规则（蓝图）：
    - 四个端点基因都必须存在于 gene_index（硬防线）
    - A 端点两基因必须在同一 chr；B 端点两基因必须在同一 chr（不自洽则丢弃）
    - 两端 chr 的 color_id 相同 且 color_id>0
    - 输出时给基因 token 前加 #HEX* 前缀（HEX 为该 color_id 对应的 color_hex）

    防爆截断（蓝图）：
    - 若候选块数 > MAX_LINKS_PER_PAIR：
      先按 chrpair total_span_bp 降序，再按输入顺序稳定截断
    """
    pair_id = f"{a}__vs__{b}"
    in_path = step02_dir / "anchors_filtered" / f"{pair_id}.anchors.simple.filtered"
    if not in_path.exists():
        raise FileNotFoundError(f"Missing anchors_filtered: {in_path}")

    anchors_in = 0
    dropped_missing = 0
    dropped_inconsistent = 0
    candidates: List[Tuple[int, int, List[str]]] = []  # (weight, order, parts_out)
    order = 0

    with open_text_auto(in_path) as fr:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                continue
            anchors_in += 1
            order += 1

            ga1 = normalize_gene_id(strip_color_prefix(parts[0]))
            ga2 = normalize_gene_id(strip_color_prefix(parts[1]))
            gb1 = normalize_gene_id(strip_color_prefix(parts[2]))
            gb2 = normalize_gene_id(strip_color_prefix(parts[3]))

            if ga1 not in gene_chr_a or ga2 not in gene_chr_a or gb1 not in gene_chr_b or gb2 not in gene_chr_b:
                dropped_missing += 1
                continue

            a_chr1 = gene_chr_a[ga1]
            a_chr2 = gene_chr_a[ga2]
            b_chr1 = gene_chr_b[gb1]
            b_chr2 = gene_chr_b[gb2]

            if a_chr1 != a_chr2 or b_chr1 != b_chr2:
                dropped_inconsistent += 1
                continue

            a_chr = a_chr1
            b_chr = b_chr1

            if a_chr not in style_a or b_chr not in style_b:
                dropped_missing += 1
                continue

            ca_id, ca_hex = style_a[a_chr]
            cb_id, cb_hex = style_b[b_chr]

            if ca_id <= 0 or cb_id <= 0 or ca_id != cb_id:
                continue

            chex = ca_hex
            if chex == "NA" or (not _HEX_COLOR_RE.match(chex)):
                continue

            out_parts = [
                f"{chex}*{ga1}",
                f"{chex}*{ga2}",
                f"{chex}*{gb1}",
                f"{chex}*{gb2}",
                parts[4],
                parts[5],
            ]

            w = pair_weight.get((a_chr, b_chr), 0)
            candidates.append((w, order, out_parts))

    n_cand = len(candidates)
    if n_cand == 0:
        logger.warn(
            f"[{pair_id}] no blocks after color-chain filter "
            f"(anchors_in_filtered={anchors_in}, dropped_missing_endpoints={dropped_missing}, dropped_inconsistent_chr={dropped_inconsistent})"
        )
        out_simple.parent.mkdir(parents=True, exist_ok=True)
        out_simple.write_text("", encoding="utf-8")
        return anchors_in, 0, 0, dropped_missing, dropped_inconsistent

    candidates.sort(key=lambda x: (-x[0], x[1]))
    if n_cand > MAX_LINKS_PER_PAIR:
        candidates = candidates[:MAX_LINKS_PER_PAIR]

    out_simple.parent.mkdir(parents=True, exist_ok=True)
    with out_simple.open("w", encoding="utf-8", newline="") as fw:
        for _w, _ord, parts_out in candidates:
            fw.write("\t".join(parts_out) + "\n")

    drawn = len(candidates)
    logger.info(
        f"[{pair_id}] anchors_in_filtered={anchors_in} candidates={n_cand} drawn={drawn} "
        f"dropped_missing_endpoints={dropped_missing} dropped_inconsistent_chr={dropped_inconsistent}"
    )
    return anchors_in, n_cand, drawn, dropped_missing, dropped_inconsistent


def write_seqids_file(track_species: List[str], seqids_species_dir: Path, out_seqids: Path, logger: Logger) -> None:
    """
    jcvi karyotype 的 seqids 文件：每行一个 track（逗号分隔 chr 列表）
    加硬防线：只允许 ChrNN / ChrNN-；禁止 ChrNNr；并校验 base chr 不重复
    """
    lines: List[str] = []
    for sid in track_species:
        p = seqids_species_dir / f"{sid}.seqids"
        if not p.exists():
            raise FileNotFoundError(f"Missing seqids for {sid}: {p}")
        line = p.read_text(encoding="utf-8", errors="replace").strip()
        tokens, _bases = validate_seqids_line(line)
        logger.info(f"[SEQIDS-OK] {sid}: n_tokens={len(tokens)} flip_tokens={sum(1 for x in tokens if x.endswith('-'))}")
        lines.append(line)
    out_seqids.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_layout_file(
    tracks: List[Dict[str, str]],
    out_layout: Path,
    bed_dir: Path,
    simple_dir: Path,
    logger: Logger
) -> None:
    """
    JCVI layout（tracks 行 + edges 行）
    tracks 行格式（jcvi help）：
      y, xstart, xend, rotation, color, label, va, bed, label_va
    edges 行：
      e, i, i+1, simplefile
    """
    rot = 0
    lines: List[str] = []

    for i, d in enumerate(tracks):
        sid = d["species_id"]
        y = float(d["y_center"])
        bar_hex = d["chr_bar_color_hex"]
        label_hex = d["label_color_hex"]
        label_text = d["display_label"]

        if not _HEX_COLOR_RE.match(bar_hex):
            logger.warn(f"[LAYOUT] {sid}: chr_bar_color_hex not #RRGGBB: {bar_hex}")
        if not _HEX_COLOR_RE.match(label_hex):
            logger.warn(f"[LAYOUT] {sid}: label_color_hex not #RRGGBB: {label_hex}")

        label = f"{label_hex}*{label_text}"
        bed = bed_dir / f"{sid}.bed"

        lines.append(
            f"{y:.6f}, {TRACK_XSTART:.2f}, {TRACK_XEND:.2f}, {rot}, {bar_hex}, {label}, {LABEL_VA}, {bed.as_posix()}, center"
        )

    for i in range(len(tracks) - 1):
        a = tracks[i]["species_id"]
        b = tracks[i + 1]["species_id"]
        simple = simple_dir / f"{a}__vs__{b}.simple"
        lines.append(f"e, {i}, {i+1}, {simple.as_posix()}")

    out_layout.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _estimate_height_inch(n_tracks: int) -> float:
    h = max(6.0, 0.66 * float(n_tracks))
    return round(h, 2)


def run_jcvi_karyotype(seqids: Path, layout: Path, fig_dir: Path, n_tracks: int, logger: Logger) -> Path:
    """
    通过 subprocess 调用 jcvi.graphics.karyotype
    - 只生成 PDF
    - 强制 --nocircles（去掉白色圆圈与染色体数字）
    - 强制 --font=Arial（默认，可由用户参数区改为 jcvi 支持的字体枚举）
    """
    fig_dir.mkdir(parents=True, exist_ok=True)

    if FONT_FAMILY not in _ALLOWED_JCVI_FONTS:
        die(logger, f"FONT_FAMILY illegal: {FONT_FAMILY}. Allowed: {', '.join(sorted(_ALLOWED_JCVI_FONTS))}")

    w = float(FIGURE_WIDTH_INCH)
    h = float(FIGURE_HEIGHT_INCH) if FIGURE_HEIGHT_INCH is not None else _estimate_height_inch(n_tracks)

    cmd = [
        JCVI_PYTHON, "-m", "jcvi.graphics.karyotype",
        str(seqids), str(layout),
        "--format=pdf",
        "--notex",
        f"--font={FONT_FAMILY}",
        f"--figsize={int(round(w))}x{int(round(h))}",
        "--nocircles",
    ]

    logger.info("[CMD] " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(fig_dir), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    out = (r.stdout or "").rstrip("\n")
    if out:
        noisy = (
            "pkg_resources is deprecated",
            "brewer2mpl/brewer2mpl.py",
            "setuptools.pypa.io",
            "slated for removal",
        )
        for ln in out.splitlines():
            if any(x in ln for x in noisy):
                continue
            logger.info(f"[JCVI] {ln}")

    if r.returncode != 0:
        raise RuntimeError(f"jcvi karyotype failed (exit={r.returncode})")

    pdf = None
    for p in sorted(fig_dir.glob("*.pdf")):
        pdf = p
        break
    if pdf is None or (not pdf.exists()):
        raise RuntimeError("jcvi did not produce any PDF in figures/")

    return pdf


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        project = Path(__file__).resolve().parents[1]
    else:
        project = Path(str(pr)).expanduser().resolve()

    step01_dir = project / "output/synteny_01_mcscan_catalog"
    step02_dir = project / "output/synteny_02_blocks_macro"
    step03_dir = project / "output/synteny_03_order_seqids"
    step04_dir = project / "output/synteny_04_layout_tracks"
    out_dir = project / "output/synteny_05_plot"
    log_dir = project / "logs"

    logger = Logger(log_dir / "synteny05_plot_karyotype.log")
    t0 = time.time()

    logger.info("========== synteny05 — plot karyotype (adjacent + defense) ==========")
    logger.info(f"PROJECT_ROOT={project}")
    logger.info(f"OUTPUT_DIR={out_dir}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"MAX_LINKS_PER_PAIR={MAX_LINKS_PER_PAIR}")
    logger.info(f"JCVI_PYTHON={JCVI_PYTHON}")
    logger.info(f"OUT_PDF={OUT_PDF}")
    logger.info(f"FONT_FAMILY={FONT_FAMILY}")
    logger.info(f"FIGURE_WIDTH_INCH={FIGURE_WIDTH_INCH}")
    logger.info(f"FIGURE_HEIGHT_INCH={FIGURE_HEIGHT_INCH}")
    logger.info("Seqids defense: allow ChrNN / ChrNN- ; forbid ChrNNr")
    logger.info("Plot defense: force --nocircles (remove white circles and chr numbers)")

    for p in (step01_dir, step02_dir, step03_dir, step04_dir):
        if not p.exists():
            die(logger, f"Missing required dir: {p}")

    layout_tracks_path = step04_dir / "layout_species_tracks.tsv"
    global_style_path = step03_dir / "global_chr_style.tsv"
    seqids_species_dir = step03_dir / "seqids_species"

    tracks = read_layout_tracks(layout_tracks_path)
    track_species = [d["species_id"] for d in tracks]
    logger.info("Track order=" + " | ".join(track_species))

    clean_dir(out_dir, CLEAN_OUTPUT, logger)
    bed_dir = out_dir / "bed_noheader"
    simple_dir = out_dir / "simple"
    layout_dir = out_dir / "layout"
    fig_dir = out_dir / "figures"
    summ_dir = out_dir / "summaries"
    for d in (bed_dir, simple_dir, layout_dir, fig_dir, summ_dir):
        d.mkdir(parents=True, exist_ok=True)

    global_style = load_global_chr_style(global_style_path)

    gene_chr: Dict[str, Dict[str, str]] = {}
    for sid in track_species:
        gene_chr[sid] = load_gene_index(step01_dir, sid)

    for sid in track_species:
        if sid not in global_style:
            die(logger, f"Missing {sid} in global_chr_style.tsv")
        prepare_bed_with_color_dups(sid, step01_dir, global_style[sid], bed_dir / f"{sid}.bed", logger)

    link_rows: List[List[str]] = []
    for i in range(len(track_species) - 1):
        a = track_species[i]
        b = track_species[i + 1]
        pair_id = f"{a}__vs__{b}"

        pair_w = load_pair_chr_weight(step02_dir, pair_id)
        out_simple = simple_dir / f"{pair_id}.simple"

        anchors_in, cand, drawn, drop_missing, drop_incon = build_simple_for_pair(
            a=a,
            b=b,
            step02_dir=step02_dir,
            gene_chr_a=gene_chr[a],
            gene_chr_b=gene_chr[b],
            style_a=global_style[a],
            style_b=global_style[b],
            pair_weight=pair_w,
            out_simple=out_simple,
            logger=logger
        )

        link_rows.append([
            pair_id,
            str(anchors_in),
            str(cand),
            str(drawn),
            str(MAX_LINKS_PER_PAIR),
            str(drop_missing),
            str(drop_incon),
        ])

    write_tsv(
        summ_dir / "step05.links.summary.tsv",
        ["pair_id", "anchors_in_filtered", "anchors_candidates", "anchors_drawn", "max_links_per_pair", "dropped_missing_endpoints", "dropped_inconsistent_chr"],
        link_rows
    )

    seqids_path = layout_dir / "seqids"
    layout_path = layout_dir / "karyotype.layout"

    try:
        write_seqids_file(track_species, seqids_species_dir, seqids_path, logger)
    except Exception as e:
        die(logger, f"Failed to build seqids file: {e}")

    write_layout_file(tracks, layout_path, bed_dir, simple_dir, logger)

    logger.info(f"Generated seqids: {seqids_path}")
    logger.info(f"Generated layout: {layout_path}")

    pdf_tmp = run_jcvi_karyotype(seqids_path, layout_path, fig_dir, n_tracks=len(track_species), logger=logger)

    pdf_out = fig_dir / OUT_PDF
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
        lg = Logger(project / "logs" / "synteny05_plot_karyotype.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

