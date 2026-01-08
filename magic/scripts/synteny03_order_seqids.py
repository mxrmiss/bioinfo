#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
synteny03_order_seqids.py
—— Step03：基于 S. constricta 的 dominant-ref 调序 + chr 级 flip（仅读 Step02 anchors.simple，不改 Step02 输出）

输入（硬位置）：
- raw_data/synteny_species_meta.tsv
- output/synteny_01_mcscan_catalog/<species>.geneorder.bed
- output/synteny_02_blocks_macro/raw_anchors/<A>__vs__<B>.anchors.simple

输出（硬位置）：
output/synteny_03_order_seqids/
  seqids_species/<species>.seqids
  chr2ref/chr2ref.tsv
  chr2ref/chr_flip.tsv
  summaries/step03.summary.tsv

核心策略：
- 以 REFERENCE_SPECIES_ID（默认 Sinonovacula_constricta）作为全图唯一参考色带坐标系
- dominant ref：每条 chr 只映射到一个 ref_chr（票数最多）
- flip：chr 级别，沿链传播（异或），不设“最小连线数门槛”
"""

# ============================================================
# 【用户参数区】（皇上只改这里；脚本不接受命令行参数）
# ============================================================

PROJECT_ROOT = None
CLEAN_OUTPUT = True

REFERENCE_SPECIES_ID = "Pecten_maximus"

META_TSV_REL = "raw_data/synteny_species_meta.tsv"
STEP01_DIR_REL = "output/synteny_01_mcscan_catalog"
STEP02_DIR_REL = "output/synteny_02_blocks_macro"
OUTPUT_DIR_REL = "output/synteny_03_order_seqids"
LOG_DIR_REL = "logs"

# flip 判向最小支持点数：<3 时默认视为同向（不翻）
MIN_POINTS_FOR_ORIENTATION = 3

# ============================================================
# 实现区（皇上勿改）
# ============================================================

import csv
import time
import shutil
import traceback
import gzip
import re
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional


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


def chr_num(chr_name: str) -> int:
    m = _CHR_NUM_RE.match(chr_name)
    if m:
        try:
            return int(m.group(1))
        except Exception:
            return 10**9
    return 10**9


def normalize_gene_id(raw: str) -> str:
    x = (raw or "").strip().strip('"').strip("'").strip()
    if not x:
        return ""
    if "*" in x and x.startswith("#") and len(x) >= 8:
        # 兼容 #RRGGBB*gene
        try:
            x = x.split("*", 1)[1]
        except Exception:
            pass
    if "|" in x:
        x = x.split("|")[-1]
    if x.startswith("rna-"):
        x = x[4:]
    return x.strip()


def median(values: List[int]) -> Optional[float]:
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 1:
        return float(s[mid])
    return (float(s[mid - 1]) + float(s[mid])) / 2.0


def corr_sign(pairs: List[Tuple[int, int]]) -> int:
    """
    最小化方向判定：看协方差符号
    返回：+1（同向） 或 -1（反向）
    """
    if len(pairs) < MIN_POINTS_FOR_ORIENTATION:
        return +1
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    n = len(xs)
    mx = sum(xs) / float(n)
    my = sum(ys) / float(n)
    cov = 0.0
    vx = 0.0
    vy = 0.0
    for x, y in zip(xs, ys):
        dx = x - mx
        dy = y - my
        cov += dx * dy
        vx += dx * dx
        vy += dy * dy
    if vx <= 0.0 or vy <= 0.0:
        return +1
    return +1 if cov >= 0.0 else -1


def load_geneorder_bed(bed_path: Path) -> Tuple[Dict[str, str], Dict[str, int], List[str]]:
    """
    返回：
      gene2chr: gene -> ChrNN
      gene2pos: gene -> midpoint (int)
      chrs:     unique chrs sorted by Chr number
    """
    if not bed_path.exists():
        raise FileNotFoundError(f"Missing geneorder bed: {bed_path}")
    gene2chr: Dict[str, str] = {}
    gene2pos: Dict[str, int] = {}
    chr_set = set()

    with bed_path.open("r", encoding="utf-8", errors="replace") as fr:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chr_, s, e, gid = parts[:4]
            gid = normalize_gene_id(gid)
            if not gid:
                continue
            try:
                ss = int(float(s))
                ee = int(float(e))
            except Exception:
                continue
            mid = (ss + ee) // 2
            gene2chr[gid] = chr_
            gene2pos[gid] = mid
            chr_set.add(chr_)

    chrs = sorted(list(chr_set), key=lambda c: (chr_num(c), c))
    if not chrs:
        raise ValueError(f"No chromosomes extracted from: {bed_path}")
    if not gene2chr:
        raise ValueError(f"No genes parsed from: {bed_path}")
    return gene2chr, gene2pos, chrs


class PairInfo:
    """
    存储相邻 pair 的 chr 票数与方向证据：
      votes[(chrA, chrB)] = count
      pospairs[(chrA, chrB)] = list of (posA, posB)
    并提供：
      bestA2B[chrA] = (chrB, support)
      bestB2A[chrB] = (chrA, support)
      sign[(chrA, chrB)] = +1/-1
      medA[(chrA, chrB)] = median posA
      medB[(chrA, chrB)] = median posB
    """
    def __init__(self):
        self.votes: Dict[Tuple[str, str], int] = {}
        self.pospairs: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
        self.bestA2B: Dict[str, Tuple[str, int]] = {}
        self.bestB2A: Dict[str, Tuple[str, int]] = {}
        self.sign: Dict[Tuple[str, str], int] = {}
        self.medA: Dict[Tuple[str, str], Optional[float]] = {}
        self.medB: Dict[Tuple[str, str], Optional[float]] = {}

    def finalize(self) -> None:
        # bestA2B
        tmpA: Dict[str, List[Tuple[str, int]]] = {}
        tmpB: Dict[str, List[Tuple[str, int]]] = {}
        for (ca, cb), v in self.votes.items():
            tmpA.setdefault(ca, []).append((cb, v))
            tmpB.setdefault(cb, []).append((ca, v))
        for ca, lst in tmpA.items():
            lst.sort(key=lambda x: (-x[1], chr_num(x[0]), x[0]))
            self.bestA2B[ca] = lst[0]
        for cb, lst in tmpB.items():
            lst.sort(key=lambda x: (-x[1], chr_num(x[0]), x[0]))
            self.bestB2A[cb] = lst[0]

        # sign + median positions
        for (ca, cb), pairs in self.pospairs.items():
            self.sign[(ca, cb)] = corr_sign(pairs)
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            self.medA[(ca, cb)] = median(xs)
            self.medB[(ca, cb)] = median(ys)


def load_pair_info(
    a: str,
    b: str,
    anchors_simple: Path,
    gene2chr_a: Dict[str, str],
    gene2pos_a: Dict[str, int],
    gene2chr_b: Dict[str, str],
    gene2pos_b: Dict[str, int],
) -> PairInfo:
    pi = PairInfo()
    if not anchors_simple.exists():
        raise FileNotFoundError(f"Missing anchors.simple: {anchors_simple}")

    with open_text_auto(anchors_simple) as fr:
        for line in fr:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 6:
                continue
            a1, a2, b1, b2, _n, _o = parts
            a1 = normalize_gene_id(a1)
            a2 = normalize_gene_id(a2)
            b1 = normalize_gene_id(b1)
            b2 = normalize_gene_id(b2)

            # 两个锚点对 (a1,b1) 与 (a2,b2) 都计票
            for ga, gb in ((a1, b1), (a2, b2)):
                if ga not in gene2chr_a or gb not in gene2chr_b:
                    continue
                ca = gene2chr_a[ga]
                cb = gene2chr_b[gb]
                pi.votes[(ca, cb)] = pi.votes.get((ca, cb), 0) + 1
                if ga in gene2pos_a and gb in gene2pos_b:
                    pa = gene2pos_a[ga]
                    pb = gene2pos_b[gb]
                    pi.pospairs.setdefault((ca, cb), []).append((pa, pb))

    pi.finalize()
    return pi


def main() -> None:
    pr = PROJECT_ROOT
    if pr is None or str(pr).strip() == "":
        project = Path(__file__).resolve().parents[1]
    else:
        project = Path(str(pr)).expanduser().resolve()

    meta_tsv = project / META_TSV_REL
    step01_dir = project / STEP01_DIR_REL
    step02_dir = project / STEP02_DIR_REL
    out_dir = project / OUTPUT_DIR_REL
    log_dir = project / LOG_DIR_REL

    logger = Logger(log_dir / "synteny03_order_seqids.log")
    t0 = time.time()

    logger.info("========== synteny03 — dominant-ref ordering + chr flip (read Step02 only) ==========")
    logger.info(f"PROJECT_ROOT={project}")
    logger.info(f"CLEAN_OUTPUT={CLEAN_OUTPUT}")
    logger.info(f"REFERENCE_SPECIES_ID={REFERENCE_SPECIES_ID}")
    logger.info(f"MIN_POINTS_FOR_ORIENTATION={MIN_POINTS_FOR_ORIENTATION}")

    if not meta_tsv.exists():
        die(logger, f"Missing meta: {meta_tsv}")
    if not step01_dir.exists():
        die(logger, f"Missing step01 dir: {step01_dir}")
    if not step02_dir.exists():
        die(logger, f"Missing step02 dir: {step02_dir}")

    species = read_meta(meta_tsv)
    if REFERENCE_SPECIES_ID not in species:
        die(logger, f"REFERENCE_SPECIES_ID not found in meta order: {REFERENCE_SPECIES_ID}")
    ref_idx = species.index(REFERENCE_SPECIES_ID)
    logger.info("META species order=" + " | ".join(species))
    logger.info(f"Reference index={ref_idx+1}/{len(species)}")

    # 读 Step01 geneorder，建立每个物种的 gene->chr/pos + chr list
    gene2chr: Dict[str, Dict[str, str]] = {}
    gene2pos: Dict[str, Dict[str, int]] = {}
    chrs: Dict[str, List[str]] = {}
    for sid in species:
        bed = step01_dir / f"{sid}.geneorder.bed"
        g2c, g2p, clist = load_geneorder_bed(bed)
        gene2chr[sid] = g2c
        gene2pos[sid] = g2p
        chrs[sid] = clist
        logger.info(f"[BED] {sid}: n_chr={len(clist)} n_gene={len(g2c)}")

    # 读 Step02 的相邻 anchors.simple，缓存每条边的 PairInfo
    raw_anchors = step02_dir / "raw_anchors"
    if not raw_anchors.exists():
        die(logger, f"Missing step02 raw_anchors dir: {raw_anchors}")

    pairinfo: Dict[Tuple[str, str], PairInfo] = {}
    for i in range(len(species) - 1):
        a = species[i]
        b = species[i + 1]
        p = raw_anchors / f"{a}__vs__{b}.anchors.simple"
        pi = load_pair_info(
            a, b, p,
            gene2chr[a], gene2pos[a],
            gene2chr[b], gene2pos[b],
        )
        pairinfo[(a, b)] = pi
        logger.info(f"[PAIR] {a}__vs__{b}: n_chr_pairs={len(pi.votes)}")

    clean_dir(out_dir, CLEAN_OUTPUT, logger)
    d_seq = out_dir / "seqids_species"
    d_map = out_dir / "chr2ref"
    d_sum = out_dir / "summaries"
    for d in (d_seq, d_map, d_sum):
        d.mkdir(parents=True, exist_ok=True)

    # 结果容器
    chr_ref: Dict[str, Dict[str, str]] = {sid: {} for sid in species}
    chr_flip: Dict[str, Dict[str, int]] = {sid: {} for sid in species}
    chr_proj: Dict[str, Dict[str, Optional[float]]] = {sid: {} for sid in species}
    chr_support: Dict[str, Dict[str, int]] = {sid: {} for sid in species}

    # 初始化参考物种
    ref = REFERENCE_SPECIES_ID
    for c in chrs[ref]:
        chr_ref[ref][c] = c
        chr_flip[ref][c] = 0
        chr_proj[ref][c] = float(chr_num(c))
        chr_support[ref][c] = 0

    # 向左传播（更靠左的物种通过右邻（更靠近 ref）继承 ref_chr 与 flip）
    for idx in range(ref_idx - 1, -1, -1):
        cur = species[idx]
        nei = species[idx + 1]  # 更靠近 ref
        pi = pairinfo[(cur, nei)]  # 文件是 cur__vs__nei (A=cur,B=nei)

        for ccur in chrs[cur]:
            # dominant mapping：cur chr -> neighbor chr
            if ccur in pi.bestA2B:
                cnei, sup = pi.bestA2B[ccur]
            else:
                cnei, sup = ("", 0)

            if cnei and cnei in chr_ref[nei] and chr_ref[nei][cnei] != "NA":
                chr_ref[cur][ccur] = chr_ref[nei][cnei]
            else:
                chr_ref[cur][ccur] = "NA"

            chr_support[cur][ccur] = sup

            # projected position：用 neighbor 坐标（B 侧）中位数做组内排序依据
            proj = None
            if cnei:
                proj = pi.medB.get((ccur, cnei), None)
            chr_proj[cur][ccur] = proj

            # flip 传播：flip_cur = flip_nei XOR (sign==-1)
            flip_nei = chr_flip[nei].get(cnei, 0) if cnei else 0
            sgn = pi.sign.get((ccur, cnei), +1) if cnei else +1
            chr_flip[cur][ccur] = flip_nei ^ (1 if sgn == -1 else 0)

        logger.info(f"[PROP-L] {cur} <= {nei}: mapped={sum(1 for c in chrs[cur] if chr_ref[cur].get(c,'NA')!='NA')}")

    # 向右传播（更靠右的物种通过左邻（更靠近 ref）继承）
    for idx in range(ref_idx + 1, len(species)):
        cur = species[idx]
        nei = species[idx - 1]  # 更靠近 ref
        pi = pairinfo[(nei, cur)]  # 文件是 nei__vs__cur (A=nei,B=cur)

        for ccur in chrs[cur]:
            # dominant mapping：cur chr (B) -> neighbor chr (A)
            if ccur in pi.bestB2A:
                cnei, sup = pi.bestB2A[ccur]
            else:
                cnei, sup = ("", 0)

            if cnei and cnei in chr_ref[nei] and chr_ref[nei][cnei] != "NA":
                chr_ref[cur][ccur] = chr_ref[nei][cnei]
            else:
                chr_ref[cur][ccur] = "NA"

            chr_support[cur][ccur] = sup

            # projected position：用 neighbor 坐标（A 侧）中位数做组内排序依据
            proj = None
            if cnei:
                proj = pi.medA.get((cnei, ccur), None)
            chr_proj[cur][ccur] = proj

            flip_nei = chr_flip[nei].get(cnei, 0) if cnei else 0
            sgn = pi.sign.get((cnei, ccur), +1) if cnei else +1
            chr_flip[cur][ccur] = flip_nei ^ (1 if sgn == -1 else 0)

        logger.info(f"[PROP-R] {cur} <= {nei}: mapped={sum(1 for c in chrs[cur] if chr_ref[cur].get(c,'NA')!='NA')}")

    # 写 seqids：按 (ref_chr_rank, projected_pos, chr_num) 排序
    def ref_rank(ref_chr: str) -> int:
        return chr_num(ref_chr) if ref_chr and ref_chr != "NA" else 10**9

    for sid in species:
        def keyfunc(c: str) -> Tuple[int, float, int, str]:
            rr = ref_rank(chr_ref[sid].get(c, "NA"))
            pj = chr_proj[sid].get(c, None)
            pjv = float(pj) if pj is not None else 1e18
            cn = chr_num(c)
            return (rr, pjv, cn, c)

        ordered = sorted(chrs[sid], key=keyfunc)
        (d_seq / f"{sid}.seqids").write_text(",".join(ordered) + "\n", encoding="utf-8")
        logger.info(f"[SEQIDS] {sid}: n_chr={len(ordered)} mapped_ref={sum(1 for c in ordered if chr_ref[sid].get(c,'NA')!='NA')} flips={sum(chr_flip[sid].get(c,0) for c in ordered)}")

    # 写 chr2ref + chr_flip
    chr2ref_rows: List[List[str]] = []
    flip_rows: List[List[str]] = []
    for sid in species:
        # 保持输出顺序与 seqids 一致
        ordered = (d_seq / f"{sid}.seqids").read_text(encoding="utf-8").strip().split(",")
        for c in ordered:
            chr2ref_rows.append([sid, c, chr_ref[sid].get(c, "NA"), str(chr_support[sid].get(c, 0))])
            flip_rows.append([sid, c, str(chr_flip[sid].get(c, 0)), str(chr_support[sid].get(c, 0))])

    write_tsv(d_map / "chr2ref.tsv", ["species_id", "chr", "ref_chr", "support"], chr2ref_rows)
    write_tsv(d_map / "chr_flip.tsv", ["species_id", "chr", "flip", "support"], flip_rows)

    # summary
    sum_rows: List[List[str]] = []
    for sid in species:
        ordered = (d_seq / f"{sid}.seqids").read_text(encoding="utf-8").strip().split(",")
        n_chr = len(ordered)
        n_map = sum(1 for c in ordered if chr_ref[sid].get(c, "NA") != "NA")
        n_flip = sum(chr_flip[sid].get(c, 0) for c in ordered)
        sum_rows.append([sid, str(n_chr), str(n_map), str(n_flip), ordered[0], ordered[-1]])
    write_tsv(d_sum / "step03.summary.tsv", ["species_id", "n_chr", "n_mapped_ref", "n_flip", "first_chr", "last_chr"], sum_rows)

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
        lg = Logger(project / LOG_DIR_REL / "synteny03_order_seqids.log")
        lg.error("Unhandled exception: " + repr(e))
        lg.error(traceback.format_exc())
        raise

