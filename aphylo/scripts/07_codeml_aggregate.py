#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
07_codeml_aggregate.py — 汇总 codeml（显式 alt/null 配对；输出到 codeml_agg_dir）

功能概述：
1) 扫描 codeml ALT/NULL 结果，解析 lnL，计算 LRT / P / Q（混合 χ²₁）。
2) 写出：
   - D_fdr_genes.tsv  —— OG × foreground 的 PSG 统计表
   - D_beb_sites.tsv  —— BEB 正选择位点表
3) 若 config.psg.enable = true，则额外：
   - 为每个前景 FG 构建 PSG cds 富集输入，符合 vNext ORA 接口形状：
       codeml_agg_dir / psg_cds_inputs / psg_<FG> /       （以 q 值显著为标准）
         ├─ test.list        # PSG cds_id 集合
         ├─ background.list  # 背景 cds_id 集合
         └─ meta.tsv         # 元信息说明
       codeml_agg_dir / psg_cds_inputs_rawp / psg_rawp_<FG> / （以原始 P 值为标准）
         ├─ test.list
         ├─ background.list
         └─ meta.tsv

注意：
- 本脚本只处理到 cds 粒度；cds → gene_id 的映射由下游脚本完成。
"""

from __future__ import annotations
import sys, logging, re, math
from pathlib import Path
from typing import Dict, Any, List, Tuple
import yaml
from collections import defaultdict

DEFAULT_CONFIG = "config.yaml"

# ===================== 通用工具 =====================

def _expand_publish_placeholders(obj, publish_dir: str):
    """递归替换字符串中的 <publish_dir> 占位符。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj

def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    """读取 config.yaml，并展开 publish_dir 相关占位符。"""
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    """确保目录存在。"""
    p.mkdir(parents=True, exist_ok=True)
    return p

def need_dir(p: Path, what: str) -> Path:
    """要求目录必须存在，否则报错。"""
    p = Path(p)
    if not p.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p

def need_file(p: Path, what: str) -> Path:
    """要求文件必须存在，否则报错。"""
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p

def get_logger(name: str, log_file: Path) -> logging.Logger:
    """配置日志：同时输出到文件和标准输出。"""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def banner(log: logging.Logger, msg: str):
    bar = "=" * 60
    log.info(bar)
    log.info(msg)
    log.info(bar)

def write_done(path: Path):
    """写一个 .done 哨兵文件，用于标记流程完成。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

# ===================== lnL / P / Q 解析 =====================

NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
FLOAT_RE = re.compile(NUM)

def last_lnl(txt: str) -> float | None:
    """
    解析 mlc 文本中的 lnL 值：
      - 优先匹配形如：lnL(ntime: 27  np: 32):  -8539.712388  +0.000000
      - 兼容简写：lnL = -8539.712388
      - 若多次出现，只取“最后一次”（视为最终优化值）
    """
    m1 = re.findall(rf"^lnL\s*\([^\n]*?\)\s*:\s*({NUM})", txt, flags=re.M)
    m2 = re.findall(rf"^lnL\s*=\s*({NUM})", txt, flags=re.M)
    m = (m1 if m1 else []) + (m2 if m2 else [])
    if not m:
        return None
    try:
        return float(m[-1])
    except Exception:
        return None

def chi2_mix_half_df1_sf(x: float) -> float:
    """
    0.5 * χ²₁ 的右尾概率（Self & Liang 1987）：
      p = 0.5 * erfc(sqrt(x/2))
    """
    if x <= 0:
        # LRT <= 0 按 0 处理，对应混合分布中 P = 0.5
        return 0.5
    return 0.5 * math.erfc(math.sqrt(x / 2.0))

def parse_beb_lines(mlc_text: str, cutoff: float = 0.95) -> List[Tuple[int, str, float]]:
    """
    解析 BEB 正选择位点表（逐行解析，不跨行）：

    支持格式示例（PAML 4.9 / 4.10）：
      Bayes Empirical Bayes (BEB) analysis ...
      Positive sites for foreground lineages Prob(w>1):
          123 K 0.999*
          234 L 0.980
      或类似不带表头的纯数据块，只要第一列是数字、第二列为氨基酸即可。

    注意：
      - 只在“BEB/Positive sites …”块内部解析；
      - 仅当第二列是氨基酸字母时才视为 BEB 位点，避免将 1 1 1.000 之类数字表误解析；
      - 仅保留 posterior ≥ cutoff 的位点。
    """
    rows: List[Tuple[int, str, float]] = []
    in_beb_block = False

    for raw in mlc_text.splitlines():
        line = raw.rstrip()

        # 进入 BEB 区（宽松匹配关键字）
        if re.search(r"Bayes Empirical Bayes|BEB analysis|Positively selected sites|Positive sites for foreground lineages", line, flags=re.I):
            in_beb_block = True
            continue

        # 离开 BEB 块：空行或开始其它小节
        if in_beb_block:
            if not line.strip():
                in_beb_block = False
                continue
            if re.match(r"\s*(The grid|Posterior on the grid|Posterior for p0-p1)", line):
                # BEB 位点段结束，后面是数值网格/后验表，不应解析为位点
                in_beb_block = False
                continue

        if not in_beb_block:
            continue

        parts = line.split()
        if len(parts) < 3:
            continue
        # 第一列必须是位点编号（纯数字）
        if not parts[0].isdigit():
            continue
        site = int(parts[0])

        aa = parts[1]
        # 第二列必须是氨基酸字母而不是数字，避免把 "1 1 1.000" 之类的行当成位点
        if not re.fullmatch(r"[A-Za-z\*]+", aa):
            continue

        # 从行尾往前找数值 token（剥离星标 *）
        post_tok = None
        for tok in reversed(parts[2:]):
            t = tok.rstrip("*")
            if FLOAT_RE.fullmatch(t):
                post_tok = t
                break
        if post_tok is None:
            continue

        try:
            post = float(post_tok)
        except Exception:
            continue

        if post >= cutoff:
            rows.append((site, aa, post))

    return rows

def bh_fdr(pvals: List[float]) -> List[float]:
    """
    Benjamini–Hochberg FDR 校正（单调回填版本）：
      - pvals: 原始 p-value 列表
      - 返回同长度的 q-value 列表
    """
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])  # 升序索引
    q = [1.0] * m
    min_q = 1.0
    for rank_from_end, i in enumerate(reversed(order), start=1):
        j = m - rank_from_end + 1  # 实际排名
        val = pvals[i] * m / j
        if val < min_q:
            min_q = val
        q[i] = min_q if min_q < 1.0 else 1.0
    return q

# ===================== PSG cds 富集输入构建 =====================

def build_psg_cds_inputs(cfg: Dict[str, Any],
                         paths: Dict[str, Any],
                         out_dir: Path,
                         rows_genes: List[List[str]],
                         qvals: List[float],
                         log: logging.Logger) -> None:
    """
    基于 codeml 聚合结果（rows_genes + qvals），构建 PSG cds 富集输入：

      1）q 值显著版（现有发布包）：
          codeml_agg_dir / psg_cds_inputs / psg_<FG> /
            ├─ test.list        # PSG cds_id 集合（cds 空间，q <= fdr_alpha）
            ├─ background.list  # 背景 cds_id 集合（cds 空间）
            └─ meta.tsv

      2）原始 P 值显著版（新增发布包）：
          codeml_agg_dir / psg_cds_inputs_rawp / psg_rawp_<FG> /
            ├─ test.list        # PSG cds_id 集合（cds 空间，P <= p_alpha）
            ├─ background.list  # 背景 cds_id 集合（cds 空间）
            └─ meta.tsv

    注意：
      - 本函数只在 cfg["psg"]["enable"] 为 True 且 rows_genes 非空时调用；
      - 使用 02 的 order/OG.order.tsv + 05 的 sets/{FG}.list；
      - 不做 cds -> gene 映射。
    """
    if not rows_genes:
        log.info("[PSG] rows_genes 为空，跳过 PSG cds 富集输入构建。")
        return

    psg_cfg = (cfg.get("psg") or {})
    if not psg_cfg.get("enable", False):
        log.info("[PSG] config.psg.enable 未开启，跳过 PSG cds 富集输入构建。")
        return

    # 解析 PSG FDR 阈值：优先 psg.fdr_alpha，其次 report.fdr_alpha，最后默认 0.05
    fdr_alpha = psg_cfg.get("fdr_alpha")
    if fdr_alpha is None:
        report_cfg = cfg.get("report") or {}
        fdr_alpha = report_cfg.get("fdr_alpha", 0.05)
        if "fdr_alpha" not in psg_cfg and "fdr_alpha" not in report_cfg:
            log.warning("[PSG] 未在 config 中找到 psg.fdr_alpha / report.fdr_alpha，使用默认 0.05。")
    try:
        fdr_alpha = float(fdr_alpha)
    except Exception:
        log.warning(f"[PSG] 无法解析 FDR 阈值（{fdr_alpha!r}），回退为 0.05。")
        fdr_alpha = 0.05

    # 解析 PSG 原始 P 值阈值：psg.p_alpha 优先，否则默认 0.01
    p_alpha = psg_cfg.get("p_alpha")
    if p_alpha is None:
        p_alpha = 0.01
        log.warning("[PSG] 未在 config 中设置 psg.p_alpha，使用默认原始 P 值阈值 0.01。")
    try:
        p_alpha = float(p_alpha)
    except Exception:
        log.warning(f"[PSG] 无法解析原始 P 值阈值（{p_alpha!r}），回退为 0.01。")
        p_alpha = 0.01

    # Step 1: 按 FG 汇总 OG 背景 & PSG 集合（分别记 q 阈值与 p 阈值）
    bg_ogs: Dict[str, set]    = defaultdict(set)
    sig_ogs_q: Dict[str, set] = defaultdict(set)
    sig_ogs_p: Dict[str, set] = defaultdict(set)

    for row, q in zip(rows_genes, qvals):
        if len(row) < 4:
            continue
        og, fg, _, p_str = row
        try:
            qv = float(q)
            pv = float(p_str)
        except Exception:
            continue
        bg_ogs[fg].add(og)
        if qv <= fdr_alpha:
            sig_ogs_q[fg].add(og)
        if pv <= p_alpha:
            sig_ogs_p[fg].add(og)

    if not bg_ogs:
        log.warning("[PSG] D_fdr_genes.tsv 中没有任何可用记录，跳过 PSG cds 富集输入构建。")
        return

    # Step 2: 前景物种集合（sets/{FG}.list）
    codeml_dir = Path(paths["codeml_dir"])
    sets_dir = codeml_dir / "sets"
    if not sets_dir.is_dir():
        log.error(f"[PSG] 缺少前景物种目录：{sets_dir}，无法构建 PSG cds 富集输入。")
        return

    fg_species: Dict[str, set] = {}
    for fg in bg_ogs.keys():
        fp = sets_dir / f"{fg}.list"
        if not fp.is_file():
            log.error(f"[PSG] 找不到前景 {fg} 的物种列表文件：{fp}，该前景将被跳过。")
            continue
        tips = [ln.strip() for ln in fp.read_text(encoding="utf-8").splitlines() if ln.strip()]
        if not tips:
            log.warning(f"[PSG] 前景 {fg} 的物种集合为空：{fp}")
        fg_species[fg] = set(tips)

    # Step 3: OG -> [(Species, cds_id)]：基于 02 产生的 order/OG.order.tsv
    bt_dir = Path(paths["bt_dir"])
    order_dir = bt_dir / "order"
    if not order_dir.is_dir():
        log.error(f"[PSG] 缺少 codon 输入的 order 目录：{order_dir}，无法构建 PSG cds 富集输入。")
        return

    og2rows: Dict[str, List[Tuple[str, str]]] = {}

    def get_og_rows(og: str) -> List[Tuple[str, str]]:
        """
        懒加载某个 OG 的 (Species, cds_id) 列表。

        修正要点：
          - 只在“第一列恰好为 'OG'”时认定为表头；
          - 真正的数据行是类似 'OG0006410' 这种，不再被误判为表头。
        """
        if og in og2rows:
            return og2rows[og]

        tsv = order_dir / f"{og}.order.tsv"
        if not tsv.is_file():
            log.warning(f"[PSG] 找不到 OG={og} 的 order.tsv：{tsv}，该 OG 将在 PSG 展开中被跳过。")
            og2rows[og] = []
            return og2rows[og]

        rows: List[Tuple[str, str]] = []
        for raw in tsv.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            # 只在首列是“字面上的表头 OG”时跳过
            if parts[0] == "OG":
                continue
            if len(parts) < 4:
                continue
            species = parts[1]
            cds_id  = parts[3]
            rows.append((species, cds_id))

        og2rows[og] = rows
        return og2rows[og]

    # Step 4: PSG cds 富集输入输出根目录（q 版 + p 版）
    cds_dir_name_q = psg_cfg.get("cds_inputs_dir") or "psg_cds_inputs"
    cds_dir_name_p = psg_cfg.get("cds_inputs_rawp_dir") or "psg_cds_inputs_rawp"
    psg_root_q = ensure_dir(out_dir / cds_dir_name_q)
    psg_root_p = ensure_dir(out_dir / cds_dir_name_p)

    fg_processed_q = 0
    fg_processed_p = 0

    for fg in sorted(bg_ogs.keys()):
        ogs_bg   = bg_ogs[fg]
        ogs_sig_q = sig_ogs_q.get(fg, set())
        ogs_sig_p = sig_ogs_p.get(fg, set())

        if fg not in fg_species:
            log.error(f"[PSG] FG={fg} 在 sets/ 目录中找不到对应 .list 文件，跳过该前景。")
            continue

        tips = fg_species[fg]
        if not tips:
            log.warning(f"[PSG] FG={fg} 的前景物种集合为空，继续构建但可能得到空的 cds 集合。")

        # 构建一个通用的背景 cds 集合（q 版/p 版共享）
        bg_cds_set: set = set()
        og2cds_cache: Dict[str, List[Tuple[str, str]]] = {}

        for og in sorted(ogs_bg):
            rows = get_og_rows(og)
            og2cds_cache[og] = rows
            if not rows:
                continue
            for species, cds_id in rows:
                if species not in tips:
                    continue
                bg_cds_set.add(cds_id)

        n_ogs_bg = len(ogs_bg)
        n_cds_bg = len(bg_cds_set)
        if n_ogs_bg == 0:
            log.warning(f"[PSG] FG={fg} 背景 OG 数为 0，跳过该前景。")
            continue

        fg_tips_str = ",".join(sorted(tips))
        head_meta = "\t".join([
            "tag","id_space","source","fg_label","fg_tips",
            "fdr_cutoff","n_ogs_bg","n_ogs_sig","n_cds_bg","n_cds_sig","note"
        ])

        # ------- 4A. q 值显著版：沿用原有发布结构 -------
        sig_cds_set_q: set = set()
        for og in sorted(ogs_sig_q):
            rows = og2cds_cache.get(og) or get_og_rows(og)
            for species, cds_id in rows:
                if species not in tips:
                    continue
                sig_cds_set_q.add(cds_id)

        n_ogs_sig_q = len(ogs_sig_q)
        n_cds_sig_q = len(sig_cds_set_q)

        tag_q = f"psg_{fg}"
        fg_dir_q = ensure_dir(psg_root_q / tag_q)

        bg_txt = "\n".join(sorted(bg_cds_set))
        if bg_txt:
            bg_txt += "\n"
        (fg_dir_q / "background.list").write_text(bg_txt, encoding="utf-8")

        test_txt_q = "\n".join(sorted(sig_cds_set_q))
        if test_txt_q:
            test_txt_q += "\n"
        (fg_dir_q / "test.list").write_text(test_txt_q, encoding="utf-8")

        note_q = "branch-site; terminals_mode; criterion=q<=fdr_alpha"
        row_meta_q = "\t".join([
            tag_q,
            "cds_id",
            "aphylo_codeml_branchsite",
            fg,
            fg_tips_str,
            f"{fdr_alpha:.6g}",
            str(n_ogs_bg),
            str(n_ogs_sig_q),
            str(n_cds_bg),
            str(n_cds_sig_q),
            note_q,
        ])
        (fg_dir_q / "meta.tsv").write_text(head_meta + "\n" + row_meta_q + "\n", encoding="utf-8")

        if n_cds_bg == 0:
            log.warning(f"[PSG] [Q] FG={fg} 无任何 cds 进入 PSG 背景（n_ogs_bg={n_ogs_bg}），background.list 为空。")
        if n_cds_sig_q == 0:
            log.warning(f"[PSG] [Q] FG={fg} 无任何显著 PSG cds（q <= {fdr_alpha}），test.list 为空。")

        log.info(
            f"[PSG] [Q] FG={fg} → tag={tag_q}, "
            f"OGs(bg/sig)={n_ogs_bg}/{n_ogs_sig_q}, CDS(bg/sig)={n_cds_bg}/{n_cds_sig_q}"
        )
        fg_processed_q += 1

        # ------- 4B. 原始 P 值显著版：新增发布结构 -------
        sig_cds_set_p: set = set()
        for og in sorted(ogs_sig_p):
            rows = og2cds_cache.get(og) or get_og_rows(og)
            for species, cds_id in rows:
                if species not in tips:
                    continue
                sig_cds_set_p.add(cds_id)

        n_ogs_sig_p = len(ogs_sig_p)
        n_cds_sig_p = len(sig_cds_set_p)

        tag_p = f"psg_rawp_{fg}"
        fg_dir_p = ensure_dir(psg_root_p / tag_p)

        # 背景依旧是同一份 bg_txt
        (fg_dir_p / "background.list").write_text(bg_txt, encoding="utf-8")

        test_txt_p = "\n".join(sorted(sig_cds_set_p))
        if test_txt_p:
            test_txt_p += "\n"
        (fg_dir_p / "test.list").write_text(test_txt_p, encoding="utf-8")

        note_p = "branch-site; terminals_mode; criterion=p<=p_alpha (raw P)"
        row_meta_p = "\t".join([
            tag_p,
            "cds_id",
            "aphylo_codeml_branchsite",
            fg,
            fg_tips_str,
            f"{p_alpha:.6g}",
            str(n_ogs_bg),
            str(n_ogs_sig_p),
            str(n_cds_bg),
            str(n_cds_sig_p),
            note_p,
        ])
        (fg_dir_p / "meta.tsv").write_text(head_meta + "\n" + row_meta_p + "\n", encoding="utf-8")

        if n_cds_sig_p == 0:
            log.info(f"[PSG] [P] FG={fg} 在原始 P 阈值 {p_alpha} 下无显著 PSG cds，test.list 为空。")
        else:
            log.info(
                f"[PSG] [P] FG={fg} → tag={tag_p}, "
                f"OGs(bg/sig)={n_ogs_bg}/{n_ogs_sig_p}, CDS(bg/sig)={n_cds_bg}/{n_cds_sig_p}"
            )
        fg_processed_p += 1

    if fg_processed_q == 0:
        log.warning("[PSG] 没有任何前景成功构建 PSG cds 富集输入（q 值版）。")
    else:
        write_done(psg_root_q / ".psg_cds_inputs.done")
        log.info(f"[PSG] PSG cds 富集输入构建完成（q 值版），前景数={fg_processed_q}，输出目录={psg_root_q}")

    if fg_processed_p == 0:
        log.warning("[PSG] 没有任何前景成功构建 PSG cds 富集输入（原始 P 值版）。")
    else:
        write_done(psg_root_p / ".psg_cds_inputs_rawp.done")
        log.info(f"[PSG] PSG cds 富集输入构建完成（原始 P 值版），前景数={fg_processed_p}，输出目录={psg_root_p}")

# ===================== 主流程 =====================

def main():
    cfg = load_config()
    paths = cfg["paths"]

    logs_dir = Path(paths["logs_dir"])
    LOG_FILE = logs_dir / "07_codeml_aggregate.log"
    log = get_logger("aphylo.07", LOG_FILE)
    banner(log, "APhylo 07 — codeml 汇总")

    raw_root = need_dir(Path(paths["codeml_dir"]) / "raw", "codeml 原始结果目录")
    out_dir  = ensure_dir(Path(paths["codeml_agg_dir"]))

    rows_genes: List[List[str]] = []
    rows_sites: List[List[str]] = []

    # 遍历 OG / FG 组合，解析 lnL，计算 LRT / P
    for og_dir in sorted(raw_root.glob("OG*")):
        og = og_dir.name
        for fg_dir in sorted(og_dir.glob("*")):
            fg = fg_dir.name
            alt = need_file(fg_dir / "alt" / "mlc.txt",  "ALT mlc")
            nul = need_file(fg_dir / "null" / "mlc.txt", "NULL mlc")

            alt_txt = alt.read_text(encoding="utf-8", errors="ignore")
            nul_txt = nul.read_text(encoding="utf-8", errors="ignore")

            la = last_lnl(alt_txt)
            ln = last_lnl(nul_txt)
            if la is None or ln is None:
                raise ValueError(f"[ERR] 无法解析 lnL：{alt} 或 {nul}")

            lrt = 2.0 * (la - ln)
            if lrt < 0:
                lrt = 0.0
            pval = chi2_mix_half_df1_sf(lrt)
            rows_genes.append([og, fg, f"{lrt:.6f}", f"{pval:.6g}"])

            # 解析 BEB 位点（只保留 post ≥ 0.95）
            for site, aa, post in parse_beb_lines(alt_txt, cutoff=0.95):
                rows_sites.append([og, fg, str(site), aa, f"{post:.3f}"])

    # 写基因层 PSG 表 + FDR
    head_g = "OG\tforeground\tLRT\tP\tQ\n"
    if rows_genes:
        pvals = [float(r[3]) for r in rows_genes]
        qvals = bh_fdr(pvals)
        lines_g = []
        for r, q in zip(rows_genes, qvals):
            lines_g.append("\t".join([r[0], r[1], r[2], r[3], f"{q:.6g}"]))
        (out_dir / "D_fdr_genes.tsv").write_text(head_g + "\n".join(lines_g) + "\n", encoding="utf-8")

        # 基于聚合结果构建 PSG cds 富集输入（可选，取决于 config.psg.enable）
        build_psg_cds_inputs(cfg, paths, out_dir, rows_genes, qvals, log)
    else:
        (out_dir / "D_fdr_genes.tsv").write_text(head_g, encoding="utf-8")

    # 写 BEB 位点表
    head_s = "OG\tforeground\tsite\taa\tpost\n"
    if rows_sites:
        (out_dir / "D_beb_sites.tsv").write_text(
            head_s + "\n".join("\t".join(x) for x in rows_sites) + "\n",
            encoding="utf-8"
        )
    else:
        (out_dir / "D_beb_sites.tsv").write_text(head_s, encoding="utf-8")

    write_done(out_dir / ".aggregate.done")
    log.info(f"[DONE] codeml 汇总完成：基因 {len(rows_genes)}；BEB 位点 {len(rows_sites)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)

