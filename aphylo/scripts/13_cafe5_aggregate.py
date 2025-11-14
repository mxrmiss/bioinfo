#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
13_cafe5_aggregate.py —— CAFE5 聚合（完全适配 12 教程增强正式版）
功能：
  ✓ primary_global 结果聚合
  ✓ large_global / large_global_e 完整合并
  ✓ primary_global_e（误差模型）并行输出
  ✓ primary_multi / primary_multi_e（多 λ）一起聚合
  ✓ 自动剔除高失败率家族（来自 flags/high_fail_ogs.list）
  ✓ 全量 / 无高失败率 两条结果线
  ✓ 误差模型 parallel 两条结果线
  ✓ branch λ 摘要（含 multi-λ 分 λ）
  ✓ BH-FDR
  ✓ 记录所有输入文件（inputs_used.tsv）
"""

from __future__ import annotations
import sys, io, logging, re
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml
import datetime

DEFAULT_CONFIG = "config.yaml"

# ============================================================
# 读取配置
# ============================================================

def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(fp: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(fp)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 配置文件不存在：{fp}")
    cfg = yaml.safe_load(p.read_text()) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

# ============================================================
# 日志
# ============================================================

def get_logger(name: str, logfile: Path):
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name)
    lg.setLevel(logging.INFO)
    lg.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s",
                            "%Y-%m-%d %H:%M:%S")

    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(logging.INFO)

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)

    lg.addHandler(fh)
    lg.addHandler(sh)

    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x):
            self.s.write(x); self.s.flush(); return len(x)

    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)

    return lg

def banner(log, text):
    bar = "=" * max(10, len(text)+2)
    log.info(bar)
    log.info(f" {text} ")
    log.info(bar)

# ============================================================
# 文件读写工具
# ============================================================

def read_table_auto(path: Path):
    txt = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not txt:
        return [], []
    delim = "\t" if "\t" in txt[0] else " "
    head = re.split(r"\s+" if delim==" " else "\t", txt[0].strip())
    body = []
    for ln in txt[1:]:
        if ln.strip():
            body.append(re.split(r"\s+" if delim==" " else "\t", ln.strip()))
    return head, body

def find_family_results_file(d: Path) -> Optional[Path]:
    """
    自动识别 family 结果文件：
      Gamma_family_results.txt
      Base_family_results.txt
      其他包含 family/result 的 txt/tsv
    """
    for name in ("Gamma_family_results.txt", "Base_family_results.txt"):
        p = d / name
        if p.is_file(): return p
    cand = [p for p in d.iterdir()
            if p.is_file()
            and p.suffix.lower() in (".txt",".tsv")
            and ("family" in p.name.lower() and "result" in p.name.lower())]
    return sorted(cand)[0] if cand else None

def take_last_block(runlog: Path) -> str:
    if not runlog.is_file(): return ""
    text = runlog.read_text(encoding="utf-8", errors="ignore")
    parts = re.split(r"^=+\s*MODEL=.*?ROUND=.*?=+\s*$",
                     text, flags=re.M)
    return parts[-1].strip() if parts else text

BRANCH_LAMBDA_PAT = re.compile(
    r"branch[:\s]+(\S+).+lambda\s*[:=]\s*([\d\.Ee+-]+)", re.I)

def parse_branch_lambda_from_block(block: str):
    out = []
    for ln in block.splitlines():
        m = BRANCH_LAMBDA_PAT.search(ln)
        if m:
            out.append((m.group(1), m.group(2)))
    return out

# ----------------------
# FDR（BH）
# ----------------------

def bh_fdr(vals: List[Optional[float]]):
    idx = [i for i,v in enumerate(vals) if (v is not None and v>=0)]
    m = len(idx)
    if m == 0:
        return [None]*len(vals)
    idx_sorted = sorted(idx, key=lambda i: vals[i])
    out = [None]*len(vals)

    min_q = 1.0
    for rank1, i in reversed(list(enumerate(idx_sorted, start=1))):
        q = vals[i]*m/rank1
        if q < min_q: min_q = q
        out[i] = min_q if min_q < 1 else 1
    return out

# ============================================================
# 主流程
# ============================================================

def main():
    cfg = load_config()
    paths = cfg["paths"]
    report_cfg = cfg.get("report", {})
    fdr_alpha = float(report_cfg.get("fdr_alpha", 0.05))

    # 日志
    logs_dir = Path(paths["logs_dir"])
    log = get_logger("aphylo.13",
                     logs_dir / "13_cafe5_aggregate.log")
    banner(log, "APhylo 13 — CAFE5 聚合（最新正式版）")

    cafe_dir = Path(paths["cafe_run_dir"])
    agg_dir = ensure_dir(Path(paths["cafe_agg_dir"]))
    done_flag = cafe_dir / ".cafe.done"
    if not done_flag.is_file():
        log.error(f"[ERR] 未发现 {done_flag} —— 12 未成功完成")
        sys.exit(2)

    models_root = cafe_dir / "models"
    if not models_root.is_dir():
        log.error("[ERR] 缺少 models/ 目录")
        sys.exit(2)

    # 容器
    fam_rows_all = []      # (model, subset, family, pval, pstr, source)
    fam_rows_all_e = []    # with_error
    branch_rows = []       # 主跑
    branch_rows_e = []     # with_error
    inputs_used = ["model\tsubdir\tfile\tmtime\tbytes\trole"]

    # ============================================================
    # 遍历所有模型（global/models/<model>）
    # ============================================================

    for mdir in sorted(models_root.glob("*")):
        if not mdir.is_dir():
            continue
        model = mdir.name
        log.info(f"[MODEL] 处理 {model}")

        # high-fail 家族表
        highfail = set()
        hf = mdir / "flags" / "high_fail_ogs.list"
        if hf.is_file():
            highfail = set(
                x.strip()
                for x in hf.read_text().splitlines()
                if x.strip() and x.strip()!="OG"
            )
            log.info(f"[{model}] high-fail OG 数量：{len(highfail)}")

        # --------------------------------------------
        # primary_global
        # --------------------------------------------
        pg = mdir / "primary_global"
        if not pg.is_dir():
            log.error(f"[ERR] {model} 缺少 primary_global/")
            sys.exit(2)

        fr_pg = find_family_results_file(pg)
        if not fr_pg:
            log.error(f"[ERR] primary_global 无 family_results 文件")
            sys.exit(2)

        hdr, body = read_table_auto(fr_pg)
        lc = [c.lower() for c in hdr]
        idx_f = next((i for i,c in enumerate(lc)
                      if c in ("family","orthogroup","og")), None)
        idx_p = next((i for i,c in enumerate(lc)
                      if c in ("p","pvalue")), None)
        if idx_f is None or idx_p is None:
            log.error(f"[ERR] 无法定位 family/p 列: {fr_pg}")
            sys.exit(2)

        for r in body:
            fam = r[idx_f] if idx_f < len(r) else ""
            pstr= r[idx_p] if idx_p < len(r) else ""
            try:
                pval = float(pstr) if pstr else None
            except:
                pval = None
            fam_rows_all.append((model, "primary", fam, pval, pstr, fr_pg.name))

        st = fr_pg.stat()
        inputs_used.append(f"{model}\tprimary_global\t{fr_pg.name}\t"
                           f"{datetime.datetime.fromtimestamp(st.st_mtime)}\t"
                           f"{st.st_size}\tfamily_results")

        # branch λ（primary）
        blk = take_last_block(pg / "run.log")
        for br, lam in parse_branch_lambda_from_block(blk):
            branch_rows.append((model, br, lam,
                                "primary_global.runlog"))

        # --------------------------------------------
        # large_global（新：必须读取）
        # --------------------------------------------
        lg = mdir / "large_global"
        if lg.is_dir():
            fr_lg = find_family_results_file(lg)
            if fr_lg:
                hdr2, body2 = read_table_auto(fr_lg)
                lc2 = [c.lower() for c in hdr2]
                jf = next((i for i,c in enumerate(lc2)
                          if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(lc2)
                          if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in body2:
                        fam = r[jf] if jf<len(r) else ""
                        pstr= r[jp] if jp<len(r) else ""
                        try:
                            pval=float(pstr) if pstr else None
                        except:
                            pval=None
                        fam_rows_all.append((model,
                                             "large",
                                             fam,
                                             pval,
                                             pstr,
                                             fr_lg.name))
                st = fr_lg.stat()
                inputs_used.append(f"{model}\tlarge_global\t{fr_lg.name}\t"
                                   f"{datetime.datetime.fromtimestamp(st.st_mtime)}\t"
                                   f"{st.st_size}\tfamily_results")

        # --------------------------------------------
        # primary_multi（多 λ）
        # --------------------------------------------
        pm = mdir / "primary_multi"
        if pm.is_dir():
            fr_pm = find_family_results_file(pm)
            if fr_pm:
                hdr3, body3 = read_table_auto(fr_pm)
                lc3 = [c.lower() for c in hdr3]
                jf = next((i for i,c in enumerate(lc3)
                          if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(lc3)
                          if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in body3:
                        fam = r[jf] if jf<len(r) else ""
                        pstr= r[jp] if jp<len(r) else ""
                        try:
                            pval=float(pstr) if pstr else None
                        except:
                            pval=None
                        fam_rows_all.append((model,
                                             "primary_multi",
                                             fam,
                                             pval,
                                             pstr,
                                             fr_pm.name))
                st = fr_pm.stat()
                inputs_used.append(f"{model}\tprimary_multi\t{fr_pm.name}\t"
                                   f"{datetime.datetime.fromtimestamp(st.st_mtime)}\t"
                                   f"{st.st_size}\tfamily_results")
            # branch λ（multi-λ）
            blk2 = take_last_block(pm / "run.log")
            for br, lam in parse_branch_lambda_from_block(blk2):
                branch_rows.append((model, br, lam,
                                    "primary_multi.runlog"))

        # --------------------------------------------
        # large_global_e / primary_global_e / primary_multi_e
        # --------------------------------------------
        # 误差模型并行集

        pge = mdir / "primary_global_e"
        if pge.is_dir():
            fr_pge = find_family_results_file(pge)
            if fr_pge:
                hdr4, body4 = read_table_auto(fr_pge)
                lc4 = [c.lower() for c in hdr4]
                jf = next((i for i,c in enumerate(lc4)
                          if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(lc4)
                          if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in body4:
                        fam = r[jf] if jf<len(r) else ""
                        pstr= r[jp] if jp<len(r) else ""
                        try:
                            pval=float(pstr) if pstr else None
                        except:
                            pval=None
                        fam_rows_all_e.append((model, "primary", fam,
                                               pval, pstr, fr_pge.name))
                st = fr_pge.stat()
                inputs_used.append(f"{model}\tprimary_global_e\t{fr_pge.name}\t"
                                   f"{datetime.datetime.fromtimestamp(st.st_mtime)}\t"
                                   f"{st.st_size}\tfamily_results_e")

            blk_e = take_last_block(pge / "run.log")
            for br, lam in parse_branch_lambda_from_block(blk_e):
                branch_rows_e.append((model, br, lam,
                                      "primary_global_e.runlog"))

        lge = mdir / "large_global_e"
        if lge.is_dir():
            fr_lge = find_family_results_file(lge)
            if fr_lge:
                hdr5, body5 = read_table_auto(fr_lge)
                lc5 = [c.lower() for c in hdr5]
                jf = next((i for i,c in enumerate(lc5)
                          if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(lc5)
                          if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in body5:
                        fam = r[jf] if jf<len(r) else ""
                        pstr= r[jp] if jp<len(r) else ""
                        try:
                            pval=float(pstr) if pstr else None
                        except:
                            pval=None
                        fam_rows_all_e.append((model, "large",
                                               fam, pval, pstr, fr_lge.name))
                st = fr_lge.stat()
                inputs_used.append(f"{model}\tlarge_global_e\t{fr_lge.name}\t"
                                   f"{datetime.datetime.fromtimestamp(st.st_mtime)}\t"
                                   f"{st.st_size}\tfamily_results_e")

        pme = mdir / "primary_multi_e"
        if pme.is_dir():
            fr_pme = find_family_results_file(pme)
            if fr_pme:
                hdr6, body6 = read_table_auto(fr_pme)
                lc6 = [c.lower() for c in hdr6]
                jf = next((i for i,c in enumerate(lc6)
                          if c in ("family","orthogroup","og")), None)
                jp = next((i for i,c in enumerate(lc6)
                          if c in ("p","pvalue")), None)
                if jf is not None and jp is not None:
                    for r in body6:
                        fam = r[jf] if jf<len(r) else ""
                        pstr= r[jp] if jp<len(r) else ""
                        try:
                            pval=float(pstr) if pstr else None
                        except:
                            pval=None
                        fam_rows_all_e.append((model, "primary_multi",
                                               fam, pval, pstr, fr_pme.name))
                st = fr_pme.stat()
                inputs_used.append(f"{model}\tprimary_multi_e\t{fr_pme.name}\t"
                                   f"{datetime.datetime.fromtimestamp(st.st_mtime)}\t"
                                   f"{st.st_size}\tfamily_results_e")

            blk3 = take_last_block(pme / "run.log")
            for br, lam in parse_branch_lambda_from_block(blk3):
                branch_rows_e.append((model, br, lam,
                                      "primary_multi_e.runlog"))

    # ============================================================
    # family 聚合（FDR）
    # ============================================================

    def write_family_outputs(rows, suffix, highfail_map):

        if not rows:
            return

        pvals   = [r[3] for r in rows]  # raw pval
        qvals   = bh_fdr(pvals)

        out_all = ["model\tsubset\tfamily\tpvalue\tqvalue\tsignif\tsource"]
        out_noh = ["model\tsubset\tfamily\tpvalue\tqvalue\tsignif\tsource"]

        # 构建每模型的 high-fail set
        hf_dict = {m:set() for m in set(r[0] for r in rows)}
        for m in highfail_map:
            hf_dict[m] = set(highfail_map[m])

        for (model, subset, fam, pval, pstr, src), qv in zip(rows, qvals):
            qstr = f"{qv:.3g}" if isinstance(qv, float) else ""
            signif = "yes" if isinstance(qv, float) and qv <= fdr_alpha else "no"
            line = "\t".join([model, subset, fam, pstr or "", qstr, signif, src])
            out_all.append(line)
            if fam not in hf_dict.get(model, set()):
                out_noh.append(line)

        (agg_dir/f"cafe_significant_families{suffix}.tsv").write_text(
            "\n".join(out_all)+"\n", encoding="utf-8"
        )
        (agg_dir/f"cafe_significant_families_no_highfail{suffix}.tsv").write_text(
            "\n".join(out_noh)+"*\n", encoding="utf-8"
        )

    # high-fail 映射
    highfail_map = {}
    for mdir in sorted(models_root.glob("*")):
        if not mdir.is_dir(): continue
        model = mdir.name
        hf = mdir/"flags"/"high_fail_ogs.list"
        if hf.is_file():
            highfail_map[model] = set(
                x.strip() for x in hf.read_text().splitlines()
                if x.strip() and x.strip()!="OG"
            )

    write_family_outputs(fam_rows_all, "", highfail_map)
    write_family_outputs(fam_rows_all_e, "_with_error", highfail_map)

    # ============================================================
    # 分支 λ 汇总
    # ============================================================

    def write_branch(rows, suffix):
        if not rows:
            return
        lines = ["model\tbranch\tlambda_estimate\tsource"]
        for model, br, lam, src in rows:
            lines.append("\t".join([model, br, lam, src]))
        (agg_dir/f"cafe_branch_summary{suffix}.tsv").write_text(
            "\n".join(lines)+"\n", encoding="utf-8"
        )

    write_branch(branch_rows, "")
    write_branch(branch_rows_e, "_with_error")

    # ============================================================
    # inputs_used
    # ============================================================

    (agg_dir/"inputs_used.tsv").write_text(
        "\n".join(inputs_used)+"\n", encoding="utf-8"
    )

    print(f"[OK] 聚合完成 → {agg_dir}")

# ============================================================

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e)+"\n")
        sys.exit(2)