#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mcmctree_postcheck.py —— MCMCTree 产物体检 + 自动联动 ESS 报告与 finetune 建议
放置：aphylo/scripts/
改动范围：仅“如何读取 config.yaml 的路径逻辑”与打印所用 config；其它逻辑未改。
"""

from pathlib import Path
import os, re, csv, subprocess, shutil, sys
from datetime import datetime

# ===================== 顶部默认路径（可被 YAML 覆盖） =====================
PREFERRED_WORK_DIR = Path("mcmctree")
FALLBACK_WORK_DIRS = [Path("results/06_cafe/mcmctree")]
QC_DIRNAME = "qc_report"

TARGET_ESS = 250
RUN_FINETUNE_SUGGEST = True

SCRIPTS_DIR = Path(__file__).resolve().parent
ESS_R_NAME = "ess_report.R"
FINETUNE_PY_NAME = "finetune_suggest.py"

_CFG_USED = None  # 仅用于提示

# -------------------- 智能定位 config.yaml --------------------
def _find_config_yaml() -> Path | None:
    envp = os.environ.get("APHYLO_CONFIG", "").strip()
    if envp:
        p = Path(envp).expanduser().resolve()
        if p.exists():
            return p
    here = Path.cwd().resolve()
    scripts_dir = Path(__file__).resolve().parent
    candidates = [
        here / "config.yaml",
        here / "aphylo" / "config.yaml",
        scripts_dir.parent / "config.yaml",  # aphylo/scripts/../config.yaml
        here.parent / "config.yaml",
    ]
    seen = set()
    for p in candidates:
        rp = p.resolve()
        if rp in seen: 
            continue
        seen.add(rp)
        if rp.exists():
            return rp
    return None

def _apply_paths_from_yaml():
    global _CFG_USED
    try:
        import yaml
    except Exception:
        return
    cfg_path = _find_config_yaml()
    if not cfg_path:
        return
    try:
        cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    except Exception:
        return
    m = cfg.get("mcmctree", {}) if isinstance(cfg, dict) else {}
    wd = m.get("work_dir")
    if isinstance(wd, str) and wd.strip():
        globals()["PREFERRED_WORK_DIR"] = Path(wd.strip())
    fbs = m.get("fallback_work_dirs")
    if isinstance(fbs, (list, tuple)):
        globals()["FALLBACK_WORK_DIRS"] = [Path(str(x)) for x in fbs if str(x).strip()]
    qn = m.get("qc_dirname")
    if isinstance(qn, str) and qn.strip():
        globals()["QC_DIRNAME"] = qn.strip()
    _CFG_USED = str(cfg_path)

# ===================== 其余函数保持原样 =====================
def pick_work_dir():
    cands = [PREFERRED_WORK_DIR] + FALLBACK_WORK_DIRS
    for d in cands:
        if (d / "out.txt").exists() or (d / "mcmc.txt").exists() or (d / "in.BV").exists():
            return d.resolve()
    return PREFERRED_WORK_DIR.resolve()

def read_text(p: Path) -> str:
    if not p.exists(): return ""
    for enc in ("utf-8","latin-1"):
        try: return p.read_text(encoding=enc, errors="ignore")
        except Exception: continue
    return ""

def human_bytes(n: int) -> str:
    u=["B","KB","MB","GB","TB"]; i=0; x=float(n)
    while x>=1024 and i<len(u)-1: x/=1024; i+=1
    return f"{x:.2f} {u[i]}"

def safe_float(x: str):
    try: return float(x)
    except Exception: return None

def parse_node_ages_from_out(out_path: Path):
    items=[]; txt = read_text(out_path)
    if not txt: return items
    regexes = [
        re.compile(r"^\s*(t[_\-\w\d]*)\s+([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*\(\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*,\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*\)", re.M),
        re.compile(r"^\s*(node\s*\d+).*?mean\s*=\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?).*?95%HPD\s*=\s*\(\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*,\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*\)", re.I|re.M),
        re.compile(r"^\s*([A-Za-z]\S*)\s+([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*\(\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*,\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*\)\s*$", re.M),
    ]
    blocks=[]
    for m in re.finditer(r"(node ages|HPD|Highest Posterior Density|95%)[\s\S]{0,4000}", txt, flags=re.I):
        blocks.append(m.group(0))
    search_texts = blocks if blocks else [txt]
    seen=set()
    for seg in search_texts:
        for rgx in regexes:
            for m in rgx.finditer(seg):
                label=m.group(1).strip()
                mean=safe_float(m.group(2)); l95=safe_float(m.group(3)); u95=safe_float(m.group(4))
                if None in (mean,l95,u95): continue
                if l95>u95: l95,u95=u95,l95
                key=(label,mean,l95,u95)
                if key in seen: continue
                seen.add(key)
                items.append({"source":"out.txt","label":label,"mean":mean,"l95":l95,"u95":u95})
    return items

def extract_tree_string(fig_txt: str):
    m = re.search(r"tree\s+\S+\s*=\s*(?:\[&R\]\s*)?(.+?);", fig_txt, flags=re.I|re.S)
    return m.group(1) if m else ""

def estimate_taxa_from_figtree(fig_path: Path):
    txt = read_text(fig_path)
    if not txt: return None
    s = extract_tree_string(txt)
    if not s: return None
    s = re.sub(r"\[&[^\]]*\]", "", s)
    tips = re.findall(r"([^\(\),:\s\[\]{}]+)\s*:", s)
    return len(tips) if tips else None

def root_hpd_from_figtree(fig_path: Path):
    txt = read_text(fig_path)
    if not txt: return None
    s = extract_tree_string(txt)
    if not s: return None
    hpds=[]
    for h in re.finditer(r"95%HPD\s*=\s*\{\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*,\s*([\-+]?\d+(?:\.\d+)?(?:[eE][\-+]?\d+)?)\s*\}", s):
        l=safe_float(h.group(1)); u=safe_float(h.group(2))
        if None in (l,u): continue
        if l>u: l,u=u,l
        hpds.append((l,u))
    if not hpds: return None
    l,u=max(hpds, key=lambda x: x[1])
    return {"source":"FigTree.tre","label":"root~approx","mean":(l+u)/2.0,"l95":l,"u95":u}

def mcmc_quick_stats(p: Path):
    if not p.exists(): return {"exists": False}
    lines=0; has_header=False
    with p.open("r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f,1):
            s=line.strip()
            if i<=5 and re.search(r"(parameter|sample|iteration|Gen)", s, flags=re.I):
                has_header=True
            if s: lines+=1
    return {"exists": True, "lines": lines, "has_header": has_header}

def main():
    _apply_paths_from_yaml()  # ← 仅路径覆盖
    project_root = Path.cwd().resolve()
    work_dir = pick_work_dir()
    report_dir = work_dir / QC_DIRNAME
    report_dir.mkdir(parents=True, exist_ok=True)

    out_txt = work_dir / "out.txt"
    mcmc_txt = work_dir / "mcmc.txt"
    fig_tree = work_dir / "FigTree.tre"
    in_bv    = work_dir / "in.BV"

    node_tsv   = report_dir / "node_ages.tsv"
    summary_md = report_dir / "summary.md"

    node_items = parse_node_ages_from_out(out_txt)
    with node_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source","label","mean","l95","u95"])
        for it in node_items:
            w.writerow([it["source"], it["label"], it["mean"], it["l95"], it["u95"]])

    # 调用 R 脚本（在项目根目录执行）
    ess_tsv     = report_dir / "ess.tsv"
    ess_summary = report_dir / "ess_summary.md"
    ess_rec     = report_dir / "ess_recommendation.txt"
    ess_r_src   = SCRIPTS_DIR / ESS_R_NAME

    ess_ok = False
    rscript = shutil.which("Rscript")
    if rscript and ess_r_src.exists():
        try:
            subprocess.run([rscript, str(ess_r_src)], check=True, cwd=str(project_root))
            ess_ok = ess_tsv.exists()
        except Exception:
            ess_ok = False

    # 可选：调用 finetune 建议脚本
    if RUN_FINETUNE_SUGGEST:
        ft_py = SCRIPTS_DIR / FINETUNE_PY_NAME
        if ft_py.exists():
            try:
                subprocess.run([sys.executable, str(ft_py)], check=True, cwd=str(project_root))
            except Exception:
                pass

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    files = []
    for p in [in_bv, out_txt, fig_tree, mcmc_txt]:
        d = {"path": str(p.relative_to(project_root)) if p.exists() else str(p),
             "exists": p.exists()}
        if p.exists():
            try: d["size"]=human_bytes(p.stat().st_size)
            except Exception: d["size"]="NA"
        files.append(d)

    lines=[]
    lines.append("# MCMCTree QC Summary\n")
    lines.append(f"- Generated at: {now}")
    lines.append(f"- Using config: `{_CFG_USED or 'N/A'}`\n")
    rel_work = str(work_dir.relative_to(project_root)) if str(work_dir).startswith(str(project_root)) else str(work_dir)
    lines.append(f"- Work dir: `{rel_work}`\n")

    lines.append("## Artifacts\n")
    for s in files:
        mark = "✅" if s["exists"] else "❌"
        size = f" ({s.get('size','')})" if s.get("size") else ""
        lines.append(f"- {mark} `{s['path']}`{size}")
    lines.append(f"- ✅ `{node_tsv.relative_to(project_root) if node_tsv.exists() else node_tsv}`")
    if ess_ok:
        lines.append(f"- ✅ `{ess_tsv.relative_to(project_root) if ess_tsv.exists() else ess_tsv}`, `{ess_summary.relative_to(project_root) if ess_summary.exists() else ess_summary}`, `{ess_rec.relative_to(project_root) if ess_rec.exists() else ess_rec}`")
    else:
        lines.append(f"- ⚠ 未生成 ESS 报告（缺 Rscript 或 `{ESS_R_NAME}`），可手动运行：`Rscript {ESS_R_NAME}`")
    lines.append("")

    lines.append("## Node Ages & 95% HPD\n")
    if node_items:
        lines.append(f"- Parsed from `out.txt`: **{len(node_items)}** nodes")
        root_like = max(node_items, key=lambda d: (d["u95"]-d["l95"]))
        lines.append(f"- Approx. root by widest HPD (out.txt): mean={root_like['mean']:.3f}, 95%HPD=({root_like['l95']:.3f}, {root_like['u95']:.3f})")
    else:
        lines.append("- Parsed from `out.txt`: 未抓取到节点年龄表")
    lines.append("")

    lines.append("## MCMC Trace\n")
    if mcmc_txt.exists():
        lines.append(f"- Lines (non-empty): {mcmc_txt.stat().st_size} bytes")
    else:
        lines.append("- 未检测到 `mcmc.txt`")
    lines.append("")

    if ess_tsv.exists() and ess_rec.exists():
        try:
            with ess_tsv.open("r", encoding="utf-8") as f:
                rows = [line.strip().split("\t") for line in f if line.strip()]
            ess_vals = []
            for r in rows[1:]:
                try: ess_vals.append(float(r[1]))
                except: pass
            min_ess = min(ess_vals) if ess_vals else None
            med_ess = sorted(ess_vals)[len(ess_vals)//2] if ess_vals else None
            rec = (report_dir / "ess_recommendation.txt").read_text(encoding="utf-8").strip()
            lines.append("## ESS Snapshot & Recommendation")
            lines.append(f"- Min ESS: **{(f'{min_ess:.2f}' if isinstance(min_ess,(int,float)) else 'NA')}**; Median ESS: **{(f'{med_ess:.2f}' if isinstance(med_ess,(int,float)) else 'NA')}**; Target: **{TARGET_ESS}**")
            lines.append(f"- Suggested `nsample` multiplier: **×{rec}**\n")
        except Exception:
            lines.append("## ESS Snapshot & Recommendation")
            lines.append("- 解析 ESS 结果时出错；请查看 `ess_summary.md` 与 `ess.tsv`\n")
    else:
        lines.append("## ESS Snapshot & Recommendation")
        lines.append("- ESS report unavailable\n")

    lines.append("## Checklist\n")
    lines.append("- [ ] 更换对齐/模型/拓扑后需重建 `in.BV` 再采样")
    lines.append("- [ ] 若 Min ESS < Target：优先提高 `nsample`；必要时微调 `finetune`")
    lines.append("- [ ] 给 CAFE5：按需要导出仅含时间分支长度的树\n")

    summary_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"[OK] Using config: {_CFG_USED or 'N/A'}")
    print(f"[OK] Report written: {summary_md}")

if __name__ == "__main__":
    main()

