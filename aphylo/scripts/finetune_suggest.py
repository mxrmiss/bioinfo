#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finetune_suggest.py —— 仅基于 ESS 的 finetune 建议（不读取 Acceptance）
放置：aphylo/scripts/
改动范围：仅“如何读取 config.yaml 的路径逻辑”与打印所用 config；其它逻辑未改。
"""

from pathlib import Path
import os, re
from datetime import datetime

# ===================== 顶部默认路径（可被 YAML 覆盖） =====================
PREFERRED_WORK_DIR = Path("mcmctree")
FALLBACK_WORK_DIRS = [Path("results/06_cafe/mcmctree")]
QC_DIRNAME = "qc_report"

# 其余参数（算法相关）保持不变
TARGET_ESS   = 250
TOP_K        = 3
MIN_FT       = 0.005
MAX_FT       = 2.0
ROUND_DIGITS = 3
EXP          = 0.5
MULT_MIN     = 0.80
MULT_MAX     = 1.25
AUTO_PATCH_CTL = False
FORCE_CTL_NAME = ""

# 记录本次实际使用的配置文件路径（仅用于提示）
_CFG_USED = None

# -------------------- 智能定位 config.yaml（只影响路径覆盖） --------------------
def _find_config_yaml() -> Path | None:
    # 1) 环境变量优先
    envp = os.environ.get("APHYLO_CONFIG", "").strip()
    if envp:
        p = Path(envp).expanduser().resolve()
        if p.exists():
            return p
    # 2) 常见候选：当前工作目录/相对路径/脚本相邻目录
    here = Path.cwd().resolve()
    scripts_dir = Path(__file__).resolve().parent
    candidates = [
        here / "config.yaml",
        here / "aphylo" / "config.yaml",
        scripts_dir.parent / "config.yaml",  # aphylo/scripts/../config.yaml
        here.parent / "config.yaml",
    ]
    seen = []
    for p in candidates:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.append(rp)
        if rp.exists():
            return rp
    return None

def _apply_paths_from_yaml():
    """仅覆盖路径相关变量；其它算法参数保持脚本默认值。"""
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

# -------------------- 其余函数保持原样（省略未变更的注释说明） --------------------
def slot_of_param(pname: str) -> int:
    p = pname.lower()
    if p.startswith("t_"):
        return 0
    if p.startswith(("mu", "sigma2")):
        return 1
    if p.startswith(("rgene", "rates")):
        return 2
    if p.startswith(("kappa", "alpha", "lambda")):
        return 4
    return 4

def read_text(p: Path) -> str:
    if not p.exists(): return ""
    for enc in ("utf-8","latin-1"):
        try: return p.read_text(encoding=enc, errors="ignore")
        except Exception: continue
    return ""

def pick_work_dir():
    for d in [PREFERRED_WORK_DIR] + FALLBACK_WORK_DIRS:
        if (d/"mcmc.txt").exists() or (d/"out.txt").exists() or (d/"in.BV").exists():
            return d.resolve()
    return PREFERRED_WORK_DIR.resolve()

def fmt(x: float) -> str:
    s = f"{x:.{ROUND_DIGITS}f}".rstrip("0").rstrip(".")
    if "." not in s: s += ".0"
    if len(s.split(".")[1]) == 1: s += "0"
    return s

def detect_ctl_path(work_dir: Path) -> Path:
    if FORCE_CTL_NAME:
        p = work_dir / FORCE_CTL_NAME
        if p.exists(): return p
    out_txt = work_dir / "out.txt"
    txt = read_text(out_txt)
    if txt:
        m = re.search(r"Reading options from\s+([^\s]+\.ctl)", txt)
        if m:
            cand = m.group(1)
            cand_path = Path(cand)
            if not cand_path.is_absolute():
                cand_path = work_dir / cand
            if cand_path.exists():
                return cand_path
    for name in ("mcmctree2.ctl", "mcmctree.ctl"):
        p = work_dir / name
        if p.exists(): return p
    ctls = sorted(work_dir.glob("*.ctl"), key=lambda p: p.stat().st_mtime if p.exists() else 0, reverse=True)
    if ctls: return ctls[0]
    return work_dir / "mcmctree.ctl"

def parse_finetune_from_ctl(ctl_text: str):
    m = re.search(r"(?im)^\s*finetune\s*=\s*([^\r\n]*)", ctl_text)
    if not m:
        return "1:", [0.10]*6
    body = m.group(1)
    for cut in ("*", "#", "!", ";"):
        if cut in body:
            body = body.split(cut, 1)[0]
    pm = re.match(r"^\s*(\d+\s*:)?\s*(.*)$", body)
    prefix = (pm.group(1) or "").strip() or "1:"
    tail = pm.group(2)
    nums = re.findall(r"\b[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\b", tail)
    arr = [float(v) for v in nums]
    if len(arr) >= 6: arr = arr[:6]
    if len(arr) == 0: arr = [0.10]*6
    while len(arr) < 6: arr.append(arr[-1])
    return prefix, arr

def replace_finetune_line(ctl_text: str, new_line: str) -> str:
    pat = re.compile(r"(?im)^\s*finetune\s*=\s*[^\r\n]*")
    if pat.search(ctl_text):
        return pat.sub(new_line, ctl_text, count=1)
    return ctl_text.rstrip() + "\n" + new_line + "\n"

def main():
    _apply_paths_from_yaml()  # ← 仅路径覆盖
    project_root = Path.cwd().resolve()
    work_dir = pick_work_dir()
    report_dir = work_dir / QC_DIRNAME
    report_dir.mkdir(parents=True, exist_ok=True)

    ess_tsv  = report_dir / "ess.tsv"
    ctl_path = detect_ctl_path(work_dir)

    MD_PATH  = report_dir / "finetune_suggestion.md"
    ONE_LINE = report_dir / "finetune_new_line.txt"

    ctl_txt = read_text(ctl_path)
    prefix, old_vals = parse_finetune_from_ctl(ctl_txt)

    tsv = read_text(ess_tsv).strip()
    new_vals = old_vals[:]
    reason = "保持原值（无 ESS 信息）"

    if tsv:
        rows = [ln.split("\t") for ln in tsv.splitlines()]
        if rows and "param" in rows[0] and "ESS" in rows[0]:
            ess_idx = rows[0].index("ESS")
            data = []
            for r in rows[1:]:
                if len(r) <= ess_idx: continue
                try:
                    e = float(r[ess_idx])
                except Exception:
                    continue
                if e > 0: data.append((r[0], e))
            data.sort(key=lambda x: x[1])
            if data and data[0][1] >= TARGET_ESS:
                reason = f"all ESS ≥ {TARGET_ESS}, keep finetune"
            else:
                worst = data[:max(1, min(TOP_K, len(data)))]
                applied = []
                for pname, ess in worst:
                    slot = slot_of_param(pname)
                    mult = (TARGET_ESS / ess) ** EXP
                    if mult < MULT_MIN: mult = MULT_MIN
                    if mult > MULT_MAX: mult = MULT_MAX
                    new_vals[slot] = max(MIN_FT, min(MAX_FT, new_vals[slot] * mult))
                    applied.append((pname, ess, slot, mult))
                if applied:
                    reason = "ESS 驱动微调：" + "; ".join(
                        [f"{p} ESS={e:.1f} → 槽{idx} ×{m:.2f}" for p, e, idx, m in applied]
                    )

    new_line = f"finetune = {prefix} " + " ".join(fmt(v) for v in new_vals[:6])

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    rel_work = str(work_dir.relative_to(project_root)) if str(work_dir).startswith(str(project_root)) else str(work_dir)
    rel_ctl  = str(ctl_path.relative_to(project_root)) if str(ctl_path).startswith(str(project_root)) else str(ctl_path)

    lines = []
    lines.append("# Finetune Recommendation (ESS-only)\n")
    lines.append(f"- Generated at: {now}")
    if _CFG_USED: lines.append(f"- Using config: `{_CFG_USED}`")
    lines.append(f"- Work dir: `{rel_work}`")
    lines.append(f"- Ctl file: `{rel_ctl}`\n")
    lines.append("## Current finetune\n")
    lines.append("`finetune = " + prefix + " " + " ".join(fmt(v) for v in old_vals) + "`\n")
    lines.append("## Recommended finetune line\n")
    lines.append("```")
    lines.append(new_line)
    lines.append("```")
    lines.append(f"\n> 依据：{reason}；规则：对 ESS 最差 Top-{TOP_K} 参数按 mult = clamp((TARGET/ESS)^{EXP}, {MULT_MIN}..{MULT_MAX}) 温和缩放到 [{MIN_FT}, {MAX_FT}]。\n")

    MD_PATH.write_text("\n".join(lines), encoding="utf-8")
    ONE_LINE.write_text(new_line + "\n", encoding="utf-8")
    print(f"[OK] Using config: {_CFG_USED or 'N/A'}")
    print(f"[OK] Wrote: {MD_PATH}")
    print(f"[OK] Wrote: {ONE_LINE}")

    if AUTO_PATCH_CTL and ctl_path.exists():
        bak = ctl_path.with_suffix(ctl_path.suffix + ".bak")
        bak.write_text(ctl_txt, encoding="utf-8")
        ctl_path.write_text(replace_finetune_line(ctl_txt, new_line), encoding="utf-8")
        print(f"[OK] Patched: {ctl_path} (backup: {bak.name})")

if __name__ == "__main__":
    main()

