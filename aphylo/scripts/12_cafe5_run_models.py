#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— CAFE5 批量模型（教程增强版·两阶段/误差/流式输出）
要点：
  - 所有参数从 config.yaml 读取（不接命令行）。
  - family.tsv 首列固定为 Desc，OG 列为 Orthogroup（大小写均可）。
  - primary 阶段：支持多轮“极端 OG 自动修正”（最多 N 轮）。
  - large 阶段：只跑一次（按皇上新指令），可固定 primary 的 λ 复算；不做极端 OG 剔除。
  - 可选误差模型（estimate/use/off），可选择只应用 primary 或 primary+large。
  - 产物目录简洁：results/06_cafe/models/<model>/{primary_global|large}/
  - 实时（流式）屏幕输出 + 写入 run.log；不自动删除旧 models（除非将 config 的 archive 开关打开）。
兼容：Python 3.8+；CAFE5 可执行名从 config.binaries.cafe5 读取，默认 'cafe5'
"""

from __future__ import annotations
import os, sys, re, io, shutil, time, datetime, subprocess
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
import yaml

# ------------------------- 基础工具 -------------------------

DEFAULT_CONFIG = "config.yaml"

def _expand_publish_placeholders(obj, pub: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", pub)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, pub) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, pub) for k,v in obj.items()}
    return obj

def load_config(path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def now_tag() -> str:
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

def tee_streamed(cmd: List[str], cwd: Path, logfile: Path) -> int:
    """
    实时把子进程 stdout/stderr 同步到屏幕与 logfile（UTF-8）
    返回进程退出码
    """
    ensure_dir(logfile.parent)
    with logfile.open("w", encoding="utf-8") as lf:
        lf.write(f"# CMD: {' '.join(cmd)}\n# CWD: {cwd}\n# START: {datetime.datetime.now()}\n\n")
        lf.flush()
        proc = subprocess.Popen(
            cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            bufsize=1, universal_newlines=True, encoding="utf-8", errors="replace"
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            # 写屏 & 写盘（流式）
            sys.stdout.write(line); sys.stdout.flush()
            lf.write(line)
        proc.wait()
        lf.write(f"\n# END: {datetime.datetime.now()}\n# RET: {proc.returncode}\n")
        lf.flush()
        return proc.returncode

def copy_file(src: Path, dst: Path):
    ensure_dir(dst.parent); shutil.copy2(src, dst)

# ------------------------- CAFE5 辅助 -------------------------

OG_RE = re.compile(r'\b(OG\d{7,})\b', re.I)  # 匹配如 OG0000123

def read_family_as_rows(path: Path) -> Tuple[List[str], List[List[str]]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines: return [], []
    header = re.split(r"\t", lines[0].rstrip("\n"))
    body = [re.split(r"\t", ln.rstrip("\n")) for ln in lines[1:] if ln.strip()]
    return header, body

def write_family_from_rows(path: Path, header: List[str], body: List[List[str]]):
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in body:
            f.write("\t".join(map(str, r)) + "\n")

def find_idx(header: List[str], oks: List[str]) -> int:
    low = [h.strip().lower() for h in header]
    for i,h in enumerate(low):
        if h in oks: return i
    return -1

def split_primary_large(family_tsv: Path, out_primary: Path, out_large: Path, copy_threshold: int) -> Tuple[int,int]:
    """
    基于任意物种计数 >= copy_threshold 的行划为 large，否则划为 primary。
    family.tsv 首列 Desc，第二列为 Orthogroup（大小写兼容）。
    """
    head, body = read_family_as_rows(family_tsv)
    if not head: raise RuntimeError(f"[ERR] 空的 family：{family_tsv}")
    idx_og = find_idx(head, ["orthogroup","og","family"])
    if idx_og < 0:
        raise RuntimeError("[ERR] family.tsv 缺少 Orthogroup/OG 列")
    # 非前两列的其余列应为物种计数列（整数）
    prim, larg = [], []
    for row in body:
        counts = []
        for x in row[2:]:
            try: counts.append(int(x))
            except: counts.append(0)
        if any(c >= copy_threshold for c in counts):
            larg.append(row)
        else:
            prim.append(row)
    write_family_from_rows(out_primary, head, prim)
    write_family_from_rows(out_large, head, larg)
    return len(prim), len(larg)

def extract_lambda_from_results(results_txt: Path) -> Optional[float]:
    if not results_txt.is_file(): return None
    txt = results_txt.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"Lambda\s*:\s*([0-9.eE+-]+)", txt)
    if m:
        try: return float(m.group(1))
        except: return None
    return None

def parse_extreme_from_runlog(logfile: Path) -> List[str]:
    """
    从 run.log 的最后一个 “Families with largest size differentials:” 块中提取 OG 列表
    兼容多轮输出；不使用 \R，逐行扫更稳健。
    """
    if not logfile.is_file(): return []
    lines = logfile.read_text(encoding="utf-8", errors="ignore").splitlines()
    locs = [i for i,ln in enumerate(lines) if "Families with largest size differentials:" in ln]
    if not locs: return []
    start = locs[-1] + 1  # 最后一个块的下一行开始
    ogs, i = [], start
    while i < len(lines):
        ln = lines[i].strip()
        if not ln: break
        m = OG_RE.search(ln)
        if m: ogs.append(m.group(1))
        else:
            # 到了非 OG 行（例如提示语/下一段），收集结束
            break
        i += 1
    return ogs

def filter_family_remove_ogs(src_tsv: Path, dst_tsv: Path, bad_ogs: List[str]) -> int:
    """在 family.tsv（首列 Desc，第2列 Orthogroup）中删除 bad_ogs 对应行。返回删除数量。"""
    if not bad_ogs: 
        copy_file(src_tsv, dst_tsv); 
        return 0
    head, body = read_family_as_rows(src_tsv)
    idx_og = find_idx(head, ["orthogroup","og","family"])
    if idx_og < 0:
        raise RuntimeError("[ERR] 无 Orthogroup 列，无法删除 OG")
    bad = set(x.upper() for x in bad_ogs)
    kept = [r for r in body if (idx_og >= len(r) or r[idx_og].upper() not in bad)]
    removed = len(body) - len(kept)
    write_family_from_rows(dst_tsv, head, kept)
    return removed

def parse_high_fail_from_runlog(logfile: Path, threshold_rate: float = 0.20) -> List[str]:
    """
    抓取 “The following families had failure rates >20% of the time:” 之后的 OG 列表（直到空行/下一段）。
    """
    if not logfile.is_file(): return []
    lines = logfile.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = -1
    key = "had failure rates >"
    for i,ln in enumerate(lines):
        if key in ln and "families" in ln:
            start = i + 1; break
    if start < 0: return []
    ogs = []
    i = start
    while i < len(lines):
        ln = lines[i].strip()
        if not ln: break
        m = OG_RE.search(ln)
        if m: ogs.append(m.group(1))
        else: break
        i += 1
    return ogs

# ------------------------- 主逻辑 -------------------------

def main():
    cfg = load_config()
    paths = cfg.get("paths", {})
    cafes = cfg.get("cafes", {})
    bins  = cfg.get("binaries", {})

    cafe_bin = str(bins.get("cafe5", "cafe5"))
    threads  = int(cafes.get("threads", 30))
    k_cycles = int(cafes.get("k_cycles", 3))
    p_alpha  = float(cafes.get("p_alpha", 0.05))
    models   = cafes.get("models", ["global"])

    # 两阶段：large
    two_stage = cafes.get("two_stage_large", {}) or {}
    large_enable = bool(two_stage.get("enable", True))
    large_copy_thr = int(two_stage.get("copy_threshold", 100))

    # 自动修正轮数（仅 primary）
    max_autofix = int(cafes.get("max_autofix_rounds", 3))

    # 误差模型
    em = cafes.get("error_model", {}) or {}
    em_enable  = bool(em.get("enable", False))
    em_mode    = str(em.get("mode", "off")).lower()   # off|estimate|use
    em_file    = em.get("file", None)
    em_apply   = str(em.get("apply_to", "primary")).lower()  # primary|both

    # 归档/清理
    archive = cafes.get("archive_on_start", False)

    # 路径
    proj_root = Path(".").resolve()
    cafe_root = ensure_dir(Path(paths.get("cafe_run_dir", "results/06_cafe")))
    input_dir = ensure_dir(cafe_root / "input")
    models_root = ensure_dir(cafe_root / "models")
    logs_dir = ensure_dir(Path(paths.get("logs_dir", "logs")))

    family_tsv = Path(paths.get("cafe_family", str(input_dir / "family.tsv")))
    utree_nwk  = Path(paths.get("cafe_tree",   str(input_dir / "utree_for_cafe.nwk")))

    print("="*60)
    print(" APhylo 12 — CAFE5（教程增强版·两阶段/误差/流式输出）")
    print("="*60)
    print(f"[IN] family: {family_tsv}")
    print(f"[IN] tree  : {utree_nwk}")
    print(f"[CFG] threads={threads}, k={k_cycles}, P={p_alpha}, models={models}")
    if large_enable:
        print(f"[CFG] two-stage large: enable, copy_threshold={large_copy_thr}")
    else:
        print(f"[CFG] two-stage large: disabled")
    if em_enable:
        print(f"[CFG] error model: mode={em_mode}, apply_to={em_apply}, file={em_file or 'N/A'}")
    else:
        print(f"[CFG] error model: off")

    # 可选归档旧模型
    if archive and models_root.exists():
        tag = now_tag()
        dst = cafe_root / f"models_archive_{tag}"
        shutil.move(str(models_root), str(dst))
        print(f"[ARCHIVE] 旧 models -> {dst}")
        ensure_dir(models_root)

    # 拆分 primary/large
    primary_tsv = input_dir / "family.primary.tsv"
    large_tsv   = input_dir / "family.large.tsv"
    prim_n, larg_n = split_primary_large(family_tsv, primary_tsv, large_tsv, large_copy_thr)
    print(f"[SPLIT] primary={prim_n} 行, large={larg_n} 行（阈值 {large_copy_thr}）")

    # 误差模型（estimate）文件落点
    estim_err_file: Optional[Path] = None

    # 遍历模型（目前多半只有 'global'）
    for model in models:
        print("\n" + "="*60)
        print(f"[MODEL={model} PRIMARY-GLOBAL] —— CAFE5 入口开始 ——")
        print("="*60)

        mroot = ensure_dir(models_root / model)
        flags_dir = ensure_dir(mroot / "flags")
        sent_dir  = ensure_dir(mroot / "sentinels")

        # ---------- PRIMARY ----------
        out_primary = ensure_dir(mroot / "primary_global")
        runlog_p = out_primary / "run.log"
        # 清空 run.log（轮转）
        if runlog_p.exists():
            bak = out_primary / f"run_{now_tag()}.log"
            shutil.move(str(runlog_p), str(bak))
        runlog_p.touch()

        # primary 输入文件（可能经历 autofix 迭代）
        cur_family = primary_tsv
        lambda_fixed: Optional[float] = None

        for rnd in range(1, max_autofix + 1):
            print(f"[PRIMARY-GLOBAL ROUND={rnd}] 入口")
            cmd = [
                cafe_bin,
                "-i", str(cur_family),
                "-t", str(utree_nwk),
                "-c", str(threads),
                "-k", str(k_cycles),
                "-P", str(p_alpha),
                "-o", str(out_primary)
            ]
            # 误差模型（仅当 enable & apply_to 包含 primary）
            em_use_this_round = (em_enable and em_mode in ("use","estimate") and em_apply in ("primary","both"))
            err_model_path_for_run: Optional[Path] = None

            if em_use_this_round:
                if em_mode == "use":
                    if not em_file or not Path(em_file).is_file():
                        print("[WARN] 指定的误差模型文件不存在，跳过误差应用。")
                    else:
                        cmd.append(f"-e{em_file}")
                elif em_mode == "estimate":
                    # 先独立估计误差模型（只需一次）
                    if estim_err_file is None:
                        print("[EM] 估计误差模型（estimate）")
                        em_out = ensure_dir(mroot / "primary_global" / "emodel_tmp")
                        em_cmd = [cafe_bin, "-i", str(cur_family), "-t", str(utree_nwk), "-e", "-o", str(em_out)]
                        _ = tee_streamed(em_cmd, cwd=proj_root, logfile=(em_out / "em_estimate.log"))
                        # CAFE5 命名可能为 Base_error_model.txt 或 Gamma_error_model.txt，做个兜底
                        cand = [em_out/"Base_error_model.txt", em_out/"Gamma_error_model.txt", em_out/"error_model.txt"]
                        for p in cand:
                            if p.is_file():
                                estim_err_file = p
                                break
                        if estim_err_file:
                            print(f"[EM] 估计完成：{estim_err_file.name}")
                        else:
                            print("[WARN] 未找到误差模型输出，继续无误差。")
                    if estim_err_file and estim_err_file.is_file():
                        cmd.append(f"-e{str(estim_err_file)}")

            # 跑 CAFE5（流式）
            ret = tee_streamed(cmd, cwd=proj_root, logfile=runlog_p)
            print(f"[PRIMARY-GLOBAL ROUND={rnd}] ret={ret}")

            # 抓 λ（用于 large 固定）
            if lambda_fixed is None:
                # Gamma_results.txt / Base_results.txt 皆可
                for nm in ("Gamma_results.txt", "Base_results.txt"):
                    p = out_primary / nm
                    lam = extract_lambda_from_results(p)
                    if lam is not None:
                        lambda_fixed = lam
                        print(f"[PRIMARY] λ = {lambda_fixed}")
                        break

            # 解析极端 OG（仅 primary 做自动修正）
            extreme = parse_extreme_from_runlog(runlog_p)
            if not extreme:
                print("[PRIMARY] 未解析到极端 OG 或已稳定，结束自动修正。")
                break
            # 生成下一轮 family
            nxt = input_dir / f"family.autofix{rnd}.tsv"
            removed = filter_family_remove_ogs(cur_family, nxt, extreme)
            # 记录被移除的 OG
            (mroot / f"autofix_removed_round{rnd}.tsv").write_text("\n".join(extreme) + "\n", encoding="utf-8")
            print(f"[PRIMARY] 解析到 {len(extreme)} 个极端 OG，实际从 family 移除 {removed} 行。")
            if removed == 0:
                print("[PRIMARY] 解析与 family 不匹配，停止自动修正。")
                break
            cur_family = nxt  # 下一轮用修正后的 family

        # 记高失败率家族
        high_fail = parse_high_fail_from_runlog(runlog_p)
        if high_fail:
            (flags_dir / "high_fail_ogs.list").write_text("\n".join(high_fail) + "\n", encoding="utf-8")
            print(f"[PRIMARY] 高失败率家族记录：{len(high_fail)}")

        # ---------- LARGE（只跑一次，不剔除极端 OG） ----------
        if large_enable and large_tsv.is_file() and larg_n > 0:
            print("\n" + "="*60)
            print(f"[MODEL={model} LARGE] —— CAFE5 入口开始 ——")
            print("="*60)

            out_large = ensure_dir(mroot / "large")
            runlog_l = out_large / "run.log"
            if runlog_l.exists():
                bak = out_large / f"run_{now_tag()}.log"
                shutil.move(str(runlog_l), str(bak))
            runlog_l.touch()

            cmdL = [
                cafe_bin,
                "-i", str(large_tsv),
                "-t", str(utree_nwk),
                "-c", str(threads),
                "-k", str(k_cycles),
                "-P", str(p_alpha),
                "-o", str(out_large)
            ]
            # 固定 λ（若 primary 成功获得）
            if lambda_fixed is not None:
                cmdL += ["-l", f"{lambda_fixed}"]

            # 误差模型应用：若 apply_to = both
            if em_enable and em_mode in ("use","estimate") and em_apply == "both":
                if em_mode == "use":
                    if em_file and Path(em_file).is_file():
                        cmdL.append(f"-e{em_file}")
                elif em_mode == "estimate":
                    if estim_err_file and estim_err_file.is_file():
                        cmdL.append(f"-e{str(estim_err_file)}")

            retL = tee_streamed(cmdL, cwd=proj_root, logfile=runlog_l)
            print(f"[LARGE] ret={retL}")

            # 记录 large 阶段的高失败率家族（额外）
            high_fail_L = parse_high_fail_from_runlog(runlog_l)
            if high_fail_L:
                (flags_dir / "high_fail_ogs_large.list").write_text("\n".join(high_fail_L) + "\n", encoding="utf-8")
                print(f"[LARGE] 高失败率家族记录：{len(high_fail_L)}")
        else:
            print("[LARGE] 跳过：未启用或无 large 行。")

        # sentinel
        (sent_dir / ".done_primary_global").write_text(now_tag()+"\n", encoding="utf-8")

    # 总完成旗标
    (cafe_root / ".cafe.done").write_text(now_tag()+"\n", encoding="utf-8")
    print("\n---- CAFÉ5 全部结束。祝皇上凯旋！----")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[FATAL] {e}\n")
        sys.exit(2)