#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— APhylo 教程增强正式版
特性：
  1) 自动剔除极端家族（largest differentials）
  2) 自动标记高失败率家族（failure rate >20%）
  3) 两阶段处理 large families（≥阈值，如 100）
  4) 支持误差模型（estimate/use）
  5) 支持 multi-λ（-y 分簇树）
  6) 输出结构稳定：primary_global / large_global / flags / sentinels
  7) 运行失败清晰提示；日志实时流式输出
"""

from __future__ import annotations
import sys, io, re, subprocess, logging
from pathlib import Path
from typing import Any, Dict, List, Tuple
import yaml

DEFAULT_CONFIG = "config.yaml"

# ============================================================
# 基础工具
# ============================================================

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p

def need_dir(p: Path, what: str):
    p = Path(p)
    if not p.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{what}: {p}")
    return p

def need_file(p: Path, what: str):
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what}: {p}")
    return p

def load_config(fp: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(fp)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{fp}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    return cfg

def get_logger(name: str, logfile: Path) -> logging.Logger:
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name)
    lg.setLevel(logging.INFO)
    lg.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s",
                            "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt); fh.setLevel(logging.INFO)

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt); sh.setLevel(logging.INFO)

    lg.addHandler(fh)
    lg.addHandler(sh)

    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x):
            self.s.write(x); self.s.flush()
            return len(x)

    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)

    return lg

def banner(log, text: str):
    bar = "=" * max(10, len(text) + 2)
    log.info(bar)
    log.info(f" {text} ")
    log.info(bar)

# ============================================================
# 解析：极端家族（largest differentials）
# ============================================================

def parse_extreme_ogs(stdout: str) -> List[str]:
    """
    提取 “Families with largest size differentials” 里所有 OG
    """
    pat = re.compile(
        r"Families with largest size differentials:\s*\n"
        r"((?:.*OG\d+:\s*\d+\s*\n)+)",
        flags=re.M
    )
    blocks = pat.findall(stdout)
    if not blocks:
        return []

    last_blk = blocks[-1]
    ogs = re.findall(r"(OG\d+):\s*(\d+)", last_blk)
    return [og for og, diff in ogs]

# ============================================================
# 解析：高失败率家族（failure >20%）
# ============================================================

def parse_high_failure_ogs(stdout: str) -> List[str]:
    """
    修复版解析器：允许前缀时间戳/INFO/模型标签。
    匹配模式：
    The following families had failure rates >20% of the time:
      OGxxxxxxx had NN failures
    """
    # 块定位
    pat = re.compile(
        r"The following families had failure rates >20% of the time:\s*\n"
        r"((?:.*OG\d+\s+had\s+\d+\s+failures.*\n)+)",
        flags=re.M
    )
    blocks = pat.findall(stdout)
    if not blocks:
        return []

    last_blk = blocks[-1]

    ogs = re.findall(r"(OG\d+)\s+had\s+(\d+)\s+failures", last_blk)
    return [og for og, num in ogs]

# ============================================================
# 自动自救：剔除极端 OG（生成新 family）
# ============================================================

def remove_ogs_from_family(family_fp: Path, remove_ogs: List[str],
                           out_fp: Path, log) -> int:
    """
    从 family.tsv 中删除给定 OG，返回删除数量
    """
    if not remove_ogs:
        return 0
    remove_set = set(remove_ogs)

    cnt_rm = 0
    with family_fp.open() as fr, out_fp.open("w") as fw:
        header = fr.readline()
        fw.write(header)
        for line in fr:
            if not line.strip():
                continue
            # Desc \t OGxxxxxx \t counts...
            og = line.strip().split("\t")[1]
            if og in remove_set:
                cnt_rm += 1
                continue
            fw.write(line)

    log.info(f"[AUTOFIX] 从 {family_fp.name} 删除 {cnt_rm} 个极端 OG -> {out_fp.name}")
    return cnt_rm

# ============================================================
# CAFE5 流式运行
# ============================================================

def run_cafe(cmd: List[str], cwd: Path, log, model_tag: str) -> str:
    log.info(f"[CMD][{model_tag}] " + " ".join(map(str, cmd)))

    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        universal_newlines=True,
        bufsize=1
    )

    full = []
    log.info(f"[{model_tag}] ---- CAFE5 输出开始 ----")
    assert proc.stdout is not None
    for line in proc.stdout:
        full.append(line)
        log.info(f"[{model_tag}] {line.rstrip()}")
    proc.wait()
    log.info(f"[{model_tag}] ---- CAFE5 输出结束 ----")

    if proc.returncode != 0:
        raise RuntimeError(f"[ERR] CAFE5 退出码 {proc.returncode}（模型 {model_tag}）")

    return "".join(full)

# ============================================================
# 主流程
# ============================================================

def main():
    cfg = load_config()
    cafe = cfg.get("cafe5", {})
    paths = cfg["paths"]
    log = get_logger("aphylo12", Path(paths["logs_dir"]) / "12_cafe5_run_models.log")
    banner(log, "APhylo 12 — CAFE5（教程增强版）")

    if not cafe.get("enable", True):
        log.info("cafe5.enable = false —— 跳过并写 .cafe.done")
        Path(paths["cafe_run_dir"], ".cafe.done").touch()
        return

    cafe_bin = cfg.get("binaries", {}).get("cafe5", "cafe5")

    # 输入目录
    inp = need_dir(Path(paths["cafe_run_dir"]) / "input", "CAFE 输入目录")
    ftsv = sorted(inp.glob("*.tsv"))
    if not ftsv:
        raise RuntimeError("[ERR] 未找到 family.tsv")
    family0 = ftsv[0]

    tree_files = sorted(inp.glob("*.nwk"))
    if not tree_files:
        raise RuntimeError("[ERR] 未找到超时钟树 .nwk")
    tree0 = tree_files[0]

    family0 = family0.resolve()
    tree0 = tree0.resolve()

    log.info(f"[IN] family: {family0}")
    log.info(f"[IN] tree:   {tree0}")

    # 输出目录
    outdir = ensure_dir(Path(paths["cafe_run_dir"]) / "models" / "global")
    flags_dir = ensure_dir(outdir / "flags")
    sent_dir = ensure_dir(outdir / "sentinels")

    # ------------------------------------------
    # 先进行 “大拷贝拆分” primary / large
    # ------------------------------------------
    ts_large = cafe.get("two_stage_large", {})
    enable_large = ts_large.get("enable", True)
    threshold = int(ts_large.get("copy_threshold", 100))

    family_primary = outdir / "family.primary.tsv"
    family_large = outdir / "family.large.tsv"

    if enable_large:
        # 读取原 family，按阈值拆分
        with family0.open() as fr, \
             family_primary.open("w") as fw_p, \
             family_large.open("w") as fw_l:

            header = fr.readline()
            fw_p.write(header)
            fw_l.write(header)

            lp = ll = 0
            for line in fr:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                counts = list(map(int, parts[2:]))

                if any(c >= threshold for c in counts):
                    fw_l.write(line); ll += 1
                else:
                    fw_p.write(line); lp += 1

        log.info(f"[SPLIT] 阈值 {threshold}: primary={lp} 行, large={ll} 行")
    else:
        # 不拆分，primary=原family，large=空
        family_primary.write_text(family0.read_text(encoding="utf-8"))
        family_large.write_text("")

    # ------------------------------------------
    # 模型参数
    # ------------------------------------------
    threads = int(cafe.get("threads", 4))
    k = cafe.get("gamma_k", None)
    pval = cafe.get("pvalue", None)
    max_round = int(cafe.get("max_autofix_rounds", 3))

    # multi-λ
    mlam = cafe.get("multi_lambda", {})
    multi_enable = mlam.get("enable", False)
    y_tree = mlam.get("y_tree", None)
    compare_global = mlam.get("compare_with_global", True)

    # 误差模型
    em = cafe.get("error_model", {})
    em_mode = em.get("mode", "off")
    em_file = em.get("file", None)
    em_apply = em.get("apply_to", ["primary"])

    # ============================================================
    # Stage 1: primary 集估计 λ（+multi-λ）
    # ============================================================

    primary_dir = ensure_dir(outdir / "primary_global")

    def build_cmd(family_fp: Path, model_tag: str,
                  fixed_lambda: float|None = None,
                  use_y: str|None = None,
                  use_e: str|None = None):

        cmd = [cafe_bin, "-i", str(family_fp), "-t", str(tree0), "-c", str(threads)]
        if fixed_lambda is not None:
            cmd += ["-l", str(fixed_lambda)]
        if k:
            cmd += ["-k", str(k)]
        if pval:
            cmd += ["-P", str(pval)]
        if use_y:
            cmd += ["-y", str(use_y)]
        if use_e:
            cmd += [f"-e{use_e}"]
        return cmd

    # ----------------------
    # 先跑 primary → λ
    # ----------------------
    family_cur = family_primary
    round_id = 1
    lambda_value = None

    while True:
        model_tag = f"PRIMARY-GLOBAL ROUND={round_id}"

        cmd = build_cmd(family_cur, model_tag)
        stdout = run_cafe(cmd, primary_dir, log, model_tag)

        # 解析极端 OG（largest differentials）
        extreme_ogs = parse_extreme_ogs(stdout)

        if extreme_ogs and round_id <= max_round:
            # 继续自救
            nxt = primary_dir / f"family.autofix{round_id}.tsv"
            remove_ogs_from_family(family_cur, extreme_ogs, nxt, log)
            family_cur = nxt
            round_id += 1
            continue

        # 不再自救，解析 λ
        # CAFE5 在 Gamma_results.txt 里有 λ
        res_file = primary_dir / "Gamma_results.txt"
        if not res_file.exists():
            raise RuntimeError("[ERR] 未找到 Gamma_results.txt")

        txt = res_file.read_text()
        m = re.search(r"Lambda:\s*([0-9eE\.\-]+)", txt)
        if not m:
            raise RuntimeError("[ERR] 未在 Gamma_results.txt 中解析 lambda")
        lambda_value = float(m.group(1))

        log.info(f"[PRIMARY] λ = {lambda_value}")
        break

    # 解析高失败率 OG
    high_fail = parse_high_failure_ogs(stdout)
    if high_fail:
        ff = flags_dir / f"high_failure_round{round_id}.tsv"
        with ff.open("w") as fw:
            fw.write("OG\n")
            for og in high_fail:
                fw.write(f"{og}\n")
        log.info(f"[HIGH-FAIL] 标记 {len(high_fail)} 个 high-failure OG -> {ff.name}")

    sent_dir.joinpath(".done_primary").touch()

    # ============================================================
    # Stage 2：误差模型（estimate/use）
    # ============================================================

    error_model_fp = None

    if em_mode == "estimate":
        log.info("[ERROR-MODEL] 开始估计误差模型（-e）")
        cmd = build_cmd(family_cur, "PRIMARY-ERROR-MODEL", use_e="")
        stdout2 = run_cafe(cmd, primary_dir, log, "PRIMARY-ERROR-MODEL")

        emfile = primary_dir / "Base_error_model.txt"
        if not emfile.exists():
            raise RuntimeError("[ERR] 未生成 Base_error_model.txt")

        error_model_fp = emfile
        sent_dir.joinpath(".done_error_model").touch()

    elif em_mode == "use":
        if not em_file:
            raise RuntimeError("[ERR] error_model.mode=use 但 file 未指定")
        error_model_fp = Path(em_file)
        need_file(error_model_fp, "误差模型文件")

    # ============================================================
    # Stage 3：Large families 用固定 λ 重跑
    # ============================================================

    large_dir = ensure_dir(outdir / "large_global")

    if enable_large and family_large.stat().st_size > 1:
        log.info("[LARGE] 大拷贝子集开始运行")

        # 若使用误差模型
        em_arg = str(error_model_fp) if (error_model_fp and "large" in em_apply) else None

        cmd = build_cmd(family_large, "LARGE-GLOBAL",
                        fixed_lambda=lambda_value,
                        use_e=em_arg)
        _ = run_cafe(cmd, large_dir, log, "LARGE-GLOBAL")
        sent_dir.joinpath(".done_large").touch()

    # ============================================================
    # Stage 4：multi-λ（可选）
    # ============================================================

    if multi_enable:
        if not y_tree:
            raise RuntimeError("[ERR] multi_lambda.enable=true 但未提供 y_tree")

        y_tree = Path(y_tree)
        need_file(y_tree, "multi-λ 分簇树")

        ml_dir = ensure_dir(outdir / "primary_multi")
        em_arg = str(error_model_fp) if (error_model_fp and "primary" in em_apply) else None

        cmd = build_cmd(family_cur, "PRIMARY-MULTI",
                        use_y=str(y_tree),
                        use_e=em_arg)
        _ = run_cafe(cmd, ml_dir, log, "PRIMARY-MULTI")

    (Path(paths["cafe_run_dir"]) / ".cafe.done").touch()
    log.info("[DONE] CAFE5 模型全部完成")

# ============================================================
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)