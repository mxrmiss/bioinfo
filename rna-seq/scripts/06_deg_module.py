#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
06_deg_module.py —— DEG 驱动（严格契约版）
遵守要点（来自“转录组计划”约定）：
  * 对比表仅认可列名：contrast, case, control（不做别名与回退）
  * 计数矩阵只认：results/05_matrix/counts/gene_counts.tsv（无回退）
  * R 端不吃命令行参数与环境变量，只读 config.yaml 与标准路径
  * 产物（每对比）：DEG_all.tsv / DEG_up.tsv / DEG_down.tsv / varTop100.list / rle_range.tsv / design.txt
  * 列名统一为 snake_case：gene_id, log2fc, lfc_se, stat, p_value, p_adjust, base_mean
  * 剔除 results/01_qc/rejects.tsv 中样本（由 R 端执行，Python 仅做存在性提示）
"""

from __future__ import annotations
import sys, csv, logging, subprocess
from pathlib import Path
from typing import Dict, Any, List

CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "data": {"samples_tsv": "data/samples.tsv", "contrasts_tsv": "data/contrasts.tsv"},
    "dirs": {
        "qc": "results/01_qc",
        "matrix": "results/05_matrix",
        "deg": "results/06_deg",
        "logs": "logs",
    },
    "binaries": {"Rscript": "Rscript"},
    "deg": {
        "lfc": 1.0,
        "fdr": 0.05,
        "use_apeglm": True,
        "independent_filter": True,
        "allow_batch": True,
    },
    "logging": {"level": "INFO", "timestamp": True},
}

def load_yaml(path: Path) -> Dict[str, Any]:
    """读取 YAML，并与默认值递归合并；缺关键项即报错。"""
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML，请先安装 pyyaml", file=sys.stderr)
        raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}

    def merge(user_d: Dict[str, Any], base_d: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(base_d)
        for k, v in user_d.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out

    cfg = merge(user, DEFAULTS)
    # 关键顶层键校验（依约定）
    for top in ["data", "dirs", "binaries", "deg"]:
        if top not in cfg:
            raise KeyError(f"配置缺少顶层键：{top}")
    return cfg

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def read_contrasts(contrasts_tsv: Path) -> List[Dict[str, str]]:
    """严格校验 contrasts.tsv 列名：contrast, case, control（不做别名）。"""
    if not contrasts_tsv.exists():
        raise FileNotFoundError(f"未找到 contrasts.tsv：{contrasts_tsv}")
    rows: List[Dict[str, str]] = []
    with open(contrasts_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        must = ["contrast", "case", "control"]
        if reader.fieldnames is None or any(col not in reader.fieldnames for col in must):
            raise ValueError(f"contrasts.tsv 必含列：{', '.join(must)}（现有列：{reader.fieldnames}）")
        for r in reader:
            rows.append({k: (r.get(k, "") or "").strip()})
    if not rows:
        raise ValueError("contrasts.tsv 为空")
    return rows

def run_cmd_stream(cmd: List[str], log_path: Path) -> int:
    """流式执行外部命令：屏幕实时打印 + 写日志到 logs/06_deg_module.log"""
    mkdir_p(log_path.parent)
    logging.info("[CMD] " + " ".join(cmd))
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write("[CMD] " + " ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert p.stdout is not None
        for line in p.stdout:
            line = line.rstrip("\n")
            print(line)
            lf.write(line + "\n")
        rc = p.wait()
        lf.write(f"[RC] {rc}\n")
    if rc != 0:
        logging.error("命令执行失败，返回码=%d（详见 %s）", rc, log_path)
    return rc

def postcheck_outputs(deg_dir: Path, labels: List[str]) -> None:
    """按契约逐对比检查 6 个产物是否存在，不存在即报错。"""
    missing: List[str] = []
    for lb in labels:
        base = deg_dir / lb
        exp = [
            base / "DEG_all.tsv",
            base / "DEG_up.tsv",
            base / "DEG_down.tsv",
            base / "varTop100.list",
            base / "rle_range.tsv",
            base / "design.txt",
        ]
        for p in exp:
            if not p.exists():
                missing.append(str(p))
    if missing:
        raise RuntimeError("06 产物缺失：\n  - " + "\n  - ".join(missing))

def main() -> None:
    # ---------- 日志 ----------
    cfg = load_yaml(Path(CONFIG_PATH))
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True)
               else "[%(levelname)s] %(message)s",
    )

    logging.info("========== 06 — DEG 驱动（严格契约） ==========")
    logging.info(f"[Info] data.samples_tsv   = {Path(cfg['data']['samples_tsv']).resolve()}")
    logging.info(f"[Info] data.contrasts_tsv = {Path(cfg['data']['contrasts_tsv']).resolve()}")
    logging.info(f"[Info] dirs.matrix        = {Path(cfg['dirs']['matrix']).resolve()}")
    logging.info(f"[Info] dirs.deg           = {Path(cfg['dirs']['deg']).resolve()}")
    logging.info(f"[Info] dirs.qc            = {Path(cfg['dirs']['qc']).resolve()}")
    logging.info(f"[Info] binaries.Rscript   = {cfg['binaries']['Rscript']}")

    # ---------- 关键文件存在性（严格路径） ----------
    counts_fp = Path(cfg["dirs"]["matrix"]) / "counts" / "gene_counts.tsv"
    if not counts_fp.exists():
        raise FileNotFoundError(f"缺少计数矩阵：{counts_fp}（请先完成 05）")
    contrasts_tsv = Path(cfg["data"]["contrasts_tsv"])
    labels = [r["contrast"] for r in read_contrasts(contrasts_tsv)]

    # ---------- 调用 R ----------
    script_dir = Path(__file__).resolve().parent
    r_fp = script_dir / "06_f_deg.R"
    if not r_fp.exists():
        raise FileNotFoundError(f"未找到 R 脚本：{r_fp}（应与本 Python 同目录）")
    mkdir_p(Path(cfg["dirs"]["deg"]))
    mkdir_p(Path(cfg["dirs"]["logs"]))

    rc = run_cmd_stream([cfg["binaries"]["Rscript"], str(r_fp)],
                        Path(cfg["dirs"]["logs"]) / "06_deg_module.log")
    if rc != 0:
        raise RuntimeError("R 侧 DEG 计算失败，请查看屏幕输出与 logs/06_deg_module.log")

    # ---------- 产物验收 ----------
    postcheck_outputs(Path(cfg["dirs"]["deg"]), labels)
    logging.info("========== 06 完成；所有对比的 6 件套已就位 ==========")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 06 执行失败：{e}", file=sys.stderr)
        sys.exit(1)

