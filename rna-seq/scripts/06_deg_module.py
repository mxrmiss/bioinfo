#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
06_deg_module.py —— DEG 驱动（严格契约版）
- 读 config.yaml（不接收命令行参数）
- 检查 05 产物是否就位（counts/gene_counts.tsv 与 tpms/gene_tpm.tsv）
- 检查 contrasts.tsv 表头：contrast, case, control
- 剔除 results/01_qc/rejects.tsv 中样本
- 调用同目录 06_f_deg.R（屏幕流式输出 + 日志）
- 事后对每个对比进行“契约产物”自检
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
        "logs": "logs"
    },
    "binaries": {"Rscript": "Rscript"},
    "deg": {
        "lfc": 1.0,
        "fdr": 0.05,
        "use_apeglm": True,
        "independent_filter": True,
        "allow_batch": True
    },
    "logging": {"level": "INFO", "timestamp": True}
}

def load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML，请安装 pyyaml", file=sys.stderr); raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}
    def merge(user_d, base_d):
        out = dict(base_d)
        for k, v in (user_d or {}).items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out
    return merge(user, DEFAULTS)

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def read_contrasts(p: Path) -> List[Dict[str, str]]:
    if not p.exists():
        raise FileNotFoundError(f"未找到 contrasts.tsv：{p}")
    rows: List[Dict[str, str]] = []
    with open(p, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        required = ["contrast", "case", "control"]
        if reader.fieldnames is None or any(c not in reader.fieldnames for c in required):
            raise ValueError(f"contrasts.tsv 必含列：{', '.join(required)}（契约）")
        for r in reader:
            rows.append({k: (r.get(k, "") or "").strip() for k in required})
    if not rows:
        raise ValueError("contrasts.tsv 为空")
    return rows

def run_cmd_stream(cmd: List[str], log_file: Path) -> int:
    mkdir_p(log_file.parent)
    logging.info("[CMD] " + " ".join(cmd))
    with open(log_file, "w", encoding="utf-8") as lf:
        lf.write("[CMD] " + " ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, bufsize=1)
        assert p.stdout is not None
        for line in p.stdout:
            line = line.rstrip("\n")
            print(line)
            lf.write(line + "\n")
        rc = p.wait()
        lf.write(f"[RC] {rc}\n")
    if rc != 0:
        logging.error("命令执行失败（返回码=%d），详见：%s", rc, log_file)
    return rc

def main() -> None:
    cfg = load_yaml(Path(CONFIG_PATH))
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True)
               else "[%(levelname)s] %(message)s",
    )

    samples_tsv   = Path(cfg["data"]["samples_tsv"])
    contrasts_tsv = Path(cfg["data"]["contrasts_tsv"])
    matrix_dir    = Path(cfg["dirs"]["matrix"])
    deg_root      = Path(cfg["dirs"]["deg"])

    # 打印关键信息
    logging.info("========== 06 — DEG 驱动（严格契约） ==========")
    logging.info("[Info] data.samples_tsv   = %s", str(samples_tsv.resolve()))
    logging.info("[Info] data.contrasts_tsv = %s", str(contrasts_tsv.resolve()))
    logging.info("[Info] dirs.matrix        = %s", str(matrix_dir.resolve()))
    logging.info("[Info] dirs.deg           = %s", str(deg_root.resolve()))
    logging.info("[Info] dirs.qc            = %s", str(Path(cfg["dirs"]["qc"]).resolve()))
    logging.info("[Info] binaries.Rscript   = %s", cfg["binaries"]["Rscript"])

    # 05 产物契约检查
    required_paths = {
        "counts": matrix_dir / "counts" / "gene_counts.tsv",
        "tpms":   matrix_dir / "tpms"   / "gene_tpm.tsv",
        "meta":   matrix_dir / "tximport_meta.tsv"
    }
    missing_names = [name for name, path in required_paths.items() if not path.exists()]
    if missing_names:
        raise SystemExit(f"[ERR] 06 执行失败：缺少以下契约产物：{', '.join(missing_names)} "
                         f"(请先完成 05)")

    # contrasts 表头契约
    contrasts = read_contrasts(contrasts_tsv)

    # 建日志目录
    mkdir_p(deg_root)
    log_path = Path(cfg["dirs"]["logs"]) / "06_deg.log"
    mkdir_p(log_path.parent)

    # 调用 R
    rscript = cfg["binaries"]["Rscript"]
    r_driver = Path(__file__).resolve().parent / "06_f_deg.R"
    if not r_driver.exists():
        raise SystemExit(f"[ERR] 未找到 R 脚本：{r_driver}")
    rc = run_cmd_stream([rscript, str(r_driver)], log_path)
    if rc != 0:
        raise SystemExit("[ERR] 06-R 侧差异分析失败，请查看日志")

    # 事后契约自检：每个对比必须有 5 个文件
    must_files = ["DEG_all.tsv", "DEG_up.tsv", "DEG_down.tsv",
                  "varTop100.list", "rle_range.tsv", "design.txt"]
    bad = []
    for row in contrasts:
        label = row["contrast"]
        deg_dir = deg_root / label
        for fn in must_files:
            if not (deg_dir / fn).exists():
                bad.append(f"{label}/{fn}")
    if bad:
        raise SystemExit("[ERR] 06 完成但契约产物缺失：\n  - " + "\n  - ".join(bad))

    logging.info("========== 06 完成：所有对比契约产物就位 ==========")

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        print(str(e), file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"[ERR] 06 执行失败：{e}", file=sys.stderr); sys.exit(1)
