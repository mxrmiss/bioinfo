#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
05_matrix_aggregate.py —— Salmon 转录本结果聚合到“基因层”（驱动脚本）
遵循契约：
  - 读 config.yaml（不接收命令行参数）
  - 读取 results/04_quant/<sample>/quant.sf 与 results/03_maps/tx2gene.clean.tsv
  - 写入：
      results/05_matrix/counts/gene_counts.tsv
      results/05_matrix/tpms/gene_tpm.tsv
      results/05_matrix/matrix_stats.tsv
      results/05_matrix/tximport_meta.tsv（供 R 侧使用）
  - 排除 results/01_qc/rejects.tsv 中样本
  - 调用同目录 05_tximport_aggregate.R，开启“屏幕流式输出+日志文件”
"""

from __future__ import annotations
import sys, csv, logging, subprocess
from pathlib import Path
from typing import Dict, Any, List

CONFIG_PATH = "config.yaml"

DEFAULTS: Dict[str, Any] = {
    "data": {"samples_tsv": "data/samples.tsv"},
    "dirs": {
        "qc": "results/01_qc",
        "ref": "ref",
        "maps": "results/03_maps",
        "quant": "results/04_quant",
        "matrix": "results/05_matrix",
        "logs": "logs"
    },
    "binaries": {"Rscript": "Rscript"},
    "tximport": {
        "counts_from_abundance": "no",
        "matrix_stats": {"report_libsize_quantiles": [0.25, 0.5, 0.75]}
    },
    "logging": {"level": "INFO", "timestamp": True}
}

def load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML，请先安装 pyyaml", file=sys.stderr)
        raise
    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}
    # 递归合并（用户配置覆盖默认）
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

def read_samples(samples_tsv: Path) -> List[Dict[str, str]]:
    if not samples_tsv.exists():
        raise FileNotFoundError(f"未找到样本清单：{samples_tsv}")
    rows: List[Dict[str, str]] = []
    with open(samples_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        required = ["sample", "group", "fastq1", "fastq2"]
        if reader.fieldnames is None or any(c not in reader.fieldnames for c in required):
            raise ValueError(f"samples.tsv 必含列：{', '.join(required)}（契约）")
        for r in reader:
            rows.append({k: (r.get(k, "") or "").strip() for k in reader.fieldnames})
    if not rows:
        raise ValueError("samples.tsv 为空")
    # sample 去重
    seen = set()
    uniq = []
    for r in rows:
        sid = r["sample"]
        if sid in seen:
            raise ValueError(f"samples.tsv 出现重复 sample：{sid}")
        seen.add(sid)
        uniq.append(r)
    return uniq

def read_rejects(rejects_tsv: Path) -> set[str]:
    if not rejects_tsv.exists():
        return set()
    rej = set()
    with open(rejects_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        if reader.fieldnames and "sample" in reader.fieldnames:
            for r in reader:
                sid = (r.get("sample") or "").strip()
                if sid:
                    rej.add(sid)
    return rej

def run_cmd_stream(cmd: List[str], log_file: Path) -> int:
    """流式执行外部命令：实时打印到屏幕，并同步落盘日志。"""
    mkdir_p(log_file.parent)
    logging.info("[CMD] " + " ".join(cmd))
    with open(log_file, "w", encoding="utf-8") as lf:
        lf.write("[CMD] " + " ".join(cmd) + "\n")
        p = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1
        )
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
    # 日志
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True) else "[%(levelname)s] %(message)s",
    )

    logging.info("========== 05 — 表达矩阵聚合（契约版） ==========")
    logging.info("[Info] data.samples_tsv = %s", str(Path(cfg["data"]["samples_tsv"]).resolve()))
    logging.info("[Info] dirs.quant       = %s", str(Path(cfg["dirs"]["quant"]).resolve()))
    logging.info("[Info] dirs.maps        = %s", str(Path(cfg["dirs"]["maps"]).resolve()))
    logging.info("[Info] dirs.matrix      = %s", str(Path(cfg["dirs"]["matrix"]).resolve()))
    logging.info("[Info] binaries.Rscript = %s", cfg["binaries"]["Rscript"])

    # 路径
    samples_tsv = Path(cfg["data"]["samples_tsv"])
    quant_root  = Path(cfg["dirs"]["quant"])
    maps_dir    = Path(cfg["dirs"]["maps"])
    matrix_dir  = Path(cfg["dirs"]["matrix"])
    qc_dir      = Path(cfg["dirs"]["qc"])
    mkdir_p(matrix_dir)
    mkdir_p(matrix_dir / "counts")
    mkdir_p(matrix_dir / "tpms")

    tx2gene = maps_dir / "tx2gene.clean.tsv"
    if not tx2gene.exists():
        raise FileNotFoundError(f"未找到 tx2gene.clean.tsv：{tx2gene}（请先完成 03）")

    # 样本 + 剔除 rejects
    samples = read_samples(samples_tsv)
    rejects = read_rejects(qc_dir / "rejects.tsv")
    kept = []
    for r in samples:
        sid = r["sample"]
        if sid in rejects:
            logging.warning("[Skip] 样本 %s 在 rejects.tsv 中，聚合时剔除", sid)
            continue
        qf = quant_root / sid / "quant.sf"
        if not qf.exists():
            raise FileNotFoundError(f"缺少 quant.sf：{qf}（请先完成 04）")
        kept.append({"sample": sid, "group": r["group"], "quant": str(qf)})
    if not kept:
        raise RuntimeError("没有可聚合的样本（04 未完成或全部被 rejects 剔除）")

    # 写 tximport 元数据（供 R 读取）
    meta_tsv = matrix_dir / "tximport_meta.tsv"
    with open(meta_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["sample", "group", "quant_sf"])
        for r in kept:
            w.writerow([r["sample"], r["group"], r["quant"]])
    logging.info("[Out] 已写 meta：%s", meta_tsv)

    # 调用同目录 R 脚本（只读 config.yaml）
    rscript = cfg["binaries"]["Rscript"]
    r_driver = Path(__file__).resolve().parent / "05_tximport_aggregate.R"
    if not r_driver.exists():
        raise FileNotFoundError(f"未找到 R 脚本：{r_driver}")
    log_path = matrix_dir / "tximport.log"
    rc = run_cmd_stream([rscript, str(r_driver)], log_path)
    if rc != 0:
        raise SystemExit("[ERR] R 侧 tximport 聚合失败，请查看 tximport.log")

    # 验收（契约产物）
    must = [
        matrix_dir / "counts" / "gene_counts.tsv",
        matrix_dir / "tpms"   / "gene_tpm.tsv",
        matrix_dir / "matrix_stats.tsv",
    ]
    for p in must:
        if not p.exists():
            raise SystemExit(f"[ERR] 期望产物缺失：{p}")

    logging.info("========== 05 完成：counts/tpms/matrix_stats 已就位 ==========")

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        print(str(e), file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"[ERR] 05 执行失败：{e}", file=sys.stderr); sys.exit(1)
