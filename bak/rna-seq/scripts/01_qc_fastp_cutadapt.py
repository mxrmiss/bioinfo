#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_qc_fastp_cutadapt.py —— 读段级质控（不出图 · 产四表 · 审稿级可追溯）
符合我们定稿契约：
  - 仅从 config.yaml 读取参数（不接收命令行参数）
  - 产出四张 QC 表：sample_qc.tsv / summary.tsv / outliers.tsv / rejects.tsv
  - 产出 fastp JSON/HTML 取证文件；可选生成 clean FASTQ（默认关闭）
  - 硬门槛（Fail-Fast）样本写入 rejects.tsv，供下游自动排除
目录与表头完全按“最终契约”：
  results/01_qc/
    ├─ sample_qc.tsv
    ├─ summary.tsv
    ├─ outliers.tsv
    ├─ rejects.tsv
    ├─ raw_json/<sample>.fastp.json
    ├─ raw_html/<sample>.fastp.html
    └─ clean/<sample>_R1.fastq.gz, <sample>_R2.fastq.gz（当启用写清洗文件时）
"""

from __future__ import annotations
import os, sys, json, csv, math, shutil, subprocess, logging, datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

# =========================
# 顶部集中参数（遵守“集中配置、无命令行”）
# =========================
CONFIG_PATH = "config.yaml"  # 若需修改，直接在此处填写或保持默认
DEFAULTS = {
    "dirs": {"qc": "results/01_qc"},
    "data": {"samples_tsv": "data/samples.tsv"},
    "binaries": {"fastp": "fastp", "cutadapt": "cutadapt"},
    "resources": {"threads": {"qc": 8}},
    "qc": {
        "thresholds": {
            "min_retention_rate": 0.50,   # 保留率 <50% → 淘汰
            "min_raw_reads": 5_000_000,   # 原始读数 <5M → 淘汰
            "min_mean_read_len": 40       # 平均读长 <40bp → 淘汰
        },
        "outlier_rules": {
            "max_adapter_percent": 30,    # 仅提示，不淘汰
            "max_dup_percent": 60
        },
        "write_clean_fastq": False,       # 是否输出 clean FASTQ（默认关闭）
        "fastp_opts": {
            "trim_poly_g": True,
            "detect_adapter": True,
            "qualified_quality_phred": 15,
            "unqualified_percent_limit": 40
        }
    },
    "logging": {"level": "INFO", "timestamp": True, "file": ""},  # 仅标准输出
}

# =========================
# 工具函数
# =========================
def load_yaml(path: Path) -> Dict[str, Any]:
    """读取 config.yaml；若缺失键则用 DEFAULTS 兜底"""
    try:
        import yaml  # 依赖：PyYAML
    except Exception as e:
        print("[ERR] 需要 PyYAML，请先安装（mamba/conda install pyyaml）", file=sys.stderr)
        raise

    if not path.exists():
        raise FileNotFoundError(f"未找到配置文件：{path}")
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    # 递归合并默认值
    def merge(d: Dict[str, Any], default: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(default)
        for k, v in d.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out
    return merge(cfg, DEFAULTS)

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def read_samples(samples_tsv: Path) -> List[Dict[str, str]]:
    """读取 data/samples.tsv；要求列：sample, group, fastq1, fastq2（fastq2 可空）"""
    if not samples_tsv.exists():
        raise FileNotFoundError(f"未找到样本清单：{samples_tsv}")
    rows: List[Dict[str, str]] = []
    with open(samples_tsv, "r", encoding="utf-8") as f:
        sniffer = csv.Sniffer()
        dialect = csv.excel_tab
        reader = csv.DictReader(f, dialect=dialect)
        required = ["sample", "group", "fastq1", "fastq2"]
        for r in reader:
            for col in required:
                if col not in r:
                    raise ValueError(f"samples.tsv 缺少必要列：{col}")
            rows.append({k: (r.get(k, "") or "").strip() for k in r})
    if not rows:
        raise ValueError("samples.tsv 为空")
    # 唯一性检查与路径标准化
    seen = set()
    for r in rows:
        sid = r["sample"]
        if sid in seen:
            raise ValueError(f"samples.tsv 中出现重复 sample：{sid}")
        seen.add(sid)
        # 路径允许相对/绝对，后续 fastp 原样使用，但先做存在性预检查
        # 仅在非空时检查存在性
        if r["fastq1"] and not Path(r["fastq1"]).exists():
            logging.warning(f"[WARN] fastq1 文件不存在：{sid} -> {r['fastq1']}")
        if r["fastq2"] and not Path(r["fastq2"]).exists():
            logging.warning(f"[WARN] fastq2 文件不存在：{sid} -> {r['fastq2']}")
    return rows

def pct1(x: Optional[float]) -> str:
    """百分比格式化为 0–100 一位小数；None → NA"""
    if x is None or math.isnan(x):
        return "NA"
    return f"{x:.1f}"

def to_int_or_na(x: Optional[int]) -> str:
    if x is None:
        return "NA"
    return str(int(x))

def write_tsv(path: Path, header: List[str], rows: List[List[str]]) -> None:
    mkdir_p(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(header)
        for row in rows:
            w.writerow(row)

def run_cmd(cmd: List[str], log_prefix: str) -> Tuple[int, str]:
    """运行外部命令，返回 (returncode, stderr/stdout 合并文本)"""
    logging.info(f"[CMD] {' '.join(cmd)}")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out_lines = []
    assert p.stdout is not None
    for line in p.stdout:
        out_lines.append(line.rstrip("\n"))
    rc = p.wait()
    out = "\n".join(out_lines)
    if rc != 0:
        logging.error(f"{log_prefix} 运行失败（返回码 {rc}）")
    return rc, out

def parse_fastp_json(jpath: Path, paired: bool) -> Dict[str, Any]:
    """解析 fastp JSON，返回标准化指标；若文件缺失则抛异常"""
    with open(jpath, "r", encoding="utf-8") as f:
        jd = json.load(f)

    # before/after 统计
    bf = jd.get("summary", {}).get("before_filtering", {})
    af = jd.get("summary", {}).get("after_filtering", {})
    dup_rate = jd.get("duplication", {}).get("rate", None)

    raw_reads = int(bf.get("total_reads") or 0)
    retained_reads = int(af.get("total_reads") or 0)

    # 百分比（0–100）
    q30 = af.get("q30_rate", None)
    q30_pct = None if q30 is None else float(q30) * 100.0
    gc = af.get("gc_content", None)
    gc_pct = None if gc is None else float(gc) * 100.0

    # 适配子修剪比例（优先按碱基数；否则按读数）
    adapter_info = jd.get("adapter_cutting", {}) or {}
    trimmed_bases = adapter_info.get("adapter_trimmed_bases", None)
    trimmed_reads = adapter_info.get("adapter_trimmed_reads", None)
    total_bases_bf = bf.get("total_bases", None)

    if trimmed_bases is not None and total_bases_bf:
        adapter_pct = float(trimmed_bases) / float(total_bases_bf) * 100.0
    elif trimmed_reads is not None and raw_reads:
        adapter_pct = float(trimmed_reads) / float(raw_reads) * 100.0
    else:
        adapter_pct = None

    # 重复率（0–100）
    dup_pct = None if dup_rate is None else float(dup_rate) * 100.0

    # 平均读长（用于硬门槛判定，但不写入表）
    if paired:
        r1 = af.get("read1_mean_length", None)
        r2 = af.get("read2_mean_length", None)
        if r1 is not None and r2 is not None:
            mean_len = (float(r1) + float(r2)) / 2.0
        elif r1 is not None:
            mean_len = float(r1)
        else:
            mean_len = None
    else:
        r1 = af.get("read1_mean_length", None)
        mean_len = None if r1 is None else float(r1)

    return {
        "raw_reads": raw_reads,
        "retained_reads": retained_reads,
        "retention_rate_pct": (retained_reads / raw_reads * 100.0) if raw_reads > 0 else 0.0,
        "q30_pct": q30_pct,
        "gc_pct": gc_pct,
        "adapter_pct": adapter_pct,
        "dup_pct": dup_pct,
        "mean_len": mean_len,
    }

# =========================
# 主流程
# =========================
def main() -> None:
    # 读取配置
    cfg = load_yaml(Path(CONFIG_PATH))
    qc_dir = Path(cfg["dirs"]["qc"])
    mkdir_p(qc_dir)
    raw_json_dir = qc_dir / "raw_json"
    raw_html_dir = qc_dir / "raw_html"
    clean_dir = qc_dir / "clean"
    mkdir_p(raw_json_dir); mkdir_p(raw_html_dir)

    # 日志初始化
    level = getattr(logging, (cfg.get("logging", {}).get("level") or "INFO").upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s" if cfg.get("logging", {}).get("timestamp", True) else "[%(levelname)s] %(message)s",
    )
    if cfg.get("logging", {}).get("file"):
        fh = logging.FileHandler(cfg["logging"]["file"])
        fh.setLevel(level)
        fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        fh.setFormatter(fmt)
        logging.getLogger().addHandler(fh)

    logging.info("========== 01 — 读段级质控（fastp） ==========")
    logging.info(f"[Info] 使用配置文件：{Path(CONFIG_PATH).resolve()}")
    logging.info(f"[Info] 输出目录：{qc_dir}")

    # 读取样本表
    samples_tsv = Path(cfg["data"]["samples_tsv"])
    samples = read_samples(samples_tsv)
    n_total = len(samples)
    logging.info(f"[Info] 样本数：{n_total}")

    # 读取阈值与开关
    thr = cfg["qc"]["thresholds"]
    outlier = cfg["qc"]["outlier_rules"]
    write_clean = bool(cfg["qc"].get("write_clean_fastq", False))

    # fastp 参数
    fastp_bin = cfg["binaries"]["fastp"]
    threads = int(cfg["resources"]["threads"]["qc"])
    fopts = cfg["qc"]["fastp_opts"]
    trim_poly_g = fopts.get("trim_poly_g", True)
    detect_adapter = fopts.get("detect_adapter", True)
    q_phred = int(fopts.get("qualified_quality_phred", 15))
    unq_pct_limit = int(fopts.get("unqualified_percent_limit", 40))

    # 产出数据容器
    sample_qc_rows: List[List[str]] = []
    outlier_rows: List[List[str]] = []
    reject_rows: List[List[str]] = []

    # 样本级循环
    for i, r in enumerate(samples, 1):
        sid = r["sample"]
        fq1 = r["fastq1"].strip()
        fq2 = r["fastq2"].strip()
        paired = bool(fq2)

        logging.info(f"[{i}/{n_total}] 处理样本：{sid}  (paired={paired})")

        # fastp 输出文件
        jpath = raw_json_dir / f"{sid}.fastp.json"
        hpath = raw_html_dir / f"{sid}.fastp.html"

        # clean FASTQ（可选）
        if write_clean:
            mkdir_p(clean_dir)
            out_r1 = clean_dir / f"{sid}_R1.fastq.gz"
            out_r2 = clean_dir / f"{sid}_R2.fastq.gz" if paired else None
        else:
            out_r1 = None
            out_r2 = None

        # 基础存在性检查
        missing = []
        if not fq1 or not Path(fq1).exists():
            missing.append("fastq1_missing")
        if paired and (not fq2 or not Path(fq2).exists()):
            missing.append("fastq2_missing")
        if missing:
            detail = ";".join(missing)
            logging.error(f"[ERR] 输入 FASTQ 缺失：{sid} -> {detail}")
            reject_rows.append([sid, "file_missing", detail])
            # 即便缺失，仍跳过 fastp 以继续其他样本
            continue

        # 构建 fastp 命令（无 cutadapt，使用 fastp 自身切接头能力）
        cmd = [fastp_bin]
        cmd += ["-i", fq1]
        if paired:
            cmd += ["-I", fq2]
        if out_r1:
            cmd += ["-o", str(out_r1)]
            if paired and out_r2:
                cmd += ["-O", str(out_r2)]
        # JSON/HTML
        cmd += ["-j", str(jpath), "-h", str(hpath)]
        # 质量参数
        cmd += ["-q", str(q_phred), "-u", str(unq_pct_limit)]
        if trim_poly_g:
            cmd += ["--trim_poly_g"]
        if detect_adapter:
            cmd += ["--detect_adapter_for_pe"] if paired else ["--adapter_fasta", ""]  # 单端时 fastp 会自动检测
        # 线程
        cmd += ["-w", str(threads)]

        # 运行 fastp
        start = datetime.datetime.now()
        rc, out = run_cmd(cmd, f"fastp::{sid}")
        elapsed = (datetime.datetime.now() - start).total_seconds()

        if rc != 0 or (not jpath.exists()):
            logging.error(f"[ERR] fastp 失败或 JSON 缺失：{sid}")
            reject_rows.append([sid, "fastp_failed", f"return_code={rc}"])
            continue

        # 解析 fastp JSON
        try:
            met = parse_fastp_json(jpath, paired=paired)
        except Exception as e:
            logging.error(f"[ERR] 解析 JSON 失败：{sid} -> {e}")
            reject_rows.append([sid, "json_parse_error", str(e)])
            continue

        # 采集指标
        raw_reads = met["raw_reads"]
        retained_reads = met["retained_reads"]
        retention_rate = met["retention_rate_pct"]
        q30_pct = met["q30_pct"]
        gc_pct = met["gc_pct"]
        adapter_pct = met["adapter_pct"]
        dup_pct = met["dup_pct"]
        mean_len = met["mean_len"]

        # 写入 sample_qc 行（百分比一律 0–100，一位小数；缺失为 NA）
        sample_qc_rows.append([
            sid,
            to_int_or_na(raw_reads),
            to_int_or_na(retained_reads),
            pct1(retention_rate if retention_rate is not None else math.nan),
            pct1(q30_pct if q30_pct is not None else math.nan),
            pct1(gc_pct if gc_pct is not None else math.nan),
            pct1(adapter_pct if adapter_pct is not None else math.nan),
            pct1(dup_pct if dup_pct is not None else math.nan),
        ])

        # 离群提示（不中断）
        if adapter_pct is not None and adapter_pct > float(outlier["max_adapter_percent"]):
            outlier_rows.append([
                sid, "adapter_percent", pct1(adapter_pct), "max_adapter_percent",
                str(outlier["max_adapter_percent"]), ""
            ])
        if dup_pct is not None and dup_pct > float(outlier["max_dup_percent"]):
            outlier_rows.append([
                sid, "dup_percent", pct1(dup_pct), "max_dup_percent",
                str(outlier["max_dup_percent"]), ""
            ])

        # 硬门槛（Fail-Fast）
        # 注意：配置阈值是比例/整数，本处与实际单位对齐
        fail_reasons = []
        if retention_rate < float(thr["min_retention_rate"]) * 100.0:
            fail_reasons.append(f"retention_rate<{thr['min_retention_rate']*100:.1f}% (obs={pct1(retention_rate)})")
        if raw_reads < int(thr["min_raw_reads"]):
            fail_reasons.append(f"raw_reads<{thr['min_raw_reads']} (obs={raw_reads})")
        if mean_len is not None and mean_len < float(thr["min_mean_read_len"]):
            fail_reasons.append(f"mean_len<{thr['min_mean_read_len']} (obs={mean_len:.1f})")

        if fail_reasons:
            reason = "fail_fast"
            detail = ";".join(fail_reasons)
            reject_rows.append([sid, reason, detail])
            logging.warning(f"[FAIL-FAST] 样本淘汰：{sid} -> {detail}")

        logging.info(f"[Done] {sid} 花费 {elapsed:.1f}s；保留率={pct1(retention_rate)}% Q30={pct1(q30_pct)}% GC={pct1(gc_pct)}%")

    # ============== 写四张表（与契约一致） ==============
    header = ["sample", "raw_reads", "retained_reads", "retention_rate", "q30", "gc_percent", "adapter_percent", "dup_percent"]
    write_tsv(qc_dir / "sample_qc.tsv", header, sample_qc_rows)
    # summary.tsv 行=样本，列同 sample_qc（与 sample_qc 结构相同，便于透传）
    write_tsv(qc_dir / "summary.tsv", header, sample_qc_rows)
    # outliers.tsv
    if not outlier_rows:
        # 仍写出表头（空表），便于下游稳定读取
        write_tsv(qc_dir / "outliers.tsv", ["sample", "metric", "value", "rule", "threshold", "note"], [])
    else:
        write_tsv(qc_dir / "outliers.tsv", ["sample", "metric", "value", "rule", "threshold", "note"], outlier_rows)
    # rejects.tsv
    if not reject_rows:
        write_tsv(qc_dir / "rejects.tsv", ["sample", "reason", "detail"], [])
    else:
        write_tsv(qc_dir / "rejects.tsv", ["sample", "reason", "detail"], reject_rows)

    # 轻量哨兵（内部使用；不纳入契约）
    with open(qc_dir / "qc.done", "w", encoding="utf-8") as f:
        f.write(datetime.datetime.now().isoformat() + "\n")

    logging.info("========== 01 质控完成；四表已产出 ==========")
    logging.info(f"[Out] {qc_dir/'sample_qc.tsv'}")
    logging.info(f"[Out] {qc_dir/'summary.tsv'}")
    logging.info(f"[Out] {qc_dir/'outliers.tsv'}")
    logging.info(f"[Out] {qc_dir/'rejects.tsv'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 01 质控执行失败：{e}", file=sys.stderr)
        sys.exit(1)
