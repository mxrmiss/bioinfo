#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
04_salmon_quant.py —— Salmon 定量（严格契约·PE）
规则（与约定一致）：
  1) 仅读取项目根的 config.yaml，不接收命令行参数；
  2) 只支持 PE；样本表必须包含列：sample, group, fastq1, fastq2；
  3) 有 clean 就用 clean（results/01_qc/clean/{sample}_R1/2.fastq.gz）；
     若缺 clean：
        - quant.require_clean_fastq=true → 直接报错退出；
        - 否则回落到 raw（samples.tsv 中的 fastq1/fastq2）；
  4) 跳过 qc/rejects.tsv 中列出的样本；
  5) 索引目录有效性：必须存在 info.json + refseq.bin + (hash.bin 或 mphf.bin)；
  6) 线程取 resources.threads.salmon；库型取 reference.salmon.libtype（通常为 "A"）；
  7) 日志与屏幕双通道：logs/04_salmon_quant.log；每个样本记录所用 FASTQ 源（clean/raw）与命令行；
  8) 若 quant.overwrite=false 且目标 quant.sf 已存在 → 跳过并在日志记录。
"""

from __future__ import annotations
import sys, csv, subprocess, shlex
from pathlib import Path
from typing import Dict, Any, List, Tuple
from datetime import datetime

# ========== 基础工具 ==========
def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log_print(msg: str, log_fp: Path, ts: bool = True) -> None:
    line = f"{now()} {msg}" if ts else msg
    print(line)
    log_fp.parent.mkdir(parents=True, exist_ok=True)
    with open(log_fp, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def load_yaml(cfg_path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except Exception:
        print("[ERR] 需要 PyYAML，请先安装：mamba/conda install pyyaml", file=sys.stderr)
        sys.exit(1)
    if not cfg_path.exists():
        print(f"[ERR] 未找到配置文件：{cfg_path}", file=sys.stderr); sys.exit(1)
    with open(cfg_path, "r", encoding="utf-8") as f:
        user = (yaml.safe_load(f) or {})

    # 默认值（用户配置覆盖默认）
    defaults: Dict[str, Any] = {
        "data": {"samples_tsv": "data/samples.tsv"},
        "dirs": {"qc": "results/01_qc", "quant": "results/04_quant", "logs": "logs"},
        "reference": {"salmon": {"index_dir": "ref/salmon_index", "libtype": "A"}},
        "resources": {"threads": {"salmon": 8}},
        "binaries": {"salmon": "salmon"},
        "quant": {"overwrite": False, "require_clean_fastq": False},
        "logging": {"timestamp": True}
    }
    # 递归合并
    def merge(u: Dict[str, Any], base: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(base)
        for k, v in u.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out

    cfg = merge(user, defaults)

    # 解析关键键位（与契约保持一致）
    samples_tsv = Path(str(cfg["data"]["samples_tsv"]))
    qc_root     = Path(str(cfg["dirs"]["qc"]))
    quant_root  = Path(str(cfg["dirs"]["quant"]))
    logs_root   = Path(str(cfg["dirs"]["logs"]))
    index_dir   = Path(str(cfg["reference"]["salmon"]["index_dir"]))
    libtype     = str(cfg["reference"]["salmon"]["libtype"])  # 按契约：libtype 位于 reference.salmon
    threads     = int(cfg["resources"]["threads"]["salmon"])
    salmon_bin  = str(cfg["binaries"]["salmon"])
    overwrite   = bool(cfg["quant"].get("overwrite", False))
    require_clean = bool(cfg["quant"].get("require_clean_fastq", False))
    use_ts = bool(cfg["logging"].get("timestamp", True))

    # 基础检查
    if not samples_tsv.exists():
        print(f"[ERR] 未找到样本表：{samples_tsv}", file=sys.stderr); sys.exit(1)
    # 索引有效：info.json + refseq.bin + (hash.bin 或 mphf.bin)
    have_any = any((index_dir / x).exists() for x in ("hash.bin", "mphf.bin"))
    if not ((index_dir / "info.json").exists() and (index_dir / "refseq.bin").exists() and have_any):
        print(f"[ERR] Salmon 索引目录无效：{index_dir}（缺少 info.json/refseq.bin 或 hash.bin|mphf.bin；请先完成 02）", file=sys.stderr)
        sys.exit(1)

    cfg["_RESOLVED"] = {
        "samples_tsv": samples_tsv, "qc_root": qc_root, "quant_root": quant_root,
        "logs_root": logs_root, "index_dir": index_dir, "libtype": libtype,
        "threads": threads, "salmon_bin": salmon_bin, "overwrite": overwrite,
        "require_clean_fastq": require_clean, "use_ts": use_ts
    }
    return cfg

def read_samples(samples_tsv: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(samples_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, dialect=csv.excel_tab)
        need = ["sample", "group", "fastq1", "fastq2"]
        if not reader.fieldnames or any(c not in reader.fieldnames for c in need):
            print(f"[ERR] samples.tsv 必含列：{', '.join(need)}（现有列：{reader.fieldnames}）", file=sys.stderr)
            sys.exit(1)
        for r in reader:
            sid = (r.get("sample") or "").strip()
            f1  = (r.get("fastq1") or "").strip()
            f2  = (r.get("fastq2") or "").strip()
            if not sid or not f1 or not f2:
                print(f"[ERR] 行缺失：sample/fastq1/fastq2 均需填写（行={r}）", file=sys.stderr); sys.exit(1)
            rows.append({"sample": sid, "group": (r.get("group") or "").strip(), "fastq1": f1, "fastq2": f2})
    # 样本 ID 去重
    seen = set()
    for r in rows:
        if r["sample"] in seen:
            print(f"[ERR] 样本重复：{r['sample']}", file=sys.stderr); sys.exit(1)
        seen.add(r["sample"])
    return rows

def read_rejects(qc_root: Path) -> set[str]:
    rej = set()
    rej_fp = qc_root / "rejects.tsv"
    if rej_fp.exists():
        with open(rej_fp, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, dialect=csv.excel_tab)
            if "sample" in (reader.fieldnames or []):
                for r in reader:
                    sid = (r.get("sample") or "").strip()
                    if sid: rej.add(sid)
    return rej

def decide_fastq(sample: str, raw_f1: str, raw_f2: str, qc_root: Path, require_clean: bool) -> Tuple[Path, Path, str]:
    """返回 (R1, R2, source)，source ∈ {'clean','raw'}；遵循“有 clean 用 clean；否则 raw/或强制报错”"""
    r1c = qc_root / "clean" / f"{sample}_R1.fastq.gz"
    r2c = qc_root / "clean" / f"{sample}_R2.fastq.gz"
    if r1c.exists() and r2c.exists():
        return r1c, r2c, "clean"
    if require_clean:
        print(f"[ERR] 要求使用 clean FASTQ，但未找到：{r1c} 或 {r2c}", file=sys.stderr); sys.exit(1)
    r1, r2 = Path(raw_f1), Path(raw_f2)
    if not r1.exists() or not r2.exists():
        print(f"[ERR] 原始 FASTQ 不存在：{r1} 或 {r2}", file=sys.stderr); sys.exit(1)
    return r1, r2, "raw"

def run_stream(cmd: List[str], log_fp: Path) -> int:
    """流式执行：屏幕实时打印 + 统一日志"""
    pretty = " ".join(shlex.quote(c) for c in cmd)
    log_print("[CMD] " + pretty, log_fp)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    assert p.stdout is not None
    with open(log_fp, "a", encoding="utf-8") as lf:
        for line in p.stdout:
            line = line.rstrip("\n")
            print(line)
            lf.write(line + "\n")
        rc = p.wait()
        lf.write(f"{now()} [RC] {rc}\n")
    if rc != 0:
        log_print(f"[ERR] 命令执行失败，返回码={rc}", log_fp)
    return rc

# ========== 主流程 ==========
def main() -> None:
    cfg = load_yaml(Path("config.yaml"))
    R = cfg["_RESOLVED"]

    logs_root: Path = R["logs_root"]; logs_root.mkdir(parents=True, exist_ok=True)
    quant_root: Path = R["quant_root"]; quant_root.mkdir(parents=True, exist_ok=True)
    qc_root: Path = R["qc_root"]
    log_fp = logs_root / "04_salmon_quant.log"

    # 电子回执（可审计）：把关键键位与取值写入日志与屏幕
    log_print("========== 04 — Salmon 定量（严格契约·PE） ==========", log_fp, R["use_ts"])
    log_print(f"[Info] data.samples_tsv            = {R['samples_tsv'].resolve()}", log_fp, R["use_ts"])
    log_print(f"[Info] reference.salmon.index_dir  = {R['index_dir'].resolve()}", log_fp, R["use_ts"])
    log_print(f"[Info] reference.salmon.libtype    = {R['libtype']}", log_fp, R["use_ts"])
    log_print(f"[Info] resources.threads.salmon    = {R['threads']}", log_fp, R["use_ts"])
    log_print(f"[Info] dirs.quant                   = {R['quant_root'].resolve()}", log_fp, R["use_ts"])
    log_print(f"[Info] dirs.logs                    = {R['logs_root'].resolve()}", log_fp, R["use_ts"])
    log_print(f"[Info] binaries.salmon             = {R['salmon_bin']}", log_fp, R["use_ts"])
    log_print(f"[Info] quant.overwrite             = {R['overwrite']}", log_fp, R["use_ts"])
    log_print(f"[Info] quant.require_clean_fastq   = {R['require_clean_fastq']}", log_fp, R["use_ts"])

    samples = read_samples(R["samples_tsv"])
    rejects = read_rejects(R["qc_root"])
    if rejects:
        log_print(f"[Info] 发现 rejects：{len(rejects)} 个样本将被剔除", log_fp, R["use_ts"])

    # 逐样本执行
    for r in samples:
        sid = r["sample"]
        if sid in rejects:
            log_print(f"[Skip] 样本 {sid} 位于 rejects.tsv，定量时剔除", log_fp, R["use_ts"])
            continue

        out_dir = R["quant_root"] / sid
        out_qsf = out_dir / "quant.sf"
        if out_qsf.exists() and not R["overwrite"]:
            log_print(f"[Keep] 已存在：{out_qsf}（overwrite=false，跳过）", log_fp, R["use_ts"])
            continue

        # 决策 clean/raw
        r1, r2, src = decide_fastq(sid, r["fastq1"], r["fastq2"], R["qc_root"], R["require_clean_fastq"])
        out_dir.mkdir(parents=True, exist_ok=True)
        log_print(f"[Use] sample={sid} | source={src} | R1={r1} | R2={r2}", log_fp, R["use_ts"])

        # 组装 Salmon 命令（契约：只读配置；线程/库型/索引从 config.yaml）
        cmd = [
            R["salmon_bin"], "quant",
            "-i", str(R["index_dir"]),
            "-l", R["libtype"],
            "-1", str(r1),
            "-2", str(r2),
            "-p", str(R["threads"]),
            "-o", str(out_dir),
            "--validateMappings",
            "--gcBias"
        ]
        rc = run_stream(cmd, log_fp)
        if rc != 0:
            print(f"[ERR] Salmon 运行失败（sample={sid}），请查看日志：{log_fp}", file=sys.stderr)
            sys.exit(1)

        # 产物快速验收
        if not out_qsf.exists() or out_qsf.stat().st_size == 0:
            print(f"[ERR] 产物缺失或为空：{out_qsf}", file=sys.stderr); sys.exit(1)
        log_print(f"[Out] {out_qsf}", log_fp, R["use_ts"])

    log_print("========== 04 完成；quant.sf 全部就位 ==========", log_fp, R["use_ts"])

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 04 执行失败：{e}", file=sys.stderr)
        sys.exit(1)

