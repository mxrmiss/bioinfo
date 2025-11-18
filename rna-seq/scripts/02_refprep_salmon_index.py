#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_refprep_salmon_index.py —— 构建 Salmon decoy-aware 索引（严格契约 + 流式输出）

契约（来自《转录组计划1.md》）：
  • 仅读 config.yaml 指定键：
      - reference.ref_fasta（必需）
      - reference.ref_gtf   （必需）
      - reference.salmon.index_dir（必需）
      - reference.salmon.kmer_len（默认 31）
      - reference.salmon.decoy_list（可选；留空则本脚本生成）
      - reference.salmon.gentrome_fa（可选；留空则本脚本生成）
      - resources.threads.salmon（必需；线程数）
  • 构建摘要写入：results/02_ref/gentrome_decoy_summary.tsv
  • 中间产物（若需生成）：ref/transcripts.fa 、（缺省）ref/decoys.txt、ref/gentrome.fa
  • 屏幕“流式输出” + 写日志 logs/02_refprep_salmon_index.log
  • 不接受命令行参数；一切以 config.yaml 为准
"""

from __future__ import annotations
import sys, os, csv, subprocess, shutil, hashlib
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

# 依赖：PyYAML
try:
    import yaml
except Exception:
    print("[ERR] 缺少 PyYAML，请先安装：mamba/conda install pyyaml", file=sys.stderr)
    sys.exit(1)

# ================= 顶部参数（其余全部从 config 读取） =================
CONFIG_PATH = "config.yaml"
DEFAULTS: Dict[str, Any] = {
    "reference": {
        "salmon": {
            "kmer_len": 31,
            # decoy_list / gentrome_fa 若在 config 留空，脚本将使用默认路径生成
            "decoy_list": "",
            "gentrome_fa": "",
            "rebuild": "if_missing"  # always / if_missing / never
        }
    },
    "resources": {
        "threads": {
            # 线程必须由 resources.threads.salmon 提供，若缺失→报错
            # 此处不设默认，避免静默回落
        }
    },
    "binaries": {
        "salmon": "salmon",
        "gffread": "gffread"
    },
    "dirs": {
        "logs": "logs"
    },
    "logging": { "timestamp": True }
}
# =====================================================================

CFG: Dict[str, Any] = {}

def now_ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log(msg: str) -> None:
    """屏幕日志（是否带时间戳由配置控制）"""
    if CFG.get("logging", {}).get("timestamp", True):
        print(f"{now_ts()} {msg}")
    else:
        print(msg)

def load_cfg(path: Path) -> Dict[str, Any]:
    if not path.exists():
        print(f"[ERR] 未找到配置文件：{path}", file=sys.stderr); sys.exit(1)
    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}

    # 浅合并：user 覆盖 defaults
    def merge(u: Dict[str, Any], base: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(base)
        for k, v in u.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(v, out[k])
            else:
                out[k] = v
        return out

    cfg = merge(user, DEFAULTS)

    # 严格取值：参考路径
    ref = cfg.get("reference", {}) or {}
    ref_fa = Path(str(ref.get("ref_fasta", "")).strip())
    ref_gtf = Path(str(ref.get("ref_gtf", "")).strip())
    sal = ref.get("salmon", {}) or {}
    index_dir = Path(str(sal.get("index_dir", "")).strip())
    if not ref_fa or not ref_fa.exists():
        print(f"[ERR] 缺少或无法读取 reference.ref_fasta：{ref_fa}", file=sys.stderr); sys.exit(1)
    if not ref_gtf or not ref_gtf.exists():
        print(f"[ERR] 缺少或无法读取 reference.ref_gtf：{ref_gtf}", file=sys.stderr); sys.exit(1)
    if not index_dir:
        print("[ERR] 缺少 reference.salmon.index_dir（必须）", file=sys.stderr); sys.exit(1)

    # 线程：严格来自 resources.threads.salmon（契约规定）
    threads_block = (cfg.get("resources", {}).get("threads", {}) or {})
    if "salmon" not in threads_block:
        print("[ERR] 缺少线程配置：resources.threads.salmon（必须）", file=sys.stderr); sys.exit(1)
    try:
        threads = int(threads_block["salmon"])
    except Exception:
        print("[ERR] resources.threads.salmon 不是有效整数", file=sys.stderr); sys.exit(1)
    if threads <= 0:
        print("[ERR] resources.threads.salmon 必须为正整数", file=sys.stderr); sys.exit(1)

    # 其它 Salmon 参数
    kmer = int(sal.get("kmer_len", 31))
    rebuild = str(sal.get("rebuild", "if_missing")).lower()
    if rebuild not in ("always", "if_missing", "never"):
        rebuild = "if_missing"

    # decoy_list 与 gentrome_fa：若为空则采用脚本默认路径
    decoy_path = Path(str(sal.get("decoy_list", "")).strip() or "ref/decoys.txt")
    gentrome_fa = Path(str(sal.get("gentrome_fa", "")).strip() or "ref/gentrome.fa")

    # 中间转录本导出路径（固定）
    transcripts_fa = Path("ref/transcripts.fa")

    # 日志与二进制
    salmon_bin = Path(str(cfg.get("binaries", {}).get("salmon", "salmon")))
    gffread_bin = Path(str(cfg.get("binaries", {}).get("gffread", "gffread")))
    log_file = Path(cfg.get("dirs", {}).get("logs", "logs")) / "02_refprep_salmon_index.log"

    # 摘要表路径（契约：results/02_ref）
    summary_tsv = Path("results/02_ref") / "gentrome_decoy_summary.tsv"

    # 把解析后的关键值塞回 cfg 便于后续使用/打印
    cfg["_RESOLVED"] = {
        "ref_fa": ref_fa, "ref_gtf": ref_gtf, "index_dir": index_dir,
        "threads": threads, "kmer": kmer, "rebuild": rebuild,
        "decoy_path": decoy_path, "gentrome_fa": gentrome_fa,
        "transcripts_fa": transcripts_fa, "salmon_bin": salmon_bin,
        "gffread_bin": gffread_bin, "log_file": log_file, "summary_tsv": summary_tsv
    }
    return cfg

def mkdir_p(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def run_cmd_stream(cmd: List[str], log_file: Path) -> int:
    """流式执行外部命令：屏幕实时打印 + 追加写日志文件"""
    mkdir_p(log_file.parent)
    log(f"[CMD] {' '.join(cmd)}")
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write(f"[{now_ts()}] [CMD] {' '.join(cmd)}\n")
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert p.stdout is not None
        for line in p.stdout:
            line = line.rstrip("\n")
            print(line)
            lf.write(line + "\n")
        rc = p.wait()
        lf.write(f"[{now_ts()}] [RC] {rc}\n")
    if rc != 0:
        log(f"[ERR] 命令执行失败，返回码={rc}（详见 {log_file}）")
    return rc

def fasta_headers(fa: Path) -> List[str]:
    """读取 FASTA 的所有头（取 > 后到第一个空白为止）"""
    names: List[str] = []
    with open(fa, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if name:
                    names.append(name)
    return names

def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(8192), b""):
            h.update(ch)
    return h.hexdigest()

def main() -> None:
    global CFG
    CFG = load_cfg(Path(CONFIG_PATH))
    R = CFG["_RESOLVED"]

    # 打印“来源键路径”以便您审核
    log("========== 02 — 构建 Salmon decoy-aware 索引 ==========")
    log(f"[Info] reference.ref_fasta      = {R['ref_fa'].resolve()}")
    log(f"[Info] reference.ref_gtf        = {R['ref_gtf'].resolve()}")
    log(f"[Info] reference.salmon.index_dir = {R['index_dir'].resolve()}")
    log(f"[Info] reference.salmon.kmer_len  = {R['kmer']}")
    log(f"[Info] resources.threads.salmon   = {R['threads']}（严格来自 resources.*）")
    log(f"[Info] decoy_list path             = {R['decoy_path']}（来自 reference.salmon.decoy_list 或默认）")
    log(f"[Info] gentrome_fa path           = {R['gentrome_fa']}（来自 reference.salmon.gentrome_fa 或默认）")
    log(f"[Info] rebuild policy             = {R['rebuild']}")

    # 准备目录
    mkdir_p(Path("ref"))
    mkdir_p(R["index_dir"])
    mkdir_p(R["summary_tsv"].parent)

    # 1) 若 gentrome_fa/decoy_list 未给或不存在，则生成
    # 1.1 提取 transcripts.fa
    need_tx = (R["rebuild"] == "always") or (not R["transcripts_fa"].exists())
    if need_tx:
        if R["transcripts_fa"].exists():
            R["transcripts_fa"].unlink()
        rc = run_cmd_stream([
            str(R["gffread_bin"]),
            "-w", str(R["transcripts_fa"]),
            "-g", str(R["ref_fa"]),
            str(R["ref_gtf"])
        ], R["log_file"])
        if rc != 0 or not R["transcripts_fa"].exists() or R["transcripts_fa"].stat().st_size == 0:
            print("[ERR] gffread 未生成有效的 transcripts.fa", file=sys.stderr); sys.exit(1)
    else:
        log(f"[Skip] 已存在：{R['transcripts_fa']}（按 rebuild 策略跳过）")

    # 1.2 decoy_list
    need_decoy = (R["rebuild"] == "always") or (not R["decoy_path"].exists())
    if need_decoy:
        if R["decoy_path"].exists():
            R["decoy_path"].unlink()
        names = fasta_headers(R["ref_fa"])
        if not names:
            print("[ERR] ref_fasta 未解析到任何序列头，无法生成 decoys.txt", file=sys.stderr); sys.exit(1)
        with open(R["decoy_path"], "w", encoding="utf-8") as f:
            for n in names:
                f.write(n + "\n")
        log(f"[Out] 生成 decoys：{R['decoy_path']}（{len(names)} 条）")
    else:
        log(f"[Skip] 已存在：{R['decoy_path']}（按 rebuild 策略跳过）")

    # 1.3 gentrome.fa
    need_gentrome = (R["rebuild"] == "always") or (not R["gentrome_fa"].exists())
    if need_gentrome:
        if R["gentrome_fa"].exists():
            R["gentrome_fa"].unlink()
        with open(R["gentrome_fa"], "wb") as out, \
             open(R["transcripts_fa"], "rb") as f1, \
             open(R["ref_fa"], "rb") as f2:
            shutil.copyfileobj(f1, out)
            shutil.copyfileobj(f2, out)
        log(f"[Out] 生成 gentrome：{R['gentrome_fa']}")
    else:
        log(f"[Skip] 已存在：{R['gentrome_fa']}（按 rebuild 策略跳过）")

    # 2) 构建索引（若需要）
    need_index = True
    if R["index_dir"].exists():
        marker = R["index_dir"] / "hash.bin"
        if R["rebuild"] == "never":
            need_index = False
        elif R["rebuild"] == "if_missing":
            need_index = not marker.exists()
        elif R["rebuild"] == "always":
            # 清空重建
            for p in R["index_dir"].glob("*"):
                if p.is_file(): p.unlink()
    if need_index:
        cmd = [
            str(R["salmon_bin"]), "index",
            "-t", str(R["gentrome_fa"]),
            "-i", str(R["index_dir"]),
            "-k", str(R["kmer"]),
            "-p", str(R["threads"]),
            "--decoys", str(R["decoy_path"])
        ]
        rc = run_cmd_stream(cmd, R["log_file"])
        if rc != 0:
            sys.exit(1)
    else:
        log(f"[Skip] 索引已存在且 rebuild={R['rebuild']}：{R['index_dir']}")

    # 3) 构建摘要（契约：results/02_ref/gentrome_decoy_summary.tsv）
    def count_fasta(fp: Path) -> int:
        c = 0
        with open(fp, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"): c += 1
        return c
    tx_n = count_fasta(R["transcripts_fa"])
    decoy_n = sum(1 for _ in open(R["decoy_path"], "r", encoding="utf-8"))
    g_size = R["gentrome_fa"].stat().st_size if R["gentrome_fa"].exists() else 0

    with open(R["summary_tsv"], "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, dialect=csv.excel_tab)
        w.writerow(["timestamp", "transcripts_fa", "transcript_count",
                    "genome_fa", "decoy_list", "decoy_count",
                    "gentrome_fa", "gentrome_bytes",
                    "index_dir", "kmer_len", "threads",
                    "ref_fasta_md5", "ref_gtf_md5"])
        w.writerow([
            now_ts(), str(R["transcripts_fa"]), tx_n,
            str(R["ref_fa"]), str(R["decoy_path"]), decoy_n,
            str(R["gentrome_fa"]), g_size,
            str(R["index_dir"]), R["kmer"], R["threads"],
            file_md5(R["ref_fa"]), file_md5(R["ref_gtf"])
        ])
    log(f"[Out] 摘要：{R['summary_tsv']}")
    log("========== 02 完成 ✅ ==========")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERR] 02 执行失败：{e}", file=sys.stderr)
        sys.exit(1)

