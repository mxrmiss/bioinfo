#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_refprep_salmon_index.py —— 构建 Salmon decoy-aware 索引（vNext 版）

契约（来自《转录组计划 vNext》）：
  • 仅读 config.yaml 指定键：
      - reference.ref_fasta（必需）
      - reference.ref_gtf   （必需）
      - reference.salmon.index_dir（必需）
      - reference.salmon.kmer_len
      - reference.salmon.decoy_list
      - reference.salmon.gentrome_fa
      - reference.salmon.rebuild
      - resources.threads.salmon
      - binaries.salmon, binaries.gffread
      - annotations.id_cleanup（统一 ID 清理开关）
      - dirs.refprep
  • 构建摘要写入：results/02_ref/gentrome_decoy_summary.tsv
  • 中间产物（若需生成）：ref/transcripts.fa 、ref/decoys.txt、ref/gentrome.fa
  • 只从 config 读取参数；不接受命令行参数。
  • transcript_id 如需修剪，只能走 annotations.id_cleanup 这一处开关。
"""

from __future__ import annotations
import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple
import datetime

DEFAULT_CONFIG = "config.yaml"

# ============================= 工具函数：加载配置 =============================

DEFAULTS: Dict[str, Any] = {
    "reference": {
        "ref_fasta": "",
        "ref_gtf": "",
        "salmon": {
            "index_dir": "ref/salmon_index",
            "kmer_len": 31,
            "decoy_list": "ref/decoys.txt",
            "gentrome_fa": "ref/gentrome.fa",
            "rebuild": "if_missing",  # always / if_missing / never
        },
    },
    "resources": {
        "threads": {
            "salmon": 1,
        }
    },
    "binaries": {
        "salmon": "salmon",
        "gffread": "gffread",
    },
    "annotations": {
        "id_cleanup": {
            "strip_prefix": False,
            "prefix": [],
            "strip_suffix": False,
            "suffix": [],
            "order": ["prefix", "suffix"],
        }
    },
    "dirs": {
        "refprep": "results/02_ref",
    },
    "logging": {
        "level": "INFO",
        "timestamp": True,
    },
}


def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml，并用 DEFAULTS 做浅层合并。"""
    try:
        import yaml
    except Exception as e:
        print("[ERR] 需要 PyYAML 支持：pip install pyyaml", file=sys.stderr)
        raise e

    if not path.exists():
        print(f"[ERR] 未找到配置文件：{path}", file=sys.stderr)
        sys.exit(1)

    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}

    def merge(base: Dict[str, Any], u: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(base)
        for k, v in u.items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = merge(out[k], v)
            else:
                out[k] = v
        return out

    return merge(DEFAULTS, user)


# ============================= ID 清理工具 =============================

def apply_id_cleanup(raw: str, policy: Dict[str, Any]) -> str:
    """
    按 annotations.id_cleanup 规则对单个 ID 进行清理。
    仅用于 transcript_id，不应用于 gene_id。
    """
    s = raw
    order = policy.get("order") or ["prefix", "suffix"]
    strip_prefix = bool(policy.get("strip_prefix"))
    strip_suffix = bool(policy.get("strip_suffix"))
    prefixes: List[str] = policy.get("prefix") or []
    suffixes: List[str] = policy.get("suffix") or []

    for step in order:
        if step == "prefix" and strip_prefix:
            for p in prefixes:
                if p and s.startswith(p):
                    s = s[len(p):]
        if step == "suffix" and strip_suffix:
            for suf in suffixes:
                if suf and s.endswith(suf):
                    s = s[:-len(suf)]
    return s


# ============================= FASTA 处理 =============================

def rewrite_fasta_headers(fa_in: Path, fa_out: Path, policy: Dict[str, Any]) -> Tuple[int, int]:
    """
    读取 FASTA，按 ID 清理规则重写头部。
    返回： (记录数, 被修改的 header 数量)
    """
    changed = 0
    total = 0
    with fa_in.open("r", encoding="utf-8") as fin, fa_out.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                total += 1
                header = line[1:].strip()
                # 一般格式：id 后跟空格描述；我们只对第一个 token 做清理
                if not header:
                    fout.write(line)
                    continue
                parts = header.split(maxsplit=1)
                tid = parts[0]
                rest = parts[1] if len(parts) > 1 else ""
                new_tid = apply_id_cleanup(tid, policy)
                if new_tid != tid:
                    changed += 1
                new_header = new_tid + ((" " + rest) if rest else "")
                fout.write(">" + new_header + "\n")
            else:
                fout.write(line)
    return total, changed


def extract_decoys_from_genome(genome_fa: Path, decoy_list: Path) -> int:
    """
    从基因组 FASTA 的 header 中提取所有序列 ID，写入 decoy_list。
    用于 Salmon decoy-aware。
    """
    ids = []
    with genome_fa.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                if not h:
                    continue
                ids.append(h.split(maxsplit=1)[0])
    ids = sorted(set(ids))
    decoy_list.parent.mkdir(parents=True, exist_ok=True)
    with decoy_list.open("w", encoding="utf-8") as out:
        for i in ids:
            out.write(i + "\n")
    return len(ids)


# ============================= 主流程 =============================

def run_cmd(cmd: List[str], log: logging.Logger) -> None:
    """运行外部命令，失败则退出。"""
    log.info("运行命令：%s", " ".join(cmd))
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        log.error("命令执行失败，退出码：%s", proc.returncode)
        sys.exit(proc.returncode)


def main() -> None:
    cfg = load_config(Path(DEFAULT_CONFIG))

    # 日志初始化
    log_level = getattr(logging, str(cfg.get("logging", {}).get("level", "INFO")).upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [02_refprep] %(levelname)s: %(message)s" if cfg["logging"].get("timestamp") else "[02_refprep] %(levelname)s: %(message)s",
    )
    log = logging.getLogger("refprep")

    ref_fa = Path(cfg["reference"]["ref_fasta"])
    ref_gtf = Path(cfg["reference"]["ref_gtf"])
    salmon_cfg = cfg["reference"]["salmon"]
    index_dir = Path(salmon_cfg["index_dir"])
    kmer_len = int(salmon_cfg.get("kmer_len", 31))

    # gentrome / decoys 路径允许为空，需回退到默认
    decoy_list = salmon_cfg.get("decoy_list") or "ref/decoys.txt"
    gentrome_fa = salmon_cfg.get("gentrome_fa") or "ref/gentrome.fa"
    decoy_list = Path(decoy_list)
    gentrome_fa = Path(gentrome_fa)

    threads = int(cfg["resources"]["threads"].get("salmon", 1))
    salmon_bin = cfg["binaries"].get("salmon", "salmon")
    gffread_bin = cfg["binaries"].get("gffread", "gffread")

    refprep_dir = Path(cfg["dirs"]["refprep"])
    refprep_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = refprep_dir / "gentrome_decoy_summary.tsv"

    # ID 清理策略（仅对 transcripts.fa 生效）
    id_policy = cfg.get("annotations", {}).get("id_cleanup", {}) or DEFAULTS["annotations"]["id_cleanup"]

    log.info("ref_fasta: %s", ref_fa)
    log.info("ref_gtf:   %s", ref_gtf)
    log.info("index_dir: %s", index_dir)
    log.info("decoy_list: %s", decoy_list)
    log.info("gentrome_fa: %s", gentrome_fa)

    if not ref_fa.exists():
        log.error("参考基因组 FASTA 不存在：%s", ref_fa)
        sys.exit(1)
    if not ref_gtf.exists():
        log.error("注释 GFF3/GTF 不存在：%s", ref_gtf)
        sys.exit(1)

    rebuild_mode = salmon_cfg.get("rebuild", "if_missing")
    index_exists = index_dir.exists() and any(index_dir.iterdir())
    log.info("索引目录状态：exists=%s, rebuild_mode=%s", index_exists, rebuild_mode)

    # 1) 生成 transcripts.fa
    transcripts_fa = Path("ref/transcripts.fa")
    transcripts_fa.parent.mkdir(parents=True, exist_ok=True)

    if not transcripts_fa.exists() or rebuild_mode == "always":
        log.info("生成转录本 FASTA：%s", transcripts_fa)
        run_cmd(
            [
                gffread_bin,
                "-g",
                str(ref_fa),
                "-w",
                str(transcripts_fa),
                str(ref_gtf),
            ],
            log,
        )
    else:
        log.info("发现已有 transcripts.fa，按 rebuild_mode=%s 保留", rebuild_mode)

    # 2) 按 ID 清理策略重写 transcripts.fa 头部（仅当策略开启）
    use_cleanup = bool(id_policy.get("strip_prefix") or id_policy.get("strip_suffix"))
    cleaned_transcripts = transcripts_fa

    total_tx = 0
    changed_tx = 0
    if use_cleanup:
        cleaned_transcripts = Path("ref/transcripts.cleaned.fa")
        log.info("启用 ID 清理策略，对 transcripts.fa 的 header 进行修剪")
        total_tx, changed_tx = rewrite_fasta_headers(transcripts_fa, cleaned_transcripts, id_policy)
        log.info("transcript 记录数：%d；其中被修改 ID 数量：%d", total_tx, changed_tx)
    else:
        # 粗略统计一下转录本数量
        t = 0
        with transcripts_fa.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    t += 1
        total_tx = t
        log.info("未启用 ID 清理；transcript 记录数：%d", total_tx)

    # 3) 生成 decoy_list
    if not decoy_list.exists() or rebuild_mode == "always":
        log.info("生成 decoy_list：%s", decoy_list)
        n_decoys = extract_decoys_from_genome(ref_fa, decoy_list)
        log.info("decoy 序列数：%d", n_decoys)
    else:
        log.info("发现已有 decoy_list，按 rebuild_mode=%s 保留", rebuild_mode)
        # 估个数
        n_decoys = sum(1 for _ in decoy_list.open("r", encoding="utf-8"))

    # 4) 生成 gentrome.fa（转录本 + 基因组）
    if not gentrome_fa.exists() or rebuild_mode == "always":
        log.info("生成 gentrome.fa：%s", gentrome_fa)
        gentrome_fa.parent.mkdir(parents=True, exist_ok=True)
        with gentrome_fa.open("w", encoding="utf-8") as out:
            with cleaned_transcripts.open("r", encoding="utf-8") as f_tx:
                for line in f_tx:
                    out.write(line)
            with ref_fa.open("r", encoding="utf-8") as f_ref:
                for line in f_ref:
                    out.write(line)
    else:
        log.info("发现已有 gentrome.fa，按 rebuild_mode=%s 保留", rebuild_mode)

    # 5) 构建 salmon 索引
    if rebuild_mode == "never" and index_exists:
        log.info("索引已存在且 rebuild_mode=never，跳过构建")
    elif rebuild_mode == "if_missing" and index_exists:
        log.info("索引已存在且 rebuild_mode=if_missing，跳过构建")
    else:
        index_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            salmon_bin,
            "index",
            "-t",
            str(gentrome_fa),
            "-d",
            str(decoy_list),
            "-i",
            str(index_dir),
            "-k",
            str(kmer_len),
            "-p",
            str(threads),
        ]
        run_cmd(cmd, log)

    # 6) 写 gentrome_decoy_summary.tsv
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    with summary_tsv.open("w", encoding="utf-8") as out:
        out.write("key\tvalue\n")
        out.write(f"ref_fasta\t{ref_fa}\n")
        out.write(f"ref_gtf\t{ref_gtf}\n")
        out.write(f"index_dir\t{index_dir}\n")
        out.write(f"kmer_len\t{kmer_len}\n")
        out.write(f"transcript_fasta\t{cleaned_transcripts}\n")
        out.write(f"n_transcripts\t{total_tx}\n")
        out.write(f"n_transcripts_id_changed\t{changed_tx}\n")
        out.write(f"genome_decoys_fasta\t{ref_fa}\n")
        out.write(f"n_decoys\t{n_decoys}\n")
        out.write(f"decoy_list\t{decoy_list}\n")
        out.write(f"gentrome_fa\t{gentrome_fa}\n")
        out.write(f"id_cleanup_enabled\t{use_cleanup}\n")
        out.write(f"timestamp\t{datetime.datetime.now().isoformat()}\n")

    log.info("gentrome/decoy 构建完成，摘要写入：%s", summary_tsv)


if __name__ == "__main__":
    main()