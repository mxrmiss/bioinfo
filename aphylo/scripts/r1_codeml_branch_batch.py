#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
r1_codeml_branch_batch.py

功能：
1）复用 aphylo 现有上游产物，批量运行 PAML/codeml 的 branch model
2）对每个 OG × FG 组合，同时运行：
   - one-ratio model  （model=0, NSsites=0）
   - two-ratio model  （model=2, NSsites=0）
3）输出到：
   results/05_sbranch_model/
     ├── runs/{OG}/{FG}/one/
     ├── runs/{OG}/{FG}/two/
     ├── logs/
     └── summary/r1_run_status.tsv

设计原则：
- 不改动你现有 03/04/05 脚本产物
- 不修改上游对齐文件
- 仅在本脚本自己的 runs 目录里复制并规范化输入副本
- 树直接复用 05_define_foregrounds.py 生成的 {FG}.nwk（其中已带 #1 标记）
- 序列 header 仍按 06 的做法，从 geneID 重写为 species，确保与树 tip 一致

运行方式：
  python scripts/r1_codeml_branch_batch.py

依赖：
- config.yaml
- results/03_codon/codon_msa/{OG}.codon.fna
- results/03_qc/keep_og.list
- results/04_codeml/sets/{FG}.nwk
- inputs.pep2cds_map
- templates/branch_model_one.ctl
- templates/branch_model_two.ctl
"""

from __future__ import annotations

import os
import re
import io
import sys
import csv
import time
import shutil
import hashlib
import datetime
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple, Any, Iterable, Optional

# ==============================
# 顶部固定参数区（皇上主要看这里）
# ==============================

PROJECT_ROOT = Path(".").resolve()
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

# 上游固定输入目录（与当前 aphylo 保持一致）
SEQ_DIR = PROJECT_ROOT / "results/03_codon/codon_msa"
QC_DIR = PROJECT_ROOT / "results/03_qc"
SETS_DIR = PROJECT_ROOT / "results/04_codeml/sets"

# 默认输出根目录（可被 config.yaml: branch_model.run_dir 覆盖）
DEFAULT_RUN_ROOT = PROJECT_ROOT / "results/05_sbranch_model"

# 日志与汇总文件
LOCK_NAME = ".codeml.lock"
FINGERPRINT_NAME = "inputs.sha1"
SUMMARY_FILENAME = "r1_run_status.tsv"
DONE_FILENAME = ".r1.done"

# 需要清理的 codeml 常见产物
CLEAN_PATTERNS = [
    "mlc.txt", "mlc",
    "rst", "rst1", "rst2", "rub", "lnf",
    "2ML.t", "2ML.dN", "2ML.dS", "2ML.out", "2ML.omega",
    "2ML.trees", "2ML.siteclass", "2ML.siteclassN",
    "*.tmp", "tmp*"
]

# ==============================
# 配置读取
# ==============================

def _expand_publish_placeholders(obj, publish_dir: str):
    """递归替换 <publish_dir> 占位符。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj


def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml；若不存在则报错。"""
    if not path.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{path}")

    try:
        import yaml
    except Exception as e:
        raise RuntimeError(f"[ERR] 读取 YAML 需要 pyyaml，但当前环境导入失败：{e}")

    cfg = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg


# ==============================
# 基础小工具
# ==============================

def ts() -> str:
    """返回当前时间字符串。"""
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def ensure_dir(p: Path) -> Path:
    """确保目录存在。"""
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha1_file(p: Path) -> str:
    """计算文件 SHA1。"""
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_fingerprint(dirpath: Path, seq_src: Path, tree_src: Path) -> None:
    """写入输入指纹，用于断点续跑与输入变化检测。"""
    fi = dirpath / FINGERPRINT_NAME
    seq_sz = seq_src.stat().st_size if seq_src.exists() else -1
    tree_sz = tree_src.stat().st_size if tree_src.exists() else -1
    seq_sha = sha1_file(seq_src) if seq_src.exists() else "NA"
    tree_sha = sha1_file(tree_src) if tree_src.exists() else "NA"

    with open(fi, "w", encoding="utf-8") as f:
        f.write(f"seq_size\t{seq_sz}\n")
        f.write(f"seq_sha1\t{seq_sha}\n")
        f.write(f"tree_size\t{tree_sz}\n")
        f.write(f"tree_sha1\t{tree_sha}\n")


def read_fingerprint(dirpath: Path) -> Dict[str, str]:
    """读取输入指纹。"""
    fi = dirpath / FINGERPRINT_NAME
    if not fi.exists():
        return {}

    out = {}
    for line in fi.read_text(encoding="utf-8").splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            out[k] = v
    return out


def inputs_changed(dirpath: Path, seq_src: Path, tree_src: Path) -> bool:
    """判断输入文件是否与上次运行时相比发生变化。"""
    fp = read_fingerprint(dirpath)
    if not fp:
        return True

    try:
        old_seq_sz = int(fp.get("seq_size", "-1"))
        old_tree_sz = int(fp.get("tree_size", "-1"))
    except Exception:
        return True

    cur_seq_sz = seq_src.stat().st_size if seq_src.exists() else -2
    cur_tree_sz = tree_src.stat().st_size if tree_src.exists() else -2

    if cur_seq_sz != old_seq_sz or cur_tree_sz != old_tree_sz:
        return True

    old_seq_sha = fp.get("seq_sha1", "")
    old_tree_sha = fp.get("tree_sha1", "")
    cur_seq_sha = sha1_file(seq_src) if seq_src.exists() else "X"
    cur_tree_sha = sha1_file(tree_src) if tree_src.exists() else "X"

    return (cur_seq_sha != old_seq_sha) or (cur_tree_sha != old_tree_sha)


def clean_dir(d: Path) -> None:
    """仅清理 codeml 产物，不删输入副本与 ctl。"""
    for pat in CLEAN_PATTERNS:
        for p in d.glob(pat):
            try:
                if p.is_file():
                    p.unlink()
                else:
                    shutil.rmtree(p, ignore_errors=True)
            except Exception as e:
                print(f"[warn] 清理失败 {p}: {e}")


def acquire_lock(d: Path) -> bool:
    """创建互斥锁，避免同一工作目录被并发写。"""
    lock = d / LOCK_NAME
    try:
        fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
        return True
    except FileExistsError:
        return False
    except Exception:
        return False


def release_lock(d: Path) -> None:
    """释放互斥锁。"""
    try:
        (d / LOCK_NAME).unlink(missing_ok=True)
    except Exception:
        pass


def run_streaming(cmd: List[str], cwd: Path, log_path: Path, tag: str) -> int:
    """
    流式执行命令：
    - 屏幕实时打印
    - 同时写入日志
    - 末尾写入退出码
    """
    ensure_dir(log_path.parent)

    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")

    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [CMD] (cwd={cwd}) $ {' '.join(cmd)}\n")
        lf.flush()

        p = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            bufsize=1,
            universal_newlines=True
        )

        assert p.stdout is not None
        for line in p.stdout:
            line = line.rstrip("\n")
            print(f"[{tag}] {line}")
            lf.write(line + "\n")
            lf.flush()

        p.wait()
        lf.write(f"[{ts()}] [EXIT] code={p.returncode}\n")
        lf.flush()
        return int(p.returncode)


# ==============================
# 结果有效性判定
# ==============================

_RE_LNL = re.compile(r"lnL\(ntime:\s*\d+", re.I)


def mlc_has_lnl(mlc: Path) -> bool:
    """判断 mlc 中是否已经有有效的 lnL 行。"""
    if not mlc.exists():
        return False

    try:
        with open(mlc, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if _RE_LNL.search(line):
                    return True
    except Exception:
        return False

    return False


# ==============================
# FASTA 与树处理
# ==============================

_TOK_TIP = re.compile(r'(?<=\(|,)\s*([^()\s:,]+)\s*(?=[:),])')


def parse_tips(nwk_text: str) -> List[str]:
    """从 Newick 文本中解析 tip 名。"""
    return _TOK_TIP.findall(nwk_text)


def strip_mark(name: str) -> str:
    """去掉 PAML tree 里的 #1 前景标记。"""
    return name.replace("#1", "")


def fasta_headers(fa: Path) -> List[str]:
    """读取 FASTA 头。"""
    out = []
    with open(fa, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                out.append(line[1:].strip().split()[0])
    return out


def normalize_key(x: str) -> str:
    """
    统一 gene/protein ID 的 key。
    规则：
    - 去两端空白
    - 只取首个空白前 token
    - 若有管道符，仅取最后一个 token
    """
    x = (x or "").strip()
    if not x:
        return ""
    x = x.split()[0]
    return x.split("|")[-1]


def iter_fasta_records(path: Path) -> Iterable[Tuple[str, str]]:
    """逐条读取 FASTA 记录。"""
    header = None
    seq_chunks: List[str] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        yield header, "".join(seq_chunks)


def write_fasta_records(path: Path, recs: Iterable[Tuple[str, str]]) -> None:
    """写出 FASTA。"""
    with open(path, "w", encoding="utf-8") as w:
        for h, s in recs:
            w.write(f">{h}\n")
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + "\n")


def normalize_fasta_inplace(fa: Path, log_path: Path) -> None:
    """
    对 FASTA 副本做规范化：
    - '.' 转 '-'
    - 小写转大写
    - U 转 T
    """
    dot2dash = 0
    lower2upper = 0
    utoT = 0
    lines: List[str] = []

    with open(fa, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                lines.append(line.rstrip("\n") + "\n")
            else:
                s = line.rstrip("\n")
                c1 = s.count(".")
                dot2dash += c1
                s = s.replace(".", "-")

                c2 = sum(1 for ch in s if ch.islower())
                lower2upper += c2
                s = s.upper()

                c3 = s.count("U")
                utoT += c3
                s = s.replace("U", "T")

                lines.append(s + "\n")

    with open(fa, "w", encoding="utf-8") as f:
        for ln in lines:
            f.write(ln)
        if not lines or not lines[-1].endswith("\n"):
            f.write("\n")

    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(
            f"[{ts()}] [NORM] dot_to_dash={dot2dash} "
            f"lower_to_upper={lower2upper} U_to_T={utoT}\n"
        )


def add_phylip_header_if_needed(src_fa: Path, dst_fa: Path) -> None:
    """
    可选：给 FASTA 文本前面加一行 “N L” 风格头。
    注意：
    - 这不是标准 FASTA 写法
    - 你的现有 branch-site 线没有这样做
    - 默认不开启，只保留这个开关以便兼容你未来某些特殊场景
    """
    recs = list(iter_fasta_records(src_fa))
    if not recs:
        raise ValueError(f"序列文件为空：{src_fa}")

    n = len(recs)
    lengths = {len(seq) for _, seq in recs}
    if len(lengths) != 1:
        raise ValueError(f"序列长度不一致，无法写 PHYLIP 头：{src_fa}")

    L = next(iter(lengths))
    text = src_fa.read_text(encoding="utf-8", errors="ignore")
    with open(dst_fa, "w", encoding="utf-8") as w:
        w.write(f"{n} {L}\n")
        w.write(text)


# ==============================
# geneID -> species 映射
# ==============================

def load_pep2cds_map(map_path: Path) -> Dict[str, Dict[str, str]]:
    """
    读取 pep2cds_resolved.tsv，生成：
      og2gid2sp[OG][geneID_key] = species
    """
    if not map_path.exists():
        raise FileNotFoundError(f"缺少映射表 inputs.pep2cds_map：{map_path}")

    txt = map_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not txt:
        raise ValueError(f"映射表为空：{map_path}")

    header = txt[0].rstrip("\n").split("\t")
    cols = {h.strip().lower(): i for i, h in enumerate(header)}

    def find_col(candidates: List[str]) -> int:
        for c in candidates:
            if c in cols:
                return cols[c]
        for k, i in cols.items():
            for c in candidates:
                if c in k:
                    return i
        return -1

    i_og = find_col(["og", "orthogroup"])
    i_sp = find_col(["species", "sp"])
    i_pep = find_col(["protein_id", "pep_id", "peptide_id", "protein", "pep", "prot", "peptide"])

    if i_og < 0 or i_sp < 0 or i_pep < 0:
        raise ValueError("pep2cds_map 表头无法识别所需列，需要至少包含 OG / species / protein(or pep) 列。")

    og2gid2sp: Dict[str, Dict[str, str]] = {}

    for line in txt[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= max(i_og, i_sp, i_pep):
            continue

        og = parts[i_og].strip()
        sp = parts[i_sp].strip()
        pep = parts[i_pep].strip()

        if not og or not sp or not pep:
            continue

        gid = normalize_key(pep)
        if not gid:
            continue

        og2gid2sp.setdefault(og, {})[gid] = sp

    if not og2gid2sp:
        raise ValueError(f"映射表解析后无有效记录：{map_path}")

    return og2gid2sp


def rewrite_codon_headers_geneid_to_species(
    seq_src: Path,
    dst: Path,
    gid2sp: Dict[str, str],
    og: str
) -> None:
    """
    将 seq_src 的 FASTA 头从 geneID 改写成 species。
    """
    src_heads = [normalize_key(h) for h in fasta_headers(seq_src)]
    src_set = set(src_heads)

    map_keys = set(gid2sp.keys())
    miss = sorted(src_set - map_keys)
    if miss:
        raise ValueError(f"{og} 映射表缺少 geneID：n={len(miss)} sample={miss[:5]}")

    sp_seen: Dict[str, str] = {}
    out_recs: List[Tuple[str, str]] = []

    for h, s in iter_fasta_records(seq_src):
        gid = normalize_key(h)
        sp = gid2sp.get(gid, "")
        if not sp:
            raise ValueError(f"{og} geneID 无法映射到 species：{gid}")
        if sp in sp_seen:
            raise ValueError(f"{og} 同一物种出现多个序列：{sp} ({sp_seen[sp]} vs {gid})")
        sp_seen[sp] = gid
        out_recs.append((sp, s))

    write_fasta_records(dst, out_recs)


# ==============================
# 路径管理
# ==============================

def get_run_root(cfg: Dict[str, Any]) -> Path:
    """获取 branch model 输出根目录。"""
    bm = cfg.get("branch_model") or {}
    run_dir = bm.get("run_dir", str(DEFAULT_RUN_ROOT))
    return PROJECT_ROOT / str(run_dir)


def get_templates(cfg: Dict[str, Any]) -> Tuple[Path, Path]:
    """获取 one/two 模板路径。"""
    bm = cfg.get("branch_model") or {}
    one = PROJECT_ROOT / str(bm.get("one_template", "templates/branch_model_one.ctl"))
    two = PROJECT_ROOT / str(bm.get("two_template", "templates/branch_model_two.ctl"))
    return one, two


def unit_paths(run_root: Path, og: str, fg: str) -> Dict[str, Path]:
    """
    返回单个 OG × FG 的 one/two 目录与文件路径。
    """
    one_dir = ensure_dir(run_root / "runs" / og / fg / "one")
    two_dir = ensure_dir(run_root / "runs" / og / fg / "two")

    return {
        "one_dir": one_dir,
        "two_dir": two_dir,

        "one_seq": one_dir / "seq.codon.fna",
        "one_tree": one_dir / "species_tree.nwk",
        "one_ctl": one_dir / "branch_model_one.ctl",
        "one_mlc": one_dir / "mlc.txt",
        "one_log": run_root / "logs" / og / fg / "one.log",

        "two_seq": two_dir / "seq.codon.fna",
        "two_tree": two_dir / "species_tree.nwk",
        "two_ctl": two_dir / "branch_model_two.ctl",
        "two_mlc": two_dir / "mlc.txt",
        "two_log": run_root / "logs" / og / fg / "two.log",
    }


# ==============================
# 输入刷新
# ==============================

def refresh_inputs(
    cfg: Dict[str, Any],
    run_root: Path,
    og: str,
    fg: str,
    tree_text: str,
    gid2sp: Dict[str, str]
) -> Tuple[Path, Path]:
    """
    将当前 OG × FG 的输入复制到 one/two 目录中，并规范化。
    返回原始 seq_src 与 tree_src，用于记录 fingerprint。
    """
    paths = unit_paths(run_root, og, fg)

    seq_src = SEQ_DIR / f"{og}.codon.fna"
    tree_src = SETS_DIR / f"{fg}.nwk"

    if not seq_src.exists():
        raise FileNotFoundError(f"缺少密码子对齐：{seq_src}")
    if not tree_src.exists():
        raise FileNotFoundError(f"缺少前景树：{tree_src}")

    # 写树副本
    paths["one_tree"].write_text(tree_text, encoding="utf-8")
    paths["two_tree"].write_text(tree_text, encoding="utf-8")

    # geneID -> species
    rewrite_codon_headers_geneid_to_species(seq_src, paths["one_seq"], gid2sp, og)
    rewrite_codon_headers_geneid_to_species(seq_src, paths["two_seq"], gid2sp, og)

    # 规范化副本序列
    ensure_dir(paths["one_log"].parent)
    ensure_dir(paths["two_log"].parent)
    normalize_fasta_inplace(paths["one_seq"], paths["one_log"])
    normalize_fasta_inplace(paths["two_seq"], paths["two_log"])

    # 可选：补 PHYLIP 头
    bm = cfg.get("branch_model") or {}
    add_phylip_header = bool(bm.get("add_phylip_header", False))
    if add_phylip_header:
        tmp_one = paths["one_seq"].with_suffix(".tmp")
        tmp_two = paths["two_seq"].with_suffix(".tmp")
        shutil.copy2(paths["one_seq"], tmp_one)
        shutil.copy2(paths["two_seq"], tmp_two)
        add_phylip_header_if_needed(tmp_one, paths["one_seq"])
        add_phylip_header_if_needed(tmp_two, paths["two_seq"])
        tmp_one.unlink(missing_ok=True)
        tmp_two.unlink(missing_ok=True)

    # 拷贝 ctl 模板
    one_tpl, two_tpl = get_templates(cfg)
    if not one_tpl.exists():
        raise FileNotFoundError(f"缺少模板：{one_tpl}")
    if not two_tpl.exists():
        raise FileNotFoundError(f"缺少模板：{two_tpl}")

    shutil.copy2(one_tpl, paths["one_ctl"])
    shutil.copy2(two_tpl, paths["two_ctl"])

    return seq_src, tree_src


# ==============================
# 是否需要重算
# ==============================

def need_recalc(
    dirpath: Path,
    mlc: Path,
    seq_src: Path,
    tree_src: Path,
    overwrite_if_inputs_changed: bool
) -> Tuple[bool, str]:
    """
    判断单个模型目录是否需要重新计算。
    """
    if not mlc.exists():
        return True, "no_mlc"

    if not mlc_has_lnl(mlc):
        return True, "invalid_result"

    if overwrite_if_inputs_changed and inputs_changed(dirpath, seq_src, tree_src):
        return True, "inputs_changed"

    return False, "reuse"


# ==============================
# 单个模型执行
# ==============================

OG2GID2SP: Dict[str, Dict[str, str]] = {}


def run_one_mode(
    cfg: Dict[str, Any],
    run_root: Path,
    og: str,
    fg: str,
    mode: str,
    codeml_bin: str
) -> Dict[str, Any]:
    """
    执行单个模式：
    - mode = one 或 two

    返回一个字典，便于最后写成汇总表。
    """
    assert mode in {"one", "two"}

    paths = unit_paths(run_root, og, fg)
    work_dir = paths["one_dir"] if mode == "one" else paths["two_dir"]
    ctl_name = "branch_model_one.ctl" if mode == "one" else "branch_model_two.ctl"
    mlc_path = paths["one_mlc"] if mode == "one" else paths["two_mlc"]
    log_path = paths["one_log"] if mode == "one" else paths["two_log"]

    start_time = time.time()

    result = {
        "OG": og,
        "foreground": fg,
        "mode": mode,
        "status": "",
        "reason": "",
        "returncode": "",
        "runtime_sec": "",
        "mlc_exists": "",
        "lnl_found": "",
    }

    if not acquire_lock(work_dir):
        result["status"] = "FAIL"
        result["reason"] = "locked"
        result["runtime_sec"] = f"{time.time() - start_time:.2f}"
        return result

    try:
        tree_src = SETS_DIR / f"{fg}.nwk"
        seq_src = SEQ_DIR / f"{og}.codon.fna"

        if not tree_src.exists():
            result["status"] = "FAIL"
            result["reason"] = f"missing_tree:{tree_src}"
            return result

        if not seq_src.exists():
            result["status"] = "FAIL"
            result["reason"] = f"missing_seq:{seq_src}"
            return result

        raw_tree = tree_src.read_text(encoding="utf-8")
        tree_text = raw_tree

        tree_tips = set(strip_mark(x) for x in parse_tips(tree_text))
        if not tree_tips:
            result["status"] = "FAIL"
            result["reason"] = "empty_tree_tips"
            return result

        gid2sp = OG2GID2SP.get(og)
        if not gid2sp:
            result["status"] = "FAIL"
            result["reason"] = "missing_gid2sp_map_for_og"
            return result

        # 一致性校验：序列映射后的物种集合必须与树一致
        fa_geneids = set(normalize_key(h) for h in fasta_headers(seq_src))
        miss_gid = sorted(fa_geneids - set(gid2sp.keys()))
        if miss_gid:
            result["status"] = "FAIL"
            result["reason"] = f"map_missing_geneid:n={len(miss_gid)} sample={miss_gid[:5]}"
            return result

        fa_species = set(gid2sp[gid] for gid in fa_geneids if gid in gid2sp)
        if fa_species != tree_tips:
            only_fa = sorted(fa_species - tree_tips)
            only_tr = sorted(tree_tips - fa_species)
            result["status"] = "FAIL"
            result["reason"] = f"set_mismatch:FA_ONLY={only_fa};TREE_ONLY={only_tr}"
            return result

        bm = cfg.get("branch_model") or {}
        overwrite_if_inputs_changed = bool(bm.get("overwrite_if_inputs_changed", True))
        allow_soft_ok = bool(bm.get("allow_soft_ok", True))

        recalc, why = need_recalc(
            dirpath=work_dir,
            mlc=mlc_path,
            seq_src=seq_src,
            tree_src=tree_src,
            overwrite_if_inputs_changed=overwrite_if_inputs_changed,
        )

        if recalc:
            clean_dir(work_dir)
            _seq_src, _tree_src = refresh_inputs(
                cfg=cfg,
                run_root=run_root,
                og=og,
                fg=fg,
                tree_text=tree_text,
                gid2sp=gid2sp
            )

            rc = run_streaming(
                cmd=[codeml_bin, ctl_name],
                cwd=work_dir,
                log_path=log_path,
                tag=f"{og}|{fg}|{mode.upper()}",
            )

            write_fingerprint(work_dir, _seq_src, _tree_src)

            has_lnl = mlc_has_lnl(mlc_path)
            result["returncode"] = str(rc)
            result["mlc_exists"] = str(mlc_path.exists())
            result["lnl_found"] = str(has_lnl)

            if rc == 0:
                result["status"] = "RUN(new)" if why == "no_mlc" else "RECALC(clean)"
                result["reason"] = why
            else:
                if allow_soft_ok and has_lnl:
                    result["status"] = "SOFT_OK"
                    result["reason"] = f"rc={rc}_but_lnl_found"
                else:
                    result["status"] = "FAIL"
                    result["reason"] = f"rc={rc}"
        else:
            result["status"] = "SKIP(resume)"
            result["reason"] = why
            result["returncode"] = ""
            result["mlc_exists"] = str(mlc_path.exists())
            result["lnl_found"] = str(mlc_has_lnl(mlc_path))

        return result

    except Exception as e:
        result["status"] = "FAIL"
        result["reason"] = f"{e.__class__.__name__}({e})"
        return result

    finally:
        result["runtime_sec"] = f"{time.time() - start_time:.2f}"
        if result["mlc_exists"] == "":
            result["mlc_exists"] = str(mlc_path.exists())
        if result["lnl_found"] == "":
            result["lnl_found"] = str(mlc_has_lnl(mlc_path))
        release_lock(work_dir)


# ==============================
# 任务与汇总
# ==============================

def write_summary_tsv(out_tsv: Path, rows: List[Dict[str, Any]]) -> None:
    """写出 r1 总状态表。"""
    ensure_dir(out_tsv.parent)

    fieldnames = [
        "OG",
        "foreground",
        "one_status",
        "two_status",
        "one_returncode",
        "two_returncode",
        "one_mlc_exists",
        "two_mlc_exists",
        "one_lnl_found",
        "two_lnl_found",
        "one_runtime_sec",
        "two_runtime_sec",
        "message",
    ]

    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def collect_status_row(
    og: str,
    fg: str,
    one_res: Dict[str, Any],
    two_res: Dict[str, Any],
) -> Dict[str, Any]:
    """把 one/two 的结果合成一行。"""
    msg_parts = []
    if one_res.get("reason"):
        msg_parts.append(f"one:{one_res['reason']}")
    if two_res.get("reason"):
        msg_parts.append(f"two:{two_res['reason']}")

    return {
        "OG": og,
        "foreground": fg,
        "one_status": one_res.get("status", ""),
        "two_status": two_res.get("status", ""),
        "one_returncode": one_res.get("returncode", ""),
        "two_returncode": two_res.get("returncode", ""),
        "one_mlc_exists": one_res.get("mlc_exists", ""),
        "two_mlc_exists": two_res.get("mlc_exists", ""),
        "one_lnl_found": one_res.get("lnl_found", ""),
        "two_lnl_found": two_res.get("lnl_found", ""),
        "one_runtime_sec": one_res.get("runtime_sec", ""),
        "two_runtime_sec": two_res.get("runtime_sec", ""),
        "message": " | ".join(msg_parts),
    }


def get_selected_ogs(cfg: Dict[str, Any]) -> List[str]:
    """获取要跑的 OG 列表。"""
    keep_list = QC_DIR / "keep_og.list"
    if not keep_list.exists():
        raise FileNotFoundError(f"[ERR] 缺少 OG 列表：{keep_list}")

    ogs = [x.strip() for x in keep_list.read_text(encoding="utf-8").splitlines() if x.strip()]
    if not ogs:
        raise ValueError("[ERR] keep_og.list 为空")

    bm = cfg.get("branch_model") or {}
    only_ogs = bm.get("only_ogs") or []
    if only_ogs:
        only_set = set(str(x).strip() for x in only_ogs if str(x).strip())
        ogs = [og for og in ogs if og in only_set]

    return ogs


def get_selected_foregrounds(cfg: Dict[str, Any]) -> List[str]:
    """获取要跑的 foreground 集合名。"""
    fgs = sorted([p.stem for p in SETS_DIR.glob("*.nwk")])
    if not fgs:
        raise FileNotFoundError("[ERR] 未发现前景树：results/04_codeml/sets/*.nwk，请先运行 05_define_foregrounds.py")

    bm = cfg.get("branch_model") or {}
    only_fgs = bm.get("only_foregrounds") or []
    if only_fgs:
        only_set = set(str(x).strip() for x in only_fgs if str(x).strip())
        fgs = [fg for fg in fgs if fg in only_set]

    if not fgs:
        raise ValueError("[ERR] 筛选后 foreground 集合为空")

    return fgs


def precheck_foreground_tip_count(cfg: Dict[str, Any], fg: str) -> None:
    """检查前景 tip 数是否达到最小要求。"""
    bm = cfg.get("branch_model") or {}
    min_fg_tips = int(bm.get("min_fg_tips", 1))

    tree_path = SETS_DIR / f"{fg}.nwk"
    text = tree_path.read_text(encoding="utf-8")
    tips = parse_tips(text)
    fg_count = sum(1 for x in tips if "#1" in x)

    if fg_count < min_fg_tips:
        raise ValueError(
            f"[ERR] foreground={fg} 的前景 tip 数为 {fg_count}，小于 branch_model.min_fg_tips={min_fg_tips}"
        )


# ==============================
# 主程序
# ==============================

def main() -> int:
    cfg = load_config(CONFIG_PATH)

    binaries = cfg.get("binaries") or {}
    codeml_bin = str(binaries.get("codeml", "codeml"))

    if shutil.which(codeml_bin) is None:
        print(f"[ERR] 未找到可执行文件：{codeml_bin}")
        return 2

    inputs = cfg.get("inputs") or {}
    pep2cds_map = inputs.get("pep2cds_map", "")
    if not pep2cds_map:
        print("[ERR] config.yaml 缺少 inputs.pep2cds_map")
        return 3

    map_path = Path(str(pep2cds_map)).expanduser()
    global OG2GID2SP
    OG2GID2SP = load_pep2cds_map(map_path)

    bm = cfg.get("branch_model") or {}
    threads_cfg = int(bm.get("threads", 1))
    cpu = max(1, (os.cpu_count() or 1))
    threads = max(1, min(threads_cfg, cpu))

    run_root = get_run_root(cfg)
    ensure_dir(run_root)
    ensure_dir(run_root / "summary")
    ensure_dir(run_root / "logs")

    ogs = get_selected_ogs(cfg)
    fgs = get_selected_foregrounds(cfg)

    # 映射表覆盖性预检
    missing_map_og = [og for og in ogs if og not in OG2GID2SP]
    if missing_map_og:
        print(f"[ERR] pep2cds_map 中缺少 OG 的映射：n={len(missing_map_og)} sample={missing_map_og[:10]}")
        return 3

    # foreground 最小 tip 数预检
    for fg in fgs:
        precheck_foreground_tip_count(cfg, fg)

    print(f"[INIT] {ts()} r1 启动  threads={threads}")
    print(f"[INIT] CODON_ROOT={SEQ_DIR}")
    print(f"[INIT] FG_ROOT={SETS_DIR}")
    print(f"[INIT] RUN_ROOT={run_root}")
    print(f"[PLAN] OG={len(ogs)} × FG={len(fgs)} × 2 (ONE/TWO)")

    all_summary_rows: List[Dict[str, Any]] = []

    total_counter = {
        "RUN(new)": 0,
        "RECALC(clean)": 0,
        "SKIP(resume)": 0,
        "SOFT_OK": 0,
        "FAIL": 0,
        "TOTAL_MODEL_RUNS": 0,
        "TOTAL_UNITS": 0,
    }

    for fg in fgs:
        print(f"\n[FG] >>> {fg}")

        # 先并行跑 one
        one_results: Dict[str, Dict[str, Any]] = {}
        with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
            futs = {
                ex.submit(run_one_mode, cfg, run_root, og, fg, "one", codeml_bin): og
                for og in ogs
            }
            for fu in as_completed(futs):
                og = futs[fu]
                try:
                    res = fu.result()
                except Exception as e:
                    res = {
                        "OG": og,
                        "foreground": fg,
                        "mode": "one",
                        "status": "FAIL",
                        "reason": f"FutureError({e})",
                        "returncode": "",
                        "runtime_sec": "",
                        "mlc_exists": "",
                        "lnl_found": "",
                    }

                one_results[og] = res
                st = res["status"]
                total_counter["TOTAL_MODEL_RUNS"] += 1
                if st in total_counter:
                    total_counter[st] += 1
                else:
                    total_counter["FAIL"] += 1

                print(f"[ONE] {og}|{fg}: {st} | {res.get('reason', '')}")

        # 再并行跑 two
        two_results: Dict[str, Dict[str, Any]] = {}
        with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
            futs = {
                ex.submit(run_one_mode, cfg, run_root, og, fg, "two", codeml_bin): og
                for og in ogs
            }
            for fu in as_completed(futs):
                og = futs[fu]
                try:
                    res = fu.result()
                except Exception as e:
                    res = {
                        "OG": og,
                        "foreground": fg,
                        "mode": "two",
                        "status": "FAIL",
                        "reason": f"FutureError({e})",
                        "returncode": "",
                        "runtime_sec": "",
                        "mlc_exists": "",
                        "lnl_found": "",
                    }

                two_results[og] = res
                st = res["status"]
                total_counter["TOTAL_MODEL_RUNS"] += 1
                if st in total_counter:
                    total_counter[st] += 1
                else:
                    total_counter["FAIL"] += 1

                print(f"[TWO] {og}|{fg}: {st} | {res.get('reason', '')}")

        # 合并 unit 级结果
        for og in ogs:
            one_res = one_results.get(og, {})
            two_res = two_results.get(og, {})
            all_summary_rows.append(collect_status_row(og, fg, one_res, two_res))
            total_counter["TOTAL_UNITS"] += 1

    summary_tsv = run_root / "summary" / SUMMARY_FILENAME
    write_summary_tsv(summary_tsv, all_summary_rows)

    done_file = run_root / DONE_FILENAME
    done_file.touch()

    print("\n[SUMMARY]")
    print(f"  RUN(new)        = {total_counter['RUN(new)']}")
    print(f"  RECALC(clean)   = {total_counter['RECALC(clean)']}")
    print(f"  SKIP(resume)    = {total_counter['SKIP(resume)']}")
    print(f"  SOFT_OK         = {total_counter['SOFT_OK']}")
    print(f"  FAIL            = {total_counter['FAIL']}")
    print(f"  TOTAL_MODEL_RUNS= {total_counter['TOTAL_MODEL_RUNS']}")
    print(f"  TOTAL_UNITS     = {total_counter['TOTAL_UNITS']}")
    print(f"  summary         = {summary_tsv}")
    print(f"  done            = {done_file}")

    return 0 if total_counter["FAIL"] == 0 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断")
        sys.exit(130)
