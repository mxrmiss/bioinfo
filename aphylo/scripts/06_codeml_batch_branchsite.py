#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_codeml_batch_branchsite.py —— 分支位点模型批跑（两段式：先全 ALT 后全 NULL）
稳定策略（与 aphylo 流水线完全对齐）：
  1) 不改写 ctl，仅复制模板：
       templates/branch_site_alt.ctl
       templates/branch_site_null.ctl
     模板中固定文件名：
       seqfile = seq.codon.fna
       treefile = species_tree.nwk
       outfile = mlc.txt
  2) 输入“副本规范化”（仅作用于 alt/null 目录中的拷贝，不碰上游源文件）：
       - 将 '.' 统一替换为 '-' 以消除 difcodonNG 的报警
       - 全部大写；U->T（核酸）
       - 确保末尾换行
  3) “最小树消毒”与集合一致性判定（忽略 #1）：
       - 去方括号注释与内部标签，仅保留叶名与 #1
       - 校验 FASTA 物种集合 == 树叶集合（忽略 #1），否则直接报错列出差集
  4) 两段式并发：
       - 对每个前景 FG：先把所有 OG 的 ALT 全部并发跑完；随后再把 NULL 全部并发跑完
       - NULL 不再以 ALT 成败作为前置条件
  5) 指纹跳过 + 清洁重算（根治 “but mlc.txt has lnL()”）：
       - 对源输入(OG.codon.fna, FG.nwk)计算 size+sha1，写入 inputs.sha1
       - 若 mlc.txt 缺失 / 指纹变更 / 结果无 lnL → 清洁后重算；否则跳过
  6) 子进程稳定：
       - OMP/BLAS 线程强制为 1
       - 流式同时写屏与写入 logs/06_branchsite/OG/FG/{alt|null}.log
       - 末行强制打印 [EXIT] code=…
  7) 软通过：
       - 若退出码 != 0 但 mlc.txt 已出现 lnL(ntime: …)，记为 SOFT_OK，不阻断流程

新增（本次修复点）：
  8) 在 06 内部完成“geneID → species”重命名（只作用于 alt/null 的 seq.codon.fna 副本）：
       - 使用 config.yaml 的 inputs.pep2cds_map（pep2cds_resolved.tsv）
       - 对每个 OG：把 codon MSA 的 header 从 geneID 改成物种名，以匹配 FG 树 tip
       - 只写副本，不修改 results/03_codon/codon_msa 的源文件

固定接口与路径（与 aphylo 其余脚本保持一致）：
  输入：
    results/03_qc/keep_og.list
    results/03_codon/codon_msa/OG*.codon.fna
    results/04_codeml/sets/<FG>.nwk        # 05_define_foregrounds.py 产物
    templates/branch_site_alt.ctl
    templates/branch_site_null.ctl
  输出（逐 OG×FG）：
    results/04_codeml/raw/OG/FG/alt/{seq.codon.fna,species_tree.nwk,branch_site_alt.ctl,mlc.txt,...}
    results/04_codeml/raw/OG/FG/null/{seq.codon.fna,species_tree.nwk,branch_site_null.ctl,mlc.txt,...}

可由 config.yaml 覆盖的键（若不存在则使用默认值）：
  binaries.codeml         # 默认 "codeml"
  codeml.threads          # 默认 1（推荐按机器核数设置）
  inputs.pep2cds_map      # 必须：geneID ↔ species 映射表（pep2cds_resolved.tsv）
"""

from __future__ import annotations
import os, re, io, sys, shutil, hashlib, datetime, subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple, Any, Iterable

# ================= 固定路径（与 aphylo 约定一致） =================
PROJECT_ROOT = Path(".").resolve()
CONFIG_PATH  = PROJECT_ROOT / "config.yaml"

SEQ_DIR      = PROJECT_ROOT / "results/03_codon/codon_msa"   # 密码子对齐源（固定）
QC_DIR       = PROJECT_ROOT / "results/03_qc"                # keep_og.list
SETS_DIR     = PROJECT_ROOT / "results/04_codeml/sets"       # <FG>.nwk（05 产物）
RAW_ROOT     = PROJECT_ROOT / "results/04_codeml/raw"        # 输出根
TPL_DIR      = PROJECT_ROOT / "templates"                    # ctl 模板
LOG_ROOT     = PROJECT_ROOT / "logs/06_branchsite"           # 日志根

ALT_CTL_NAME  = "branch_site_alt.ctl"
NULL_CTL_NAME = "branch_site_null.ctl"

LOCK_NAME        = ".codeml.lock"    # 互斥锁
FINGERPRINT_NAME = "inputs.sha1"     # 指纹文件

# 清洁模式（仅清产物，不动输入拷贝与 ctl）
CLEAN_PATTERNS = [
    "mlc.txt", "rst", "rst1", "rst2", "rub", "lnf",
    "2ML.t", "2ML.dN", "2ML.dS", "2ML.out", "2ML.omega",
    "2ML.trees", "2ML.siteclass", "2ML.siteclassN", "*.tmp", "tmp*"
]

# ======== 配置加载 ========
def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str):  return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(path: Path) -> Dict[str, Any]:
    """读取 config.yaml（若不存在则返回空 dict；允许 <publish_dir> 占位符）"""
    if not path.exists(): return {}
    try:
        import yaml
    except Exception:
        return {}
    try:
        cfg = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        pub = cfg.get("publish_dir")
        if pub:
            cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
        return cfg
    except Exception:
        return {}

# ======== 小工具 ========
def ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def sha1_file(p: Path) -> str:
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()

def write_fingerprint(dirpath: Path, seq_src: Path, tree_src: Path):
    fi = dirpath / FINGERPRINT_NAME
    seq_sz  = seq_src.stat().st_size if seq_src.exists() else -1
    tree_sz = tree_src.stat().st_size if tree_src.exists() else -1
    seq_sha = sha1_file(seq_src) if seq_src.exists() else "NA"
    tree_sha= sha1_file(tree_src) if tree_src.exists() else "NA"
    with open(fi, "w", encoding="utf-8") as f:
        f.write(f"seq_size\t{seq_sz}\nseq_sha1\t{seq_sha}\n")
        f.write(f"tree_size\t{tree_sz}\ntree_sha1\t{tree_sha}\n")

def read_fingerprint(dirpath: Path) -> Dict[str,str]:
    fi = dirpath / FINGERPRINT_NAME
    if not fi.exists(): return {}
    out={}
    for line in fi.read_text(encoding="utf-8").splitlines():
        if "\t" in line:
            k,v = line.split("\t",1); out[k]=v
    return out

def inputs_changed(dirpath: Path, seq_src: Path, tree_src: Path) -> bool:
    fp = read_fingerprint(dirpath)
    if not fp: return True
    try:
        old_seq_sz=int(fp.get("seq_size","-1")); old_tree_sz=int(fp.get("tree_size","-1"))
    except Exception:
        return True
    cur_seq_sz = seq_src.stat().st_size if seq_src.exists() else -2
    cur_tree_sz= tree_src.stat().st_size if tree_src.exists() else -2
    if cur_seq_sz!=old_seq_sz or cur_tree_sz!=old_tree_sz: return True
    old_seq_sha = fp.get("seq_sha1","")
    old_tree_sha= fp.get("tree_sha1","")
    cur_seq_sha = sha1_file(seq_src) if seq_src.exists() else "X"
    cur_tree_sha= sha1_file(tree_src) if tree_src.exists() else "X"
    return (cur_seq_sha!=old_seq_sha) or (cur_tree_sha!=old_tree_sha)

def clean_dir(d: Path):
    for pat in CLEAN_PATTERNS:
        for p in d.glob(pat):
            try:
                p.unlink() if p.is_file() else shutil.rmtree(p, ignore_errors=True)
            except Exception as e:
                print(f"[warn] 清理失败 {p}: {e}")

def acquire_lock(d: Path)->bool:
    lock = d / LOCK_NAME
    try:
        fd = os.open(str(lock), os.O_CREAT|os.O_EXCL|os.O_WRONLY)
        os.close(fd); return True
    except FileExistsError:
        return False
    except Exception:
        return False

def release_lock(d: Path):
    try: (d/LOCK_NAME).unlink(missing_ok=True)
    except Exception: pass

def run_streaming(cmd: List[str], cwd: Path, log_path: Path, tag: str)->int:
    """流式执行命令：同时写屏与写日志；末行强制写 [EXIT] code=…"""
    ensure_dir(log_path.parent)
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [CMD] (cwd={cwd}) $ {' '.join(cmd)}\n"); lf.flush()
        p = subprocess.Popen(
            cmd, cwd=str(cwd),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, env=env, bufsize=1, universal_newlines=True
        )
        assert p.stdout is not None
        for line in p.stdout:
            line=line.rstrip("\n")
            print(f"[{tag}] {line}")
            lf.write(line+"\n"); lf.flush()
        p.wait()
        lf.write(f"[{ts()}] [EXIT] code={p.returncode}\n"); lf.flush()
        return int(p.returncode)

# ======== MLC 判定 ========
_RE_LNL = re.compile(r"lnL\(ntime:\s*\d+", re.I)
def mlc_has_lnl(mlc: Path)->bool:
    if not mlc.exists(): return False
    try:
        with open(mlc,"r",encoding="utf-8",errors="ignore") as f:
            for line in f:
                if _RE_LNL.search(line): return True
    except Exception:
        return False
    return False

# ======== 树 & FASTA 处理 ========
_TOK_TIP = re.compile(r'(?<=\(|,)\s*([^()\s:,;]+)\s*(?=[:),;])')

def sanitize_foreground_tree(raw_text: str) -> str:
    """最小树消毒：去方括号注释、内部标签，压缩空白，补分号"""
    s = raw_text.strip()
    s = re.sub(r'\[.*?\]', '', s)
    s = re.sub(r'(\))\s*[^():,;\s]+\s*(?=[:),;])', r'\1', s)
    s = re.sub(r'\s+', ' ', s).strip()
    if not s.endswith(';'): s += ';'
    if s.count('(') != s.count(')'):
        raise ValueError("物种树括号不匹配，无法用于 codeml")
    return s

def parse_tips(nwk_text: str) -> List[str]:
    return _TOK_TIP.findall(nwk_text)

def strip_mark(name: str) -> str:
    return name.replace("#1","")

def fasta_headers(fa: Path) -> List[str]:
    out=[]
    with open(fa,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                out.append(line[1:].strip().split()[0])
    return out

def normalize_key(x: str) -> str:
    x = (x or "").strip()
    if not x: return ""
    x = x.split()[0]
    return x.split("|")[-1]

def iter_fasta_records(path: Path) -> Iterable[Tuple[str, str]]:
    h=None; s=[]
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(s)
                h=line[1:].strip().split()[0]
                s=[]
            else:
                s.append(line)
    if h is not None:
        yield h, "".join(s)

def write_fasta_records(path: Path, recs: Iterable[Tuple[str,str]]):
    with open(path, "w", encoding="utf-8") as w:
        for h, s in recs:
            w.write(f">{h}\n")
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + "\n")

# ======== 副本规范化（仅作用于 alt/null 目录内的 seq.codon.fna 副本） ========
def normalize_fasta_inplace(fa: Path, log_path: Path):
    """将 '.'→'-'、u→t、全大写；保证末尾换行；统计修改量并写入日志"""
    dot2dash = lower2upper = utoT = 0
    lines=[]
    with open(fa,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                lines.append(line.rstrip("\n")+"\n")
            else:
                s=line.rstrip("\n")
                c1=s.count(".");   dot2dash += c1; s=s.replace(".","-")
                c2=sum(1 for ch in s if ch.islower()); lower2upper += c2; s=s.upper()
                c3=s.count("U");   utoT += c3; s=s.replace("U","T")
                lines.append(s+"\n")
    with open(fa,"w",encoding="utf-8") as f:
        for ln in lines: f.write(ln)
        if not lines or not lines[-1].endswith("\n"): f.write("\n")
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [NORM] dot_to_dash={dot2dash} lower_to_upper={lower2upper} U_to_T={utoT}\n")

# ======== geneID → species 映射表读取 ========
def load_pep2cds_map(map_path: Path) -> Dict[str, Dict[str, str]]:
    """
    读取 pep2cds_resolved.tsv，生成：
      og2gid2sp[OG][geneID_key] = Species_name
    其中 geneID_key 采用 normalize_key（取最后一个 '|' 后的 token）。
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
        raise ValueError(
            "pep2cds_map 表头无法识别所需列。需要至少包含：OG / species / protein(or pep) 列。"
        )

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

def rewrite_codon_headers_geneid_to_species(seq_src: Path, dst: Path, gid2sp: Dict[str, str], og: str):
    """
    将 seq_src（header=geneID）改写为 dst（header=species），不修改序列本体。
    """
    src_heads = [normalize_key(h) for h in fasta_headers(seq_src)]
    src_set = set(src_heads)

    map_keys = set(gid2sp.keys())
    miss = sorted(src_set - map_keys)
    if miss:
        raise ValueError(f"{og} 映射表缺少 geneID：n={len(miss)} sample={miss[:5]}")

    # 防止一个物种被写两次（严格 SCO 理论上不会发生，但这里防御）
    sp_seen = {}
    out_recs = []
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

# ======== 路径拼装 ========
def unit_paths(og: str, fg: str) -> Dict[str, Path]:
    unit = RAW_ROOT / og / fg
    alt_d = ensure_dir(unit / "alt")
    nul_d = ensure_dir(unit / "null")
    paths = dict(
        alt_dir = alt_d,
        nul_dir = nul_d,
        alt_seq = alt_d / "seq.codon.fna",
        alt_tree= alt_d / "species_tree.nwk",
        alt_ctl = alt_d / ALT_CTL_NAME,
        alt_mlc = alt_d / "mlc.txt",
        alt_log = LOG_ROOT / og / fg / "alt.log",
        nul_seq = nul_d / "seq.codon.fna",
        nul_tree= nul_d / "species_tree.nwk",
        nul_ctl = nul_d / NULL_CTL_NAME,
        nul_mlc = nul_d / "mlc.txt",
        nul_log = LOG_ROOT / og / fg / "null.log",
    )
    return paths

# ======== 复制 + 消毒/规范化（仅对副本操作） ========
def refresh_inputs(og: str, fg: str, clean_tree: str, gid2sp: Dict[str, str]) -> Tuple[Path, Path]:
    """将源输入拷贝到 alt/null，并对副本进行规范化；返回 (seq_src, tree_src)"""
    seq_src  = SEQ_DIR / f"{og}.codon.fna"
    tree_src = SETS_DIR / f"{fg}.nwk"
    if not seq_src.exists():  raise FileNotFoundError(f"缺少密码子对齐：{seq_src}")
    if not tree_src.exists(): raise FileNotFoundError(f"缺少前景树：{tree_src}")

    paths = unit_paths(og, fg)

    # 写入树副本（已经过消毒）
    paths["alt_tree"].write_text(clean_tree, encoding="utf-8")
    paths["nul_tree"].write_text(clean_tree, encoding="utf-8")

    # 生成序列副本：geneID → species（只写到 alt/null）
    rewrite_codon_headers_geneid_to_species(seq_src, paths["alt_seq"], gid2sp, og)
    rewrite_codon_headers_geneid_to_species(seq_src, paths["nul_seq"], gid2sp, og)

    # 对副本做规范化，并在各自日志中记录统计
    ensure_dir(paths["alt_log"].parent); ensure_dir(paths["nul_log"].parent)
    normalize_fasta_inplace(paths["alt_seq"], paths["alt_log"])
    normalize_fasta_inplace(paths["nul_seq"], paths["nul_log"])

    # 复制 ctl 模板
    alt_tpl = TPL_DIR / ALT_CTL_NAME
    nul_tpl = TPL_DIR / NULL_CTL_NAME
    if not alt_tpl.exists(): raise FileNotFoundError(f"缺少模板：{alt_tpl}")
    if not nul_tpl.exists(): raise FileNotFoundError(f"缺少模板：{nul_tpl}")
    shutil.copy2(alt_tpl, paths["alt_ctl"])
    shutil.copy2(nul_tpl, paths["nul_ctl"])

    return seq_src, tree_src

# ======== 判定是否需要重算 ========
def need_recalc(dirpath: Path, mlc: Path, is_alt: bool, seq_src: Path, tree_src: Path) -> Tuple[bool, str]:
    if not mlc.exists():
        return True, "no_mlc"
    if inputs_changed(dirpath, seq_src, tree_src):
        return True, "inputs_changed"
    if not mlc_has_lnl(mlc):
        return True, "invalid_result"
    return False, "reuse"

# ======== 单个 ALT/NULL 子任务 ========
OG2GID2SP: Dict[str, Dict[str, str]] = {}

def run_one(og: str, fg: str, mode: str, codeml_bin: str) -> Tuple[str, str, str]:
    """
    mode: "alt" or "null"
    返回: (og, fg, 状态字符串)
    状态字符串：RUN(new)/RECALC(clean)/SKIP(resume)/SOFT_OK(rc!=0且有lnL)/FAIL:原因
    """
    paths = unit_paths(og, fg)
    work_dir = paths["alt_dir"] if mode=="alt" else paths["nul_dir"]
    ctl_name = ALT_CTL_NAME if mode=="alt" else NULL_CTL_NAME
    mlc_path = paths["alt_mlc"] if mode=="alt" else paths["nul_mlc"]
    log_path = paths["alt_log"] if mode=="alt" else paths["nul_log"]

    if not acquire_lock(work_dir):
        return og, fg, f"FAIL:{mode}_locked"

    try:
        tree_src = SETS_DIR / f"{fg}.nwk"
        if not tree_src.exists(): return og, fg, f"FAIL:missing_tree({tree_src})"
        raw_tree = tree_src.read_text(encoding="utf-8")
        clean_tree = sanitize_foreground_tree(raw_tree)

        tree_tips = set(strip_mark(x) for x in parse_tips(clean_tree))
        if not tree_tips:
            return og, fg, "FAIL:empty_tree_tips"

        seq_src = SEQ_DIR / f"{og}.codon.fna"
        if not seq_src.exists(): return og, fg, f"FAIL:missing_seq({seq_src})"

        gid2sp = OG2GID2SP.get(og)
        if not gid2sp:
            return og, fg, "FAIL:missing_gid2sp_map_for_og"

        # 用映射表得到该 OG 的物种集合，并与树 tip 集合一致性校验
        fa_geneids = set(normalize_key(h) for h in fasta_headers(seq_src))
        miss_gid = sorted(fa_geneids - set(gid2sp.keys()))
        if miss_gid:
            return og, fg, f"FAIL:map_missing_geneid(n={len(miss_gid)} sample={miss_gid[:5]})"

        fa_species = set(gid2sp[gid] for gid in fa_geneids if gid in gid2sp)
        if fa_species != tree_tips:
            only_fa = sorted(fa_species - tree_tips)
            only_tr = sorted(tree_tips - fa_species)
            return og, fg, f"FAIL:set_mismatch(FA_ONLY={only_fa},TREE_ONLY={only_tr})"

        recalc, why = need_recalc(work_dir, mlc_path, mode=="alt", seq_src, tree_src)

        if recalc:
            clean_dir(work_dir)
            _seq_src, _tree_src = refresh_inputs(og, fg, clean_tree, gid2sp)
            rc = run_streaming([codeml_bin, ctl_name], work_dir, log_path, f"{og}|{fg}|{mode.upper()}")
            write_fingerprint(work_dir, _seq_src, _tree_src)
            if rc != 0 and mlc_has_lnl(mlc_path):
                return og, fg, "SOFT_OK"
            elif rc != 0 and not mlc_has_lnl(mlc_path):
                return og, fg, f"FAIL:{mode}_rc={rc}"
            else:
                return og, fg, "RUN(new)" if why=="no_mlc" else "RECALC(clean)"
        else:
            return og, fg, "SKIP(resume)"

    except Exception as e:
        return og, fg, f"FAIL:{e.__class__.__name__}({e})"
    finally:
        release_lock(work_dir)

# ======== 两段式执行：先 ALT 后 NULL ========
def run_pass(fg: str, mode: str, ogs: List[str], threads: int, codeml_bin: str, totals: Dict[str,int]):
    print(f"\n[FG] >>> {fg}  [{mode.upper()}-PASS]")
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        futs = {ex.submit(run_one, og, fg, mode, codeml_bin):(og, fg) for og in ogs}
        for fu in as_completed(futs):
            og, _fg = futs[fu]
            try:
                og2, fg2, st = fu.result()
            except Exception as e:
                print(f"[FAIL] {og}|{fg}|{mode.upper()}: {e}")
                totals["FAIL"]+=1; totals["TOTAL"]+=1
                continue
            totals["TOTAL"]+=1
            if st.startswith("FAIL"):
                totals["FAIL"]+=1
                print(f"[FAIL] {og}|{fg}|{mode.upper()}: {st}")
            elif st=="SOFT_OK":
                totals["SOFT_OK"]+=1
                print(f"[SOFT] {og}|{fg}|{mode.upper()}: rc!=0 but mlc has lnL()")
            else:
                totals[st]+=1
                print(f"[OK]   {og}|{fg}|{mode.upper()}: {st}")

# ======== 主流程 ========
def main()->int:
    cfg = load_config(CONFIG_PATH)
    codeml_bin = (cfg.get("binaries") or {}).get("codeml", "codeml")
    threads_cfg = int((cfg.get("codeml") or {}).get("threads", 1))
    cpu = max(1, (os.cpu_count() or 1))
    threads = max(1, min(threads_cfg, cpu))

    if shutil.which(codeml_bin) is None:
        print(f"[ERR] 未找到可执行文件：{codeml_bin}")
        return 2

    # 读取映射表（必须）
    inputs = cfg.get("inputs") or {}
    pep2cds_map = inputs.get("pep2cds_map", "")
    if not pep2cds_map:
        print("[ERR] config.yaml 缺少 inputs.pep2cds_map（pep2cds_resolved.tsv）")
        return 3
    map_path = Path(str(pep2cds_map)).expanduser()
    global OG2GID2SP
    OG2GID2SP = load_pep2cds_map(map_path)

    keep_list = QC_DIR / "keep_og.list"
    if not keep_list.exists():
        print(f"[ERR] 缺少 OG 列表：{keep_list}")
        return 3
    ogs = [x.strip() for x in keep_list.read_text(encoding="utf-8").splitlines() if x.strip()]
    if not ogs:
        print("[ERR] keep_og.list 为空"); return 3

    fgs = sorted([p.stem for p in SETS_DIR.glob("*.nwk")])
    if not fgs:
        print("[ERR] 未发现前景树：results/04_codeml/sets/*.nwk，请先运行 05_define_foregrounds.py")
        return 3

    # 映射表覆盖性预检（只检查 keep_og.list 里出现的 OG）
    missing_map_og = [og for og in ogs if og not in OG2GID2SP]
    if missing_map_og:
        print(f"[ERR] pep2cds_map 中缺少 OG 的映射：n={len(missing_map_og)} sample={missing_map_og[:10]}")
        return 3

    ensure_dir(RAW_ROOT); ensure_dir(LOG_ROOT)

    print(f"[INIT] {ts()} 06 启动  two-pass=ALT-then-NULL  threads={threads}")
    print(f"[INIT] CODON_ROOT={SEQ_DIR}")
    print(f"[PLAN] OG={len(ogs)} × FG={len(fgs)} × 2 (ALT/NULL)")

    totals = {
        "RUN(new)":0, "RECALC(clean)":0, "SKIP(resume)":0,
        "SOFT_OK":0, "FAIL":0, "TOTAL":0
    }

    for fg in fgs:
        run_pass(fg, "alt", ogs, threads, codeml_bin, totals)
        run_pass(fg, "null", ogs, threads, codeml_bin, totals)

    print("\n[SUMMARY]")
    for k in ["RUN(new)","RECALC(clean)","SKIP(resume)","SOFT_OK","FAIL","TOTAL"]:
        print(f"  {k:14s}= {totals[k]}")

    return 0 if totals["FAIL"]==0 else 1

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断"); sys.exit(130)

