#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_build_concat_clean_phy.py —— 由 codon MSA 目录拼接生成 concat.clean.phy（供 MCMCTree 使用）

输入契约（全部来自 config.yaml:mcmctree）：
  - codon_dir    : 03/04 阶段的 codon 对齐目录（文件名形如 OGXXXX.codon.fna/fa/fasta；表头=物种名）
  - keep_list    : 可选；若存在，仅保留该清单中的 OG（首列为 OGID）
  - codon_pos    : "12" 或 "123"（默认 "12"）
  - outdir       : 输出目录（通常=results/06_cafe/mcmctree）
  - tree / nobl_trees : 两行树（首行 N 1，第二行 newick）；仅用于**确定叶序**（不要求无分支长度）
产物：
  - <outdir>/concat.clean.phy    —— PHYLIP：多物种 × 拼接后的核酸序列（供 MCMCTree 的 seqfile）
  - logs/09_build_concat_clean_phy.log

设计要点：
  - 仅取 .trees 的**第二行 newick**，通过 stdin 喂给 gotree，调用：`gotree labels -i -`（不再使用 --type）
  - 对齐字符集正规化：非 A/C/G/T 统一置为 N（保持长度，避免 PAML 行为不确定）
  - 只纳入“覆盖全部物种”的 OG；可选 keep_og.list 过滤
  - PHYLIP 名称列固定宽度（默认 30）+ 至少两个空格；兼容 PAML 4.10.9
  - 失败不写“完成”，日志标注 [FAIL]；成功路径才写“完成”
"""

from __future__ import annotations
from pathlib import Path
from datetime import datetime
import subprocess, re, sys

# ===================== 顶部参数（可被 config.yaml:mcmctree 覆盖） =====================
CONFIG_PATH          = Path("config.yaml")
CODON_DIR_DEFAULT    = Path("results/03_codon/codon_msa")
KEEP_LIST_DEFAULT    = Path("results/03_qc/keep_og.list")
CODON_POS_DEFAULT    = "12"   # 允许 "12" 或 "123"
OUTDIR_DEFAULT       = Path("results/06_cafe/mcmctree")

# 兼容 08 的命名：优先读取 config 中 mcmctree.tree / nobl_trees；若均未设置则从候选路径中自动探测
NOBL_TREES_CANDIDATES = [
    Path("results/06_cafe/mcmctree/species_calib.trees"),       # 两行，有分支长度
    Path("results/06_cafe/mcmctree/species_calib.nobl.trees"),  # 老命名（两行，无分支长度）
]

LOG_PATH             = Path("logs/09_build_concat_clean_phy.log")
PHYLIP_NAME_WIDTH    = 30      # 名称列固定宽度，保证与 PAML/PHYLIP 的兼容性
# =====================================================================================

# ----------------------------- YAML 读取（ruamel 优先） -----------------------------
def load_yaml_dict(p: Path):
    if not p.exists(): return {}
    try:
        from ruamel.yaml import YAML
        y = YAML(typ='rt'); y.preserve_quotes = True
        return y.load(p.read_text(encoding="utf-8")) or {}
    except Exception:
        try:
            import yaml
            return yaml.safe_load(p.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}

# -------------------------------- 简易日志器 --------------------------------
class Logger:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.f = path.open("w", encoding="utf-8", newline="\n")
        self.ok = False  # 成功标记
        self.write(f"[{datetime.now()}] 09_build_concat_clean_phy 启动")
    def write(self, s: str):
        print(s)
        self.f.write(s + ("\n" if not s.endswith("\n") else ""))
        self.f.flush()
    def success(self):
        self.ok = True
    def close(self):
        self.write(f"[{datetime.now()}] {'完成' if self.ok else 'FAIL'}")
        self.f.close()

# --------------------------------- 工具函数 ---------------------------------
def require_gotree():
    import shutil
    if shutil.which("gotree") is None:
        raise SystemExit("[ERR] 未找到 gotree；请安装：mamba install -y -c bioconda gotree")

def read_two_line_newick_second(trees_path: Path) -> str:
    """读取两行 .trees，仅返回第二行 newick 字符串（去首尾空白；确保以 ; 结尾）。"""
    lines = [ln.strip() for ln in trees_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(lines) < 2:
        raise SystemExit(f"[ERR] {trees_path} 不是“两行树”文件（首行 N 1，第二行 newick）")
    nwk = lines[1].strip()
    if not nwk.endswith(";"):
        nwk += ";"
    return nwk

def ordered_species_from_trees(trees_path: Path) -> list[str]:
    """
    仅把第二行 newick 通过 stdin 喂给 gotree，稳定提取**按树顺序**的叶名。
    使用 `gotree labels -i -`（tips 为默认；部分版本也可显式 `--tips`）。
    """
    nwk = read_two_line_newick_second(trees_path)
    p = subprocess.run(
        ["gotree", "labels", "-i", "-"],  # 兼容你环境的 labels 接口；不使用 --type
        input=nwk, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if p.returncode != 0:
        err = (p.stderr or "").strip()
        raise SystemExit(f"[ERR] gotree labels 失败：{err or '（无stderr）'}")
    names = [ln.strip() for ln in p.stdout.splitlines() if ln.strip()]
    if not names:
        raise SystemExit("[ERR] 无法从 .trees 第二行解析到叶名，请检查 newick 是否正常")
    return names

def read_fasta(path: Path) -> dict[str, str]:
    """读取 FASTA：键=header 第一段（物种名），值=序列（大写，去空白）"""
    d = {}; name=None; buf=[]
    with path.open("r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf).upper().replace(" ","")
                name = ln[1:].strip().split()[0]; buf=[]
            else:
                buf.append(ln.strip())
    if name is not None:
        d[name] = "".join(buf).upper().replace(" ","")
    return d

def sanitize_dna(seq: str) -> str:
    """非 A/C/G/T 统一置为 N（保持长度）"""
    return re.sub(r'[^ACGT]', 'N', seq.upper())

def pick_pos(seq: str, mode: str) -> str:
    """选择密码子位点；'12' 取前两位，'123' 保留全部。"""
    if mode == "123":
        return seq
    out=[]; n=len(seq)
    # 以 3 为步长安全切片，末端不足 3 的尾巴丢弃（不影响阅读框）
    for i in range(0, n, 3):
        cod = seq[i:i+3]
        if len(cod) == 3:
            out.append(cod[:2])
    return "".join(out)

def list_codon_files(codon_dir: Path) -> list[Path]:
    """枚举 codon MSA 文件，支持 .fna/.fa/.fasta；按 OGID 去重。"""
    files=[]
    for pat in ("OG*.codon.fna","OG*.codon.fa","OG*.codon.fasta"):
        files += list(codon_dir.glob(pat))
    uniq={}
    for p in sorted(files):
        ogid = p.stem.split(".")[0]  # OG0001234.codon.fna -> OG0001234
        if ogid not in uniq:
            uniq[ogid] = p
    return [uniq[k] for k in sorted(uniq.keys())]

def phylip_write_line(name: str, seq: str, width: int = PHYLIP_NAME_WIDTH) -> str:
    """
    生成一行 PHYLIP：固定宽度名字列 + 至少两个空格 + 序列。
    - 过长截断到 width；过短右侧填充空格到 width
    - 始终追加两个空格作为分隔
    """
    safe = (name[:width]).ljust(width, ' ')
    return f"{safe}  {seq}\n"

# ----------------------------------- 主流程 -----------------------------------
def main():
    log = Logger(LOG_PATH)
    try:
        # 读取配置
        cfg = load_yaml_dict(CONFIG_PATH)
        m = (cfg.get("mcmctree") or {})

        CODON_DIR = Path(m.get("codon_dir", CODON_DIR_DEFAULT))
        KEEP_LIST = Path(m.get("keep_list", KEEP_LIST_DEFAULT))
        CODON_POS = str(m.get("codon_pos", CODON_POS_DEFAULT))
        OUTDIR    = Path(m.get("outdir", OUTDIR_DEFAULT))

        # 树路径解析：优先 mcmctree.tree / nobl_trees；否则从候选路径中自动探测
        tree_from_cfg = m.get("tree") or m.get("nobl_trees")
        if tree_from_cfg:
            NOBL_TREES = Path(tree_from_cfg)
        else:
            cand = [p for p in NOBL_TREES_CANDIDATES if p.exists()]
            if not cand:
                NOBL_TREES = NOBL_TREES_CANDIDATES[0]  # 给出一个合理的默认提示路径
            else:
                NOBL_TREES = cand[0]

        # 路径存在性与参数校验
        if not CODON_DIR.exists():
            raise SystemExit(f"[ERR] 未找到 codon MSA 目录：{CODON_DIR}")
        if not NOBL_TREES.exists():
            raise SystemExit(f"[ERR] 未找到两行树文件（首行 N 1，第二行 newick）：{NOBL_TREES}")
        if CODON_POS not in {"12","123"}:
            raise SystemExit(f"[ERR] codon_pos 应为 '12' 或 '123'，当前={CODON_POS}")

        OUTDIR.mkdir(parents=True, exist_ok=True)
        LOG_PATH.parent.mkdir(parents=True, exist_ok=True)

        log.write("[PATH] codon_dir   = " + str(CODON_DIR))
        log.write("[PATH] keep_list   = " + (str(KEEP_LIST) if KEEP_LIST.exists() else "(无)"))
        log.write("[PATH] trees(two-line) = " + str(NOBL_TREES))
        log.write("[PATH] outdir      = " + str(OUTDIR))

        # gotree 依赖
        require_gotree()

        # 从 .trees 第二行提取**保序**的叶名
        species = ordered_species_from_trees(NOBL_TREES)
        log.write(f"[OK] 叶序（{len(species)}）：")
        for i in range(0, len(species), 5):
            log.write("  - " + ", ".join(species[i:i+5]))

        # 枚举 OG 文件
        files = list_codon_files(CODON_DIR)
        if not files:
            raise SystemExit(f"[ERR] {CODON_DIR} 中未找到 OG*.codon.fna/fa/fasta")

        # 可选 keep_og.list 过滤
        if KEEP_LIST.exists():
            keep = {ln.strip().split()[0] for ln in KEEP_LIST.read_text(encoding="utf-8").splitlines() if ln.strip()}
            before = len(files)
            files = [p for p in files if p.stem.split(".")[0] in keep]
            log.write(f"[INFO] keep_og.list 过滤：{before} -> {len(files)}")
            if not files:
                raise SystemExit("[ERR] keep_og.list 过滤后无 OG，请检查 03/04 输出或清单")

        species_set = set(species)
        per_og: dict[str, dict[str,str]] = {}
        union_species = set()

        pep2cds_path = None
        try:
            I = (cfg.get("inputs") or {})
            pep2cds_raw = I.get("pep2cds_map")
            if isinstance(pep2cds_raw, str):
                publish_dir = cfg.get("publish_dir")
                if isinstance(publish_dir, str) and "<publish_dir>" in pep2cds_raw:
                    pep2cds_raw = pep2cds_raw.replace("<publish_dir>", publish_dir)
                pep2cds_path = Path(pep2cds_raw)
        except Exception:
            pep2cds_path = None

        og2gid2sp = None

        def _norm_id(x: str) -> str:
            x = (x or "").strip().split()[0]
            if "|" in x:
                x = x.split("|")[-1]
            return x

        def _load_og2gid2sp(map_path: Path) -> dict[str, dict[str, str]]:
            import csv
            with map_path.open("r", encoding="utf-8") as f:
                r = csv.reader(f, delimiter="\t")
                header = next(r, None)
                if not header:
                    raise SystemExit(f"[ERR] pep2cds_resolved.tsv 为空：{map_path}")
                h = [c.strip().lower() for c in header]

                def _idx(cands):
                    for c in cands:
                        if c in h:
                            return h.index(c)
                    return None

                i_og  = _idx(["og","orthogroup","family"])
                i_sp  = _idx(["species","taxon","sp"])
                i_pep = _idx(["protein_id","pep_id","protein","pep"])
                if i_og is None or i_sp is None or i_pep is None:
                    raise SystemExit("[ERR] pep2cds_resolved.tsv 缺少必要列（OG/species/protein_id）")

                out: dict[str, dict[str, str]] = {}
                for row in r:
                    if (not row) or (len(row) <= max(i_og, i_sp, i_pep)):
                        continue
                    og = row[i_og].strip()
                    sp = row[i_sp].strip()
                    pep = _norm_id(row[i_pep])
                    if not og or not sp or not pep:
                        continue
                    d = out.setdefault(og, {})
                    if pep in d and d[pep] != sp:
                        raise SystemExit(f"[ERR] pep2cds_resolved.tsv 冲突：OG={og} {pep} -> {d[pep]} vs {sp}")
                    d[pep] = sp
            return out

        # 读取并挑选“全覆盖”的 OG
        for p in files:
            msa = read_fasta(p)
            if not msa:
                continue

            if not species_set.issubset(msa.keys()):
                if pep2cds_path is None or (not pep2cds_path.exists()):
                    raise SystemExit("[ERR] codon MSA 表头不是物种名，但 inputs.pep2cds_map 未指向有效的 pep2cds_resolved.tsv")
                if og2gid2sp is None:
                    og2gid2sp = _load_og2gid2sp(pep2cds_path)
                ogid_tmp = p.stem.split(".")[0]
                gid2sp = og2gid2sp.get(ogid_tmp) or {}
                if not gid2sp:
                    raise SystemExit(f"[ERR] pep2cds_resolved.tsv 中找不到该 OG：{ogid_tmp}")
                msa2: dict[str, str] = {}
                for h, seq in msa.items():
                    k = _norm_id(h)
                    sp = gid2sp.get(k)
                    if not sp:
                        raise SystemExit(f"[ERR] OG={ogid_tmp} FASTA header 无法映射到 species：{h}")
                    if sp in msa2:
                        raise SystemExit(f"[ERR] OG={ogid_tmp} 映射后同一 species 出现多条序列：{sp}")
                    msa2[sp] = seq
                msa = msa2

            union_species.update(msa.keys())
            if species_set.issubset(msa.keys()):
                ogid = p.stem.split(".")[0]
                # 保序采样 + 位点选择（此处先不 sanitize，统一在拼接后处理）
                per_og[ogid] = {sp: pick_pos(msa[sp], CODON_POS) for sp in species}

        if not per_og:
            miss = sorted(species_set - union_species)
            log.write("[DIAG] 没有任何 OG 覆盖全部物种；排查建议：")
            log.write(" - 检查 03/04 阶段 fasta 表头是否为“仅物种名”（不应包含 | 或 geneID）")
            log.write(" - 是否漏跑 03/04；或 keep_og.list 过严（首列必须是 OGID）")
            log.write(" - 完全缺失的物种： " + (', '.join(miss) if miss else "（无明显缺失）"))
            raise SystemExit("[ERR] 无可用 OG，无法拼接 PHYLIP")

        # 拼接
        concat: dict[str, list[str]] = {sp: [] for sp in species}
        ogids_sorted = sorted(per_og.keys())
        log.write(f"[OK] 纳入 OG 数量：{len(ogids_sorted)}")
        if ogids_sorted:
            show_head = ", ".join(ogids_sorted[:3])
            show_tail = ", ".join(ogids_sorted[-3:]) if len(ogids_sorted) > 3 else ""
            log.write(f"[INFO] 例：首尾 OG = {show_head}" + (f" ... {show_tail}" if show_tail else ""))

        for og in ogids_sorted:
            msa = per_og[og]
            for sp in species:
                concat[sp].append(msa[sp])

        # 序列拼接 + 字符集正规化
        for sp in species:
            concat[sp] = sanitize_dna("".join(concat[sp]))

        # 一致性自检：各物种拼接后的长度必须一致
        lengths = {sp: len(concat[sp]) for sp in species}
        if len(set(lengths.values())) != 1:
            bad = "\n  ".join(f"{sp}: {L}" for sp, L in lengths.items())
            raise SystemExit("[ERR] 拼接后长度不一致，请排查上游对齐或 codon_pos 选取：\n  " + bad)

        aln_len = next(iter(lengths.values()))
        out_phy = OUTDIR / "concat.clean.phy"
        with out_phy.open("w", encoding="utf-8") as f:
            f.write(f"{len(species)} {aln_len}\n")
            for sp in species:
                f.write(phylip_write_line(sp, concat[sp], PHYLIP_NAME_WIDTH))

        log.write(f"[OK] 写入 {out_phy.resolve()} —— {len(species)} taxa × {aln_len} nt（codon_pos={CODON_POS}）")
        log.write("     已满足 10_mcmctree_run_and_publish.py 对 seqfile 的输入契约。")
        log.write("[DONE] 09 完成")
        log.success()
    finally:
        log.close()

if __name__ == "__main__":
    main()

