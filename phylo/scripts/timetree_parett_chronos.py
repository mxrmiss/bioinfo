#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
timetree_parett_chronos.py（“一脚本流 · TimeTree→treePL”）

目标：不再手工准备 chronos_calibrations.pruned.tsv；只要有 IQ-TREE 的 species_tree.nwk，
本脚本就会自动：
  1) 解析叶名→学名→NCBI 物种级 taxid（严格，不做近亲替代）
  2) 为每个内部节点挑选“两叶代表”并从 TimeTree 获取分化时间（pairwise 优先，mrca 两叶兜底）
  3) 自动确定根年龄（全叶 MRCA；或外群 vs 其余；或 ROOT_AGE_MANUAL）
  4) 生成 treePL 控制文件并拟合，输出超时树 + 节点年龄表

输出：
  - data/mcmctree/species_names_clean.nwk
  - results/mcmctree/node_ages.tsv
  - results/mcmctree/chronos_calibrations.pruned.tsv（仅两叶校准，便于审阅留痕）
  - results/mcmctree/treepl.ctl（可复跑留档）
  - logs/treepl.log（treePL 运行日志）
  - 缓存：results/mcmctree/cache_timetree_api/

依赖：
  - python: requests, ete3
  - 外部：treePL（conda: mamba install -y -c bioconda treepl）

参数规范（皇上长期偏好）：全部集中在脚本顶部；无需命令行参数。
"""

import os, re, json, time, csv, math, subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import requests
from ete3 import Tree

# ======================
# 参数区（可按需修改）
# ======================
TREE_IN      = "results/trees/species_tree.nwk"         # IQ-TREE 拓扑输入
OUT_TREE     = "data/mcmctree/species_names_clean.nwk"  # 输出超时树
NODE_AGES    = "results/mcmctree/node_ages.tsv"         # 节点年龄表
CALIB_TSV    = "results/mcmctree/chronos_calibrations.pruned.tsv"  # 仅两叶校准留痕
CTL_PATH     = "results/mcmctree/treepl.ctl"            # treePL 控制文件
LOG_PATH     = "logs/treepl.log"                        # treePL 日志
CACHE_DIR    = "results/mcmctree/cache_timetree_api"    # 缓存目录（JSON）

# 外群设置：优先使用 config.yaml 的外群（由 Snakefile 传参不再需要；这里手填或保持默认）
OUTGROUPS: List[str] = ["Nematostella_vectensis"]       # 可填多个叶名，精确匹配树上叶名（含下划线）

# 根年龄手工兜底（单位 Myr）：留空则尝试从 TimeTree 自动获取；仍失败时报错提醒
ROOT_AGE_MANUAL: Optional[float] = None    # 例：ROOT_AGE_MANUAL = 715.483

# HTTP 设置
TT_BASE   = "http://timetree.temple.edu/api"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TIMEOUT  = 12
SLEEP    = 0.4   # 轻微限速，减少 429

# treePL 调参（通常不改）
SMOOTH_START = 0.1
SMOOTH_END   = 1000
SMOOTH_STEPS = 8
THOROUGH     = 1

# 结束后是否自动清理缓存（True=移动到 results/mcmctree/cache_timetree_api.done；False=保留）
CLEAN_CACHE_AFTER = False

# ======================
# 工具函数
# ======================
def log(msg: str):
    print(time.strftime("[%H:%M:%S] "), msg, flush=True)

def ensure_dirs():
    Path(OUT_TREE).parent.mkdir(parents=True, exist_ok=True)
    Path(NODE_AGES).parent.mkdir(parents=True, exist_ok=True)
    Path(CTL_PATH).parent.mkdir(parents=True, exist_ok=True)
    Path(LOG_PATH).parent.mkdir(parents=True, exist_ok=True)
    Path(CACHE_DIR).mkdir(parents=True, exist_ok=True)

def leaf_to_sciname(leaf: str) -> str:
    # IQ-TREE 叶名形如 "Genus_species"；转为 "Genus species"
    return leaf.replace("_", " ").strip()

def cache_path(name: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.+-]", "_", name)
    return Path(CACHE_DIR) / safe

def http_get_json(url: str) -> Optional[dict]:
    time.sleep(SLEEP)
    r = requests.get(url, timeout=TIMEOUT, headers={"Accept": "application/json"})
    if r.status_code != 200:
        return None
    try:
        js = r.json()
        return js
    except Exception:
        return None

def ncbi_taxid_for_species(sci: str) -> Optional[int]:
    """
    严格解析到“物种级” taxid：先 SCIN 精确搜，再 ALL，再 esummary 判 rank=species。
    返回 int taxid 或 None。缓存两级：
      - esearch_SCIN_*.json / esearch_ALL_*.json
      - esummary_{tid}.json
    """
    nm = sci.replace(" ", "_")
    # 1) SCIN
    p1 = cache_path(f"ncbi_esearch_SCIN_{nm}.json")
    if p1.exists():
        js = json.loads(p1.read_text())
    else:
        url = f"{NCBI_EUTILS}/esearch.fcgi?db=taxonomy&retmode=json&term={requests.utils.quote(sci)}%5BSCIN%5D"
        log(f"[HTTP GET] {url}")
        js = http_get_json(url)
        if js: p1.write_text(json.dumps(js), encoding="utf-8")
    if js:
        ids = js.get("esearchresult", {}).get("idlist", [])
        for tid in ids:
            p2 = cache_path(f"ncbi_esummary_{tid}.json")
            if p2.exists():
                js2 = json.loads(p2.read_text())
            else:
                url = f"{NCBI_EUTILS}/esummary.fcgi?db=taxonomy&retmode=json&id={tid}"
                log(f"[HTTP GET] {url}")
                js2 = http_get_json(url)
                if js2: p2.write_text(json.dumps(js2), encoding="utf-8")
            if js2:
                rec = js2.get("result", {}).get(str(tid), {})
                rank = (rec.get("rank") or rec.get("Rank") or "").lower()
                if rank == "species":
                    return int(tid)

    # 2) ALL
    p3 = cache_path(f"ncbi_esearch_ALL_{nm}.json")
    if p3.exists():
        js = json.loads(p3.read_text())
    else:
        url = f"{NCBI_EUTILS}/esearch.fcgi?db=taxonomy&retmode=json&term={requests.utils.quote(sci)}"
        log(f"[HTTP GET] {url}")
        js = http_get_json(url)
        if js: p3.write_text(json.dumps(js), encoding="utf-8")
    if js:
        ids = js.get("esearchresult", {}).get("idlist", [])
        for tid in ids:
            p2 = cache_path(f"ncbi_esummary_{tid}.json")
            if p2.exists():
                js2 = json.loads(p2.read_text())
            else:
                url = f"{NCBI_EUTILS}/esummary.fcgi?db=taxonomy&retmode=json&id={tid}"
                log(f"[HTTP GET] {url}")
                js2 = http_get_json(url)
                if js2: p2.write_text(json.dumps(js2), encoding="utf-8")
            if js2:
                rec = js2.get("result", {}).get(str(tid), {})
                rank = (rec.get("rank") or rec.get("Rank") or "").lower()
                if rank == "species":
                    return int(tid)
    return None

def timetree_pair_age(tid1: int, tid2: int) -> Optional[float]:
    """
    TimeTree pairwise（JSON）拿分化时间；失败则 mrca 两叶兜底。
    返回 Myr 或 None。
    """
    # pairwise
    url = f"{TT_BASE}/pairwise/{tid1}/{tid2}/json"
    log(f"[HTTP GET] {url}")
    js = http_get_json(url)
    if js:
        # 多种字段容忍：precomputed_age 或 time 或 studies.precomputed_age
        age = None
        if isinstance(js, dict):
            age = js.get("precomputed_age") or js.get("time")
            if age is None and "studies" in js and isinstance(js["studies"], dict):
                age = js["studies"].get("precomputed_age") or js["studies"].get("time")
        if age is not None:
            try:
                return float(age)
            except Exception:
                pass
    # mrca 两叶
    url = f"{TT_BASE}/mrca/id/{tid1}+{tid2}/json"
    log(f"[HTTP GET] {url}")
    js = http_get_json(url)
    if js:
        age = None
        if isinstance(js, dict):
            age = js.get("precomputed_age") or js.get("time")
        if age is not None:
            try:
                return float(age)
            except Exception:
                pass
    return None

def choose_two_representatives(node) -> Tuple[str, str]:
    """
    从 node 的两个子支各挑一个叶名（叶名保持原始下划线形式）
    """
    ch = node.get_children()
    left = ch[0].get_leaves()
    right = ch[1].get_leaves()
    return left[0].name, right[0].name

def build_two_tip_calibrations(tr: Tree, leaf2tid: Dict[str,int]) -> List[Tuple[Tuple[str,str], float]]:
    """
    遍历所有二分内部节点：为每个节点挑两叶，打到 TimeTree 拿年龄。
    仅保留命中成功的两叶校准。
    """
    calibs = []
    i = 0
    for node in tr.traverse("postorder"):
        if node.is_leaf(): continue
        if len(node.get_children()) != 2: continue
        a, b = choose_two_representatives(node)
        if a not in leaf2tid or b not in leaf2tid:
            continue
        tid1, tid2 = leaf2tid[a], leaf2tid[b]
        i += 1
        age = timetree_pair_age(tid1, tid2)
        if age is not None:
            calibs.append(((a, b), age))
    return calibs

def dedup_pairs(pairs: List[Tuple[Tuple[str,str], float]]) -> List[Tuple[Tuple[str,str], float]]:
    seen = set()
    out = []
    for (a,b),age in pairs:
        key = tuple(sorted((a,b)))
        if key in seen: continue
        seen.add(key)
        out.append(((a,b), age))
    return out

def write_pruned_tsv(pairs: List[Tuple[Tuple[str,str], float]], path: str):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tips","age_mya","source"])
        for (a,b),age in pairs:
            w.writerow([f"{a},{b}", f"{age:.6f}", "PairwiseOrMRCA2"])

def root_age_from_timetree_all(tids: List[int]) -> Optional[float]:
    url = f"{TT_BASE}/mrca/id/{'+'.join(str(t) for t in tids)}/json"
    log(f"[HTTP GET] {url}")
    js = http_get_json(url)
    if js:
        age = js.get("precomputed_age") or js.get("time")
        if age is not None:
            try:
                return float(age)
            except Exception:
                pass
    return None

def root_age_from_outgroup_vs_rest(out_tids: List[int], other_tids: List[int]) -> Optional[float]:
    # 多外群时：取外群之一与其余融合；为稳妥，挑第一个外群
    if not out_tids or not other_tids: return None
    url = f"{TT_BASE}/mrca/id/{out_tids[0]}+{'+'.join(str(t) for t in other_tids)}/json"
    log(f"[HTTP GET] {url}")
    js = http_get_json(url)
    if js:
        age = js.get("precomputed_age") or js.get("time")
        if age is not None:
            try:
                return float(age)
            except Exception:
                pass
    return None

def write_treepl_ctl(treefile: str, outtree: str,
                     pairs: List[Tuple[Tuple[str,str], float]],
                     root_age: float, ctl_path: str):
    lines=[]
    lines.append(f"treefile = {treefile}")
    lines.append(f"outfile = {outtree}")
    lines.append(f"smoothcv = {SMOOTH_START} {SMOOTH_END} {SMOOTH_STEPS}")
    if THOROUGH:
        lines.append("thorough")
    # 所有两叶校准：mrca + fixage
    for i,((a,b),age) in enumerate(pairs, 1):
        nm = f"N{i}"
        lines.append(f"mrca {nm} {a} {b}")
        lines.append(f"fixage {nm} = {age:.6f}")
    # 根固定：沿用第一对命名 ROOT 定位（仅用于定义名字，不改变那对的 fixage）
    a0,b0 = pairs[0][0]
    lines.append(f"mrca ROOT {a0} {b0}")
    lines.append(f"fixage ROOT = {root_age:.6f}")
    lines.append("prime")
    lines.append("autocrossvalidation")
    lines.append("optimize")
    Path(ctl_path).write_text("\n".join(lines)+"\n", encoding="utf-8")

def run_treepl(ctl: str, log_path: str):
    with open(log_path, "w") as lg:
        proc = subprocess.run(["treePL", ctl], stdout=lg, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"[ERR] treePL 运行失败；详见 {log_path}")

def export_node_ages(tree_path: str, out_tsv: str, root_age: float):
    """
    读回 ultrametric 树，计算每个内部节点年龄：
      age(node) = root_age - dist(root -> node)
    叶节点也给出（等于其“现今”= 0 Myr），便于校核。
    """
    tr = Tree(tree_path, format=1)
    root = tr.get_tree_root()
    # root-to-leaf 距离（应相等≈root_age），用于比例校验
    dists = [tr.get_distance(root, lf) for lf in tr.get_leaves()]
    mean_rt = sum(dists)/len(dists) if dists else 0.0
    scale = (root_age/mean_rt) if mean_rt>0 else 1.0

    rows=[]
    nid=0
    for n in tr.traverse("postorder"):
        dist_root = tr.get_distance(root, n) * scale
        age = max(root_age - dist_root, 0.0)
        if n.is_leaf():
            tips = n.name
        else:
            # 记录此节点覆盖的叶集合（逗号拼接，便于比对）
            tips = ",".join(sorted([lf.name for lf in n.get_leaves()]))
        rows.append((tips, f"{age:.6f}", "0", "treePL"))
        nid+=1
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tips","age","fossil","source"])
        for r in rows:
            w.writerow(r)

def maybe_cleanup_cache():
    if not CLEAN_CACHE_AFTER:
        return
    dst = Path(str(CACHE_DIR) + ".done")
    if dst.exists():
        for p in dst.iterdir():
            pass
    Path(CACHE_DIR).rename(dst)

# ======================
# 主流程
# ======================
def main():
    ensure_dirs()
    if not Path(TREE_IN).exists():
        raise SystemExit(f"[ERR] 物种树不存在：{TREE_IN}")
    tr = Tree(TREE_IN, format=1)
    leaves = [lf.name for lf in tr.get_leaves()]
    log(f"[Stage 0] 读取拓扑： {TREE_IN}")
    log(f"[INFO] 叶子数 = {len(leaves)}")
    log(f"[INFO] 物种列表： {', '.join(leaves)}")

    # 1) 叶名 → 学名 → NCBI 物种级 taxid（严格）
    log("[Stage 1] 解析学名 → NCBI taxid")
    leaf2tid: Dict[str,int] = {}
    for i, leaf in enumerate(leaves, 1):
        sci = leaf_to_sciname(leaf)
        tid = ncbi_taxid_for_species(sci)
        if tid:
            log(f"[Leaf {i}/{len(leaves)}] {leaf} -> taxid={tid}")
            leaf2tid[leaf] = tid
        else:
            log(f"[MISS] {leaf} 未解析到物种级 taxid（严格模式）")
    if len(leaf2tid) < len(leaves):
        log(f"[WARN] 有 {len(leaves)-len(leaf2tid)} 个叶未命中 taxid，将在相应节点跳过校准对。")

    # 2) 为每个内部节点挑“两叶代表”，抓 TimeTree 年龄（pairwise→mrca 两叶）
    log("[Stage 2] 枚举内部节点并抓取 TimeTree 年龄（两叶优先 Pairwise）")
    pairs = build_two_tip_calibrations(tr, leaf2tid)
    pairs = dedup_pairs(pairs)
    log(f"[STAT] 全量两叶校准对 = {len(pairs)}")
    if not pairs:
        raise SystemExit("[ERR] 未获取到任何两叶校准；请检查网络/TimeTree/物种名解析。")

    # 3) 根年龄：全叶 MRCA；若失败则 外群 vs 其余；仍失败 → ROOT_AGE_MANUAL
    tids_all = [leaf2tid[x] for x in leaves if x in leaf2tid]
    root_age = None
    if len(tids_all) == len(leaves):
        root_age = root_age_from_timetree_all(tids_all)
    if root_age is None:
        outs = [x for x in OUTGROUPS if x in leaf2tid]
        rest = [x for x in leaves if x not in OUTGROUPS and x in leaf2tid]
        out_tids = [leaf2tid[x] for x in outs]
        other_tids = [leaf2tid[x] for x in rest]
        if out_tids and other_tids:
            root_age = root_age_from_outgroup_vs_rest(out_tids, other_tids)
    if root_age is None and ROOT_AGE_MANUAL is not None:
        root_age = float(ROOT_AGE_MANUAL)
    if root_age is None:
        raise SystemExit("[ERR] 根年龄获取失败：请在脚本顶部 ROOT_AGE_MANUAL 填写一个值（如 715.483）")
    log(f"[STAT] 根年龄 = {root_age:.6f} Myr")

    # 4) 写出两叶 TSV（留痕审阅）+ treePL 控制文件 + 运行
    write_pruned_tsv(pairs, CALIB_TSV)
    write_treepl_ctl(TREE_IN, OUT_TREE, pairs, root_age, CTL_PATH)
    log(f"[Stage 3] 运行 treePL：控制文件 {CTL_PATH}")
    run_treepl(CTL_PATH, LOG_PATH)
    log(f"[OK] treePL 完成；输出树：{OUT_TREE}")

    # 5) 导出节点年龄表（统一口径）
    export_node_ages(OUT_TREE, NODE_AGES, root_age)
    log(f"[DONE] 节点年龄表：{NODE_AGES}")

    # 6) 可选清理缓存
    maybe_cleanup_cache()

if __name__ == "__main__":
    main()

