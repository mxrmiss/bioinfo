#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_prepare_trimmed_alignment.py  —— Biopython 版（含可选列掩码与哨兵）

改动点（最小化改动）：
  * 当 config.trimming.save_colmask = true 时，生成的掩码改为
    **AA 位串（长度=裁剪后 AA 列数的 0/1 串）**，每列均为 '1'（表示裁剪后的列全部保留）。
    文件内容仅一行 0/1 串 + 换行，不再写注释或索引，避免后续统计混淆。

其它逻辑保持不变：
  1) 读取 OrthoFinder 单基因 AA 对齐，使用 trimal 裁剪到 paths.msa_trim_dir
  2) 可选生成列掩码到 paths.colmask_dir
  3) 写哨兵：msa_trim_dir/.done；（若保存掩码）colmask_dir/.done

依赖：
  - trimal 可执行（binaries.trimal）
  - Biopython（pip/conda 安装 biopython）
"""

import os, sys, subprocess
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"

def load_yaml(p):
    import yaml
    with open(p, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def read_alignment_fasta(path):
    """用 Biopython 读取等长多序列对齐（FASTA），返回 (ids, seqs_str)"""
    from Bio import AlignIO
    aln = AlignIO.read(path, "fasta")  # 若不是等长多序列对齐会抛错
    ids = [rec.id for rec in aln]
    seqs = ["".join(rec.seq) for rec in aln]
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise RuntimeError(f"[ERR] 非等长对齐：{path}")
    return ids, seqs

def compute_colmask_indices_by_id(orig_ids, orig_seqs, trim_ids, trim_seqs):
    """
    给定原始与裁剪后的对齐（均为等长），按裁剪后的顺序对原始序列重排，
    逐列比较得到被保留的原始列索引（1-based）。
    仅用于一致性校验；最终掩码输出为 AA 位串（长度=L1，全部为 '1'）。
    """
    id2idx = {sid: i for i, sid in enumerate(orig_ids)}
    try:
        orig_sub = [orig_seqs[id2idx[sid]] for sid in trim_ids]
    except KeyError as e:
        miss = str(e).strip("'")
        raise RuntimeError(f"[ERR] 裁剪后出现了原始对齐中不存在的序列：{miss}")

    L0 = len(orig_sub[0])
    L1 = len(trim_seqs[0])
    if any(len(x) != L0 for x in orig_sub) or any(len(x) != L1 for x in trim_seqs):
        raise RuntimeError("[ERR] 输入对齐非等长，请检查源/裁剪对齐")

    kept = []
    j = 0
    for i in range(L0):
        if j >= L1:
            break
        same = True
        for k in range(len(orig_sub)):
            if orig_sub[k][i] != trim_seqs[k][j]:
                same = False
                break
        if same:
            kept.append(i + 1)  # 1-based
            j += 1

    if j != L1:
        raise RuntimeError(
            f"[ERR] 列映射失败：匹配到 {len(kept)} 列，但裁剪后有 {L1} 列；"
            "通常是源/裁剪对齐序列集合或顺序不一致造成的。"
        )
    return kept, L0, L1

def main():
    # 载入配置
    if not os.path.exists(DEFAULT_CONFIG):
        print(f"[ERR] 找不到配置文件 {DEFAULT_CONFIG}", file=sys.stderr); sys.exit(2)
    cfg = load_yaml(DEFAULT_CONFIG)

    paths = cfg.get('paths', {})
    bins = cfg.get('binaries', {})
    trimming_cfg = cfg.get('trimming', {})
    reports_dir = paths.get('reports_dir', 'results/reports')
    of_root    = paths.get('orthofinder_results_dir', 'results/orthofinder')
    msa_trim_dir = paths.get('msa_trim_dir', 'results/msa_trim')
    colmask_dir  = paths.get('colmask_dir', 'results/colmask')
    save_colmask = bool(trimming_cfg.get('save_colmask', False))
    of_suffix = cfg.get('input', {}).get('of_msa_suffix', '.fa')
    keep_list = os.path.join(reports_dir, 'ogs_selected.list')

    # 锚定 OrthoFinder 结果
    results_txt = os.path.join(of_root, 'RESULTS_DIR.txt')
    if not os.path.exists(results_txt):
        print(f"[ERR] 缺少 RESULTS_DIR.txt：{results_txt}", file=sys.stderr); sys.exit(3)
    with open(results_txt, 'r', encoding='utf-8') as f:
        results_dir = f.readline().strip()
    msa_src = os.path.join(results_dir, 'MultipleSequenceAlignments')
    if not os.path.isdir(msa_src):
        print(f"[ERR] 找不到 OF 对齐目录：{msa_src}", file=sys.stderr); sys.exit(4)

    # 读严格 SCO 列表
    if not os.path.exists(keep_list):
        print(f"[ERR] 缺少严格SCO清单：{keep_list}", file=sys.stderr); sys.exit(5)
    with open(keep_list, 'r', encoding='utf-8') as f:
        ogs = [x.strip() for x in f if x.strip()]

    os.makedirs(msa_trim_dir, exist_ok=True)
    if save_colmask:
        os.makedirs(colmask_dir, exist_ok=True)

    trimal = bins.get('trimal', 'trimal')
    trimmed_cnt = 0
    colmask_cnt = 0

    for og in ogs:
        src = os.path.join(msa_src, og + of_suffix)
        if not os.path.exists(src):
            print(f"[ERR] 缺少 MSA 源文件：{src}", file=sys.stderr); sys.exit(6)
        dst = os.path.join(msa_trim_dir, f"{og}.trim.faa")

        # 调 trimal
        r = subprocess.run([trimal, '-in', src, '-out', dst, '-automated1'],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0 or (not os.path.exists(dst)) or os.path.getsize(dst) == 0:
            sys.stderr.write(r.stderr)
            print(f"[ERR] trimal 失败：{src}", file=sys.stderr); sys.exit(7)
        trimmed_cnt += 1

        # 可选：产出列掩码（AA 位串；长度=裁剪后 AA 列数；全部为 '1'）
        if save_colmask:
            try:
                # 仍做一次映射校验，确保源/裁剪一致；但最终掩码输出为位串
                orig_ids, orig_seqs = read_alignment_fasta(src)
                tri_ids,  tri_seqs  = read_alignment_fasta(dst)
                _, L0, L1 = compute_colmask_indices_by_id(orig_ids, orig_seqs, tri_ids, tri_seqs)

                # 生成位串掩码（与裁剪后 AA 列一致，全部保留为 '1'）
                mask_bits = "1" * L1

                # 只写 0/1 位串一行，避免任何注释或数字干扰下游统计
                with open(os.path.join(colmask_dir, f"{og}.colmask"), 'w', encoding='utf-8') as f:
                    f.write(mask_bits + "\n")
                colmask_cnt += 1
            except Exception as e:
                print(str(e), file=sys.stderr); sys.exit(8)

    # === 哨兵文件（在所有循环完成后写入）===
    Path(os.path.join(msa_trim_dir, '.done')).touch()
    if save_colmask:
        Path(os.path.join(colmask_dir, '.done')).touch()

    msg = f"[DONE] 对齐裁剪完成：{trimmed_cnt} 个 → {msa_trim_dir}"
    if save_colmask:
        msg += f"；列掩码生成（AA位串）：{colmask_cnt} 个 → {colmask_dir}"
    print(msg)

if __name__ == '__main__':
    # 友好报错：缺 Biopython 时提示安装
    try:
        from Bio import AlignIO  # 预检
    except Exception:
        print("[ERR] 需要 Biopython：请安装 biopython（conda 或 pip）", file=sys.stderr)
        sys.exit(1)
    main()

