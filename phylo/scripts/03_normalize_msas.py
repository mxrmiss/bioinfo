#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
表头规范化（零猜测）：
- 输入：paths.msa_trim_dir 下的 OG*.trim.faa
- 要求：输入表头必须已带分隔符“Species|SeqID”，否则报错（不再“推断物种名”）
- 可选：species.alias_tsv 或 alias_map 将左侧 Species 统一为标准名
- 输出：paths.normalize_dir 下同名文件；并写 reports/header_map.tsv；生成 .norm.done
"""
import os, sys
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"


class Fasta:
    @staticmethod
    def read(path):
        recs = []
        name = None
        seqs = []
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    if name is not None:
                        recs.append((name, ''.join(seqs)))
                    name = line[1:].strip()
                    seqs = []
                else:
                    seqs.append(line.strip())
        if name is not None:
            recs.append((name, ''.join(seqs)))
        return recs

    @staticmethod
    def write(path, recs):
        with open(path, 'w', encoding='utf-8') as f:
            for h, s in recs:
                f.write('>' + h + '\n')
                for i in range(0, len(s), 60):
                    f.write(s[i:i+60] + '\n')


def load_yaml(p):
    import yaml
    with open(p, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def get_config():
    cfg = os.environ.get('PHYLO_CONFIG', DEFAULT_CONFIG)
    if not os.path.exists(cfg):
        print(f"[ERR] 未找到配置文件：{cfg}", file=sys.stderr)
        sys.exit(2)
    return load_yaml(cfg)


def load_alias(cfg):
    sp = cfg.get('species', {})
    if sp.get('alias_tsv'):
        mp = {}
        with open(sp['alias_tsv'], 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                a, b = line.split('\t')[:2]
                mp[a] = b
        return mp
    mp = sp.get('alias_map', None)
    if not mp:
        return None
    return mp


def canon_species(sp_raw: str) -> str:
    """
    规范化物种名以去除重复前缀：
    - 若形如 A_B_A_B（用下划线分隔，且左右两半完全一致），折叠为 A_B
    - 递归处理，直到不再可折叠
    """
    s = sp_raw.strip()
    # 最多递归 3 次，防御极端重复
    for _ in range(3):
        toks = s.split('_')
        if len(toks) % 2 != 0 or len(toks) == 0:
            break
        half = len(toks) // 2
        left = '_'.join(toks[:half])
        right = '_'.join(toks[half:])
        if left == right and left:  # 左右两半完全一致
            s = left
        else:
            break
    return s


def main():
    cfg = get_config()
    P = cfg['paths']
    reports = P['reports_dir']
    os.makedirs(P['normalize_dir'], exist_ok=True)
    os.makedirs(reports, exist_ok=True)

    alias = load_alias(cfg)

    in_dir = P['msa_trim_dir']
    out_dir = P['normalize_dir']

    files = [x for x in sorted(os.listdir(in_dir)) if x.endswith('.trim.faa') and x.startswith('OG')]
    if not files:
        print(f"[ERR] 未在 {in_dir} 发现 OG*.trim.faa", file=sys.stderr)
        sys.exit(3)

    hmap_path = os.path.join(reports, 'header_map.tsv')
    with open(hmap_path, 'w', encoding='utf-8') as hmap:
        hmap.write('OG\told_header\tnew_header\tSpecies\tSeqID\n')

        converted = 0
        for fn in files:
            og = fn.split('.', 1)[0]
            recs = Fasta.read(os.path.join(in_dir, fn))
            out_recs = []
            for hdr, seq in recs:
                if '|' not in hdr:
                    print(f"[ERR] 表头缺少分隔符'|'：{hdr}（零猜测模式不再自动推断物种名）", file=sys.stderr)
                    sys.exit(4)
                sp_raw, sid = hdr.split('|', 1)

                # ① 去除 OF 带来的重复物种前缀（如 A_B_A_B）
                sp_norm = canon_species(sp_raw)

                # ② 可选：别名映射成标准名
                if alias:
                    sp_norm = alias.get(sp_norm, sp_norm)

                new_h = f"{sp_norm}|{sid}"
                out_recs.append((new_h, seq))
                hmap.write(f"{og}\t{hdr}\t{new_h}\t{sp_norm}\t{sid}\n")

            Fasta.write(os.path.join(out_dir, fn), out_recs)
            converted += 1

    if converted == 0:
        print("[ERR] 无任何文件完成规范化", file=sys.stderr)
        sys.exit(5)

    Path(os.path.join(out_dir, '.norm.done')).touch()
    print(f"[DONE] 表头规范化完成：{converted} 个 → {out_dir}")


if __name__ == '__main__':
    main()
