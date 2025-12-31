#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基于 filelist 的超矩阵拼接（零猜测）：
- 只读 reports/sco_filelist.txt（绝对路径）
- 阈值：min_align_len / min_taxa_per_locus / top_n_partitions（按长度排序取前 N）
- 输出：supermatrix.faa + partitions.txt + .supermatrix.done
"""
import os, sys
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"

class Fasta:
    @staticmethod
    def read(path):
        recs=[]; name=None; seqs=[]
        with open(path,'r',encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    if name is not None:
                        recs.append((name, ''.join(seqs)))
                    name=line[1:].strip(); seqs=[]
                else:
                    seqs.append(line.strip())
        if name is not None:
            recs.append((name, ''.join(seqs)))
        return recs
    @staticmethod
    def write(path, recs):
        with open(path,'w',encoding='utf-8') as f:
            for h,s in recs:
                f.write('>'+h+'\n')
                for i in range(0,len(s),60):
                    f.write(s[i:i+60]+'\n')

def load_yaml(p):
    import yaml
    with open(p,'r',encoding='utf-8') as f:
        return yaml.safe_load(f)

def get_config():
    cfg=os.environ.get('PHYLO_CONFIG', DEFAULT_CONFIG)
    if not os.path.exists(cfg):
        print(f"[ERR] 未找到配置文件：{cfg}", file=sys.stderr); sys.exit(2)
    return load_yaml(cfg)

def main():
    cfg=get_config(); P=cfg['paths']
    C=cfg['species_tree']['concat']
    rep=P['reports_dir']; sup=P['supermatrix_dir']
    os.makedirs(sup, exist_ok=True)

    filelist=os.path.join(rep,'sco_filelist.txt')
    if not os.path.exists(filelist):
        print(f"[ERR] 缺少 filelist：{filelist}", file=sys.stderr); sys.exit(3)

    min_len=int(C.get('min_align_len',0))
    min_taxa=int(C.get('min_taxa_per_locus',0))
    topN=int(C.get('top_n_partitions',0))

    loci=[]  # (og, recs)
    species=set()

    with open(filelist,'r',encoding='utf-8') as f:
        for line in f:
            p=line.strip()
            if not p: continue
            og=os.path.basename(p).split('.',1)[0]
            recs=Fasta.read(p)
            if not recs: continue
            L=len(recs[0][1])
            if min_len and L<min_len: continue
            taxa=set(h for h,_ in recs)
            if min_taxa and len(taxa)<min_taxa: continue
            loci.append((og, recs))
            species |= taxa

    if not loci:
        print("[ERR] 过滤后无可用分区", file=sys.stderr); sys.exit(4)

    loci.sort(key=lambda x: len(x[1][0][1]), reverse=True)
    if topN and len(loci)>topN:
        loci=loci[:topN]

    species=sorted(species)
    concat={sp: [] for sp in species}
    parts=[]  # (og, start, end)

    offset=0
    for og, recs in loci:
        m = {h:s for h,s in recs}
        L = len(recs[0][1])
        for sp in species:
            seq = m.get(sp, '-'*L)
            concat[sp].append(seq)
        parts.append((og, offset+1, offset+L))
        offset += L

    out_faa=os.path.join(sup,'supermatrix.faa')
    Fasta.write(out_faa, [(sp, ''.join(concat[sp])) for sp in species])

    with open(os.path.join(sup,'partitions.txt'),'w',encoding='utf-8') as g:
        for og,a,b in parts:
            g.write(f"{og}\t{a}\t{b}\n")

    Path(os.path.join(sup,'.supermatrix.done')).touch()
    print(f"[DONE] 超矩阵完成：{len(loci)} 分区，物种 {len(species)}，长度 {offset} → {out_faa}")

if __name__=='__main__':
    main()
