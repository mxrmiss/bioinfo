#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
覆盖度矩阵报告（零猜测）：
- 输入：paths.supermatrix_dir/supermatrix.faa
- 输出：paths.reports_dir/matrix.tsv + .matrix.done
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
    sup_dir=P['supermatrix_dir']; rep=P['reports_dir']
    os.makedirs(rep, exist_ok=True)
    sup=os.path.join(sup_dir,'supermatrix.faa')
    if not os.path.exists(sup):
        print(f"[ERR] 缺少 supermatrix：{sup}", file=sys.stderr); sys.exit(3)

    recs=Fasta.read(sup)
    if not recs:
        print("[ERR] 超矩阵为空", file=sys.stderr); sys.exit(4)

    L=len(recs[0][1])
    out=os.path.join(rep,'matrix.tsv')
    with open(out,'w',encoding='utf-8') as f:
        f.write('species\ttotal_sites\tnongap_sites\toccupancy_percent\n')
        for h,s in recs:
            nongap=sum(1 for c in s if c!='-')
            occ= (nongap*100.0/L) if L>0 else 0.0
            f.write(f"{h}\t{L}\t{nongap}\t{occ:.2f}\n")

    Path(os.path.join(rep,'.matrix.done')).touch()
    print(f"[DONE] 覆盖度矩阵：{out}")

if __name__=='__main__':
    main()
