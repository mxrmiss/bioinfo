#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
仅保留物种名作为表头（供拼接/树构建/aphylo）：
- 输入：paths.collapse_dir 的 OG*.trim.faa（Species|SeqID）
- 输出：paths.species_collapse_dir 的 OG*.trim.faa（仅 >Species）
- 写 .done 哨兵
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
    in_dir=P['collapse_dir']; out_dir=P['species_collapse_dir']
    os.makedirs(out_dir, exist_ok=True)

    files=[x for x in sorted(os.listdir(in_dir)) if x.endswith('.trim.faa') and x.startswith('OG')]
    if not files:
        print(f"[ERR] 未在 {in_dir} 发现 OG*.trim.faa", file=sys.stderr); sys.exit(3)

    n=0
    for fn in files:
        recs=Fasta.read(os.path.join(in_dir,fn))
        out=[]
        for h,s in recs:
            if '|' not in h:
                print(f"[ERR] 规范化前置未达成（缺 '|'): {h}", file=sys.stderr); sys.exit(4)
            sp=h.split('|',1)[0]
            out.append((sp, s))
        Fasta.write(os.path.join(out_dir,fn), out)
        n+=1

    if n==0:
        print("[ERR] species-only 产物为空", file=sys.stderr); sys.exit(5)

    Path(os.path.join(out_dir,'.done')).touch()
    print(f"[DONE] species-only 完成：{n} 个 → {out_dir}")

if __name__=='__main__':
    main()
