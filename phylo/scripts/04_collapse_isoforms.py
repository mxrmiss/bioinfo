#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按物种折叠同一 OG 的多个 isoform：
- 输入：paths.normalize_dir 中的 OG*.trim.faa（表头须为 Species|SeqID）
- 规则：每物种仅保留最长序列；输出到 paths.collapse_dir
- 附件：results/reports/og_species_protein.tsv（OG,Species,SeqID,length）
- 写 .done 哨兵
"""
import os, sys
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"

class Fasta:
    @staticmethod
    def read(path):
        recs = []; name=None; seqs=[]
        with open(path, 'r', encoding='utf-8') as f:
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
        with open(path, 'w', encoding='utf-8') as f:
            for h,s in recs:
                f.write('>'+h+'\n')
                for i in range(0,len(s),60):
                    f.write(s[i:i+60]+'\n')

def load_yaml(p):
    import yaml
    with open(p,'r',encoding='utf-8') as f:
        return yaml.safe_load(f)

def get_config():
    cfg = os.environ.get('PHYLO_CONFIG', DEFAULT_CONFIG)
    if not os.path.exists(cfg):
        print(f"[ERR] 未找到配置文件：{cfg}", file=sys.stderr); sys.exit(2)
    return load_yaml(cfg)

def main():
    cfg=get_config(); P=cfg['paths']
    in_dir=P['normalize_dir']; out_dir=P['collapse_dir']; rep=P['reports_dir']
    os.makedirs(out_dir, exist_ok=True); os.makedirs(rep, exist_ok=True)

    files=[x for x in sorted(os.listdir(in_dir)) if x.endswith('.trim.faa') and x.startswith('OG')]
    if not files:
        print(f"[ERR] 未在 {in_dir} 发现 OG*.trim.faa", file=sys.stderr); sys.exit(3)

    tab=open(os.path.join(rep,'og_species_protein.tsv'),'w',encoding='utf-8')
    tab.write('OG\tSpecies\tSeqID\tlength\n')

    ok=0
    for fn in files:
        og=fn.split('.',1)[0]
        recs=Fasta.read(os.path.join(in_dir,fn))
        best={}  # sp -> (sid, seq)
        for h,s in recs:
            if '|' not in h:
                print(f"[ERR] 规范化前置未达成（缺 '|'): {h}", file=sys.stderr); sys.exit(4)
            sp,sid=h.split('|',1)
            if sp not in best or len(s)>len(best[sp][1]):
                best[sp]=(sid,s)
        out=[(f"{sp}|{sid}", seq) for sp,(sid,seq) in sorted(best.items())]
        Fasta.write(os.path.join(out_dir,fn), out)
        for sp,(sid,seq) in sorted(best.items()):
            tab.write(f"{og}\t{sp}\t{sid}\t{len(seq)}\n")
        ok+=1

    tab.close()
    if ok==0:
        print("[ERR] 折叠后无产物", file=sys.stderr); sys.exit(5)

    Path(os.path.join(out_dir,'.done')).touch()
    print(f"[DONE] 折叠完成：{ok} 个 → {out_dir}")

if __name__=='__main__':
    main()
