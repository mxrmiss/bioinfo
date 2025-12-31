#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
严格按白名单筛选 MSA，生成唯一的 filelist（绝对路径）：
- 输入：paths.species_collapse_dir 与 reports/ogs_selected.list
- 输出：reports/sco_filelist.txt（绝对路径，每行一个 OG*.trim.faa），以及 .sco_filelist.done
- 不再创建软链，不再扫描候选目录
"""
import os, sys
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"

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
    spc=P['species_collapse_dir']; rep=P['reports_dir']
    keep=os.path.join(rep,'ogs_selected.list')
    outl=os.path.join(rep,'sco_filelist.txt')
    os.makedirs(rep, exist_ok=True)

    if not os.path.isdir(spc):
        print(f"[ERR] species-only 目录不存在：{spc}", file=sys.stderr); sys.exit(3)
    if not os.path.exists(keep):
        print(f"[ERR] 缺少严格 SCO 列表：{keep}", file=sys.stderr); sys.exit(4)

    n_hit=0
    with open(keep,'r',encoding='utf-8') as fin, open(outl,'w',encoding='utf-8') as fout:
        for line in fin:
            og=line.strip()
            if not og: continue
            p=os.path.join(spc, og+'.trim.faa')
            if not os.path.exists(p):
                print(f"[WARN] 丢失 OG 文件（跳过）：{p}")
                continue
            fout.write(os.path.abspath(p)+'\n')
            n_hit+=1

    if n_hit==0:
        print("[ERR] sco_filelist.txt 为空（OG 命名或输入目录不匹配）", file=sys.stderr); sys.exit(5)

    Path(os.path.join(rep,'.sco_filelist.done')).touch()
    print(f"[DONE] 命中 {n_hit} 个 OG → {outl}")

if __name__=='__main__':
    main()
