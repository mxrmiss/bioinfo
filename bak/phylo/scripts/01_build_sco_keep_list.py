#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
构建严格 SCO 白名单（零猜测）：
- 只读 paths.orthofinder_results_dir/RESULTS_DIR.txt 锚点
- 优先使用 Orthogroups_SingleCopyOrthologues.txt；
  若缺失，则从 Orthogroups.tsv 精确筛选“每个物种恰一条”的 OG
- 输出：ogs_selected.list、ogs_policy.json、.ogs_selected.done
"""
import os, sys, json
from pathlib import Path

DEFAULT_CONFIG = "config.yaml"

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

def read_results_anchor(out_dir:str)->str:
    anchor = os.path.join(out_dir, 'RESULTS_DIR.txt')
    if not os.path.exists(anchor):
        print(f"[ERR] 未找到 RESULTS_DIR.txt：{anchor}", file=sys.stderr)
        sys.exit(3)
    return open(anchor, 'r', encoding='utf-8').read().strip()

def parse_singlecopy_txt(p):
    ogs = []
    with open(p, 'r', encoding='utf-8') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            og = line.split(':',1)[0].strip()
            if og: ogs.append(og)
    return ogs

def parse_sco_from_tsv(p):
    ogs = []
    with open(p, 'r', encoding='utf-8') as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            arr = line.rstrip('\n').split('\t')
            og = arr[0]
            genes = arr[1:]
            ok = True
            for cell in genes:
                if cell == '':
                    ok = False; break
                parts = [x for x in cell.split(',') if x]
                if len(parts) != 1:
                    ok = False; break
            if ok:
                ogs.append(og)
    return ogs

def main():
    cfg = get_config()
    P = cfg['paths']
    reports = P['reports_dir']
    os.makedirs(reports, exist_ok=True)

    of_out = P['orthofinder_results_dir']
    R = read_results_anchor(of_out)

    sc_txt = os.path.join(R, 'Orthogroups', 'Orthogroups_SingleCopyOrthologues.txt')
    og_tsv = os.path.join(R, 'Orthogroups', 'Orthogroups.tsv')

    ogs = []
    policy = {"source": None}
    if os.path.exists(sc_txt):
        ogs = parse_singlecopy_txt(sc_txt)
        policy["source"] = "Orthogroups_SingleCopyOrthologues.txt"
    elif os.path.exists(og_tsv):
        ogs = parse_sco_from_tsv(og_tsv)
        policy["source"] = "Orthogroups.tsv (strict single-copy by columns)"
    else:
        print(f"[ERR] 缺失 OrthoGroups 文件：{sc_txt} / {og_tsv}", file=sys.stderr)
        sys.exit(4)

    if not ogs:
        print("[ERR] 未筛出任何严格 SCO", file=sys.stderr)
        sys.exit(5)

    out_list = os.path.join(reports, 'ogs_selected.list')
    with open(out_list, 'w', encoding='utf-8') as f:
        for og in ogs:
            f.write(og + '\n')
    with open(os.path.join(reports, 'ogs_policy.json'), 'w', encoding='utf-8') as f:
        json.dump(policy, f, ensure_ascii=False, indent=2)

    Path(os.path.join(reports, '.ogs_selected.done')).touch()
    print(f"[DONE] 严格SCO选择完成：{len(ogs)} 个OG → {out_list}")

if __name__ == '__main__':
    main()
