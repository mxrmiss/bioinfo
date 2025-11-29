#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
helper_generate_pep2cds.py (V6 自动表头版)
功能：
  扫描 data/proteomes 和 data/cds 目录，自动生成映射关系。
  
核心修正：
  1. 在写入 pep2cds.tsv 之前，强制清洗掉 'gnl|' 前缀。
  2. 自动检查输出文件：如果文件不存在或为空，自动写入表头 (species\tprotein_id\tcds_id)。
  3. 自动创建父目录 (data/maps)，防止目录缺失报错。
"""

import os, sys
from pathlib import Path

# ================= 配置区 =================
PEP_DIR = "data/proteomes"
CDS_DIR = "data/cds"
MAP_FILE = "data/maps/pep2cds.tsv"

# 支持的后缀
EXTS = (".fa", ".faa", ".fna", ".fasta")
# =========================================

def clean_id_string(header_str):
    """
    清洗 ID 字符串，移除数据库前缀
    例如：'gnl|WGS:AMQN|CAPTEDRAFT_mRNA210255' -> 'WGS:AMQN|CAPTEDRAFT_mRNA210255'
    """
    # 1. 如果以 gnl| 开头，去掉它
    if header_str.startswith("gnl|"):
        header_str = header_str[4:]
    
    # 2. 如果被 rename 过了 (Species|gnl|...)，中间夹着 gnl|，也去掉
    if "|gnl|" in header_str:
        header_str = header_str.replace("|gnl|", "|")
        
    return header_str

def get_headers(fasta_path):
    """
    读取 FASTA，返回原始 ID 列表
    """
    headers = []
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                # 取第一个空格前的部分作为原始 ID
                raw_id = line.strip()[1:].split()[0]
                headers.append(raw_id)
    return headers

def load_existing_species(map_file):
    """读取已有的物种，避免重复生成"""
    existing_species = set()
    if not os.path.exists(map_file):
        return existing_species
        
    with open(map_file, 'r', encoding='utf-8') as f:
        for line in f:
            # 跳过表头和空行
            if line.startswith("species") or not line.strip(): continue
            parts = line.split('\t')
            if len(parts) > 0:
                existing_species.add(parts[0])
    return existing_species

def main():
    if not os.path.exists(PEP_DIR) or not os.path.exists(CDS_DIR):
        print(f"[ERR] 目录不存在：{PEP_DIR} 或 {CDS_DIR}")
        sys.exit(1)

    # 1. 加载已有记录
    recorded_species = load_existing_species(MAP_FILE)
    print(f"[INFO] 已有记录的物种：{len(recorded_species)} 个")

    # 2. 扫描蛋白目录
    pep_files = [f for f in os.listdir(PEP_DIR) if f.endswith(EXTS)]
    
    new_entries = []
    
    print(f"[INFO] 开始扫描外部数据...")
    for p_file in pep_files:
        # 推断物种名
        stem = Path(p_file).stem
        while Path(stem).suffix: stem = Path(stem).stem
        species = stem
        
        if species in recorded_species:
            continue

        # 寻找对应的 CDS 文件
        c_file = None
        for ext in EXTS:
            cand = os.path.join(CDS_DIR, species + ext)
            if os.path.exists(cand):
                c_file = cand
                break
        
        if not c_file:
            print(f"[WARN] {species}: 找到蛋白文件但未找到同名 CDS 文件，跳过。")
            continue

        # 读取原始 ID (不做任何修改)
        p_path = os.path.join(PEP_DIR, p_file)
        pep_ids_raw = get_headers(p_path)
        
        # 将 CDS 原始 ID 存为集合
        cds_lookup = {} # { cleaned_id : original_id }
        for cid in get_headers(c_file):
            cleaned = clean_id_string(cid)
            cds_lookup[cleaned] = cid

        # 匹配
        matched_count = 0
        
        for pid in pep_ids_raw:
            # 1. 对蛋白质 ID 进行清洗 (作为匹配键)
            pid_clean = clean_id_string(pid)
            
            # 剥离 Species| 前缀逻辑
            if '|' in pid_clean and not pid_clean.startswith('WGS:'): 
                try:
                    if pid_clean.startswith(species + "|"):
                         pid_clean = pid_clean.split('|', 1)[1]
                except:
                    pass

            # 2. 尝试匹配
            if pid_clean in cds_lookup:
                final_id = pid_clean 
                # 写入：species \t 清洗后ID \t 清洗后ID
                new_entries.append(f"{species}\t{final_id}\t{final_id}")
                matched_count += 1
        
        if matched_count == 0:
            print(f"[WARN] {species}: 未找到任何 ID 匹配。")
        else:
            print(f"[ADD] {species}: 新增 {matched_count} 条映射")

    # 3. 写入文件 (核心修改部分)
    if new_entries:
        # 确保输出目录存在
        Path(MAP_FILE).parent.mkdir(parents=True, exist_ok=True)
        
        # 判断是否需要写表头：文件不存在 OR 文件为空
        need_header = False
        if not os.path.exists(MAP_FILE):
            need_header = True
        elif os.path.getsize(MAP_FILE) == 0:
            need_header = True
            
        with open(MAP_FILE, 'a', encoding='utf-8') as f:
            if need_header:
                print(f"[INFO] 创建新文件并写入表头：{MAP_FILE}")
                f.write("species\tprotein_id\tcds_id\n")
            
            for line in new_entries:
                f.write(line + "\n")
                
        print(f"[DONE] 成功追加 {len(new_entries)} 条记录到 {MAP_FILE}")
    else:
        print("[INFO] 没有新的数据需要添加。")

if __name__ == "__main__":
    main()
