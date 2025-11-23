#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
小工具：用 tx2gene.clean.tsv 把第三方提供的 cds/transcript 列表映射成 gene_id

使用方式（示例）：
  1. 确保本脚本与 tx2gene.clean.tsv 在同一目录下。
  2. 把第三方给你的 cds/transcript 列表文件路径，填到 INPUT_CDS_FILE 变量中。
  3. 根据实际情况修改：输入文件是否有表头、ID 在第几列等参数。
  4. 运行：python map_cds_to_gene.py
  5. 输出：
       - 映射表：OUTPUT_MAPPED_FILE
       - 去重的 gene_id 列表：OUTPUT_GENE_LIST_FILE
"""

import csv
from pathlib import Path

# ===================== 参数区（皇上只要改这里就行） =====================

# 第三方提供的 cds/transcript 列表文件
# 可以是：
#   - 一列：每行一个 cds_id
#   - 多列：ID 在第 ID_COLUMN 列
INPUT_CDS_FILE = "third_party_cds.list"   # ← 皇上改成自己的文件名

# tx2gene.clean.tsv 文件名（与本脚本同目录）
TX2GENE_FILE = "tx2gene.clean.tsv"

# 输入文件是否有表头
INPUT_HAS_HEADER = False  # 若第一行为表头则改为 True

# cds/transcript ID 在第几列（从 1 开始计数）
ID_COLUMN = 1

# 输出文件：映射后的完整表（=原始内容 + 最后一列 gene_id）
OUTPUT_MAPPED_FILE = "cds_mapped_to_gene.tsv"

# 输出文件：仅包含去重后的 gene_id，一列，以便直接用作 gene 集合
OUTPUT_GENE_LIST_FILE = "cds_mapped_gene.list"

# tx2gene.clean.tsv 是否有表头（推荐 True）
TX2GENE_HAS_HEADER = True

# ===================== 下面不用改了 =====================


def load_tx2gene(tx2gene_path: Path, has_header: bool = True) -> dict:
    """
    读取 tx2gene.clean.tsv，建立 transcript_id -> gene_id 映射字典。
    默认认为前两列分别为 transcript_id, gene_id。
    """
    mapping = {}
    total_rows = 0
    duplicated_keys = 0

    if not tx2gene_path.exists():
        raise FileNotFoundError(f"找不到映射文件：{tx2gene_path}")

    with tx2gene_path.open("r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        if has_header:
            header = next(reader, None)
            if header is None or len(header) < 2:
                raise ValueError("tx2gene.clean.tsv 表头不合法，前两列应为 transcript_id, gene_id")
        for row in reader:
            if not row:
                continue
            if len(row) < 2:
                continue
            total_rows += 1
            tid = row[0].strip()
            gid = row[1].strip()
            if not tid:
                continue
            if tid in mapping and mapping[tid] != gid:
                # 同一个 transcript_id 对应不同 gene_id，记录一下，但以第一次为准
                duplicated_keys += 1
            else:
                mapping[tid] = gid

    print(f"[INFO] 已从 {tx2gene_path} 读入 {total_rows} 行，得到 {len(mapping)} 个唯一 transcript_id")
    if duplicated_keys > 0:
        print(f"[WARN] 发现 {duplicated_keys} 个 transcript_id 对应多个 gene_id，仅保留首次出现的映射")

    return mapping


def map_cds_to_gene(
    input_file: Path,
    tx2gene_map: dict,
    id_column: int = 1,
    has_header: bool = False,
    output_mapped: Path = Path("cds_mapped_to_gene.tsv"),
    output_gene_list: Path = Path("cds_mapped_gene.list"),
) -> None:
    """
    对输入文件进行逐行映射：
      - 保留原始所有列；
      - 在最后追加一列 gene_id（若找不到映射则留空）；
      - 同时收集所有成功映射到的 gene_id，输出去重列表。
    """
    if not input_file.exists():
        raise FileNotFoundError(f"找不到输入文件：{input_file}")

    # ID_COLUMN 从 1 开始计数，转成 0-based 索引
    col_idx = id_column - 1
    if col_idx < 0:
        raise ValueError("ID_COLUMN 必须 >= 1")

    mapped_gene_ids = set()
    total_lines = 0
    mapped_lines = 0
    unmapped_lines = 0

    with input_file.open("r", encoding="utf-8", errors="replace") as fin, \
            output_mapped.open("w", encoding="utf-8", newline="") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")

        # 处理表头
        if has_header:
            header = next(reader, None)
            if header is None:
                raise ValueError("输入文件声明有表头，但内容为空")
            # 追加 gene_id 列
            new_header = list(header) + ["gene_id"]
            writer.writerow(new_header)
        else:
            header = None

        # 逐行处理
        for row in reader:
            if not row:
                continue
            total_lines += 1
            # 防止这一行列数比声明的列数少
            if col_idx >= len(row):
                # 这一行在指定列没有 ID，视为 unmapped
                writer.writerow(list(row) + [""])
                unmapped_lines += 1
                continue

            cds_id = row[col_idx].strip()
            gene_id = tx2gene_map.get(cds_id, "")

            if gene_id:
                mapped_lines += 1
                mapped_gene_ids.add(gene_id)
            else:
                unmapped_lines += 1

            writer.writerow(list(row) + [gene_id])

    # 写出去重后的 gene_id 列表
    with output_gene_list.open("w", encoding="utf-8") as fgene:
        for gid in sorted(mapped_gene_ids):
            fgene.write(gid + "\n")

    print(f"[INFO] 输入文件行数（不含表头）：{total_lines}")
    print(f"[INFO] 成功映射行数：{mapped_lines}")
    print(f"[INFO] 未找到映射行数：{unmapped_lines}")
    print(f"[INFO] 输出映射表：{output_mapped}")
    print(f"[INFO] 输出去重 gene_id 列表：{output_gene_list}")
    if unmapped_lines > 0:
        print("[WARN] 存在未映射的 cds/transcript，请检查 ID 口径是否与 tx2gene.clean.tsv 一致")


def main():
    # 当前工作目录
    cwd = Path(".").resolve()

    tx2gene_path = cwd / TX2GENE_FILE
    input_path = cwd / INPUT_CDS_FILE
    output_mapped_path = cwd / OUTPUT_MAPPED_FILE
    output_gene_list_path = cwd / OUTPUT_GENE_LIST_FILE

    print(f"[INIT] 当前工作目录：{cwd}")
    print(f"[INIT] 使用 tx2gene 文件：{tx2gene_path}")
    print(f"[INIT] 输入 cds 列表文件：{input_path}")

    tx2gene_map = load_tx2gene(tx2gene_path, has_header=TX2GENE_HAS_HEADER)

    map_cds_to_gene(
        input_file=input_path,
        tx2gene_map=tx2gene_map,
        id_column=ID_COLUMN,
        has_header=INPUT_HAS_HEADER,
        output_mapped=output_mapped_path,
        output_gene_list=output_gene_list_path,
    )

    print("[DONE] 映射完成。")


if __name__ == "__main__":
    main()

