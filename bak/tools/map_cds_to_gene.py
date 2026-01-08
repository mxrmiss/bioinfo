#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
小工具：用 tx2gene.clean.tsv 把第三方提供的 cds/transcript 列表映射成 gene_id

使用方式（示例）：
- 与脚本同路径下要有results/03_maps/tx2gene.clean.tsv

  1. 确保本脚本与 tx2gene.clean.tsv 在同一目录下（脚本目录）。
  2. 把第三方给你的 cds/transcript 列表文件路径，填到 INPUT_CDS_FILE 变量中，
     或在命令行中显式传入一个或多个输入文件。
  3. 根据实际情况修改：输入文件是否有表头、ID 在第几列等参数。
  4. 运行方式：
       - 老模式（使用参数区）：
           python map_cds_to_gene.py
       - 单文件（命令行）：
           python map_cds_to_gene.py file1
       - 单文件 + 备份（命令行）：
           python map_cds_to_gene.py file1 -c
       - 多文件（命令行）：
           python map_cds_to_gene.py file1 file2 file3
       - 多文件 + 备份（命令行）：
           python map_cds_to_gene.py file1 file2 file3 -c
  5. 输出规则：
       - 产物文件名的前缀一律与输入文件前缀（去掉扩展名后的部分）相同。
       - 老模式：
           * 不修改原输入文件；
           * 映射表：脚本目录 / <前缀>.mapped.tsv
           * gene 列表：脚本目录 / <前缀>.list
       - 命令行模式：
           * 若输入文件与脚本在同一目录：
               - 映射表：覆盖原输入文件（可选 .bak 备份）
               - gene 列表：脚本目录 / <前缀>.list
           * 若输入文件与脚本不在同一目录：
               - 原输入文件保持不变；
               - 映射表：脚本目录 / <前缀>.mapped.tsv
               - gene 列表：脚本目录 / <前缀>.list
"""

import csv
import sys
import shutil
from pathlib import Path

# ===================== 参数区（皇上只要改这里就行） =====================

# 第三方提供的 cds/transcript 列表文件
# 可以是：
#   - 一列：每行一个 cds_id
#   - 多列：ID 在第 ID_COLUMN 列
INPUT_CDS_FILE = "third_party_cds.list"   # ← 皇上改成自己的文件名（支持相对/绝对路径）

# tx2gene.clean.tsv 文件名（与本脚本同目录）
TX2GENE_FILE = "tx2gene.clean.tsv"

# 输入文件是否有表头
INPUT_HAS_HEADER = False  # 若第一行为表头则改为 True

# cds/transcript ID 在第几列（从 1 开始计数）
ID_COLUMN = 1

# 下面两个名字仅作占位，实际输出文件名会根据输入文件前缀自动生成：
#   <前缀>.mapped.tsv, <前缀>.list
OUTPUT_MAPPED_FILE = "cds_mapped_to_gene.tsv"
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


def _parse_cli_args():
    """
    简单解析命令行参数：
      - 无参数：返回 (None, False) —— 使用参数区 INPUT_CDS_FILE（老模式）
      - 有参数：
          python map_cds_to_gene.py file1 [file2 ... fileN]
          python map_cds_to_gene.py file1 [file2 ... fileN] -c
        其中：
          * 最后一个参数若为 -c，则表示对所有文件进行备份
    返回：
      ( [file1, file2, ...] 或 None, backup_flag_bool )
    """
    args = sys.argv[1:]
    if not args:
        return None, False

    backup = False
    if args[-1] == "-c":
        backup = True
        file_args = args[:-1]
    else:
        file_args = args

    # 不允许在中间出现 -c
    if any(a == "-c" for a in file_args):
        print("用法: python map_cds_to_gene.py file1 [file2 ...] [-c]", file=sys.stderr)
        sys.exit(1)

    if not file_args:
        print("用法: python map_cds_to_gene.py file1 [file2 ...] [-c]", file=sys.stderr)
        sys.exit(1)

    return file_args, backup


def main():
    # 脚本所在目录（tx2gene 与所有产物统一放在这里）
    script_dir = Path(__file__).resolve().parent
    # 当前工作目录（运行 python 命令的地方）
    cwd = Path(".").resolve()

    cli_files, backup_flag = _parse_cli_args()

    # tx2gene 固定从脚本目录读取
    tx2gene_path = script_dir / TX2GENE_FILE

    print(f"[INIT] 当前工作目录：{cwd}")
    print(f"[INIT] 脚本所在目录：{script_dir}")
    print(f"[INIT] 使用 tx2gene 文件：{tx2gene_path}")

    tx2gene_map = load_tx2gene(tx2gene_path, has_header=TX2GENE_HAS_HEADER)

    if cli_files is not None:
        # —— 命令行模式：一个或多个显式输入 —— 
        input_paths = []
        for name in cli_files:
            p = Path(name)
            if p.is_absolute():
                input_paths.append(p.resolve())
            else:
                input_paths.append((cwd / p).resolve())

        # -c：对所有输入文件备份 .bak（在各自目录下）
        if backup_flag:
            for p in input_paths:
                bak_path = p.with_name(p.name + ".bak")
                shutil.copy2(p, bak_path)
                print(f"[INIT] 备份原文件为：{bak_path}")

        for p in input_paths:
            same_dir = (p.parent == script_dir)
            prefix = p.stem

            print(f"[INIT] 正在处理输入 cds 列表文件：{p}")
            print(f"[INIT] 输入文件与脚本同目录？ {'是' if same_dir else '否'}")
            print(f"[INIT] 当前前缀：{prefix}")

            # 映射表输出路径：
            if same_dir:
                # 同一目录：覆盖原文件
                output_mapped_path = p
            else:
                # 不同目录：保留原文件，映射表写在脚本目录
                output_mapped_path = script_dir / f"{prefix}.mapped.tsv"

            # 基因列表一律写在脚本目录，前缀一致，后缀为 .list
            output_gene_list_path = script_dir / f"{prefix}.list"

            map_cds_to_gene(
                input_file=p,
                tx2gene_map=tx2gene_map,
                id_column=ID_COLUMN,
                has_header=INPUT_HAS_HEADER,
                output_mapped=output_mapped_path,
                output_gene_list=output_gene_list_path,
            )
    else:
        # —— 老模式：使用参数区 INPUT_CDS_FILE —— 
        input_path = Path(INPUT_CDS_FILE)
        if not input_path.is_absolute():
            input_path = (cwd / input_path).resolve()

        prefix = input_path.stem
        output_mapped_path = script_dir / f"{prefix}.mapped.tsv"
        output_gene_list_path = script_dir / f"{prefix}.list"

        print(f"[INIT] 输入 cds 列表文件：{input_path}")
        print(f"[INIT] 映射表输出到脚本目录：{output_mapped_path}")
        print(f"[INIT] gene 列表输出到脚本目录：{output_gene_list_path}")
        print(f"[INIT] 当前前缀：{prefix}")

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

