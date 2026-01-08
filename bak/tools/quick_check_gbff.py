#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quick_check_gbff.py
用于快速检查一个 NCBI 的 .gbff 文件是否“有希望”被可靠地转换为 GFF3。

检查内容：
  1. 能否被 Biopython 正常解析；
  2. 是否包含常见的注释特征（gene / CDS 等）；
  3. 简单统计记录数、特征数量，并给出文字结论。

使用方式（重要说明）：
  * 按皇上的特殊规定：本脚本“唯一允许”的命令行参数就是 gbff 文件名；
  * 其它参数（如最小记录数、需要检查的特征类型）必须在脚本顶部手动修改。

示例：
    python quick_check_gbff.py GCA_977019795.1_xbRudPhil1.hap2.1_genomic.gbff
"""

# =========================
# 1. 参数区（手动填写）
# =========================

# 认为“有希望转换为 GFF3”所需的最小记录数（通常 1 就够，整条基因组也可能就 1 条 record）
MIN_RECORDS = 1

# 在 gbff 中至少需要看到的特征类型（常见为 gene、CDS）
REQUIRED_FEATURE_TYPES = ["gene", "CDS"]


# =========================
# 2. 正文逻辑
# =========================

import os
import sys
import argparse

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write(
        "【错误】未找到 Biopython（Bio 模块）。\n"
        "请先在当前环境安装 Biopython，例如：\n"
        "  mamba install -c conda-forge biopython\n"
        "或：\n"
        "  conda install -c conda-forge biopython\n"
    )
    sys.exit(1)


def quick_check_gbff(path: str):
    """
    对指定 gbff 文件进行快速体检，打印总结信息。
    """
    if not os.path.exists(path):
        print(f"【错误】找不到 gbff 文件：{path}")
        return

    print("==========================================")
    print("  GBFF 文件快速体检报告")
    print("==========================================")
    print(f"文件路径：{os.path.abspath(path)}")

    record_count = 0
    feature_type_counts = {}  # 各类 feature 的计数
    required_found = {ftype: 0 for ftype in REQUIRED_FEATURE_TYPES}

    try:
        # 逐条读取 GenBank 记录
        for record in SeqIO.parse(path, "genbank"):
            record_count += 1
            for feat in getattr(record, "features", []):
                ftype = getattr(feat, "type", None)
                if not ftype:
                    continue
                feature_type_counts[ftype] = feature_type_counts.get(ftype, 0) + 1
                if ftype in required_found:
                    required_found[ftype] += 1
    except Exception as e:
        print("")
        print("【结果】❌ 解析失败")
        print("原因可能是：")
        print("  - 文件并不是标准的 GenBank/GBFF 格式；")
        print("  - 文件损坏或下载不完整；")
        print("  - Biopython 版本过旧，无法处理当前格式。")
        print(f"具体异常信息：{repr(e)}")
        return

    print("")
    print("【一、基本情况】")
    print(f"  记录（record）数量：{record_count}")

    if record_count == 0:
        print("")
        print("【结果】❌ 未读到任何记录，文件可能不是有效的 GBFF。")
        return

    print("")
    print("【二、特征类型统计（feature type）】")
    # 按出现次数降序输出前若干类特征
    for ftype, cnt in sorted(feature_type_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  - {ftype}: {cnt} 条")

    print("")
    print("【三、关键特征检查】")
    all_required_ok = True
    for ftype in REQUIRED_FEATURE_TYPES:
        cnt = required_found.get(ftype, 0)
        if cnt > 0:
            print(f"  ✔ 找到 {cnt} 条 '{ftype}' 特征：OK")
        else:
            print(f"  ✘ 未找到任何 '{ftype}' 特征：可能会影响后续 GFF3 转换和基因级分析")
            all_required_ok = False

    print("")
    print("【四、综合结论】")
    if record_count < MIN_RECORDS:
        print(
            f"  ❌ 记录数仅为 {record_count} (< {MIN_RECORDS})，"
            "不符合预期的基因组 GBFF 文件结构，建议重新检查文件来源。"
        )
    elif not all_required_ok:
        print(
            "  ⚠ 虽然文件可以被解析，但缺少部分关键特征类型（如 gene/CDS），\n"
            "    说明该 GBFF 可能只包含部分注释或非常简化的标注。\n"
            "    若后续仍要转换为 GFF3，建议先确认：\n"
            "      - 该物种是否有更完整的注释版本；\n"
            "      - 或者是否只需要序列层面的信息。"
        )
        print("")
        print("  【结论】该 GBFF 文件“可以技术上转换为 GFF3”，但注释信息不完整，")
        print("          不一定适合作为完整功能注释 / 共线性分析的唯一来源。")
    else:
        print(
            "  ✅ 文件可被正常解析，且包含常见的 gene/CDS 等注释特征。\n"
            "    从技术角度看，这个 GBFF 文件“适合转换为 GFF3”，\n"
            "    后续可以编写/使用转换脚本（例如基于 Biopython 的自定义脚本）\n"
            "    生成下游流程所需的 GFF3 注释文件。"
        )


def parse_args() -> str:
    """
    解析命令行参数。

    按皇上的特别要求：
      - 这里只允许传入“一个位置参数”：gbff 文件名/路径；
      - 其它阈值（如 MIN_RECORDS、REQUIRED_FEATURE_TYPES）必须在脚本顶部手动修改。
    """
    parser = argparse.ArgumentParser(
        description=(
            "快速检查指定的 .gbff 文件是否适合转换为 GFF3。"
        )
    )
    parser.add_argument(
        "gbff_file",
        help="需要检查的 GBFF 文件路径（唯一允许的命令行参数）",
    )
    args = parser.parse_args()
    return args.gbff_file


if __name__ == "__main__":
    gbff_path = parse_args()
    quick_check_gbff(gbff_path)

