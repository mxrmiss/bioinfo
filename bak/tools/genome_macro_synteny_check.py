#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genome_macro_synteny_check.py
用于快速评估一个基因组是否适合作为“宏观染色体共线性图”的材料。

使用方式（命令行）：
    python genome_macro_synteny_check.py your_genome.fa

本脚本会：
  1. 统计每条序列长度，区分“明显特别长的序列”和“零碎序列”；
  2. 计算长序列数量、长度占比、N50 等指标；
  3. 在当前目录生成一个文本报告：<基因组前缀>.macro_synteny_report.txt；
  4. 给出一句结论：推荐 / 勉强可用 / 不推荐 用于宏观染色体共线性图。
"""

import argparse
import gzip
import os
from typing import List, Tuple

# ==============================
#  参数区（可根据需要手动修改）
# ==============================

# 定义“明显特别长的序列”的长度阈值（单位：bp）
# 一般宏观染色体级别的 scaffold / pseudochromosome 通常在 10 Mb 以上
LONG_LEN_THRESHOLD = 10_000_000  # 皇上可根据物种大小略微调整

# 长序列覆盖全基因组的最低比例：
#   比如 >= 0.8 表示：>=80% 的碱基都集中在这些长 scaffold 上
MIN_LONG_COVERAGE_STRICT = 0.80  # 达到则判定为“推荐用于宏观共线性图”
MIN_LONG_COVERAGE_LOOSE = 0.60   # 介于宽松阈值和严格阈值之间则为“勉强可用”

# 期望的长序列条数范围（可选，仅用于给出提醒，不作为硬性条件）
# 比如 2n=38 的贝类，单倍体染色体条数约 19 条，因此可以设置 (10, 30)
EXPECTED_LONG_COUNT_MIN = 10
EXPECTED_LONG_COUNT_MAX = 30

# 报告文件名后缀（会自动接在基因组文件名前缀后）
REPORT_SUFFIX = ".macro_synteny_report.txt"


def read_fasta_lengths(path: str) -> List[int]:
    """
    从 FASTA / FASTA.GZ 文件中读取所有序列长度。

    参数
    ----
    path : str
        基因组 FASTA 文件路径。

    返回
    ----
    lengths : List[int]
        每条序列的长度列表。
    """
    # 根据后缀名判断是否使用 gzip 打开
    if path.endswith(".gz"):
        fh = gzip.open(path, "rt")
    else:
        fh = open(path, "r")

    lengths: List[int] = []
    current_len = 0

    with fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                # 遇到新的序列头，把上一条的长度先存进去
                if current_len > 0:
                    lengths.append(current_len)
                    current_len = 0
            else:
                current_len += len(line.strip())
        # 文件结束后，把最后一条序列长度也存进去
        if current_len > 0:
            lengths.append(current_len)

    return lengths


def calc_n50(lengths: List[int]) -> int:
    """
    计算 N50（基于序列长度列表）。

    参数
    ----
    lengths : List[int]
        序列长度列表（单位：bp）。

    返回
    ----
    n50 : int
        N50 值（单位：bp），若列表为空则返回 0。
    """
    if not lengths:
        return 0

    lengths_sorted = sorted(lengths, reverse=True)
    total = sum(lengths_sorted)
    running = 0
    for L in lengths_sorted:
        running += L
        if running >= total * 0.5:
            return L
    return lengths_sorted[-1]


def classify_genome(
    total_len: int,
    n_long: int,
    len_long: int,
    n_short: int,
    len_short: int,
) -> Tuple[str, str]:
    """
    根据统计结果给出宏观染色体共线性适用性的结论。

    返回
    ----
    level : str
        结论等级： "A_推荐", "B_勉强可用", "C_不推荐"。
    message : str
        中文解释说明。
    """
    if total_len == 0:
        return (
            "C_不推荐",
            "基因组总长度为 0，可能 FASTA 文件为空或格式有误，无法用于任何分析。",
        )

    long_cov = len_long / total_len

    # 是否处于期望的长序列条数范围内
    count_in_range = (
        EXPECTED_LONG_COUNT_MIN <= n_long <= EXPECTED_LONG_COUNT_MAX
    )

    if long_cov >= MIN_LONG_COVERAGE_STRICT and count_in_range:
        level = "A_推荐"
        msg = (
            "长序列覆盖度较高（%.1f%%），且长序列条数处于合理范围内，"
            "可以作为“pseudochromosome 级”材料，用于绘制宏观染色体共线性图。"
        ) % (long_cov * 100.0)
    elif long_cov >= MIN_LONG_COVERAGE_STRICT:
        level = "A_推荐"
        msg = (
            "长序列覆盖度较高（%.1f%%），可以视为“准染色体级”基因组，"
            "适合作为宏观染色体共线性图的主材料；"
            "但长序列条数（%d 条）略偏离预期范围（%d–%d 条），"
            "建议在方法部分说明染色体条数的不确定性。"
        ) % (long_cov * 100.0, n_long, EXPECTED_LONG_COUNT_MIN, EXPECTED_LONG_COUNT_MAX)
    elif long_cov >= MIN_LONG_COVERAGE_LOOSE:
        level = "B_勉强可用"
        msg = (
            "长序列覆盖度为 %.1f%%，说明大部分序列已经被锚定到较长 scaffold 上，"
            "但仍有较多碎片序列存在。建议在绘制宏观染色体共线性图时："
            "只使用最长的若干条 scaffold 作为“伪染色体”，"
            "并在正文/方法中说明组装仍然存在一定碎片化。"
        ) % (long_cov * 100.0)
    else:
        level = "C_不推荐"
        msg = (
            "长序列覆盖度仅为 %.1f%%，大部分碱基分布在较短的碎片 scaffold 中，"
            "更适合用于基因/局部共线性分析，而不适合绘制顶层的宏观染色体共线性图。"
        ) % (long_cov * 100.0)

    return level, msg


def format_bp(x: int) -> str:
    """
    将碱基数格式化为带单位的字符串（Mb 精度）。
    """
    if x >= 1_000_000_000:
        return f"{x/1_000_000_000:.2f} Gb"
    elif x >= 1_000_000:
        return f"{x/1_000_000:.2f} Mb"
    elif x >= 1_000:
        return f"{x/1_000:.2f} Kb"
    else:
        return f"{x} bp"


def main():
    # ===========================================================
    #  解析命令行参数
    #  重要说明（按皇上要求的“唯一例外”）：
    #  —— 本脚本中，只有“基因组 FASTA 文件名”允许通过命令行传入；
    #  —— 所有阈值（长序列长度、覆盖比例等）都必须在脚本顶部手动修改。
    # ===========================================================
    parser = argparse.ArgumentParser(
        description=(
            "统计基因组中长序列与碎片序列的数量与长度，"
            "评估是否适合作为宏观染色体共线性图的材料。"
        )
    )
    parser.add_argument(
        "genome_fasta",
        help="基因组 FASTA 文件路径（可为 .fa/.fasta 或 .gz 压缩文件）",
    )
    args = parser.parse_args()

    genome_path = args.genome_fasta
    if not os.path.exists(genome_path):
        raise FileNotFoundError(f"找不到基因组文件：{genome_path}")

    # 读取序列长度
    lengths = read_fasta_lengths(genome_path)
    if not lengths:
        raise RuntimeError("未在 FASTA 中解析到任何序列，请检查文件格式。")

    total_len = sum(lengths)
    n_total = len(lengths)

    # 按长度降序排序，方便查看
    lengths_sorted = sorted(lengths, reverse=True)

    # 统计长序列 / 碎片序列
    long_lengths = [L for L in lengths_sorted if L >= LONG_LEN_THRESHOLD]
    short_lengths = [L for L in lengths_sorted if L < LONG_LEN_THRESHOLD]

    n_long = len(long_lengths)
    len_long = sum(long_lengths)

    n_short = len(short_lengths)
    len_short = sum(short_lengths)

    n50 = calc_n50(lengths_sorted)

    # 给出适用性结论
    level, message = classify_genome(
        total_len=total_len,
        n_long=n_long,
        len_long=len_long,
        n_short=n_short,
        len_short=len_short,
    )

    # 生成报告文件名（与基因组在同一目录下）
    base_name = os.path.basename(genome_path)
    prefix = os.path.splitext(os.path.splitext(base_name)[0])[0]  # 兼容 .fa.gz
    report_name = prefix + REPORT_SUFFIX
    report_path = os.path.join(os.path.dirname(genome_path), report_name)

    # 组装报告内容
    lines = []
    lines.append("=== 宏观染色体共线性适用性评估报告 ===")
    lines.append(f"基因组文件：{base_name}")
    lines.append("")
    lines.append("【一、整体统计】")
    lines.append(f"  序列总条数：{n_total}")
    lines.append(f"  基因组总长度：{total_len} bp（约 {format_bp(total_len)}）")
    lines.append(f"  N50：{n50} bp（约 {format_bp(n50)}）")
    lines.append("")
    lines.append("【二、按长度阈值划分】")
    lines.append(f"  长序列判定阈值：>= {LONG_LEN_THRESHOLD} bp（约 {format_bp(LONG_LEN_THRESHOLD)}）")
    lines.append("")
    lines.append("  1. 明显特别长的序列（近似染色体 / 大 scaffold）：")
    lines.append(f"     条数：{n_long}")
    lines.append(f"     总长度：{len_long} bp（约 {format_bp(len_long)}）")
    lines.append(f"     占全基因组比例：{len_long/total_len*100:.2f}%")
    lines.append("")
    lines.append("  2. 零碎序列（小 scaffold / contig）：")
    lines.append(f"     条数：{n_short}")
    lines.append(f"     总长度：{len_short} bp（约 {format_bp(len_short)}）")
    lines.append(f"     占全基因组比例：{len_short/total_len*100:.2f}%")
    lines.append("")
    lines.append("【三、适用性结论】")
    lines.append(f"  等级：{level}")
    lines.append(f"  说明：{message}")
    lines.append("")
    lines.append("【四、使用建议】")
    lines.append(
        "  - 若等级为 A_推荐：可以将这些长序列视为 pseudochromosomes，"
        " 直接用于绘制类似“跨物种染色体带状共线性图”的宏观图。"
    )
    lines.append(
        "  - 若等级为 B_勉强可用：建议仅选取最长的若干条 scaffold 参与宏观图，"
        " 其余零碎部分用于局部共线性或基因家族分析。"
    )
    lines.append(
        "  - 若等级为 C_不推荐：不适合作为染色体级共线性图材料，"
        " 更适合作为基因尺度或局部 synteny 的参考。"
    )

    report_text = "\n".join(lines)

    # 写入报告文件
    with open(report_path, "w", encoding="utf-8") as out_f:
        out_f.write(report_text)

    # 同时在屏幕上打印简要信息，方便皇上即时查看
    print(report_text)
    print(f"\n报告已写入：{report_path}")


if __name__ == "__main__":
    main()

