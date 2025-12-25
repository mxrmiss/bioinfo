#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_cafe_ids_tree.py

作用：
- 读取 CAFE5 的 Base_report.cafe（或 *_report.cafe）
- 提取 "# IDs of nodes:" 这一行的 Newick（带 <node_id>）
- 画出一棵简单的树：叶子显示 "species<id>"，内部节点显示 "<id>"
- 输出 PNG + PDF，方便你核对节点编号

使用方式：
- 把本脚本放到包含 Base_report.cafe 的目录里
- 修改下面的参数（如有需要）
- 运行后会在输出路径生成图片
"""

# =========================
# 皇上手动填写的参数区（不走命令行参数）
# =========================
INPUT_CAFE_REPORT = "Base_report.cafe"

OUTPUT_PNG = "cafe_node_ids_tree.png"
OUTPUT_PDF = "cafe_node_ids_tree.pdf"

FIG_WIDTH_INCH = 14
FIG_HEIGHT_INCH = 10
DPI = 300

FONT_SIZE = 10  # 树比较大时可以调小，比如 7 或 6
SHOW_BRANCH_LENGTHS = False  # 只为了看节点编号，默认不显示分支长度

# =========================
# 正式逻辑区
# =========================
import re
import sys
from io import StringIO

def die(msg: str, code: int = 1) -> None:
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(code)

def read_ids_newick_from_report(path: str) -> str:
    """
    从 Base_report.cafe 中提取 "# IDs of nodes:" 后面的 Newick 串
    """
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                # 兼容前导空格
                if line.lstrip().startswith("# IDs of nodes:"):
                    newick = line.split(":", 1)[1].strip()
                    if not newick.endswith(";"):
                        newick += ";"
                    return newick
    except FileNotFoundError:
        die(f"[ERROR] 找不到输入文件：{path}")
    die(f"[ERROR] 在 {path} 中没有找到 '# IDs of nodes:' 这一行。")

def sanitize_newick_for_parser(newick: str):
    """
    Bio.Phylo 对名字里包含 < > 有时会解析不稳定（不同版本差异）。
    所以这里先把:
      - 物种名<1>  -> 物种名__ID1
      - 内部节点 <31> -> __ID31
    解析后再把显示名还原成你想看的形式。
    """
    # 把叶子 "name<123>" 改成 "name__ID123"
    newick2 = re.sub(r"([A-Za-z0-9_.-]+)<(\d+)>", r"\1__ID\2", newick)
    # 把内部节点 ")<123>" 改成 ")__ID123"
    newick2 = re.sub(r"\)<(\d+)>", r")__ID\1", newick2)
    return newick2

def restore_display_name(name: str) -> str:
    """
    将解析后的 name（例如 Anadara__ID1 或 __ID31）
    还原成显示用的 "Anadara<1>" 或 "<31>"
    """
    if name is None:
        return None
    # 内部节点
    if name.startswith("__ID"):
        return "<" + name.replace("__ID", "") + ">"
    # 叶节点
    if "__ID" in name:
        sp, nid = name.split("__ID", 1)
        return f"{sp}<{nid}>"
    return name

def main():
    # 依赖检查
    try:
        from Bio import Phylo
    except Exception:
        die("[ERROR] 缺少依赖：biopython（模块 Bio）。请先在你的环境里安装 biopython。")

    try:
        import matplotlib
        matplotlib.use("Agg")  # 强制无界面后端，服务器/ssh 也能画
        import matplotlib.pyplot as plt
    except Exception:
        die("[ERROR] 缺少依赖：matplotlib。请先在你的环境里安装 matplotlib。")

    raw_newick = read_ids_newick_from_report(INPUT_CAFE_REPORT)
    safe_newick = sanitize_newick_for_parser(raw_newick)

    # 读取树
    handle = StringIO(safe_newick)
    tree = Phylo.read(handle, "newick")

    # 还原显示名：叶子和内部节点都显示 <id>
    for clade in tree.find_clades(order="preorder"):
        if clade.name:
            clade.name = restore_display_name(clade.name)
        # 你这条 IDs tree 没有分支长度（CAFE 这行通常不带长度）
        # 这里不强行补，保证“只看节点”

    # 作图
    fig = plt.figure(figsize=(FIG_WIDTH_INCH, FIG_HEIGHT_INCH))
    ax = fig.add_subplot(1, 1, 1)

    # Bio.Phylo.draw 会自动画树；branch_labels 控制是否显示长度
    branch_labels = (lambda c: c.branch_length) if SHOW_BRANCH_LENGTHS else None
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        branch_labels=branch_labels
    )

    # 调整字体大小（Bio.Phylo.draw 画出来的文本是 matplotlib 的 Text 对象）
    for text in ax.texts:
        text.set_fontsize(FONT_SIZE)

    ax.set_title("CAFE5 node-ID tree (from # IDs of nodes:)", fontsize=FONT_SIZE + 2)
    fig.tight_layout()

    fig.savefig(OUTPUT_PNG, dpi=DPI)
    fig.savefig(OUTPUT_PDF)

    sys.stderr.write(f"[OK] 已输出：{OUTPUT_PNG}\n")
    sys.stderr.write(f"[OK] 已输出：{OUTPUT_PDF}\n")

if __name__ == "__main__":
    main()

