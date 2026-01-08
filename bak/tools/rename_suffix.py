#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ============================================================
# 批量修改文件后缀名脚本（皇上御用版）
# ------------------------------------------------------------
# 功能：
#   将脚本所在文件夹中的所有文件的后缀名批量改为指定的新后缀。
# 使用方法：
#   直接执行本脚本即可（无需命令行参数）。
# ============================================================

import os

# ==========【参数区：皇上可自行修改】==========
# 想要修改成的新后缀名（不带点号）
NEW_EXT = "faa"

# 是否递归修改子文件夹内的文件
RECURSIVE = False

# 是否保留原扩展名（如 True，则会在原扩展名后追加新的扩展名）
KEEP_OLD_EXT = False
# 例如：example.csv → example.csv.txt
# ===========================================


def rename_files_in_dir(folder):
    """在指定文件夹中批量修改文件后缀"""
    this_script = os.path.basename(__file__)

    for root, dirs, files in os.walk(folder):
        for filename in files:
            # 跳过脚本自身
            if filename == this_script:
                continue

            old_path = os.path.join(root, filename)
            name, ext = os.path.splitext(filename)

            # 生成新文件名
            if KEEP_OLD_EXT and ext:
                new_name = f"{filename}.{NEW_EXT}"
            else:
                new_name = f"{name}.{NEW_EXT}"

            new_path = os.path.join(root, new_name)

            # 如果目标文件已存在则跳过
            if os.path.exists(new_path):
                print(f"⚠️ 跳过：{new_path} 已存在")
                continue

            os.rename(old_path, new_path)
            print(f"✅ 已改名：{filename} → {new_name}")

        if not RECURSIVE:
            break


if __name__ == "__main__":
    folder_path = os.path.dirname(os.path.abspath(__file__))
    print(f"📂 当前工作目录：{folder_path}")
    print(f"🎯 目标后缀：.{NEW_EXT}\n")
    rename_files_in_dir(folder_path)
    print("\n✨ 全部处理完毕！皇上可安坐无忧。")
