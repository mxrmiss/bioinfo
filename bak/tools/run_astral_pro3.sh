#!/usr/bin/env bash
set -euo pipefail

############################################
# ========== 用户自定义参数区域 ==========
############################################

# astral-pro3 可执行文件路径
# 如果你是通过 aster 安装的，且 astral-pro3 在 PATH 中，保持默认即可
ASTRAL_PRO3_BIN="astral-pro3"

# 线程数
THREADS=20

# 是否指定外群（留空则不定根）
OUTGROUP="Lottia_gigantea"

# 输入文件（已存在）
GENE_TREES="Resolved_Gene_Trees.txt"

# 输出前缀
OUT_PREFIX="astralpro3"

############################################
# ========== 不要修改下面的内容 ==========
############################################

echo "[INFO] extracting gene leaves with '|' ..."
tr '(),;' '\n' < "${GENE_TREES}" | sed 's/:.*//' | grep '|' | sort -u > genes.withpipe.unique.txt

echo "[INFO] building gene-to-species mapping ..."
perl -ne '
chomp;
$g=$_;
($sp)=split(/\|/,$g,2);
$sp =~ s/^(.+)_\1$/$1/;
print "$g\t$sp\n";
' genes.withpipe.unique.txt > gene2species.map

echo "[INFO] running ASTRAL-Pro3 ..."
if [[ -n "${OUTGROUP}" ]]; then
    ${ASTRAL_PRO3_BIN} \
        -t ${THREADS} \
        --root ${OUTGROUP} \
        -a gene2species.map \
        -i ${GENE_TREES} \
        -o ${OUT_PREFIX}.stree.tre \
        2> >(tee ${OUT_PREFIX}.log >&2)
else
    ${ASTRAL_PRO3_BIN} \
        -t ${THREADS} \
        -a gene2species.map \
        -i ${GENE_TREES} \
        -o ${OUT_PREFIX}.stree.tre \
        2> >(tee ${OUT_PREFIX}.log >&2)
fi

echo "[INFO] done."
