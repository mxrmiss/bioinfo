#!/usr/bin/env bash
set -euo pipefail

# 基因→GO 映射表路径（可以按需改）
GENE2GO="results/07_annot/gene2go.tsv"

# 简单检查
if [ ! -f "$GENE2GO" ]; then
  echo "ERROR: gene2go mapping not found at: $GENE2GO" >&2
  exit 1
fi

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 gene_list1 [gene_list2 ...]" >&2
  exit 1
fi

# 解析脚本所在目录，用来判断“是否和输入在同一目录”
case "$0" in
  /*)
    SCRIPT_FULL="$0"
    ;;
  *)
    SCRIPT_FULL="$(pwd)/$0"
    ;;
esac
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_FULL")" && pwd)"

# 先把 gene2go 去掉表头并按 gene_id 排序，只做一次
TMP_GENE2GO="$(mktemp)"
trap 'rm -f "$TMP_GENE2GO"' EXIT

LC_ALL=C tail -n +2 "$GENE2GO" | LC_ALL=C sort -k1,1 > "$TMP_GENE2GO"

# 逐个输入文件处理
for INPUT in "$@"; do
  if [ ! -f "$INPUT" ]; then
    echo "WARNING: input not found, skip: $INPUT" >&2
    continue
  fi

  IN_DIR="$(cd "$(dirname "$INPUT")" && pwd)"
  IN_BASE="$(basename "$INPUT")"

  # 取前缀（去掉最后一个扩展名）
  PREFIX="${IN_BASE%.*}"
  OUT_FILE="${IN_DIR}/${PREFIX}.tsv"

  # 如果脚本和输入在同一目录，且输入后缀为 .tsv，则先备份
  if [[ "$IN_DIR" == "$SCRIPT_DIR" && "$IN_BASE" == *.tsv ]]; then
    cp "$INPUT" "${INPUT}.bak"
    echo "Backup created: ${INPUT}.bak" >&2
  fi

  echo "Processing: $INPUT -> $OUT_FILE" >&2

  # 排序输入基因列表，然后与已排序 gene2go 做 join
  join -t $'\t' -1 1 -2 1 \
    <(LC_ALL=C sort -k1,1 "$INPUT") \
    "$TMP_GENE2GO" \
    | awk 'BEGIN{OFS="\t"; print "gene_id","GO_id"} {print $1,$2}' \
    > "$OUT_FILE"
done

