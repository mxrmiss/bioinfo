#!/usr/bin/env bash
# =============================================================================
# 模块 C 扩展工具：构建 Salmon 可识别的 geneMap 文件（脚本内设参版）
# 生成两列（transcript_id   gene_id）
#
# 功能：
#   自动判断 GFF3 / GTF 格式，从注释文件中提取转录本与基因对应关系。
#   结果文件将保存为 ref/tx2gene.geneMap.tsv
#
# 作者：小茹（为皇上定制）
# =============================================================================
set -euo pipefail

# ----------------------------- 皇上在此设置参数 -----------------------------
ANNOT="ref/Sinonovacula_constricta_genome.gff3"   # 注释文件路径（.gff3 或 .gtf）
OUT="ref/tx2gene.geneMap.tsv"                     # 输出文件路径
# -----------------------------------------------------------------------------

mkdir -p "$(dirname "$OUT")"

echo "🔍 检测注释文件格式..."
# 判断是 GFF3 还是 GTF
FORMAT="gff3"
if grep -q "transcript_id" "$ANNOT"; then
  FORMAT="gtf"
fi
echo "🧠 格式识别结果：$FORMAT"

TMP=$(mktemp)

# ------------------------- 解析逻辑 -------------------------
if [[ "$FORMAT" == "gff3" ]]; then
  echo "🧩 从 GFF3 提取 ID / Parent..."
  awk -F'\t' '
    $3=="mRNA" || $3=="transcript" {
      match($9,/ID=([^;]+)/,a);
      match($9,/Parent=([^;]+)/,b);
      if(a[1]!="" && b[1]!="") print a[1]"\t"b[1];
    }' "$ANNOT" > "$TMP"

elif [[ "$FORMAT" == "gtf" ]]; then
  echo "🧩 从 GTF 提取 transcript_id / gene_id..."
  awk -F'\t' '
    $3=="transcript" || $3=="mRNA" {
      match($9,/transcript_id "([^"]+)"/,a);
      match($9,/gene_id "([^"]+)"/,b);
      if(a[1]!="" && b[1]!="") print a[1]"\t"b[1];
    }' "$ANNOT" > "$TMP"
fi

# 去重 + 检查
sort -u "$TMP" > "$OUT"
rm "$TMP"

LINES=$(wc -l < "$OUT" | tr -d ' ')
echo "✅ 已生成：$OUT"
echo "📊 共 ${LINES} 条转录本-基因映射记录"
echo "📁 示例："
head -5 "$OUT"

# 完整性检查
LINES=$(wc -l < "$OUT" | tr -d '[:space:]')
if (( LINES < 100 )); then
  echo "⚠️ 警告：生成的映射条目少于100，可能解析不完全。"
else
  echo "🎉 geneMap 构建完成，可直接用于 Salmon --geneMap 参数（共 ${LINES} 条）。"
fi

