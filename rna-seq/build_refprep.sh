#!/usr/bin/env bash
# =============================================================================
# build_refprep.sh —— 参考准备（可选仅 EVM）+ Salmon decoy-aware 索引
# 作者：小茹（为皇上定制）
#
# 功能概述：
#   1) （可选）过滤注释：当 EVM_ONLY=true 时，仅保留 EVM 体系条目
#   2) 用 gffread 从 GTF/GFF3 + 基因组提取转录本序列 transcripts.fa
#   3) 生成 decoys.txt（来自基因组序列头）
#   4) 生成 gentrome.fa = transcripts.fa + genome.fa
#   5) 构建 decoy-aware Salmon 索引：results/refprep/salmon_index
#
# 依赖：
#   gffread（Bioconda: mamba install -c bioconda gffread）
#   salmon  （Bioconda: mamba install -c bioconda salmon）
# =============================================================================
set -euo pipefail

# ----------------------------- 脚本内参数（皇上在此处修改） -----------------------------
GTF="ref/Sinonovacula_constricta_genome.gff3"   # 可为 .gtf 或 .gff3
FA="ref/Sinonovacula_constricta_genome.fa"      # 基因组
OUT="results/refprep"                           # 输出目录
THREADS=10
EVM_ONLY=false                                  # 是否启用“仅 EVM”严格模式（true/false）

# ----------------------------- 依赖检查 -----------------------------
command -v gffread >/dev/null || { echo "❌ 缺少 gffread"; exit 1; }
command -v salmon  >/dev/null || { echo "❌ 缺少 salmon";  exit 1; }

[[ -f "$GTF" ]] || { echo "❌ 找不到 GTF/GFF3: $GTF"; exit 1; }
[[ -f "$FA"  ]] || { echo "❌ 找不到基因组: $FA";   exit 1; }

# ----------------------------- 目录与路径 -----------------------------
mkdir -p "$OUT"
TMP="$OUT/_tmp_v3"; mkdir -p "$TMP"

EVM_GTF="$OUT/annotation.evm_only.gtf"   # 当 EVM_ONLY=true 时生成
USE_GTF="$GTF"                           # 实际用于 gffread 的注释路径（默认原始 GTF/GFF3）
TRANSCRIPTS="$OUT/transcripts.fa"
DECOYS="$OUT/decoys.txt"
GENTROME="$OUT/gentrome.fa"
INDEX_DIR="$OUT/salmon_index"

# ----------------------------- 可选：EVM 严格过滤 -----------------------------
if $EVM_ONLY; then
  echo "🧽 启用严格模式：仅保留 EVM 体系（EVM_ONLY=true）..."
  # 说明：
  #  - 同时兼容 GTF 与 GFF3 的属性字段；
  #  - 原逻辑只匹配 evm* 容易漏掉皇上的 Sco… ID。
  #  - 【最小改动】在不删原判定的前提下，新增两条“并集条件”：
  #      A) 第二列 source 为 EVM 也保留；
  #      B) gene_id/transcript_id/ID/Parent 以 Sco 开头也保留。
  awk -F'\t' '
    BEGIN{IGNORECASE=1}
    function get_attr_val(attr, key,   r, s) {
      # GTF: key "value";
      r=key "[[:space:]]+\"[^\"]+\""
      if (match(attr, r)) { s=substr(attr, RSTART, RLENGTH); gsub(key "[[:space:]]+\"|\"", "", s); return s }
      # GFF3: key=value;
      r=key "=[^;]+"
      if (match(attr, r)) { s=substr(attr, RSTART, RLENGTH); sub("^" key "=", "", s); return s }
      return ""
    }
    {
      if ($0 ~ /^#/) { print; next }
      attr=$9
      gid = get_attr_val(attr, "gene_id")
      tid = get_attr_val(attr, "transcript_id")
      id  = get_attr_val(attr, "ID")
      par = get_attr_val(attr, "Parent")

      gl = tolower(gid); tl = tolower(tid); il = tolower(id); pl = tolower(par)

      is_evm = 0
      # 原有：匹配 evm* 或 URL 编码 evm%2
      if (gl ~ /^evm[._]/ || index(gl,"evm%2")>0) is_evm=1
      if (tl ~ /^evm[._]/ || index(tl,"evm%2")>0) is_evm=1
      if (il ~ /^evm[._]/ || index(il,"evm%2")>0) is_evm=1
      if (pl ~ /^evm[._]/ || index(pl,"evm%2")>0) is_evm=1

      # 【新增并集条件 · 最小改动】
      # A) 第二列为 EVM 的，直接保留（与 ID 文本无关）
      if ($2 == "EVM") is_evm=1
      # B) 皇上物种的 ID/Parent 前缀是 Sco…，也要保留
      if (gid ~ /^Sco/ || tid ~ /^Sco/ || id ~ /^Sco/ || par ~ /^Sco/) is_evm=1

      if (is_evm) print $0
    }
  ' "$GTF" | sort -u > "$EVM_GTF"

  # 质量闸门（沿用原逻辑，不改动）
  if grep -q $'\tStringTie\t' "$EVM_GTF" || grep -q 'transcript_id[[:space:]]*"novel\.' "$EVM_GTF"; then
    echo "❌ 过滤后仍检测到 StringTie/novel 条目，请检查输入或规则。示例："
    (grep -n $'\tStringTie\t' "$EVM_GTF" || true) | head -3
    (grep -n 'transcript_id[[:space:]]*"novel\.' "$EVM_GTF" || true) | head -3
    exit 2
  fi

  echo "✅ EVM 过滤完成：$(wc -l < "$EVM_GTF") 行 → $EVM_GTF"
  USE_GTF="$EVM_GTF"
else
  echo "✅ 默认模式：不过滤 GTF/GFF3，直接使用原始注释 → $GTF"
fi

# ----------------------------- 提取转录本序列 -----------------------------
echo "🧵 gffread 提取转录本序列 ..."
gffread "$USE_GTF" -g "$FA" -w "$TRANSCRIPTS"

# ----------------------------- 自适应体检（替代固定 5 万阈值） -----------------------------
TXN_NUM=$(grep -c '^>' "$TRANSCRIPTS" || true)
echo "📊 transcripts.fa 条目数：$TXN_NUM"

# 从当前注释中估计期望的 mRNA/transcript 数（自适应基准）
EXPECTED_MRNA=$(awk '$3=="mRNA"||$3=="transcript"{n++}END{print n+0}' "$USE_GTF")
EXPECTED_MRNA=${EXPECTED_MRNA:-0}
echo "📐 注释中的 mRNA/transcript 数：$EXPECTED_MRNA"

# 先保留硬闸：明显异常直接拦截（防止假索引）
if [[ "${TXN_NUM:-0}" -lt 1000 ]]; then
  echo "❌ 异常：转录本数量 < 1000，疑似筛选规则或输入路径错误。已停止以防生成无效索引。"
  exit 3
fi

# 自适应提示：按比例而不是固定 50,000
if [[ "$EXPECTED_MRNA" -ge 5000 ]]; then
  RATIO=$(awk -v a="$TXN_NUM" -v b="$EXPECTED_MRNA" 'BEGIN{if(b<1)print 0; else printf "%.3f", a/b}')
  if awk 'BEGIN{exit !('"$RATIO"' < 0.60)}'; then
    echo "⚠️ 注意：transcripts.fa / 注释mRNA 比例 = $RATIO (<0.60)。可能过滤过严或提取规则不匹配，请确认。"
  elif awk 'BEGIN{exit !('"$RATIO"' < 0.90)}'; then
    echo "ℹ️ 提示：比例 = $RATIO（略低于期望，若非刻意过滤可再核对）。"
  else
    echo "✅ 自检通过：transcripts.fa 与注释规模一致（比例 $RATIO）。"
  fi
else
  echo "ℹ️ 提示：注释规模较小（mRNA=$EXPECTED_MRNA），跳过比例告警，仅保留硬闸检查。"
fi

# ----------------------------- 生成 decoys/gentrome -----------------------------
echo "📜 生成 decoys.txt ..."
grep "^>" "$FA" | sed "s/^>//; s/ .*//" > "$DECOYS"

echo "🧷 生成 gentrome.fa ..."
cat "$TRANSCRIPTS" "$FA" > "$GENTROME"

# ----------------------------- Salmon 索引 -----------------------------
echo "⚒️ 构建 Salmon 索引（threads=$THREADS） ..."
rm -rf "$INDEX_DIR"
mkdir -p "$INDEX_DIR"
# 说明：-t 使用 gentrome，-d 提供 decoys 列表；为 decoy-aware 工作流的标准做法
salmon index -t "$GENTROME" -d "$DECOYS" -p "$THREADS" -i "$INDEX_DIR"

# ----------------------------- 收尾与自检 -----------------------------
rm -rf "$TMP"
echo "🎉 完成！产物位于：$OUT"
du -h --max-depth=1 "$OUT" || true

echo "🔎 索引自检（应看到非空 decoys 列表）："
grep -n '"decoys"' "$INDEX_DIR/info.json" || true

echo "✅ 使用方法提示：定量时指向该索引（示例）"
echo "  salmon quant -i $INDEX_DIR -l A -1 sample_R1.fq.gz -2 sample_R2.fq.gz -p $THREADS --validateMappings --gcBias --seqBias -o results/quant/xxx"

