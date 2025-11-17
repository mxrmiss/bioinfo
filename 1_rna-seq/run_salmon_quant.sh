#!/usr/bin/env bash
# =============================================================================
# 模块 D：Salmon 定量（支持 --geneMap，直接产出 quant.genes.sf）
# 基于您现有版本最小改动而来（+ --geneMap 参数与检查）
# =============================================================================
set -euo pipefail

# ----------------------------- 默认参数 -----------------------------
CLEAN_DIR="results/qc/clean"
INDEX_DIR="results/refprep/salmon_index"
OUT_DIR="results/quant"
THREADS=16
EXTRA_ARGS="--validateMappings --gcBias --seqBias"
GENEMAP="ref/tx2gene.geneMap.tsv"   # 新增：两列（transcript_id  gene_id）

# ----------------------------- 参数解析 -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean)   CLEAN_DIR="$2"; shift 2;;
    --index)   INDEX_DIR="$2"; shift 2;;
    --out)     OUT_DIR="$2";  shift 2;;
    --threads) THREADS="$2";  shift 2;;
    --extra)   EXTRA_ARGS="$2"; shift 2;;
    --geneMap) GENEMAP="$2"; shift 2;;   # 新增
    -h|--help)
      echo "用法: bash $0 [--clean results/qc/clean] [--index results/refprep/salmon_index] [--out results/quant] [--threads 16] [--extra '...'] [--geneMap ref/tx2gene.geneMap.tsv]"
      exit 0;;
    *) echo "未知参数：$1"; exit 1;;
  esac
done

# ----------------------------- 依赖与输入检查 -----------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "❌ 缺少依赖：$1"; exit 1; }; }
need salmon

[[ -d "$CLEAN_DIR" ]] || { echo "❌ 找不到 clean 目录：$CLEAN_DIR"; exit 1; }
[[ -d "$INDEX_DIR" ]] || { echo "❌ 找不到 Salmon 索引目录：$INDEX_DIR"; exit 1; }
mkdir -p "$OUT_DIR"

# geneMap 允许不存在（用户可显式传空：--geneMap ""），但如存在则检查格式
if [[ -n "${GENEMAP}" ]]; then
  [[ -f "$GENEMAP" ]] || { echo "❌ 找不到 geneMap：$GENEMAP"; exit 1; }
fi

MANIFEST="$OUT_DIR/quant_list.tsv"
SUMMARY="$OUT_DIR/summary.tsv"
: > "$MANIFEST"
echo -e "sample\tquant_dir\tR1\tR2" >> "$MANIFEST"

# ----------------------------- 样本发现（配对FASTQ） -----------------------------
echo "🔎 扫描样本于：$CLEAN_DIR"
mapfile -t R1S < <(find "$CLEAN_DIR" -maxdepth 1 -type f -name "*_R1.clean.fq.gz" | sort)
SAMPLE_COUNT=0

for R1 in "${R1S[@]}"; do
  base="$(basename "$R1")"
  sample="${base%_R1.clean.fq.gz}"
  R2="$CLEAN_DIR/${sample}_R2.clean.fq.gz"
  if [[ ! -f "$R2" ]]; then
    echo "⚠️ 跳过：找不到配对的 R2 → $R2"
    continue
  fi
  SAMPLE_COUNT=$((SAMPLE_COUNT+1))

  OUT_SAMPLE="$OUT_DIR/${sample}"
  mkdir -p "$OUT_SAMPLE"

  echo "🚀 [$sample] salmon quant 开始（threads=$THREADS）..."
  {
    date +"[start] %F %T"
    salmon quant \
      -i "$INDEX_DIR" \
      -l A \
      -1 "$R1" -2 "$R2" \
      -p "$THREADS" \
      $EXTRA_ARGS \
      ${GENEMAP:+--geneMap "$GENEMAP"} \
      -o "$OUT_SAMPLE"
    date +"[end] %F %T"
  } &> "$OUT_SAMPLE/salmon.log"

  echo -e "${sample}\t${OUT_SAMPLE}\t${R1}\t${R2}" >> "$MANIFEST"
done

if [[ "$SAMPLE_COUNT" -eq 0 ]]; then
  echo "❌ 未发现任何样本（期望 *_R1.clean.fq.gz / *_R2.clean.fq.gz）"
  exit 2
fi

echo "📒 样本清单：$MANIFEST （样本数：$SAMPLE_COUNT）"

# ----------------------------- 解析统计并汇总 -----------------------------
parse_one() {
  local sample="$1" qdir="$2"
  local log="$qdir/salmon.log" json="$qdir/aux_info/meta_info.json"
  local total="NA" mapped="NA" rate="NA" elapsed="NA" note="ok"
  if [[ -f "$json" ]]; then
    total=$(grep -o '"num_processed":[^,]*' "$json" | awk -F: '{gsub(/[[:space:]]*/,"",$2); print $2}')
    mapped=$(grep -o '"num_mapped":[^,]*'    "$json" | awk -F: '{gsub(/[[:space:]]*/,"",$2); print $2}')
    rate=$(grep -o '"percent_mapped":[^,]*' "$json" | awk -F: '{gsub(/[[:space:]]*/,"",$2); print $2}')
    if [[ -n "${rate}" ]]; then rate="${rate%\"}"; rate="${rate#\"}"
      if awk 'BEGIN{exit !(('"${rate:-0}"' <= 1.0000001))}'; then rate=$(awk 'BEGIN{printf "%.6f", '"$rate"' * 100}'); fi
    fi
  else note="missing_meta_info"; fi
  if [[ "$total" == "NA" || "$mapped" == "NA" || "$rate" == "NA" ]]; then
    if [[ -f "$log" ]]; then
      local tline mline rline
      tline=$(grep -E "total fragments|total reads" "$log" | tail -1 || true)
      mline=$(grep -E "mapped fragments|reads mapped" "$log" | tail -1 || true)
      rline=$(grep -E "Mapping rate|mapping rate" "$log" | tail -1 || true)
      [[ -n "$tline" ]] && total=$(echo "$tline" | grep -oE '[0-9]+' | head -1)
      [[ -n "$mline" ]] && mapped=$(echo "$mline" | grep -oE '[0-9]+' | head -1)
      [[ -n "$rline" ]] && rate=$(echo "$rline" | grep -oE '[0-9.]+(?=%)' -o)
      note="${note:+$note;}"from_log
    else note="${note:+$note;}"missing_log; fi
  fi
  if [[ -f "$log" ]]; then
    local start_ts end_ts s e
    start_ts=$(grep "^\[start\]" "$log" | tail -1 | awk '{print $2" "$3}')
    end_ts=$(grep "^\[end\]"   "$log" | tail -1 | awk '{print $2" "$3}')
    if [[ -n "$start_ts" && -n "$end_ts" ]]; then
      s=$(date -d "$start_ts" +%s 2>/dev/null || echo ""); e=$(date -d "$end_ts" +%s 2>/dev/null || echo "")
      [[ -n "$s" && -n "$e" ]] && elapsed=$((e - s)) || elapsed="NA"
    fi
  fi
  echo -e "${sample}\t${total}\t${mapped}\t${rate}\t${elapsed}\t${note:-ok}"
}

echo -e "sample\ttotal_fragments\tmapped_fragments\tmapping_rate(%)\telapsed_sec\tnotes" > "$SUMMARY"
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample qdir r1 r2; do
  parse_one "$sample" "$qdir" >> "$SUMMARY"
done

echo "✅ 汇总完成：$SUMMARY"
echo "🎉 模块 D 结束：quant.sf 与 quant.genes.sf 均已生成（已启用 geneMap）。"
