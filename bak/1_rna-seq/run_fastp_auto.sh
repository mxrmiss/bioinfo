#!/usr/bin/env bash
# =====================================================================
# 公司口径（仅保留 clean FASTQ + 合并报告）
# - Clean FASTQ 放 results/qc/clean
# - 合并报告 HTML/JSON 放 results/qc/report
# - 中间文件与单独 raw/clean 报告自动删除
# =====================================================================

set -euo pipefail

THREADS=${THREADS:-8}
INDIR=${1:-data}
CLEANDIR="results/qc/clean"
REPORTDIR="results/qc/report"
TMPDIR="results/qc/_tmp"
mkdir -p "$CLEANDIR" "$REPORTDIR" "$TMPDIR"

# Illumina TruSeq（如公司提供 adapters.fa，请替换）
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

shopt -s nullglob
declare -A seen

for R1 in "$INDIR"/*_1.fq.gz "$INDIR"/*_R1.fq.gz "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq.gz \
          "$INDIR"/*_1_001.fastq.gz "$INDIR"/*_R1_001.fastq.gz "$INDIR"/*_1_001.fq.gz "$INDIR"/*_R1_001.fq.gz; do
  [ -e "$R1" ] || continue
  base=$(basename "$R1")

  sample=${base}
  sample=${sample%_R1_001.fastq.gz}; sample=${sample%_R1.fastq.gz}
  sample=${sample%_1_001.fastq.gz};   sample=${sample%_1.fastq.gz}
  sample=${sample%_R1_001.fq.gz};     sample=${sample%_R1.fq.gz}
  sample=${sample%_1_001.fq.gz};      sample=${sample%_1.fq.gz}

  [[ -n "${seen[$sample]:-}" ]] && continue
  seen[$sample]=1

  case "$base" in
    *_R1_001.fastq.gz) R2="$INDIR/${sample}_R2_001.fastq.gz" ;;
    *_R1.fastq.gz)     R2="$INDIR/${sample}_R2.fastq.gz" ;;
    *_1_001.fastq.gz)  R2="$INDIR/${sample}_2_001.fastq.gz" ;;
    *_1.fastq.gz)      R2="$INDIR/${sample}_2.fastq.gz" ;;
    *_R1_001.fq.gz)    R2="$INDIR/${sample}_R2_001.fq.gz" ;;
    *_R1.fq.gz)        R2="$INDIR/${sample}_R2.fq.gz" ;;
    *_1_001.fq.gz)     R2="$INDIR/${sample}_2_001.fq.gz" ;;
    *_1.fq.gz)         R2="$INDIR/${sample}_2.fq.gz" ;;
    *) echo "无法识别 R1/R2 模式：$base" >&2; continue ;;
  esac
  if [[ ! -f "$R2" ]]; then
    echo "⚠️ 未找到配对文件：$R2，跳过 $sample"
    continue
  fi

  # ① fastp 原始报告（写 TMPDIR）
  echo "📄 fastp 原始数据报告：$sample"
  fastp -i "$R1" -I "$R2" \
    -o "$TMPDIR/${sample}.raw.R1.tmp.fq.gz" \
    -O "$TMPDIR/${sample}.raw.R2.tmp.fq.gz" \
    --thread "$THREADS" \
    --disable_adapter_trimming \
    --disable_length_filtering \
    --disable_quality_filtering \
    --html "$TMPDIR/${sample}_raw.fastp.html" \
    --json "$TMPDIR/${sample}_raw.fastp.json" \
    >/dev/null

  # ② cutadapt 丢接头
  echo "🔎 cutadapt 丢接头对：$sample"
  cutadapt -a "$ADAPTER_R1" -A "$ADAPTER_R2" -O 5 -e 0.1 \
    --pair-filter=any --discard-trimmed -j "$THREADS" \
    -o "$TMPDIR/${sample}.noAdapter.R1.fq.gz" \
    -p "$TMPDIR/${sample}.noAdapter.R2.fq.gz" \
    "$R1" "$R2" > "$TMPDIR/${sample}.cutadapt.log"

  # ③ fastp 清洗（仅 N/Q 过滤；不剪长度）—— clean FASTQ + clean 报告（TMPDIR）
  echo "🚀 fastp 清洗：$sample"
  fastp -i "$TMPDIR/${sample}.noAdapter.R1.fq.gz" \
    -I "$TMPDIR/${sample}.noAdapter.R2.fq.gz" \
    -o "$CLEANDIR/${sample}_R1.clean.fq.gz" \
    -O "$CLEANDIR/${sample}_R2.clean.fq.gz" \
    --disable_adapter_trimming \
    --disable_length_filtering \
    --n_base_limit 0 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 50 \
    --thread "$THREADS" \
    --html "$TMPDIR/${sample}_clean.fastp.html" \
    --json "$TMPDIR/${sample}_clean.fastp.json" \
    >/dev/null

  # ④ 合并报告：输出到 REPORTDIR
  echo "🧩 合并报告：$sample"
  python3 - "$TMPDIR/${sample}_raw.fastp.json" "$TMPDIR/${sample}_clean.fastp.json" "$REPORTDIR/${sample}.fastp.merged.json" "$REPORTDIR/${sample}.fastp.merged.html" "$sample" <<'PY'
import json, sys, html
raw_json, clean_json, out_json, out_html, samp = sys.argv[1:]
r = json.load(open(raw_json))
c = json.load(open(clean_json))
rb = r["summary"]["before_filtering"]
ca = c["summary"]["after_filtering"]
merged = {"sample": samp, "raw": rb, "clean": ca}
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(merged, f, ensure_ascii=False, indent=2)

def row(s):
    return f"""
    <tr><td>total reads</td><td>{s['total_reads']:,}</td></tr>
    <tr><td>total bases</td><td>{s['total_bases']/1e9:.6f} G</td></tr>
    <tr><td>Q20 bases</td><td>{s['q20_bases']/1e9:.6f} G ({s['q20_rate']*100:.6f}%)</td></tr>
    <tr><td>Q30 bases</td><td>{s['q30_bases']/1e9:.6f} G ({s['q30_rate']*100:.6f}%)</td></tr>
    <tr><td>GC content</td><td>{s['gc_content']*100:.6f}%</td></tr>
    """

html_text = f"""<!doctype html>
<meta charset="utf-8">
<title>{html.escape(samp)} • merged fastp report</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;line-height:1.5;margin:20px}}
h2{{margin-top:22px}} table{{border-collapse:collapse;min-width:380px}}
td,th{{border:1px solid #ccc;padding:6px 10px}} th{{background:#f5f5f5;text-align:left}}
.wrap{{display:flex;gap:24px;flex-wrap:wrap}} .note{{color:#666;margin:8px 0 18px}}
</style>
<h1>Sample: {html.escape(samp)}</h1>
<p class="note">合并视图：左侧 <b>原始 (raw)</b>；右侧 <b>clean</b>。仅保留本页与同名 JSON。</p>
<div class="wrap">
  <div><h2>Raw</h2><table><tbody>{row(rb)}</tbody></table></div>
  <div><h2>Clean</h2><table><tbody>{row(ca)}</tbody></table></div>
</div>
"""
open(out_html, "w", encoding="utf-8").write(html_text)
PY

  # ⑤ 清理该样本的临时文件
  rm -f "$TMPDIR/${sample}"*
  echo "✅ 完成并清理：$sample"
  echo "--------------------------------------"
done

echo "🎉 全部完成！"
echo "   • Clean FASTQ 在: $CLEANDIR"
echo "   • 报告在: $REPORTDIR"
