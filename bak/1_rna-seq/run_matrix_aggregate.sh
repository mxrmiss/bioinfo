#!/usr/bin/env bash
set -euo pipefail

# ========= 皇上在此设定参数 =========
GTF="ref/Sinonovacula_constricta_genome.gff3"   # 可为 .gff3 或 .gtf
QUANT_DIR="results/quant"                       # 定量结果目录
OUT_DIR="results/matrix"                        # 输出目录
CFA="no"                                        # no | scaledTPM | lengthScaledTPM
# ===================================

[[ -f "$GTF" ]] || { echo "❌ 找不到注释：$GTF"; exit 2; }
[[ -d "$QUANT_DIR" ]] || { echo "❌ 找不到定量目录：$QUANT_DIR"; exit 2; }
mkdir -p "$OUT_DIR/logs"

# ========= 收集 quant.sf 与“转录本白名单” =========
QSF_LIST="$OUT_DIR/logs/quant_sf.list"; : > "$QSF_LIST"
if [[ -f "$QUANT_DIR/quant_list.tsv" ]]; then
  tail -n +2 "$QUANT_DIR/quant_list.tsv" | awk -F'\t' '{print $2"/quant.sf"}' \
    | while read -r q; do [[ -f "$q" ]] && echo "$q" >> "$QSF_LIST"; done
else
  find "$QUANT_DIR" -mindepth 2 -maxdepth 2 -type f -name quant.sf | sort > "$QSF_LIST"
fi
NQ=$(wc -l < "$QSF_LIST" | tr -d ' ')
[[ $NQ -gt 0 ]] || { echo "❌ 未发现 quant.sf"; exit 2; }
echo "🔎 发现 quant.sf 数量：$NQ"

TX_USED="$OUT_DIR/logs/transcripts.used.list"; : > "$TX_USED"
while read -r qsf; do awk 'NR>1{print $1}' "$qsf"; done < "$QSF_LIST" | sort -u > "$TX_USED"
NTX=$(wc -l < "$TX_USED" | tr -d ' ')
echo "🧾 quant 中实际转录本种类：$NTX"

# ========= 判定注释格式 =========
IS_GFF3=$(awk -F'\t' 'NR<=2000 && $0!~/^#/ && NF>=9 {if(index($9,"=")>0){print "YES"; exit}}' "$GTF")
echo "🧠 解析模式："$([[ "${IS_GFF3:-}" == "YES" ]] && echo "GFF3" || echo "GTF")

# ========= 构建 name→gene 与 id→gene、以及 gene 显示名 =========
NAME2GENE="$OUT_DIR/logs/_name2gene.tsv"
ID2GENE="$OUT_DIR/logs/_id2gene.tsv"
GENE_LABEL="$OUT_DIR/logs/_gene_label.tsv"
: > "$NAME2GENE"; : > "$ID2GENE"; : > "$GENE_LABEL"

if [[ "${IS_GFF3:-}" == "YES" ]]; then
  # gene 显示名
  awk -F'\t' '
    $0!~/^#/ && tolower($3)=="gene" {
      a=$9; gid=""; gname="";
      if (match(a, /(^|;)ID=[^;]+/))   { gid=substr(a,RSTART,RLENGTH); sub(/(^|;)ID=/,"",gid);   sub(/;.*/,"",gid) }
      if (match(a, /(^|;)Name=[^;]+/)) { gname=substr(a,RSTART,RLENGTH); sub(/(^|;)Name=/,"",gname); sub(/;.*/,"",gname) }
      print gid "\t" (gname!=""?gname:gid)
  }' "$GTF" | sort -u > "$GENE_LABEL"

  # mRNA：Name→Parent
  awk -F'\t' '
    $0!~/^#/ && tolower($3) ~ /^(transcript|mrna)$/ {
      a=$9; name=""; par="";
      if (match(a, /(^|;)Name=[^;]+/))   { name=substr(a,RSTART,RLENGTH);   sub(/(^|;)Name=/,"",name); sub(/;.*/,"",name) }
      if (match(a, /(^|;)Parent=[^;]+/)) { par=substr(a,RSTART,RLENGTH);    sub(/(^|;)Parent=/,"",par); sub(/;.*/,"",par) }
      if (name!="" && par!="") print name "\t" par
  }' "$GTF" | sort -u > "$NAME2GENE"

  # mRNA：ID→Parent
  awk -F'\t' '
    $0!~/^#/ && tolower($3) ~ /^(transcript|mrna)$/ {
      a=$9; id=""; par="";
      if (match(a, /(^|;)ID=[^;]+/))     { id=substr(a,RSTART,RLENGTH);     sub(/(^|;)ID=/,"",id);     sub(/;.*/,"",id) }
      if (match(a, /(^|;)Parent=[^;]+/)) { par=substr(a,RSTART,RLENGTH);    sub(/(^|;)Parent=/,"",par); sub(/;.*/,"",par) }
      if (id!="" && par!="") print id "\t" par
  }' "$GTF" | sort -u > "$ID2GENE"
else
  # GTF gene 显示名
  awk -F'\t' '
    $0!~/^#/ && $3=="gene" {
      a=$9; gid=""; gname="";
      if (match(a,/gene_id[[:space:]]+"[^"]+"/))   { gid=substr(a,RSTART,RLENGTH);   gsub(/gene_id[[:space:]]+"|"/,"",gid) }
      if (match(a,/gene_name[[:space:]]+"[^"]+"/)) { gname=substr(a,RSTART,RLENGTH); gsub(/gene_name[[:space:]]+"|"/,"",gname) }
      print gid "\t" (gname!=""?gname:gid)
  }' "$GTF" | sort -u > "$GENE_LABEL"

  # transcript_name → gene_id
  awk -F'\t' '
    $0!~/^#/ && ($3=="transcript" || $3=="mRNA") {
      a=$9; tname=""; gid="";
      if (match(a,/transcript_name[[:space:]]+"[^"]+"/)) { tname=substr(a,RSTART,RLENGTH); gsub(/transcript_name[[:space:]]+"|"/,"",tname) }
      if (match(a,/gene_id[[:space:]]+"[^"]+"/))         { gid=substr(a,RSTART,RLENGTH);   gsub(/gene_id[[:space:]]+"|"/,"",gid) }
      if (tname!="" && gid!="") print tname "\t" gid
  }' "$GTF" | sort -u > "$NAME2GENE"

  # transcript_id → gene_id
  awk -F'\t' '
    $0!~/^#/ && ($3=="transcript" || $3=="mRNA") {
      a=$9; tid=""; gid="";
      if (match(a,/transcript_id[[:space:]]+"[^"]+"/))   { tid=substr(a,RSTART,RLENGTH);   gsub(/transcript_id[[:space:]]+"|"/,"",tid) }
      if (match(a,/gene_id[[:space:]]+"[^"]+"/))         { gid=substr(a,RSTART,RLENGTH);   gsub(/gene_id[[:space:]]+"|"/,"",gid) }
      if (tid!="" && gid!="") print tid "\t" gid
  }' "$GTF" | sort -u > "$ID2GENE"
fi

# ========= 生成 tx2gene.final.tsv：第一列 = quant 名（Name优先→ID→核心号兜底） =========
TX_USED="$OUT_DIR/logs/transcripts.used.list"
TX2GENE="$OUT_DIR/logs/tx2gene.final.tsv"
echo -e "transcript_id\tgene_id" > "$TX2GENE"

gawk -v TX="$TX_USED" -v N2G="$NAME2GENE" -v I2G="$ID2GENE" -v GL="$GENE_LABEL" -v OUT="$TX2GENE" '
  function core(s, t){ t=s; gsub(/\.[0-9]+$/,"",t); sub(/\|.*/,"",t);
                       if (match(t,/g[0-9]+/)) return substr(t,RSTART,RLENGTH); return "" }
  BEGIN{
    FS=OFS="\t"
    while((getline < GL)>0){ gid=$1; gl=$2; gid2label[gid]=gl } close(GL)
    while((getline < N2G)>0){ nm=$1; pg=$2; name2gene[nm]=pg } close(N2G)
    while((getline < I2G)>0){ id=$1; pg=$2; id2gene[id]=pg }   close(I2G)

    while((getline < TX)>0){
      t=$0; if(t=="") continue;
      gid=""
      if      (t in name2gene) gid=name2gene[t]
      else if (t in id2gene)   gid=id2gene[t]
      else {
        c=core(t)
        if (c!=""){
          for (k in name2gene){ if (core(k)==c){ gid=name2gene[k]; break } }
          if (gid=="") for (k in id2gene){ if (core(k)==c){ gid=id2gene[k]; break } }
        }
      }
      label = (gid in gid2label ? gid2label[gid] : (gid!=""?gid:t))
      print t, label >> OUT
    }
    close(TX)
  }'

# ========= 校验与统计（原逻辑保留） =========
TX2_FIRST="$OUT_DIR/logs/_tx2gene.first.list"
tail -n +2 "$OUT_DIR/logs/tx2gene.final.tsv" | cut -f1 | sort -u > "$TX2_FIRST"

SAMPLE="$OUT_DIR/logs/_tx_sample.list"
NS=$(wc -l < "$TX_USED" | tr -d ' ')
if [[ "$NS" -le 50 ]]; then
  cp "$TX_USED" "$SAMPLE"
else
  shuf -n 50 "$TX_USED" > "$SAMPLE"
fi

HIT=$(grep -Fxf "$SAMPLE" "$TX2_FIRST" | wc -l | tr -d ' ')
TOT=$(wc -l < "$SAMPLE" | tr -d ' ')
OKRATE=$(awk -v h="$HIT" -v t="$TOT" 'BEGIN{ if(t==0) print 0; else printf "%.2f", 100*h/t }')

echo "🧪 校验抽样命中：$HIT / $TOT (${OKRATE}%)"
if awk -v r="$OKRATE" 'BEGIN{ exit !(r<95) }'; then
  echo "⚠️  抽样命中率低于 95%，请留意注释/索引是否同源（继续运行）。"
fi

# ========= 覆盖率（原逻辑保留） =========
COVER="$OUT_DIR/logs/tx_coverage.tsv"
echo -e "sample\ttx_in_quant\ttx_in_tx2gene\ttx_mapped\ttx_coverage(%)" > "$COVER"
TX_IN_MAP=$(tail -n +2 "$TX2GENE" | cut -f1 | sort -u | wc -l)
while read -r qsf; do
  s=$(basename "$(dirname "$qsf")")
  qn=$(awk 'NR>1{print $1}' "$qsf" | sort -u | wc -l)
  mn=$(join -1 1 -2 1 <(awk 'NR>1{print $1}' "$qsf" | sort -u) <(tail -n +2 "$TX2GENE" | cut -f1 | sort -u) | wc -l)
  cov=$(awk -v m="$mn" -v q="$qn" 'BEGIN{printf (q? "%.2f":"0.00"), 100*m/q}')
  echo -e "$s\t$qn\t$TX_IN_MAP\t$mn\t$cov" >> "$COVER"
done < "$QSF_LIST"
echo "📊 已写出覆盖率表：$COVER"

# ========= 调用 R =========
Rscript "$(dirname "$0")/tximport_aggregate.R" \
  --quant_list "$QSF_LIST" \
  --tx2gene "$OUT_DIR/logs/tx2gene.final.tsv" \
  --outdir "$OUT_DIR" \
  --countsFromAbundance "$CFA"

