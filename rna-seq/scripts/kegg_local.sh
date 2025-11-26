#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# kegg_local.sh — 本地构建 KEGG 词典（遵循约定：表头与落点）
# 产物：
#   results/07_annot/kegg/term2gene.tsv  （pathway_id\tgene_id）
#   results/07_annot/kegg/term2name.tsv  （pathway_id\tterm_name）
# 依赖：
#   ref/annotations/ko_to_pathway.tsv
#   ref/annotations/kegg_pathway.tsv
#   results/07_annot/gene2ko.tsv  （两列：gene_id\tko 或 gene_id\tko_id）
# =============================================================================

REF="ref/annotations"
OUTDIR="results/07_annot/kegg"
RAW_PATH="$REF/kegg_offline/ko_to_pathway.raw"
RAW_MOD="$REF/kegg_offline/ko_to_module.raw"
MANUAL="$REF/ko_to_pathway.tsv.manual"

KO2PATH_TSV="$REF/ko_to_pathway.tsv"
PATHNAME_TSV="$REF/kegg_pathway.tsv"
GENE2KO="results/07_annot/gene2ko.tsv"
T2G_PW="$OUTDIR/term2gene.tsv"
T2N_PW="$OUTDIR/term2name.tsv"
LOG="$REF/.kegg_build.log"

mkdir -p "$REF/kegg_offline" "$REF/kegg_legacy" "$OUTDIR"
: > "$LOG"
log(){ echo "$*" | tee -a "$LOG"; }

# 规范化 RAW -> TSV（两列：Kxxxxx \t koYYYYY）
norm_ko_path(){
  awk -vOFS='\t' '{
    gsub(/\r/,"");
    k=$1; p=$2;
    gsub(/^ko:/,"",k); k=toupper(k);
    gsub(/^path:/,"",p); sub(/^map/,"ko",p);
    if (k ~ /^K[0-9]{5}$/ && p ~ /^ko[0-9]{5}$/) print k,p;
  }' "$1" | LC_ALL=C sort -u
}

# A) ko_to_pathway.tsv
if [[ -s "$MANUAL" ]]; then
  log "[1] 使用手工基准映射：$MANUAL"
  cp -f "$MANUAL" "$KO2PATH_TSV"
else
  if [[ ! -s "$RAW_PATH" ]]; then
    log "[1] 下载 raw: KO→Pathway"
    curl -sSL --retry 3 --retry-delay 2 "https://rest.kegg.jp/link/pathway/ko" -o "$RAW_PATH"
  else
    log "[1] 使用本地 raw: $RAW_PATH"
  fi
  log "[1] 规范化为 $KO2PATH_TSV"
  norm_ko_path "$RAW_PATH" > "$KO2PATH_TSV"
fi
log "[✔] ko_to_pathway.tsv: pairs=$(wc -l < "$KO2PATH_TSV")  KO=$(cut -f1 "$KO2PATH_TSV"|LC_ALL=C sort -u|wc -l)  pathways=$(cut -f2 "$KO2PATH_TSV"|LC_ALL=C sort -u|wc -l)"

# B) pathway_names.tsv → 确保带表头（关键修复：按 TAB 分列，保留完整描述）
if [[ ! -s "$PATHNAME_TSV" ]]; then
  log "[2] 获取 pathway_names.tsv"
  if curl -sSL --max-time 60 "https://rest.kegg.jp/list/pathway/ko" \
    | awk -vFS='\t' -vOFS='\t' '{
        gsub(/\r/,"");
        gsub(/^path:/,"",$1);
        sub(/^map/,"ko",$1);
        if ($1 ~ /^ko[0-9]{5}$/) print $1,$2
      }' \
    | LC_ALL=C sort -u > "$PATHNAME_TSV"; then
    :
  else
    # 在线失败兜底：ID=Name
    cut -f2 "$KO2PATH_TSV" | LC_ALL=C sort -u | awk -vOFS='\t' '{print $1,$1}' > "$PATHNAME_TSV"
  fi
fi

# 补表头（pathway_id\tname）
if [[ "$(head -n1 "$PATHNAME_TSV")" != $'pathway_id\tname' ]]; then
  { echo -e "pathway_id\tname"; cat "$PATHNAME_TSV"; } > "$REF/.kegg_pathway.with_header.tsv"
  mv -f "$REF/.kegg_pathway.with_header.tsv" "$PATHNAME_TSV"
fi
log "[✔] pathway_names.tsv: $(wc -l < "$PATHNAME_TSV") 行"

# C) gene→KO 预处理（CRLF + 自动列定位 + KO 规范化；不改基因ID）
if [[ ! -s "$GENE2KO" ]]; then
  log "[❌] 缺少 $GENE2KO"
  exit 2
fi

# 统一去 CRLF（安全幂等）
sed -i 's/\r$//' "$GENE2KO"

# 生成 .gene2ko.clean
awk '
  BEGIN{FS="\t"; OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      t=tolower($i)
      if(t=="gene_id"||t=="gene"||t=="id") gc=i
      if(t=="ko_id" ||t=="ko")             kc=i
    }
    if(!gc) gc=1; if(!kc) kc=2; next
  }
  {
    g=$gc
    n=split($kc,arr,/[;, ]+/)
    for(i=1;i<=n;i++){
      x=arr[i]
      sub(/^ko:/,"",x)
      gsub(/[ \t]/,"",x)
      x=toupper(x)
      if(x ~ /^K[0-9]{5}$/ && g!="") print g,x
    }
  }
' "$GENE2KO" | LC_ALL=C sort -u > "$REF/.gene2ko.clean"

log "[✔] gene2ko.clean: pairs=$(wc -l < "$REF/.gene2ko.clean")  genes=$(cut -f1 "$REF/.gene2ko.clean"|LC_ALL=C sort -u|wc -l)  KO=$(cut -f2 "$REF/.gene2ko.clean"|LC_ALL=C sort -u|wc -l)"

# D) join：Pathway→gene
tmp_t2g="$REF/.term2gene_kegg_pathway.tmp"
join -t $'\t' -1 1 -2 2 \
  <(LC_ALL=C sort -k1,1 "$KO2PATH_TSV") \
  <(LC_ALL=C sort -k2,2 "$REF/.gene2ko.clean") \
| awk -vOFS='\t' '{print $2,$3}' | LC_ALL=C sort -u > "$tmp_t2g"

{ echo -e "pathway_id\tgene_id"; cat "$tmp_t2g"; } > "$T2G_PW"
rm -f "$tmp_t2g"

if [[ $(wc -l < "$T2G_PW") -le 1 ]]; then
  log "[❌] term2gene.tsv 为空"
  exit 3
else
  log "[✔] term2gene.tsv: $(wc -l < "$T2G_PW") 行"
fi

# E) 名称字典落点（保持流水线约定：pathway_id\tterm_name）
{ echo -e "pathway_id\tterm_name"; tail -n +2 "$PATHNAME_TSV"; } > "$T2N_PW"
log "[✔] term2name.tsv: $(wc -l < "$T2N_PW") 行"

