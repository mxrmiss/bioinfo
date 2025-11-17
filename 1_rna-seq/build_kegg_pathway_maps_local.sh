#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# build_kegg_pathway_maps_local.sh  v9.2 — 自愈构建器
# - 手工基准映射优先：ref/ko_to_pathway.tsv.manual 若存在，直接用，不覆盖
# - 若无手工表：若缺 RAW 自动下载，再规范化为 ko_to_pathway.tsv
# - 生成：term2gene_kegg_pathway.tsv / term2name_kegg_pathway.tsv
# - 可选：KEGG-MODULE（term2gene_kegg_module.tsv / term2name_kegg_module.tsv）
# - 规范化严格且不丢行：只在正确列做去前缀/大小写/正则校验，多 KO 拆分，CR/LANG 兼容
# =============================================================================

# ==== 开关 ====
ENABLE_MODULE=${ENABLE_MODULE:-1}   # 1=生成 MODULE，0=跳过

# ==== 路径 ====
REF="ref"
RAW_PATH="$REF/kegg_offline/ko_to_pathway.raw"
RAW_MOD="$REF/kegg_offline/ko_to_module.raw"
MANUAL="$REF/ko_to_pathway.tsv.manual"

KO2PATH_TSV="$REF/ko_to_pathway.tsv"
PATHNAME_TSV="$REF/pathway_names.tsv"

T2G_PW="$REF/term2gene_kegg_pathway.tsv"
T2N_PW="$REF/term2name_kegg_pathway.tsv"

T2G_MOD="$REF/term2gene_kegg_module.tsv"
T2N_MOD="$REF/term2name_kegg_module.tsv"

GENE2KO="$REF/kegg_gene2ko.tsv"
LOG="$REF/.kegg_build.log"

mkdir -p "$REF/kegg_offline"
: > "$LOG"
log(){ echo "$*" | tee -a "$LOG"; }

# ==== 函数：规范化 RAW -> TSV（两列：Kxxxxx \t koYYYYY）====
norm_ko_path(){
  awk -vOFS='\t' '{
    gsub(/\r/,"");
    k=$1; p=$2;
    gsub(/^ko:/,"",k); k=toupper(k);                  # 仅 KO 列
    gsub(/^path:/,"",p); sub(/^map/,"ko",p);         # 仅 Path 列
    if (k ~ /^K[0-9]{5}$/ && p ~ /^ko[0-9]{5}$/) print k,p;
  }' "$1" | LC_ALL=C sort -u
}

# ==== A) 准备 ko_to_pathway.tsv ====
if [[ -s "$MANUAL" ]]; then
  log "[1] 检测到手工基准映射：$MANUAL  （直接使用，不覆盖）"
  cp -f "$MANUAL" "$KO2PATH_TSV"
else
  if [[ ! -s "$RAW_PATH" ]]; then
    log "[1] 自动下载 raw: KO→Pathway"
    curl -sSL --retry 3 --retry-delay 2 "https://rest.kegg.jp/link/pathway/ko" -o "$RAW_PATH"
  else
    log "[1] 使用本地 raw: $RAW_PATH"
  fi
  log "[1] 规范化为 $KO2PATH_TSV"
  norm_ko_path "$RAW_PATH" > "$KO2PATH_TSV"
fi

pairs=$(wc -l < "$KO2PATH_TSV")
kos=$(cut -f1 "$KO2PATH_TSV" | LC_ALL=C sort -u | wc -l)
paths=$(cut -f2 "$KO2PATH_TSV" | LC_ALL=C sort -u | wc -l)
log "[✔] ko_to_pathway.tsv: pairs=$pairs  KO=$kos  pathways=$paths"

# ==== B) pathway_names.tsv ====
if [[ ! -s "$PATHNAME_TSV" ]]; then
  log "[2] 获取 pathway_names.tsv"
  if curl -sSL --max-time 60 "https://rest.kegg.jp/list/pathway/ko" \
    | awk -vOFS='\t' '{gsub(/\r/,""); gsub(/^path:/,"",$1); sub(/^map/,"ko",$1); if ($1 ~ /^ko[0-9]{5}$/) print $1,$2}' \
    | LC_ALL=C sort -u > "$PATHNAME_TSV"; then
    :
  fi
  if [[ ! -s "$PATHNAME_TSV" ]]; then
    log "[W] 在线失败，使用 ID=Name 兜底"
    cut -f2 "$KO2PATH_TSV" | LC_ALL=C sort -u | awk -vOFS='\t' '{print $1,$1}' > "$PATHNAME_TSV"
  fi
fi
log "[✔] pathway_names.tsv: $(wc -l < "$PATHNAME_TSV") 行"

# ==== C) gene→KO 预处理：任意分隔符→TAB、多 KO 拆分、清洗 ====
if [[ ! -s "$GENE2KO" ]]; then
  log "[致命] 缺少 $GENE2KO"
  exit 1
fi
log "[3] 预处理 gene2ko"
awk '
  BEGIN{FS="[\t, ]+"; OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){low=tolower($i); if(low=="gene_id"||low=="gene"||low=="id") gc=i;
                       if(low=="ko_id"||low=="ko") kc=i}
    if(!gc) gc=1; if(!kc) kc=2; next
  }
  {
    g=$gc; sub(/\|.*/,"",g); sub(/\.[0-9]+$/,"",g);
    n=split($kc,arr,/[;, ]+/);
    for(i=1;i<=n;i++){
      x=arr[i]; gsub(/^ko:/,"",x); x=toupper(x);
      if (x ~ /^K[0-9]{5}$/ && g!="") print g,x;
    }
  }
' "$GENE2KO" | LC_ALL=C sort -u > "$REF/.gene2ko.clean"

g2k_pairs=$(wc -l < "$REF/.gene2ko.clean")
g2k_genes=$(cut -f1 "$REF/.gene2ko.clean" | LC_ALL=C sort -u | wc -l)
g2k_kos=$(cut -f2 "$REF/.gene2ko.clean"  | LC_ALL=C sort -u | wc -l)
log "[✔] gene2ko.clean: pairs=$g2k_pairs  genes=$g2k_genes  KO=$g2k_kos"

# ==== D) join：Pathway→gene ====
log "[4] 连接 KO→Pathway × gene→KO → Pathway→gene"
join -t $'\t' -1 1 -2 2 \
  <(LC_ALL=C sort -k1,1 "$KO2PATH_TSV") \
  <(LC_ALL=C sort -k2,2 "$REF/.gene2ko.clean") \
| awk -vOFS='\t' '{print $2,$3}' | LC_ALL=C sort -u > "$T2G_PW"

t2g_lines=$(wc -l < "$T2G_PW")
if [[ $t2g_lines -eq 0 ]]; then
  log "[❌] term2gene_kegg_pathway.tsv 为空，打印未命中 KO 样例："
  comm -23 \
    <(cut -f2 "$REF/.gene2ko.clean" | LC_ALL=C sort -u) \
    <(cut -f1 "$KO2PATH_TSV"       | LC_ALL=C sort -u) | head -n 10 | sed 's/^/[KO 未映射] /' | tee -a "$LOG"
else
  log "[✔] term2gene_kegg_pathway.tsv: $t2g_lines 行"
fi
cp -f "$PATHNAME_TSV" "$T2N_PW"
log "[✔] term2name_kegg_pathway.tsv: $(wc -l < "$T2N_PW") 行"

# ==== E) 可选：MODULE ====
if [[ "$ENABLE_MODULE" -eq 1 ]]; then
  log "[5] 生成 KEGG-MODULE（可用于补充分析）"
  if [[ ! -s "$RAW_MOD" ]]; then
    log "[5] 自动下载 raw: KO→MODULE"
    curl -sSL --retry 3 --retry-delay 2 "https://rest.kegg.jp/link/module/ko" -o "$RAW_MOD" || true
  fi
  if [[ -s "$RAW_MOD" ]]; then
    awk -vOFS='\t' '{
      gsub(/\r/,"");
      k=$1; m=$2;
      gsub(/^ko:/,"",k); k=toupper(k);
      gsub(/^md:/,"",m);
      if (k ~ /^K[0-9]{5}$/ && m ~ /^M[0-9]+$/) print k,m;
    }' "$RAW_MOD" | LC_ALL=C sort -u > "$REF/ko_to_module.tsv"

    join -t $'\t' -1 1 -2 2 \
      <(LC_ALL=C sort -k1,1 "$REF/ko_to_module.tsv") \
      <(LC_ALL=C sort -k2,2 "$REF/.gene2ko.clean") \
    | awk -vOFS='\t' '{print $2,$3}' | LC_ALL=C sort -u > "$T2G_MOD"

    cut -f2 "$REF/ko_to_module.tsv" | LC_ALL=C sort -u \
      | awk -vOFS='\t' '{print $1,$1}' > "$T2N_MOD"

    log "[✔] term2gene_kegg_module.tsv: $(wc -l < "$T2G_MOD") 行"
    log "[✔] term2name_kegg_module.tsv: $(wc -l < "$T2N_MOD") 行"
  else
    log "[W] KO→MODULE 原始文件缺失，跳过 MODULE"
  fi
fi

log "[完成] 构建日志：$LOG"
