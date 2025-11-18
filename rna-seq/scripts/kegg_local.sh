#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# kegg_local.sh  v9.2 — 自愈构建器（最小化口径统一：表头与落点）
# - 手工基准映射优先：ref/ko_to_pathway.tsv.manual 若存在，直接用，不覆盖
# - 若无手工表：若缺 RAW 自动下载，再规范化为 ko_to_pathway.tsv
# - 生成：results/07_annot/kegg/term2gene.tsv / results/07_annot/kegg/term2name.tsv
# - 可选：KEGG-MODULE（term2gene_kegg_module.tsv / term2name_kegg_module.tsv）
# - 统一表头：
#     term2name.tsv                 -> pathway_id\tterm_name
#     term2gene.tsv                  -> pathway_id\tgene_id
# =============================================================================

# ==== 开关 ====
ENABLE_MODULE=1   # 1=生成 MODULE，0=跳过

# ==== 路径 ====
REF="ref/annotations"
OUTDIR="results/07_annot/kegg"
RAW_PATH="$REF/kegg_offline/ko_to_pathway.raw"
RAW_MOD="$REF/kegg_offline/ko_to_module.raw"
MANUAL="$REF/ko_to_pathway.tsv.manual"

KO2PATH_TSV="$REF/ko_to_pathway.tsv"
PATHNAME_TSV="$REF/kegg_pathway.tsv"

# —— 改动：输入路径改为 results/07_annot —— 
GENE2KO="results/07_annot/gene2ko.tsv"  # [修改] 输入路径改为正确的路径
T2G_PW="$OUTDIR/term2gene.tsv"  # 输出路径改为 results/07_annot/kegg
T2N_PW="$OUTDIR/term2name.tsv"  # 输出路径改为 results/07_annot/kegg

T2G_MOD="$REF/kegg_legacy/term2gene_kegg_module.tsv"
T2N_MOD="$REF/kegg_legacy/term2name_kegg_module.tsv"

LOG="$REF/.kegg_build.log"

mkdir -p "$REF/kegg_offline" "$REF/kegg_legacy" "$OUTDIR"
: > "$LOG"
log(){ echo "$*" | tee -a "$LOG"; }

# ==== 函数：规范化 RAW -> TSV（两列：Kxxxxx \t koYYYYY）====
norm_ko_path(){
  awk -vOFS='\t' '{
    gsub(/\r/,"");
    k=$1; p=$2;
    gsub(/^ko:/,"",k); k=toupper(k);
    gsub(/^path:/,"",p); sub(/^map/,"ko",p);
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
  else
    log "[W] 在线获取失败，将用 ko_to_pathway.tsv 的路径集合构造名称占位"
    cut -f2 "$KO2PATH_TSV" | LC_ALL=C sort -u | awk -vOFS='\t' '{print $1,$1}' > "$PATHNAME_TSV"
  fi
fi

# —— 改动①：确保 kegg_pathway.tsv 带表头 —— 
if [[ -s "$PATHNAME_TSV" ]]; then
  if [[ "$(head -n1 "$PATHNAME_TSV")" != $'pathway_id\tname' ]]; then
    { echo -e "pathway_id\tname"; cat "$PATHNAME_TSV"; } > "$REF/.kegg_pathway.with_header.tsv"
    mv -f "$REF/.kegg_pathway.with_header.tsv" "$PATHNAME_TSV"
  fi
fi
log "[✔] pathway_names.tsv: $(wc -l < "$PATHNAME_TSV") 行"

# ==== C) gene→KO 预处理：任意分隔符→TAB、多 KO 拆分、清洗 ====
if [[ ! -s "$GENE2KO" ]]; then
  log "[3] 未发现 gene2ko.tsv：$GENE2KO"
  log "[❌] 缺少 gene2ko.tsv，无法展开 Pathway→gene。请先准备 results/07_annot/gene2ko.tsv"
  exit 2
fi

# —— 修正：严格解析 gene2ko.tsv，确保正确生成 gene2ko.clean 文件 —— 
awk -vOFS='\t' '
  NR==1 { next }  # 忽略表头
  NF>=2 {
    g=$1; kc=$2;
    gsub(/\r/,"",g); gsub(/^[ \t]+|[ \t]+$/,"",g);
    if (g=="") next;
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
tmp_t2g="$REF/.term2gene_kegg_pathway.tmp"
join -t $'\t' -1 1 -2 2 \
  <(LC_ALL=C sort -k1,1 "$KO2PATH_TSV") \
  <(LC_ALL=C sort -k2,2 "$REF/.gene2ko.clean") \
| awk -vOFS='\t' '{print $2,$3}' | LC_ALL=C sort -u > "$tmp_t2g"

# —— 改动②：确保 term2gene 表头与落点一致 —— 
{ echo -e "pathway_id\tgene_id"; cat "$tmp_t2g"; } > "$T2G_PW"
rm -f "$tmp_t2g"

t2g_lines=$(wc -l < "$T2G_PW")
if [[ $t2g_lines -le 1 ]]; then
  log "[❌] term2gene.tsv 为空，打印未命中 KO 样例："
  comm -23 \
    <(cut -f2 "$REF/.gene2ko.clean" | LC_ALL=C sort -u) \
    <(cut -f1 "$KO2PATH_TSV"       | LC_ALL=C sort -u) | head -n 10 | sed 's/^/[KO 未映射] /' | tee -a "$LOG"
else
  log "[✔] term2gene.tsv: $t2g_lines 行"
fi

# —— 改动③：名称字典落到约定目录，并改成 term_name 表头 —— 
{ echo -e "pathway_id\tterm_name"; tail -n +2 "$PATHNAME_TSV"; } > "$T2N_PW"
log "[✔] term2name.tsv: $(wc -l < "$T2N_PW") 行"

# ==== E) 可选：MODULE（原逻辑不动，仍落在 legacy 目录） ====
if [[ "$ENABLE_MODULE" -eq 1 ]]; then
  log "[5] 生成 KEGG-MODULE（可用于补充分析）"
  if [[ -s "$RAW_MOD" ]]; then
    if [[ ! -s "$REF/ko_to_module.tsv" ]]; then
      log "[5] 规范化 RAW: KO→MODULE"
      awk -vOFS='\t' '{
        gsub(/\r/,"");
        k=$1; m=$2;
        gsub(/^ko:/,"",k); k=toupper(k);
        gsub(/^md:/,"",m); m=toupper(m);
        if (k ~ /^K[0-9]{5}$/ && m ~ /^M[0-9]{5}$/) print k,m;
      }' "$RAW_MOD" | LC_ALL=C sort -u > "$REF/ko_to_module.tsv"
    fi

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