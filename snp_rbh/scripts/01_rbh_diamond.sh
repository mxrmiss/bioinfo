#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 01_rbh_diamond.sh
# 功能：
#   1) 两物种蛋白 DIAMOND blastp：A→B、B→A（TSV带表头，真正Tab分隔）
#   2) 过滤（evalue + qcovhsp/scovhsp 双向覆盖度）（TSV带表头）
#   3) 每个 query 取 best hit（按 bitscore 最大）（TSV带表头）
#   4) 取 RBH 得到 1:1 pairs（TSV带表头，第一列永远是参考物种ID）
#
# 特点：
#   - 物种缩写自动从文件名提取：genus 首字母 + species 前两字母（全小写）
#     例：Sinonovacula_constricta → sco；Sinonovacula_rivularis → sri
#   - 所有目录与路径均使用相对路径（脚本自动 cd 到项目根目录）
# ==============================================================================

# ----------------------------- 用户自定义区（皇上可改） -----------------------------
THREADS=16

EVALUE="1e-20"
MIN_QCOV=70
MIN_SCOV=70
MAX_TARGET_SEQS=200
SENS_MODE="--more-sensitive"

# 物种A（推荐把参考物种放A，读起来更直观；但REFERENCE也可选B）
A_PROTEOME="data/proteomes/Sinonovacula_constricta.faa"
A_CDS="data/cds/Sinonovacula_constricta.fna"
A_GFF="data/annotation/Sinonovacula_constricta.gff3.gz"

# 物种B
B_PROTEOME="data/proteomes/Sinonovacula_rivularis.faa"
B_CDS="data/cds/Sinonovacula_rivularis.fna"
B_GFF="data/annotation/Sinonovacula_rivularis.gff.gz"

# 参考物种：填 "A" 或 "B"
REFERENCE="B"
# ------------------------------------------------------------------------------

# 自动定位项目根目录（保证相对路径可用）
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

mkdir -p results/step1_diamond results/step2_rbh logs

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }
}

need_cmd diamond
need_cmd awk
need_cmd sort
need_cmd head
need_cmd printf
need_cmd basename
need_cmd tr

# 从文件名推断物种缩写：genus首字母 + species前两字母（全小写）
# 期望 basename 形如 Genus_species.xxx(.gz)
infer_code_from_path() {
  local p="$1"
  local bn
  bn="$(basename "$p")"
  bn="${bn%.gz}"
  bn="${bn%.gff3}"
  bn="${bn%.gff}"
  bn="${bn%.faa}"
  bn="${bn%.fna}"
  bn="${bn%.fa}"
  bn="${bn%.fasta}"

  local genus species
  genus="$(echo "$bn" | awk -F'_' '{print $1}')"
  species="$(echo "$bn" | awk -F'_' '{print $2}')"

  if [ -z "${genus}" ] || [ -z "${species}" ]; then
    echo "ERROR: cannot infer species code from filename: ${bn} (need Genus_species.*)" >&2
    exit 1
  fi

  local g1 s2
  g1="$(echo "$genus" | awk '{print substr($0,1,1)}' | tr 'A-Z' 'a-z')"
  s2="$(echo "$species" | awk '{print substr($0,1,2)}' | tr 'A-Z' 'a-z')"
  echo "${g1}${s2}"
}

A_CODE="$(infer_code_from_path "${A_PROTEOME}")"
B_CODE="$(infer_code_from_path "${B_PROTEOME}")"

# 输出文件命名使用缩写
A_DB="results/step1_diamond/${A_CODE}.dmnd"
B_DB="results/step1_diamond/${B_CODE}.dmnd"

A_TO_B_RAW="results/step1_diamond/${A_CODE}_to_${B_CODE}.raw.tsv"
B_TO_A_RAW="results/step1_diamond/${B_CODE}_to_${A_CODE}.raw.tsv"

A_TO_B_FILT="results/step1_diamond/${A_CODE}_to_${B_CODE}.filt.tsv"
B_TO_A_FILT="results/step1_diamond/${B_CODE}_to_${A_CODE}.filt.tsv"

A_BEST="results/step2_rbh/${A_CODE}_best_to_${B_CODE}.tsv"
B_BEST="results/step2_rbh/${B_CODE}_best_to_${A_CODE}.tsv"

A2B_MAP="results/step2_rbh/${A_CODE}2${B_CODE}.map.tsv"
B2A_MAP="results/step2_rbh/${B_CODE}2${A_CODE}.map.tsv"

RBH_PAIRS="results/step2_rbh/rbh_pairs.tsv"

LOG_A2B="logs/diamond.${A_CODE}_to_${B_CODE}.log"
LOG_B2A="logs/diamond.${B_CODE}_to_${A_CODE}.log"

echo "[INFO] A_CODE=${A_CODE}  B_CODE=${B_CODE}  REFERENCE=${REFERENCE}"
echo "[INFO] EVALUE=${EVALUE}  MIN_QCOV=${MIN_QCOV}  MIN_SCOV=${MIN_SCOV}  THREADS=${THREADS}"

# 1) 建库
diamond makedb --in "${A_PROTEOME}" --db "${A_DB}"
diamond makedb --in "${B_PROTEOME}" --db "${B_DB}"

# 2) blastp（TSV带表头；注意这里写入的是真正Tab，不会出现字面 \t）
printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovhsp\tscovhsp\n" > "${A_TO_B_RAW}"
diamond blastp --query "${A_PROTEOME}" --db "${B_DB}" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp \
  --evalue "${EVALUE}" --max-target-seqs "${MAX_TARGET_SEQS}" --threads "${THREADS}" ${SENS_MODE} \
  >> "${A_TO_B_RAW}" 2> "${LOG_A2B}"

printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovhsp\tscovhsp\n" > "${B_TO_A_RAW}"
diamond blastp --query "${B_PROTEOME}" --db "${A_DB}" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp \
  --evalue "${EVALUE}" --max-target-seqs "${MAX_TARGET_SEQS}" --threads "${THREADS}" ${SENS_MODE} \
  >> "${B_TO_A_RAW}" 2> "${LOG_B2A}"

# 3) 覆盖度过滤（带表头）
head -n 1 "${A_TO_B_RAW}" > "${A_TO_B_FILT}"
awk -F'\t' -v q="${MIN_QCOV}" -v s="${MIN_SCOV}" 'NR==1{next} ($15+0)>=q && ($16+0)>=s' "${A_TO_B_RAW}" >> "${A_TO_B_FILT}"

head -n 1 "${B_TO_A_RAW}" > "${B_TO_A_FILT}"
awk -F'\t' -v q="${MIN_QCOV}" -v s="${MIN_SCOV}" 'NR==1{next} ($15+0)>=q && ($16+0)>=s' "${B_TO_A_RAW}" >> "${B_TO_A_FILT}"

# 4) best hit（按 bitscore 最大）（带表头）
printf "qseqid\tsseqid\tbitscore\tevalue\tqcovhsp\tscovhsp\n" > "${A_BEST}"
awk -F'\t' 'NR==1{next}{print $1"\t"$2"\t"$12"\t"$11"\t"$15"\t"$16}' "${A_TO_B_FILT}" \
  | sort -k1,1 -k3,3gr \
  | awk -F'\t' '!seen[$1]++{print}' >> "${A_BEST}"

printf "qseqid\tsseqid\tbitscore\tevalue\tqcovhsp\tscovhsp\n" > "${B_BEST}"
awk -F'\t' 'NR==1{next}{print $1"\t"$2"\t"$12"\t"$11"\t"$15"\t"$16}' "${B_TO_A_FILT}" \
  | sort -k1,1 -k3,3gr \
  | awk -F'\t' '!seen[$1]++{print}' >> "${B_BEST}"

# 5) map（带表头）
printf "%s_id\t%s_id\n" "${A_CODE}" "${B_CODE}" > "${A2B_MAP}"
awk -F'\t' 'NR==1{next}{print $1"\t"$2}' "${A_BEST}" | sort -u >> "${A2B_MAP}"

printf "%s_id\t%s_id\n" "${B_CODE}" "${A_CODE}" > "${B2A_MAP}"
awk -F'\t' 'NR==1{next}{print $1"\t"$2}' "${B_BEST}" | sort -u >> "${B2A_MAP}"

# 6) RBH（带表头；第一列永远是参考物种）
tmp_kv="$(mktemp)"
trap 'rm -f "${tmp_kv}"' EXIT

# kv: B_id \t A_id
awk -F'\t' 'NR==1{next}{print $1"\t"$2}' "${B2A_MAP}" > "${tmp_kv}"

if [ "${REFERENCE}" = "A" ]; then
  printf "%s_id\t%s_id\n" "${A_CODE}" "${B_CODE}" > "${RBH_PAIRS}"
  awk -F'\t' 'NR==1{next}{print $1"\t"$2}' "${A2B_MAP}" \
    | awk -F'\t' 'NR==FNR{rev[$1]=$2;next}{a=$1;b=$2; if(rev[b]==a) print a"\t"b}' "${tmp_kv}" - \
    | sort -u >> "${RBH_PAIRS}"
elif [ "${REFERENCE}" = "B" ]; then
  printf "%s_id\t%s_id\n" "${B_CODE}" "${A_CODE}" > "${RBH_PAIRS}"
  awk -F'\t' 'NR==1{next}{print $1"\t"$2}' "${A2B_MAP}" \
    | awk -F'\t' 'NR==FNR{rev[$1]=$2;next}{a=$1;b=$2; if(rev[b]==a) print a"\t"b}' "${tmp_kv}" - \
    | awk -F'\t' '{print $2"\t"$1}' \
    | sort -u >> "${RBH_PAIRS}"
else
  echo "ERROR: REFERENCE must be A or B" >&2
  exit 1
fi

echo "[DONE] Outputs:"
echo "  ${A_TO_B_RAW}"
echo "  ${B_TO_A_RAW}"
echo "  ${A_TO_B_FILT}"
echo "  ${B_TO_A_FILT}"
echo "  ${A_BEST}"
echo "  ${B_BEST}"
echo "  ${A2B_MAP}"
echo "  ${B2A_MAP}"
echo "  ${RBH_PAIRS}"
