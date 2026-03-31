#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# 可调参数区
# ============================================================

MAX_JOBS=2
THREADS_QUAST=4
THREADS_BUSCO=12
BUSCO_LINEAGE="mollusca_odb12"

# BUSCO 数据库固定目录
BUSCO_DOWNLOAD_DIR="busco_downloads/lineages"

# 0 = 必须使用本地已存在数据库，不允许自动下载
# 1 = 若本地不存在，则先自动下载一次
ALLOW_BUSCO_DOWNLOAD_IF_MISSING=0

# BUSCO 专用 Java 内存
# 用于规避 bbtools / AssemblyStats2 的 Java heap space 问题
BUSCO_JAVA_XMS="1g"
BUSCO_JAVA_XMX="8g"

# BUSCO 输入 FASTA 统一换行宽度
BUSCO_WRAP_WIDTH=60

# 1 = 每个物种 BUSCO 结束后删除临时目录
# 0 = 保留临时目录，便于排错
CLEANUP_BUSCO_TMP=1

# 0 = 若 BUSCO 已完成则跳过
# 1 = 强制重跑 BUSCO
FORCE_RERUN_BUSCO=0

# ============================================================
# 固定目录
# ============================================================

GENOMIC_DIR="data/genomic"
ANNOTATION_DIR="data/annotation"
CDS_DIR="data/cds"
PROTEOME_DIR="data/proteomes"

QC_DIR="results/qc"
QUAST_DIR="${QC_DIR}/quast"
BUSCO_DIR="${QC_DIR}/busco"
AGAT_DIR="${QC_DIR}/agat"
SEQKIT_DIR="${QC_DIR}/seqkit"
SUPP_DIR="results/supp_table"

# BUSCO 临时目录根目录
BUSCO_TMP_ROOT="${BUSCO_DIR}/_tmp"

mkdir -p \
  "${QUAST_DIR}" \
  "${BUSCO_DIR}" \
  "${AGAT_DIR}" \
  "${SEQKIT_DIR}" \
  "${SUPP_DIR}" \
  "${BUSCO_DOWNLOAD_DIR}" \
  "${BUSCO_TMP_ROOT}"

# ============================================================
# 基础检查
# ============================================================

command -v quast >/dev/null 2>&1 || { echo "[ERROR] quast not found"; exit 1; }
command -v busco >/dev/null 2>&1 || { echo "[ERROR] busco not found"; exit 1; }
command -v agat_sp_statistics.pl >/dev/null 2>&1 || { echo "[ERROR] agat_sp_statistics.pl not found"; exit 1; }
command -v seqkit >/dev/null 2>&1 || { echo "[ERROR] seqkit not found"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "[ERROR] python not found"; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo "[ERROR] gzip not found"; exit 1; }
command -v find >/dev/null 2>&1 || { echo "[ERROR] find not found"; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "[ERROR] awk not found"; exit 1; }
command -v grep >/dev/null 2>&1 || { echo "[ERROR] grep not found"; exit 1; }
command -v sed >/dev/null 2>&1 || { echo "[ERROR] sed not found"; exit 1; }

# ============================================================
# 日志函数
# ============================================================

log_info() {
  echo "[INFO] $*"
}

log_warn() {
  echo "[WARN] $*" >&2
}

log_error() {
  echo "[ERROR] $*" >&2
}

# ============================================================
# 物种识别
# ============================================================

mapfile -t species_list < <(
  find "${GENOMIC_DIR}" -maxdepth 1 -type f \
  | sed 's#^.*/##' \
  | grep -E '\.(fa|fna|fasta)(\.gz)?$' \
  | sed -E 's/\.(fa|fna|fasta)(\.gz)?$//' \
  | sort -u
)

if [ "${#species_list[@]}" -eq 0 ]; then
  log_error "no species detected from ${GENOMIC_DIR}"
  exit 1
fi

log_info "detected species count = ${#species_list[@]}"
printf '%s\n' "${species_list[@]}"

# ============================================================
# 文件定位函数
# ============================================================

find_unique_file() {
  local species="$1"
  local label="$2"
  local search_dir="$3"
  shift 3

  local files=()
  while IFS= read -r line; do
    [ -n "$line" ] && files+=("$line")
  done < <(find "${search_dir}" -maxdepth 1 -type f "$@")

  if [ "${#files[@]}" -eq 0 ]; then
    log_error "no ${label} file found for ${species}"
    return 1
  fi

  if [ "${#files[@]}" -gt 1 ]; then
    log_error "multiple ${label} files found for ${species}"
    printf '%s\n' "${files[@]}" >&2
    return 1
  fi

  printf '%s\n' "${files[0]}"
}

find_genome_file() {
  local species="$1"
  find_unique_file "${species}" "genome" "${GENOMIC_DIR}" \
    \( -name "${species}.fa" -o -name "${species}.fa.gz" -o -name "${species}.fna" -o -name "${species}.fna.gz" -o -name "${species}.fasta" -o -name "${species}.fasta.gz" \)
}

find_gff_file() {
  local species="$1"
  find_unique_file "${species}" "annotation" "${ANNOTATION_DIR}" \
    \( -name "${species}.gff" -o -name "${species}.gff.gz" -o -name "${species}.gff3" -o -name "${species}.gff3.gz" \)
}

find_cds_file() {
  local species="$1"
  find_unique_file "${species}" "cds" "${CDS_DIR}" \
    \( -name "${species}.fa" -o -name "${species}.fa.gz" -o -name "${species}.fna" -o -name "${species}.fna.gz" -o -name "${species}.fasta" -o -name "${species}.fasta.gz" \)
}

find_protein_file() {
  local species="$1"
  find_unique_file "${species}" "proteome" "${PROTEOME_DIR}" \
    \( -name "${species}.faa" -o -name "${species}.faa.gz" -o -name "${species}.fa" -o -name "${species}.fa.gz" -o -name "${species}.fasta" -o -name "${species}.fasta.gz" \)
}

# ============================================================
# BUSCO 数据集检查
# ============================================================

detect_busco_lineage_dir() {
  local root="$1"
  local lineage="$2"

  local candidates=(
    "${root}/lineages/${lineage}"
    "${root}/${lineage}"
    "${root}/busco_downloads/lineages/${lineage}"
    "${root}/busco_downloads/${lineage}"
  )

  local x
  for x in "${candidates[@]}"; do
    if [ -d "${x}" ]; then
      printf '%s\n' "${x}"
      return 0
    fi
  done

  return 1
}

prepare_busco_dataset() {
  local lineage="$1"
  local download_dir="$2"

  local lineage_dir=""
  if lineage_dir=$(detect_busco_lineage_dir "${download_dir}" "${lineage}"); then
    log_info "BUSCO lineage found locally: ${lineage_dir}" >&2
    printf '%s\n' "${lineage_dir}"
    return 0
  fi

  if [ "${ALLOW_BUSCO_DOWNLOAD_IF_MISSING}" -eq 1 ]; then
    log_info "BUSCO lineage not found locally, downloading once: ${lineage}" >&2
    busco --download "${lineage}" --download_path "${download_dir}"
    lineage_dir=$(detect_busco_lineage_dir "${download_dir}" "${lineage}") || {
      log_error "BUSCO lineage download finished but dataset directory still not found: ${lineage}"
      exit 1
    }
    log_info "BUSCO lineage prepared: ${lineage_dir}" >&2
    printf '%s\n' "${lineage_dir}"
    return 0
  fi

  log_error "BUSCO lineage not found locally: ${lineage}"
  log_error "expected under: ${download_dir}"
  log_error "please download it first, or set ALLOW_BUSCO_DOWNLOAD_IF_MISSING=1"
  exit 1
}

BUSCO_LINEAGE_PATH=$(prepare_busco_dataset "${BUSCO_LINEAGE}" "${BUSCO_DOWNLOAD_DIR}")
log_info "BUSCO_LINEAGE_PATH=${BUSCO_LINEAGE_PATH}"

# ============================================================
# 完成状态判定
# ============================================================

is_done_quast() {
  local sp="$1"
  [ -s "${QUAST_DIR}/${sp}/report.tsv" ]
}

find_busco_summary_file_for_species() {
  local sp="$1"
  local d="${BUSCO_DIR}/${sp}"

  [ -d "${d}" ] || return 1

  local f=""
  f=$(find "${d}" -maxdepth 3 -type f -name 'short_summary*.txt' | head -n 1 || true)

  if [ -n "${f}" ] && [ -s "${f}" ]; then
    printf '%s\n' "${f}"
    return 0
  fi

  return 1
}

is_done_busco() {
  local sp="$1"

  if [ "${FORCE_RERUN_BUSCO}" -eq 1 ]; then
    return 1
  fi

  local summary_file=""
  if ! summary_file=$(find_busco_summary_file_for_species "${sp}"); then
    return 1
  fi

  local run_dir=""
  run_dir=$(find "${BUSCO_DIR}/${sp}" -maxdepth 1 -mindepth 1 -type d -name 'run_*' | head -n 1 || true)

  if [ -n "${run_dir}" ] && [ -d "${run_dir}" ]; then
    return 0
  fi

  return 1
}

is_done_agat() {
  local sp="$1"
  [ -s "${AGAT_DIR}/${sp}/stats.txt" ]
}

is_done_seqkit_genome() {
  local sp="$1"
  [ -s "${SEQKIT_DIR}/${sp}.genome.tsv" ]
}

is_done_seqkit_protein() {
  local sp="$1"
  [ -s "${SEQKIT_DIR}/${sp}.protein.tsv" ]
}

# ============================================================
# BUSCO 输入预处理
# ============================================================

cleanup_busco_tmp_dir() {
  local tmp_dir="$1"

  if [ "${CLEANUP_BUSCO_TMP}" -eq 1 ] && [ -n "${tmp_dir}" ] && [ -d "${tmp_dir}" ]; then
    rm -rf "${tmp_dir}"
  fi
}

validate_fasta_file() {
  local fasta="$1"

  [ -s "${fasta}" ] || {
    log_error "FASTA file is empty: ${fasta}"
    return 1
  }

  grep -m 1 '^>' "${fasta}" >/dev/null 2>&1 || {
    log_error "FASTA file does not look valid (no header found): ${fasta}"
    return 1
  }
}

prepare_busco_input_fasta() {
  local sp="$1"
  local genome="$2"
  local prep_dir="$3"

  mkdir -p "${prep_dir}"

  local wrapped_fa="${prep_dir}/${sp}.busco.wrap.fa"

  if [[ "${genome}" =~ \.gz$ ]]; then
    gzip -dc "${genome}" | seqkit seq -w "${BUSCO_WRAP_WIDTH}" > "${wrapped_fa}"
  else
    seqkit seq -w "${BUSCO_WRAP_WIDTH}" "${genome}" > "${wrapped_fa}"
  fi

  validate_fasta_file "${wrapped_fa}"
  printf '%s\n' "${wrapped_fa}"
}

# ============================================================
# 核心流程：单物种
# ============================================================

run_one_species() {
  local sp="$1"

  local genome=""
  local gff=""
  local cds=""
  local protein=""

  local busco_species_dir="${BUSCO_DIR}/${sp}"
  local busco_tmp_dir="${BUSCO_TMP_ROOT}/${sp}"
  local busco_prep_dir="${busco_tmp_dir}/prep"
  local busco_tmp_logs_dir="${busco_tmp_dir}/logs"

  mkdir -p "${QUAST_DIR}/${sp}" "${AGAT_DIR}/${sp}"

  log_info "start ${sp}"

  genome=$(find_genome_file "${sp}")
  gff=$(find_gff_file "${sp}")
  cds=$(find_cds_file "${sp}")
  protein=$(find_protein_file "${sp}")

  # ----------------------------
  # QUAST
  # ----------------------------
  if is_done_quast "${sp}"; then
    log_info "${sp} | QUAST skipped"
  else
    log_info "${sp} | QUAST"
    quast "${genome}" --threads "${THREADS_QUAST}" -o "${QUAST_DIR}/${sp}"
  fi

  # ----------------------------
  # BUSCO
  # ----------------------------
  if is_done_busco "${sp}"; then
    log_info "${sp} | BUSCO skipped"
  else
    log_info "${sp} | BUSCO"

    rm -rf "${busco_species_dir}"
    rm -rf "${busco_tmp_dir}"
    mkdir -p "${busco_prep_dir}" "${busco_tmp_logs_dir}"

    local busco_input_fa=""
    busco_input_fa=$(prepare_busco_input_fasta "${sp}" "${genome}" "${busco_prep_dir}")
    validate_fasta_file "${busco_input_fa}"

    local busco_stdout="${busco_tmp_logs_dir}/busco.stdout.log"
    local busco_stderr="${busco_tmp_logs_dir}/busco.stderr.log"
    local busco_cmdlog="${busco_tmp_logs_dir}/busco.command.txt"

    {
      echo "species=${sp}"
      echo "original_genome=${genome}"
      echo "busco_input=${busco_input_fa}"
      echo "busco_lineage=${BUSCO_LINEAGE_PATH}"
      echo "threads_busco=${THREADS_BUSCO}"
      echo "java_xms=${BUSCO_JAVA_XMS}"
      echo "java_xmx=${BUSCO_JAVA_XMX}"
    } > "${busco_cmdlog}"

    if ! env _JAVA_OPTIONS="-Xms${BUSCO_JAVA_XMS} -Xmx${BUSCO_JAVA_XMX}" \
      busco \
        -i "${busco_input_fa}" \
        -o "${sp}" \
        -m genome \
        -l "${BUSCO_LINEAGE_PATH}" \
        --cpu "${THREADS_BUSCO}" \
        --out_path "${BUSCO_DIR}" \
        --download_path "${BUSCO_DOWNLOAD_DIR}" \
        --offline \
        -f \
        > >(tee "${busco_stdout}") \
        2> >(tee "${busco_stderr}" >&2); then

      log_error "${sp} | BUSCO failed"
      log_error "${sp} | BUSCO stdout log: ${busco_stdout}"
      log_error "${sp} | BUSCO stderr log: ${busco_stderr}"
      log_error "${sp} | last 20 lines of stderr:"
      tail -n 20 "${busco_stderr}" >&2 || true

      cleanup_busco_tmp_dir "${busco_tmp_dir}"
      return 1
    fi

    if ! is_done_busco "${sp}"; then
      log_error "${sp} | BUSCO command exited but final summary file not found"
      log_error "${sp} | BUSCO stdout log: ${busco_stdout}"
      log_error "${sp} | BUSCO stderr log: ${busco_stderr}"

      cleanup_busco_tmp_dir "${busco_tmp_dir}"
      return 1
    fi

    mkdir -p "${busco_species_dir}/logs"
    cp -f "${busco_stdout}" "${busco_species_dir}/logs/" || true
    cp -f "${busco_stderr}" "${busco_species_dir}/logs/" || true
    cp -f "${busco_cmdlog}" "${busco_species_dir}/logs/" || true

    cleanup_busco_tmp_dir "${busco_tmp_dir}"
  fi

  # ----------------------------
  # AGAT
  # ----------------------------
  if is_done_agat "${sp}"; then
    log_info "${sp} | AGAT skipped"
  else
    log_info "${sp} | AGAT"
    agat_sp_statistics.pl --gff "${gff}" > "${AGAT_DIR}/${sp}/stats.txt" 2>&1
  fi

  # ----------------------------
  # SeqKit genome
  # ----------------------------
  if is_done_seqkit_genome "${sp}"; then
    log_info "${sp} | SeqKit genome skipped"
  else
    log_info "${sp} | SeqKit genome"
    seqkit stats -Ta "${genome}" > "${SEQKIT_DIR}/${sp}.genome.tsv"
  fi

  # ----------------------------
  # SeqKit protein
  # ----------------------------
  if is_done_seqkit_protein "${sp}"; then
    log_info "${sp} | SeqKit protein skipped"
  else
    log_info "${sp} | SeqKit protein"
    seqkit stats -Ta "${protein}" > "${SEQKIT_DIR}/${sp}.protein.tsv"
  fi

  log_info "done ${sp}"
}

# ============================================================
# 并发调度
# ============================================================

running=0
failed=0

for sp in "${species_list[@]}"; do
  run_one_species "${sp}" &
  running=$((running + 1))

  if [ "${running}" -ge "${MAX_JOBS}" ]; then
    if ! wait -n; then
      failed=$((failed + 1))
    fi
    running=$((running - 1))
  fi
done

while [ "${running}" -gt 0 ]; do
  if ! wait -n; then
    failed=$((failed + 1))
  fi
  running=$((running - 1))
done

if [ "${failed}" -gt 0 ]; then
  log_error "one or more species failed: failed_jobs=${failed}"
  exit 1
fi

log_info "all per-species QC finished"

# ============================================================
# 生成补充表
# ============================================================

if [ -f "scripts/genome_supp_table.py" ]; then
  log_info "build supplementary genome table"
  python scripts/genome_supp_table.py
else
  log_warn "scripts/genome_supp_table.py not found, skip supplementary table generation"
fi

log_info "all tasks finished"
