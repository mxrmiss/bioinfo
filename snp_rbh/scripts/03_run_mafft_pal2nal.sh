#!/usr/bin/env bash
set -euo pipefail

# ----------------------------- 用户自定义区（皇上可改） -----------------------------
TOTAL_THREADS=25

JOBS=16
MAFFT_THREADS=1

PROGRESS_EVERY=100

PAL2NAL="scripts/pal2nal.pl"

PROT_DIR="results/step3_align_prot/pairs_faa"
CDS_DIR="results/step4_codon/pairs_fna"

OUT_AA_DIR="results/step3_align_prot/mafft_aa"
OUT_CODON_DIR="results/step4_codon/codon_aln"

LOG_DIR="logs"
# ------------------------------------------------------------------------------

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

mkdir -p "${OUT_AA_DIR}" "${OUT_CODON_DIR}" "${LOG_DIR}"

command -v mafft >/dev/null 2>&1 || { echo "ERROR: missing mafft" >&2; exit 1; }
command -v perl  >/dev/null 2>&1 || { echo "ERROR: missing perl" >&2; exit 1; }
[ -s "${PAL2NAL}" ] || { echo "ERROR: pal2nal not found: ${PAL2NAL}" >&2; exit 1; }

n_total=$(ls -1 "${PROT_DIR}"/*.faa 2>/dev/null | wc -l | awk '{print $1}')
if [ "${n_total}" -eq 0 ]; then
  echo "ERROR: no pair faa found in ${PROT_DIR}" >&2
  exit 1
fi

if [ "${JOBS}" -lt 1 ]; then JOBS=1; fi
if [ "${MAFFT_THREADS}" -lt 1 ]; then MAFFT_THREADS=1; fi

if [ $((JOBS * MAFFT_THREADS)) -gt "${TOTAL_THREADS}" ]; then
  JOBS=$(( TOTAL_THREADS / MAFFT_THREADS ))
  if [ "${JOBS}" -lt 1 ]; then JOBS=1; fi
fi

MAFFT_HAS_THREAD=0
if mafft --help 2>&1 | grep -q -- '--thread'; then
  MAFFT_HAS_THREAD=1
fi

echo "[03_run_mafft_pal2nal] total_pairs=${n_total}"
echo "[03_run_mafft_pal2nal] TOTAL_THREADS=${TOTAL_THREADS} JOBS=${JOBS} MAFFT_THREADS=${MAFFT_THREADS} (mafft_has_thread=${MAFFT_HAS_THREAD})"

run_one() {
  local faa="$1"
  local base
  base=$(basename "${faa}" .faa)

  local aa_out="${OUT_AA_DIR}/${base}.aa.aln"
  local codon_out="${OUT_CODON_DIR}/${base}.codon.fa"
  local mafft_log="${LOG_DIR}/mafft.${base}.log"
  local pal2nal_log="${LOG_DIR}/pal2nal.${base}.log"

  if [ "${MAFFT_HAS_THREAD}" -eq 1 ]; then
    mafft --auto --thread "${MAFFT_THREADS}" "${faa}" > "${aa_out}" 2> "${mafft_log}"
  else
    mafft --auto "${faa}" > "${aa_out}" 2> "${mafft_log}"
  fi

  perl "${PAL2NAL}" "${aa_out}" "${CDS_DIR}/${base}.fna" -output fasta > "${codon_out}" 2> "${pal2nal_log}"
}

export -f run_one
export OUT_AA_DIR OUT_CODON_DIR LOG_DIR CDS_DIR PAL2NAL MAFFT_THREADS MAFFT_HAS_THREAD

tmp_list=$(mktemp)
trap 'rm -f "${tmp_list}"' EXIT

ls -1 "${PROT_DIR}"/*.faa > "${tmp_list}"

tmp_done=$(mktemp)
trap 'rm -f "${tmp_list}" "${tmp_done}"' EXIT

touch "${tmp_done}"

worker() {
  while IFS= read -r faa; do
    run_one "${faa}"
    echo 1 >> "${tmp_done}"
    local done
    done=$(wc -l < "${tmp_done}" | awk '{print $1}')
    if [ $((done % PROGRESS_EVERY)) -eq 0 ] || [ "${done}" -eq "${n_total}" ]; then
      echo "[03_run_mafft_pal2nal] ${done}/${n_total} done"
    fi
  done
}

export -f worker
export tmp_done n_total PROGRESS_EVERY

if command -v xargs >/dev/null 2>&1; then
  cat "${tmp_list}" | xargs -I{} -P "${JOBS}" bash -c 'run_one "$1"; echo 1 >> "'"${tmp_done}"'"; done=$(wc -l < "'"${tmp_done}"'" | awk "{print \$1}"); if [ $((done % '"${PROGRESS_EVERY}"')) -eq 0 ] || [ "${done}" -eq "'"${n_total}"'" ]; then echo "[03_run_mafft_pal2nal] ${done}/'"${n_total}"' done"; fi' _ {}
else
  while IFS= read -r faa; do
    run_one "${faa}"
  done < "${tmp_list}"
  echo "[03_run_mafft_pal2nal] ${n_total}/${n_total} done"
fi

echo "[03_run_mafft_pal2nal] DONE"
