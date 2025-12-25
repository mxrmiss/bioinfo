#!/usr/bin/env bash
set -euo pipefail

############################################
# 13_run_genespace.sh
# —— 在 phylo 项目中一键运行 GENESPACE（25 线程）
#
# 关键修复：
# - 用 <<'RS' 传入 R 代码，避免 bash 展开 gpar$synteny 里的 $synteny
# - 用环境变量把 WD / THREADS / path2mcscanx 传给 R
############################################

# =========================
# 皇上参数区（只改这里）
# =========================

THREADS=25

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WD="${PROJECT_ROOT}/results/genespace"

CONDA_ENV_NAME="gspace"
CONDA_BASE="${HOME}/miniconda3"
MCSCANX_BIN="${CONDA_BASE}/envs/${CONDA_ENV_NAME}/bin"

AUTO_FIX_SPECIESIDS_SUFFIX="yes"
SUFFIX_TO_STRIP_REGEX='\.fa[[:space:]]*$'

LOG_DIR="${WD}/logs"

# =========================
# 工具函数
# =========================

log() {
  local msg="$1"
  echo "[INFO] ${msg}" | tee -a "${LOG_FILE}"
}

err() {
  local msg="$1"
  echo "[ERR] ${msg}" | tee -a "${LOG_FILE}" >&2
}

pick_latest_results_dir() {
  local parent="$1"
  local latest=""
  local d=""
  shopt -s nullglob
  for d in "${parent}"/Results_*; do
    [[ -d "$d" ]] || continue
    if [[ -z "$latest" ]]; then
      latest="$d"
    else
      if [[ "$d" -nt "$latest" ]]; then
        latest="$d"
      fi
    fi
  done
  shopt -u nullglob
  echo "$latest"
}

is_usable_orthofinder_run() {
  local results_dir="$1"
  [[ -d "$results_dir" ]] || return 1
  [[ -f "${results_dir}/WorkingDirectory/SpeciesIDs.txt" ]] || return 1
  [[ -f "${results_dir}/Orthogroups/Orthogroups.tsv" || -f "${results_dir}/Orthogroups/Orthogroups.tsv.gz" ]] || return 1
  return 0
}

fix_speciesids_suffix_if_needed() {
  local speciesids="$1"
  [[ -f "$speciesids" ]] || return 0

  if grep -qE "${SUFFIX_TO_STRIP_REGEX}" "$speciesids"; then
    if [[ "${AUTO_FIX_SPECIESIDS_SUFFIX}" != "yes" ]]; then
      err "检测到 SpeciesIDs.txt 中存在 .fa 后缀，但 AUTO_FIX_SPECIESIDS_SUFFIX=no：${speciesids}"
      return 1
    fi

    local ts
    ts="$(date +%Y%m%d_%H%M%S)"
    cp -a "$speciesids" "${speciesids}.bak_${ts}"

    sed -E -i 's/\.fa([[:space:]]*)$/\1/' "$speciesids"

    log "已自动修复 SpeciesIDs.txt：删除末尾 .fa 后缀（备份：${speciesids}.bak_${ts}）"
  else
    log "SpeciesIDs.txt 未检测到末尾 .fa 后缀：${speciesids}"
  fi
}

# =========================
# 初始化与基础检查
# =========================

mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/genespace.$(date +%Y%m%d_%H%M%S).log"

log "PROJECT_ROOT=${PROJECT_ROOT}"
log "WD=${WD}"
log "THREADS=${THREADS}"
log "CONDA_BASE=${CONDA_BASE}"
log "CONDA_ENV_NAME=${CONDA_ENV_NAME}"
log "MCSCANX_BIN=${MCSCANX_BIN}"
log "AUTO_FIX_SPECIESIDS_SUFFIX=${AUTO_FIX_SPECIESIDS_SUFFIX}"

if [[ ! -d "${WD}" ]]; then
  err "wd 不存在：${WD}"
  exit 2
fi

if [[ ! -d "${WD}/peptide" || ! -d "${WD}/bed" ]]; then
  err "缺少 genespace 输入目录：${WD}/peptide 或 ${WD}/bed"
  exit 3
fi

RSCRIPT_BIN="$(command -v Rscript || true)"
if [[ -z "${RSCRIPT_BIN}" ]]; then
  err "找不到系统 Rscript"
  exit 4
fi
log "System Rscript=${RSCRIPT_BIN}"

if [[ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
  err "找不到 conda.sh：${CONDA_BASE}/etc/profile.d/conda.sh"
  exit 5
fi

# =========================
# 激活 conda 环境（提供 orthofinder/diamond/MCScanX_h）
# =========================

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

log "which orthofinder: $(command -v orthofinder || true)"
log "which diamond: $(command -v diamond || true)"
log "which MCScanX_h: $(command -v MCScanX_h || true)"

if [[ -z "$(command -v orthofinder || true)" ]]; then
  err "环境中未找到 orthofinder"
  exit 6
fi
if [[ -z "$(command -v MCScanX_h || true)" ]]; then
  err "环境中未找到 MCScanX_h"
  exit 7
fi

export PATH="${MCSCANX_BIN}:${PATH}"

export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="1"
export MKL_NUM_THREADS="1"
export NUMEXPR_NUM_THREADS="1"

# =========================
# 自动修复 OrthoFinder SpeciesIDs.txt（仅在已有可复用结果时）
# =========================

OF_PARENT="${WD}/orthofinder"
if [[ -d "${OF_PARENT}" ]]; then
  OF_RESULTS_DIR="$(pick_latest_results_dir "${OF_PARENT}")"
  if [[ -n "${OF_RESULTS_DIR}" ]] && is_usable_orthofinder_run "${OF_RESULTS_DIR}"; then
    log "检测到可复用 OrthoFinder 结果：${OF_RESULTS_DIR}"
    fix_speciesids_suffix_if_needed "${OF_RESULTS_DIR}/WorkingDirectory/SpeciesIDs.txt"
  else
    log "未检测到可复用的 OrthoFinder Results_*（GENESPACE 可能会选择重跑）"
  fi
else
  log "未发现 ${OF_PARENT}（GENESPACE 将自行创建/运行 OrthoFinder）"
fi

# =========================
# 运行 GENESPACE（用系统 Rscript）
# =========================

log "Running GENESPACE..."

export GS_WD="${WD}"
export GS_MCSCANX_BIN="${MCSCANX_BIN}"
export GS_THREADS="${THREADS}"

"${RSCRIPT_BIN}" - <<'RS' 2>&1 | tee -a "${LOG_FILE}"
library(GENESPACE)

wd <- Sys.getenv("GS_WD")
path2mcscanx <- paste0(Sys.getenv("GS_MCSCANX_BIN"), "/")
threads <- as.integer(Sys.getenv("GS_THREADS"))

init_formals <- names(formals(GENESPACE::init_genespace))
init_args <- list(wd = wd, path2mcscanx = path2mcscanx)

if ("nCores" %in% init_formals) {
  init_args[["nCores"]] <- threads
} else if ("nc" %in% init_formals) {
  init_args[["nc"]] <- threads
} else if ("cores" %in% init_formals) {
  init_args[["cores"]] <- threads
}

gpar <- do.call(GENESPACE::init_genespace, init_args)

if (!"synteny" %in% names(gpar) || is.null(gpar$synteny)) {
  gpar$synteny <- list()
}

if (!"SpeciesIDs" %in% names(gpar$synteny) || is.null(gpar$synteny$SpeciesIDs) || length(gpar$synteny$SpeciesIDs) == 0) {
  gpar$synteny$SpeciesIDs <- NA
}

if (!"hogs" %in% names(gpar$synteny) || is.null(gpar$synteny$hogs) || length(gpar$synteny$hogs) == 0) {
  gpar$synteny$hogs <- NA
}

run_formals <- names(formals(GENESPACE::run_genespace))
run_args <- list(gsParam = gpar)

if ("nCores" %in% run_formals) {
  run_args[["nCores"]] <- threads
} else if ("nc" %in% run_formals) {
  run_args[["nc"]] <- threads
} else if ("cores" %in% run_formals) {
  run_args[["cores"]] <- threads
}

out <- do.call(GENESPACE::run_genespace, run_args)
saveRDS(out, file = file.path(wd, "genespace_run_output.rds"))
RS

log "DONE. Log saved to: ${LOG_FILE}"

