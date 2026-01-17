#!/usr/bin/env bash
set -euo pipefail

# =========================
# 参数区（按需改这里就行）
# =========================

# 1) ascp 可执行文件路径（你服务器上已有）
ASCP_BIN="/biostack/tools/utils/aspera-connect-4.2.6/bin/ascp"

# 2) Aspera 私钥（EBI public key，一般就是这个文件）
ASCP_KEY="/biostack/tools/utils/aspera-connect-4.2.6/etc/asperaweb_id_dsa.openssh"

# 3) EBI SRA Aspera 服务器信息（fastq 走这个最省事）
ASCP_USER="era-fasp"
ASCP_HOST="fasp.sra.ebi.ac.uk"
ASCP_PORT="33001"

# 4) 输入：一行一个 ftp 链接（ftp://ftp.sra.ebi.ac.uk/vol1/fastq/...）
URLS_FILE="urls.txt"

# 5) 输出目录（fastq.gz 会下载到这里）
OUT_DIR="."

# 6) 并发数（你要的：同时 3 个）
PARALLEL_JOBS=3

# 7) 限速（按需改；比如 0 表示不设上限；300m 表示 300 Mbps）
RATE_LIMIT="300m"

# 8) 断点续传级别：-k 1（推荐；更强校验可用 2/3，但会更慢）
RESUME_LEVEL=1

# 9) 是否显示速度/进度：1=显示；0=不显示（会写日志）
SHOW_PROGRESS=1

# 10) 是否跳过“已完整下载”的文件：1=跳过；0=每次都尝试（不推荐）
SKIP_COMPLETE=1

# 11) 对已存在的 .fastq.gz 做 gzip 完整性检查：1=检查（更稳）；0=不检查（更快）
CHECK_GZIP=1

# 12) 日志目录（SHOW_PROGRESS=0 时更有用）
LOG_DIR="logs_ascp"

# =========================
# 主流程（一般不用改）
# =========================

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [[ ! -s "$URLS_FILE" ]]; then
  echo "[ERROR] 找不到或为空：$URLS_FILE"
  exit 1
fi

if [[ ! -x "$ASCP_BIN" ]]; then
  echo "[ERROR] ASCP_BIN 不可执行：$ASCP_BIN"
  exit 1
fi

if [[ ! -f "$ASCP_KEY" ]]; then
  echo "[ERROR] ASCP_KEY 不存在：$ASCP_KEY"
  exit 1
fi

# 把 ftp 链接转为 ascp 远端路径：去掉域名，只保留 /vol1/fastq/...
PATHS_FILE="$(mktemp)"
trap 'rm -f "$PATHS_FILE"' EXIT

grep -v '^[[:space:]]*$' "$URLS_FILE" \
  | grep -v '^[[:space:]]*#' \
  | sed 's%^ftp://ftp\.sra\.ebi\.ac\.uk%%' \
  > "$PATHS_FILE"

download_one() {
  local remote_path="$1"

  # 远端路径必须以 / 开头（例如 /vol1/fastq/SRR...）
  if [[ "$remote_path" != /* ]]; then
    echo "[WARN] 跳过异常行（不是以 / 开头）：$remote_path"
    return 0
  fi

  local fname
  fname="$(basename "$remote_path")"
  local out_file="$OUT_DIR/$fname"
  local aspx_file="${out_file}.aspx"
  local log_file="$LOG_DIR/${fname}.log"

  # 已完整文件：直接跳过
  if [[ "$SKIP_COMPLETE" -eq 1 && -s "$out_file" && ! -f "$aspx_file" ]]; then
    if [[ "$CHECK_GZIP" -eq 1 ]]; then
      if gzip -t "$out_file" >/dev/null 2>&1; then
        echo "[SKIP] ok (gzip pass): $fname"
        return 0
      else
        echo "[REDO] gzip fail, re-download: $fname"
        rm -f "$out_file"
      fi
    else
      echo "[SKIP] exists: $fname"
      return 0
    fi
  fi

  # 组装 ascp 参数：
  # -T：禁用传输加密换吞吐 ；-k：断点续传；-l：限速；-P：SSH 端口
  # 注意：不要加 -q，否则会关闭进度/速度显示
  local -a cmd
  cmd=(
    "$ASCP_BIN"
    -T
    -k "$RESUME_LEVEL"
    -l "$RATE_LIMIT"
    -P "$ASCP_PORT"
    -i "$ASCP_KEY"
    "${ASCP_USER}@${ASCP_HOST}:${remote_path}"
    "$OUT_DIR"
  )

  if [[ "$SHOW_PROGRESS" -eq 1 ]]; then
    echo "[DL] $fname"
    "${cmd[@]}"
  else
    echo "[DL] $fname  (log: $log_file)"
    "${cmd[@]}" >"$log_file" 2>&1
  fi
}

export -f download_one
export ASCP_BIN ASCP_KEY ASCP_USER ASCP_HOST ASCP_PORT OUT_DIR RATE_LIMIT RESUME_LEVEL SHOW_PROGRESS SKIP_COMPLETE CHECK_GZIP LOG_DIR

echo "[INFO] paths: $(wc -l < "$PATHS_FILE")"
echo "[INFO] parallel: $PARALLEL_JOBS"
echo "[INFO] progress: $SHOW_PROGRESS (提示：ascp -q 会关闭进度显示)"

/usr/bin/env bash -c 'cat "$0" | xargs -r -n 1 -P "'"$PARALLEL_JOBS"'" bash -c "download_one \"\$1\"" _' "$PATHS_FILE"

echo "[DONE] all tasks finished."

# ========SRR urls auto download=======
# need srr.txt 
#rm -f urls.txt
#while read -r acc; do
  #curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
  #| awk -F'\t' 'NR==2{print $2}' \
  #| tr ';' '\n' \
  #| sed 's%^ftp\.sra\.ebi\.ac\.uk%ftp://ftp.sra.ebi.ac.uk%g'
#done < srr.txt | sed '/^ftp:\/\//!d' >> urls.txt

