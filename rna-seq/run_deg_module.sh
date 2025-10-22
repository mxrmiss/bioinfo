#!/usr/bin/env bash
set -euo pipefail

#############################################
# F 模块运行参数（仅改 sh；R 主体 f_deg.R 不动）
#############################################
MATRIX_DIR="results/matrix"

# 样本与对比设计
SAMPLES_TSV="samples.tsv"       # 必含列：sample, group（可选 batch）
CONTRASTS_TSV="contrasts.tsv"   # 两列：case, control

# 输出目录
OUT_DEG="results/deg"
OUT_PCA="results/pca"
LOG_DIR="results/logs"

# 统计阈值（与 f_deg.R 顶部 CFG 字段保持一致）
PADJ_CUTOFF="0.05"
LFC_CUTOFF="1"
MIN_COUNT="10"
TRANSFORM="vst"    # vst 或 rlog
USE_BATCH="false"  # 是否启用 batch（samples.tsv 需有 batch 列）

#############################################
# 运行前准备
#############################################
mkdir -p "${OUT_DEG}" "${OUT_PCA}" "${LOG_DIR}"

[[ -f "${SAMPLES_TSV}"   ]] || { echo "[F][致命] 缺少 ${SAMPLES_TSV}" >&2; exit 1; }
[[ -f "${CONTRASTS_TSV}" ]] || { echo "[F][致命] 缺少 ${CONTRASTS_TSV}" >&2; exit 1; }
[[ -d "${MATRIX_DIR}"    ]] || { echo "[F][致命] 缺少目录 ${MATRIX_DIR}" >&2; exit 1; }

echo "[F] 探测 ${MATRIX_DIR} 下的候选矩阵文件..."
ls -lh "${MATRIX_DIR}" || true

# 优先使用 txi RDS（若存在）
TXI_RDS="$(ls -1 "${MATRIX_DIR}"/txi*lengthScaledTPM*.rds 2>/dev/null | head -n1 || true)"

# 常见的 counts 命名
CANDIDATES=(
  "${MATRIX_DIR}/matrix_gene_counts.tsv"
  "${MATRIX_DIR}/counts.tsv"
  "${MATRIX_DIR}/gene_counts.tsv"
)

COUNTS_TSV=""
for f in "${CANDIDATES[@]}"; do
  if [[ -f "$f" ]]; then
    COUNTS_TSV="$f"
    break
  fi
done

if [[ -n "${TXI_RDS}" && -f "${TXI_RDS}" ]]; then
  MODE="txi"
  echo "[F] 使用 tximport RDS 对象: ${TXI_RDS}"
elif [[ -n "${COUNTS_TSV}" ]]; then
  MODE="counts"
  echo "[F] 使用 counts 矩阵: ${COUNTS_TSV}"
else
  echo "[F][致命] 未在 ${MATRIX_DIR} 找到以下任一文件：" >&2
  echo "        - matrix_gene_counts.tsv" >&2
  echo "        - counts.tsv" >&2
  echo "        - gene_counts.tsv" >&2
  echo "      或 txi*lengthScaledTPM*.rds" >&2
  exit 1
fi

# 友好预警：估算每组最小生物重复数（k），k<2 提醒（不拦截）
if command -v awk >/dev/null 2>&1; then
  MIN_REP="$(awk -F '\t' 'NR>1{cnt[$2]++} END{min=1e9; for(g in cnt){if(cnt[g]<min)min=cnt[g]} if(min==1e9)min=0; print min}' "${SAMPLES_TSV}")"
  if [[ "${MIN_REP}" -lt 2 ]]; then
    echo "[F][提示] 检测到每组最小生物重复数 k=${MIN_REP}。DESeq2 可能无法估计离散度；"
    echo "         如遇停止，可先仅生成 PCA/热图或补充生物重复后再跑差异分析。"
  fi
fi

#############################################
# 生成临时 runner（不改 f_deg.R 本体，仅在内存替换 CFG）
#############################################
export MODE TXI_RDS COUNTS_TSV SAMPLES_TSV CONTRASTS_TSV OUT_DEG OUT_PCA PADJ_CUTOFF LFC_CUTOFF MIN_COUNT TRANSFORM USE_BATCH
MODE="${MODE}"
TXI_RDS="${TXI_RDS:-}"
COUNTS_TSV="${COUNTS_TSV:-}"

RUNNER="${LOG_DIR}/f_runner.R"

# 用单引号 heredoc，避免 bash 展开；在 R 内读取 Sys.getenv，规避转义问题
cat > "${RUNNER}" <<'RSCRIPT'
getenv <- function(k, def="") {
  v <- Sys.getenv(k, unset = def)
  # 统一反斜杠为正斜杠，避免 R 字符串转义问题（如 \s、\t）
  v <- gsub("\\\\", "/", v)
  # 转义双引号，安全注入
  v <- gsub('"', '\\"', v, fixed = TRUE)
  v
}

suppressWarnings({
  code <- readLines("f_deg.R", warn = FALSE)
})

# 定位 CFG <- list( ... ) 段
start <- grep("^\\s*CFG\\s*<-\\s*list\\(", code)
if (length(start) != 1) stop("[runner] 无法唯一定位 CFG 起始行。")

cnt <- 0; close <- NA_integer_
for (k in start:length(code)) {
  opens  <- gregexpr("\\(", code[k])[[1]]; opens  <- if (opens[1] == -1) 0 else length(opens)
  closes <- gregexpr("\\)", code[k])[[1]]; closes <- if (closes[1] == -1) 0 else length(closes)
  cnt <- cnt + opens - closes
  if (cnt == 0) { close <- k; break }
}
if (is.na(close)) stop("[runner] 未找到 CFG 结束括号。")

mode        <- getenv("MODE")
txi_rds     <- getenv("TXI_RDS")
counts_tsv  <- getenv("COUNTS_TSV")
samples     <- getenv("SAMPLES_TSV")
contrasts   <- getenv("CONTRASTS_TSV")
outdeg      <- getenv("OUT_DEG")
outpca      <- getenv("OUT_PCA")
padj        <- Sys.getenv("PADJ_CUTOFF", unset = "0.05")
lfc         <- Sys.getenv("LFC_CUTOFF",  unset = "1")
minCount    <- Sys.getenv("MIN_COUNT",   unset = "10")
transform   <- getenv("TRANSFORM")
use_batch   <- Sys.getenv("USE_BATCH",   unset = "false")

# —— 修正点：将布尔开关规范成 R 识别的 TRUE/FALSE —— #
use_batch_tok <- if (tolower(use_batch) %in% c("1","true","t","yes","y")) "TRUE" else "FALSE"

new_cfg <- c(
  "CFG <- list(",
  sprintf('  mode        = "%s",', mode),
  sprintf('  txi_rds     = "%s",', txi_rds),
  sprintf('  counts_tsv  = "%s",', counts_tsv),
  sprintf('  samples     = "%s",', samples),
  sprintf('  contrasts   = "%s",', contrasts),
  sprintf('  outdeg      = "%s",', outdeg),
  sprintf('  outpca      = "%s",', outpca),
  sprintf('  padj        = %s,',  padj),
  sprintf('  lfc         = %s,',  lfc),
  sprintf('  minCount    = %s,',  minCount),
  sprintf('  transform   = "%s",', transform),
  sprintf('  use_batch   = %s',   use_batch_tok),
  ")"
)

code2 <- c(code[1:(start-1)], new_cfg, code[(close+1):length(code)])

# 打印当前 CFG（供日志核对）
cat("[runner] 当前 CFG 已注入：\n")
cat(paste(new_cfg, collapse = "\n"), "\n\n")

# 执行修改后的脚本主体
eval(parse(text = code2), envir = .GlobalEnv)
RSCRIPT

#############################################
# 执行 runner 并记录日志
#############################################
echo "[F] 调用 R 执行 f_deg.R（由 runner 注入 CFG）..."
Rscript "${RUNNER}" | tee "${LOG_DIR}/f_module_session.txt"
echo "[F] 运行完毕。日志见：${LOG_DIR}/f_module_session.txt"
