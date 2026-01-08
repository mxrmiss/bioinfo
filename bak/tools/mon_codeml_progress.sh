#!/usr/bin/env bash
ROOT="."
QC_KEEP="$ROOT/results/03_qc/keep_og.list"
SEQ_DIR="$ROOT/results/03_codon/codon_msa"
SETS_DIR="$ROOT/results/04_codeml/sets"
RAW_DIR="$ROOT/results/04_codeml/raw"
LOG_FILE="$ROOT/logs/06_codeml_batch_branchsite.log"
INTERVAL=10
# 是否持续刷新:1 yes; 2 no
LOOP=0

set -euo pipefail

calc() {
  if [ -s "$QC_KEEP" ]; then
    OG_N=$(grep -vc '^[[:space:]]*$' "$QC_KEEP")
  else
    OG_N=$(find "$SEQ_DIR" -maxdepth 1 -type f -name 'OG*.codon.fna' 2>/dev/null | wc -l)
  fi
  FG_N=$(ls "$SETS_DIR"/*.nwk 2>/dev/null | wc -l)
  NEED=$((OG_N*FG_N))

  ALT_LIST=$(grep -R "lnL(ntime:" -l "$RAW_DIR"/*/*/alt/mlc.txt 2>/dev/null | sed 's#/alt/mlc.txt##' | LC_ALL=C sort -u || true)
  NULL_LIST=$(grep -R "lnL(ntime:" -l "$RAW_DIR"/*/*/null/mlc.txt 2>/dev/null | sed 's#/null/mlc.txt##' | LC_ALL=C sort -u || true)

  ALT_DONE=$(printf "%s\n" "$ALT_LIST" | sed '/^$/d' | wc -l)
  NULL_DONE=$(printf "%s\n" "$NULL_LIST" | sed '/^$/d' | wc -l)
  PAIR_DONE=$(comm -12 <(printf "%s\n" "$ALT_LIST") <(printf "%s\n" "$NULL_LIST") | sed '/^$/d' | wc -l)

  TOTAL_ENUM=$(find "$RAW_DIR" -mindepth 2 -maxdepth 2 -type d 2>/dev/null | LC_ALL=C sort -u | wc -l)
  PLAN=$(grep -o 'CMD\[[^]]*/alt\]' "$LOG_FILE" 2>/dev/null | sed 's/^CMD\[//;s/\/alt\]$//' | LC_ALL=C sort -u | wc -l)

  RUNNING=$(pgrep -af codeml 2>/dev/null | wc -l)
  RECENT=$(find "$RAW_DIR" -type f -name mlc.txt -mmin -1 2>/dev/null | wc -l)
  DONE_FLAG="RUNNING"; [ -f "$ROOT/results/04_codeml/.codeml.done" ] && DONE_FLAG="DONE"

  PROG_PAIRS=$(awk -v a="$PAIR_DONE" -v n="$NEED" 'BEGIN{if(n>0)printf("%.2f",100*a/n);else printf("0.00")}')
  SUM=$((ALT_DONE+NULL_DONE))
  TARGET=$((2*NEED))
  PROG_SUM=$(awk -v s="$SUM" -v t="$TARGET" 'BEGIN{if(t>0)printf("%.2f",100*s/t);else printf("0.00")}')

  printf "og=%s fg=%s need_pairs=%s plan_pairs=%s enum_pairs=%s alt_done=%s null_done=%s pair_done=%s progress_pairs=%s%% progress_sum=%s%% running=%s recent_writes=%s status=%s\n" \
    "$OG_N" "$FG_N" "$NEED" "$PLAN" "$TOTAL_ENUM" "$ALT_DONE" "$NULL_DONE" "$PAIR_DONE" "$PROG_PAIRS" "$PROG_SUM" "$RUNNING" "$RECENT" "$DONE_FLAG"
}

if [ "$LOOP" -eq 1 ]; then
  while true; do
    clear
    date +"%F %T"
    calc
    sleep "$INTERVAL"
  done
else
  calc
fi

