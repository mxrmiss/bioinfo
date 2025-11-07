#!/usr/bin/env bash
set -euo pipefail

USER_NAME="${USER:-$(id -un)}"
PATTERNS=(
  "snakemake"
  "python .*snakemake"
  "orthofinder"
  "/diamond\b|^diamond\b"   # diamond blastp/makedb 等
  "\biqtree2?\b"            # iqtree 或 iqtree2
  "\bFastTree(MP)?\b|\bfasttree\b"
  "\bmafft\b"
  "\bastral\b"
)

echo "==> 终止 Snakemake 及其进程组（含子进程）"
for spid in $(pgrep -u "$USER_NAME" -f "snakemake" || true); do
  PGID=$(ps -o pgid= -p "$spid" | tr -d ' ')
  if [[ -n "$PGID" ]]; then
    echo "  -> kill 进程组 PGID=$PGID (leader PID=$spid)"
    # 先优雅终止，再强杀
    kill -TERM -"$PGID" 2>/dev/null || true
  fi
done
sleep 1

echo "==> 逐类结束常见工具进程（优雅终止）"
for pat in "${PATTERNS[@]}"; do
  pkill -u "$USER_NAME" -f -e "$pat" 2>/dev/null || true
done
sleep 2

echo "==> 仍存活则强制结束（KILL）"
for pat in "${PATTERNS[@]}"; do
  pkill -9 -u "$USER_NAME" -f -e "$pat" 2>/dev/null || true
done
sleep 1

echo "==> 剩余相关进程（若有）："
ps -o pid,ppid,pgid,stat,etime,cmd -u "$USER_NAME" \
  | egrep -E "snakemake|orthofinder|diamond|iqtree2?|FastTree|fasttree|mafft|astral" || true

echo "==> 若仍有遗留 diamond，可用以下诊断："
echo "  pstree -aps <PID>"
echo "  ps -o pid,ppid,pgid,stat,etime,cmd -p <PID>"
echo "  lsof -p <PID> | head"

