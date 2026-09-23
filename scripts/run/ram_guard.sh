#!/usr/bin/env bash
CEIL_GB=${CEIL_GB:-36}
WARN_GB=${WARN_GB:-30}
LOG=/12TBDrive1/mega_pep_bench/runs/ram_guard.log
PRIORITY=(
  "pepscorer"        # feature calc: restartable per-pose, cheapest to lose
  "inference.py"     # RAPiDock: resumable per complex
  "score_fresh.py"   # DockQ: fully cached per complex
)
echo "[$(date '+%F %T')] guard up: ceiling ${CEIL_GB}G, warn ${WARN_GB}G" >> "$LOG"
while true; do
  U=$(free -g | awk '/Mem:/{print $3}')
  if [ "${U:-0}" -ge "$CEIL_GB" ]; then
    echo "[$(date '+%F %T')] RAM ${U}G >= ${CEIL_GB}G" >> "$LOG"
    for pat in "${PRIORITY[@]}"; do
      if pgrep -f "$pat" >/dev/null 2>&1; then
        echo "  killing '$pat'" >> "$LOG"
        pkill -f "$pat"
        break                      # one at a time, then re-measure
      fi
    done
    sleep 20
  elif [ "${U:-0}" -ge "$WARN_GB" ]; then
    echo "[$(date '+%F %T')] warn: RAM ${U}G" >> "$LOG"
    sleep 30
  fi
  sleep 10
done
