#!/usr/bin/env bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
B=$REPO/runs/track_c/boltz2_affinity
IN=${1:-$B/all_yaml}
OUT=${2:-$B/out_final}
LOG=$B/run_final.log

FLOOR_GB=${FLOOR_GB:-8}
mkdir -p "$OUT"

conda run -n boltz --no-capture-output boltz predict "$IN" \
  --out_dir "$OUT" \
  --accelerator gpu --devices 1 \
  --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 \
  --sampling_steps_affinity 200 --diffusion_samples_affinity 5 \
  --output_format mmcif \
  --num_workers 2 --preprocessing-threads 4 \
  >>"$LOG" 2>&1 &
BOLTZ=$!
echo "[$(date '+%F %T')] boltz pid=$BOLTZ in=$IN out=$OUT floor=${FLOOR_GB}G" | tee -a "$LOG"

while kill -0 $BOLTZ 2>/dev/null; do
  AVAIL=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
  if [ "$AVAIL" -lt "$FLOOR_GB" ]; then
    echo "[$(date '+%F %T')] ABORT: MemAvailable ${AVAIL}G < ${FLOOR_GB}G" | tee -a "$LOG"
    pkill -TERM -P $BOLTZ; kill -TERM $BOLTZ; sleep 20; kill -KILL $BOLTZ 2>/dev/null
    exit 3
  fi
  sleep 30
done
wait $BOLTZ; rc=$?
echo "[$(date '+%F %T')] boltz exit=$rc  predictions=$(find "$OUT" -name 'affinity_*.json' | wc -l)" | tee -a "$LOG"
exit $rc
