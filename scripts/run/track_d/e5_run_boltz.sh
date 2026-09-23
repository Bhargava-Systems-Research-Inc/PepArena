#!/usr/bin/env bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
E5=$REPO/runs/track_d/e5_zeroshot
IN=${1:-$E5/yaml}
OUT=${2:-$E5/out}
LOG=$E5/run.log
FLOOR_GB=${FLOOR_GB:-10}
mkdir -p "$OUT"

conda run -n boltz --no-capture-output boltz predict "$IN" \
  --out_dir "$OUT" \
  --accelerator gpu --devices 1 \
  --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 \
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
echo "[$(date '+%F %T')] boltz exit=$rc  confidences=$(find "$OUT" -name 'confidence_*.json' | wc -l)" | tee -a "$LOG"
exit $rc
