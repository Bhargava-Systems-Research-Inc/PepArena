#!/bin/bash
cd /12TBDrive1/mega_pep_bench || exit 1
IB=/14TBDrive/6TBDrive1_backup/ImmunoBench/ESMFold2Slow_Bench
export BIOHUB_API_TOKEN=$(grep BIOHUB_API_TOKEN "$IB/.env" | cut -d= -f2)
LOG=runs/track_a/logs/esmfold2_daily.log
echo "[$(date '+%F %T')] esmfold2 daily start ($(find runs/track_a/predictions/esmfold2 -name '*_model_0.cif'|wc -l)/334 done)" >> "$LOG"

timeout 9000 ~/anaconda3/envs/esmfold2_api/bin/python scripts/run/track_a/run_esmfold2_api.py \
  --input-csv runs/track_a/esmfold2_inputs_canonical.csv \
  --out-dir runs/track_a/predictions/esmfold2 \
  --summary-csv runs/track_a/predictions/esmfold2/summary.csv \
  --model esmfold2-2026-05 --num-loops 20 --num-sampling-steps 100 \
  --max-new 90 >> "$LOG" 2>&1

~/anaconda3/bin/python scripts/run/track_a/score_dockq.py --model esmfold2 >> "$LOG" 2>&1
~/anaconda3/bin/python scripts/run/track_a/build_leaderboard.py >> "$LOG" 2>&1
echo "[$(date '+%F %T')] esmfold2 daily done ($(find runs/track_a/predictions/esmfold2 -name '*_model_0.cif'|wc -l)/334 done)" >> "$LOG"
