#!/usr/bin/env bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
OUT=$REPO/runs/track_b/pepscorer
LOG=$OUT/run.log
VENV=$REPO/methods/external/pepscorer_venv
JOBS=${JOBS:-2}
mkdir -p "$OUT"
say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

"$VENV/bin/python" - <<'PY' > "$OUT/targets.txt"
import pandas as pd, os
m = pd.read_csv('/12TBDrive1/mega_pep_bench/runs/track_b/scoring/manifest.csv')
ind = m[(m.peptide_len >= 3) & (m.peptide_len <= 10)]
for p in ind.pdb:
    if os.path.isdir(f'/12TBDrive1/mega_pep_bench/data/external/MM_PBGBSA-CP/decoy/{p}'):
        print(p)
PY
N=$(wc -l < "$OUT/targets.txt")
say "start: $N in-domain complexes, $JOBS in parallel"

xargs -a "$OUT/targets.txt" -n1 -P "$JOBS" "$REPO/scripts/run/track_b/pepscorer_one.sh" >>"$LOG" 2>&1

say "predicted: $(ls $OUT/*/predicted.csv 2>/dev/null | wc -l)/$N"
/home/yash/anaconda3/envs/pepgym-curate/bin/python "$REPO/scripts/run/track_b/collect_pepscorer.py" >>"$LOG" 2>&1
say "DONE -> $REPO/runs/track_b/pepscorer_leaderboard.json"
