#!/bin/bash
cd /12TBDrive1/mega_pep_bench
PYP="$HOME/anaconda3/bin/python"
log(){ echo "[$(date '+%F %T')] $*"; }
log "RFAA chain: waiting for HelixFold3 to finish (frees GPU)..."
while pgrep -f run_helixfold3.py >/dev/null; do sleep 120; done
log "HelixFold3 done -> starting RFAA (full GPU)"
$PYP scripts/run/track_a/run_rfaa.py --workers 4 > runs/track_a/logs/rfaa_canonical.log 2>&1
log "RFAA done -> scoring"
$PYP scripts/run/track_a/score_dockq.py --model rfaa > runs/track_a/logs/score_rfaa.log 2>&1
$PYP scripts/run/track_a/build_leaderboard.py > runs/track_a/logs/leaderboard.log 2>&1
log "RFAA scored: $(tail -2 runs/track_a/logs/score_rfaa.log|head -1)"
