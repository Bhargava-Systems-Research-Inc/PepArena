#!/bin/bash
set -u
cd /12TBDrive1/mega_pep_bench
PYP="$HOME/anaconda3/bin/python"                  # pandas env for scoring + orchestration
CHAIPY="$HOME/anaconda3/envs/chai/bin/python"
log() { echo "[$(date '+%F %T')] $*"; }

log "waiting for Chai-1 to finish..."
while pgrep -f run_chai.py >/dev/null; do sleep 120; done
log "Chai-1 finished -> scoring"
$PYP scripts/run/track_a/score_dockq.py --model chai1 > runs/track_a/logs/score_chai.log 2>&1
log "Chai-1: $(tail -2 runs/track_a/logs/score_chai.log | head -1)"
$PYP scripts/run/track_a/build_leaderboard.py > runs/track_a/logs/leaderboard.log 2>&1

log "HighFold3 (canonical, WITH MSA)"
$PYP scripts/run/track_a/run_highfold3.py > runs/track_a/logs/highfold3_canonical.log 2>&1
log "HighFold3 finished -> scoring"
$PYP scripts/run/track_a/score_dockq.py --model highfold3 > runs/track_a/logs/score_highfold3.log 2>&1
log "HighFold3: $(tail -2 runs/track_a/logs/score_highfold3.log | head -1)"
$PYP scripts/run/track_a/build_leaderboard.py > runs/track_a/logs/leaderboard.log 2>&1
log "overnight chain complete"
