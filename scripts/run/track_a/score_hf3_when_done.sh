#!/bin/bash
cd /12TBDrive1/mega_pep_bench
while pgrep -f run_helixfold3.py >/dev/null; do sleep 60; done
~/anaconda3/bin/python scripts/run/track_a/score_dockq.py --model helixfold3 > runs/track_a/logs/score_helixfold3.log 2>&1
~/anaconda3/bin/python scripts/run/track_a/build_leaderboard.py > runs/track_a/logs/leaderboard.log 2>&1
echo "[$(date '+%F %T')] HelixFold3 scored: $(tail -2 runs/track_a/logs/score_helixfold3.log|head -1)" >> runs/track_a/logs/hf3_autoscore.log
