#!/bin/bash
set -uo pipefail
R=/12TBDrive1/mega_pep_bench
export FRESH=$R/runs/track_b/fresh
export PATH=$HOME/anaconda3/bin:$PATH
PY="conda run -n pepgym-curate python"
cd "$R/scripts/run/track_b" || exit 1

step () { echo "=== $1 $(date +%H:%M)"; shift; "$@"; echo "    rc=$?"; }

for t in haddock3 adcp unidock rapidock; do
  step "score $t" $PY score_fresh.py --tool "$t" --workers 3
done

for f in /tmp/adcp_shard*.log; do [ -e "$f" ] && cp "$f" "$FRESH/adcp_$(basename "$f" .log).log"; done
for f in /tmp/unidock_cap15_s*.log /tmp/unidock_band*.log; do
  [ -e "$f" ] && cp "$f" "$FRESH/unidock_$(basename "$f" .log).log"; done
step runtime $PY extract_runtime.py

step "unidock sensitivity" $PY unidock_sensitivity.py

step "rapidock coverage" ~/anaconda3/bin/python rapidock_coverage.py
step "merge rapidock into coverage" ~/anaconda3/bin/python merge_rapidock_coverage.py
step "concurrency effect" ~/anaconda3/bin/python concurrency_effect.py
step summary ~/anaconda3/bin/python tb_summary.py
step gate ~/anaconda3/bin/python worklist_gate.py

step "numbers the manuscript quotes" ~/anaconda3/bin/python report_numbers.py

echo "TRACKB_FINISH_DONE $(date +%H:%M)"
