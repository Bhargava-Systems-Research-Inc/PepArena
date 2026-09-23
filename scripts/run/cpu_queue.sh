#!/usr/bin/env bash
set -uo pipefail
source ~/anaconda3/etc/profile.d/conda.sh
REPO=/12TBDrive1/mega_pep_bench
NWORK=${NWORK:-6}
say() { echo "[$(date '+%F %T')] $*"; }

wait_for_ram() {
  for _ in $(seq 1 60); do
    a=$(free -g | awk '/^Mem:/{print $7}')
    [ "${a:-0}" -ge 8 ] && return 0
    say "only ${a}G RAM available, waiting 60s"; sleep 60
  done
  say "RAM never recovered; skipping stage"; return 1
}

say "=== queue start (workers=$NWORK, gmx alive=$(pgrep -xc gmx_mpi)) ==="

say "--- stage 1/4: receptor sequences (network-bound) ---"
conda activate pepgym-curate
( cd $REPO/scripts/run/track_c && timeout 3600 python fetch_receptor_seqs.py ) \
  && say "stage 1 ok" || say "stage 1 FAILED (non-fatal; only blocks the GPU affinity run)"

say "--- stage 2/4: Track B scoring power, DockQ over all decoys ---"
if wait_for_ram; then
  cd $REPO/scripts/run/track_b
  tail -n +2 $REPO/runs/track_b/scoring/manifest.csv | cut -d, -f1 \
    | xargs -P $NWORK -I{} python -c "import score_docking as s; s.score_complex('{}')" 
  python score_docking.py > $REPO/runs/track_b/scoring/scoring_power.log 2>&1 \
    && say "stage 2 ok -> runs/track_b/scoring/leaderboard.json" || say "stage 2 summary FAILED"
fi

say "--- stage 3/4: MM-GBSA on the 50 cyclic-Kd complexes ---"
if wait_for_ram; then
  mkdir -p $REPO/runs/track_c/mmgbsa
  ls $REPO/data/external/MM_PBGBSA-CP/dataset1 \
    | xargs -P $NWORK -I{} bash $REPO/scripts/run/track_c/mmgbsa_one.sh {}
  say "stage 3 done -> $(wc -l < $REPO/runs/track_c/mmgbsa/results.tsv 2>/dev/null || echo 0) rows"
fi

say "--- stage 4/4: HADDOCK3 fresh docking (long pole) ---"
if wait_for_ram; then
  export NCORES_BUSY=$NWORK NCORES_FREE=24
  bash $REPO/scripts/run/track_b/run_haddock3.sh
  say "stage 4 done"
fi
say "=== queue complete ==="
