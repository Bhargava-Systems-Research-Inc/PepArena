#!/usr/bin/env bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
EXT=$REPO/methods/external/RAPiDock
FRESH=${FRESH:-$REPO/runs/track_b/fresh}
OUT=$FRESH/rapidock
LOG=$OUT/run.log
VENV=$REPO/methods/external/rapidock_venv39
BASEPY=$HOME/anaconda3/envs/dockq/bin/python3
CEIL_GB=${CEIL_GB:-34}
mkdir -p "$OUT"
say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
stamp() { touch "$OUT/.$1.done"; }
done_p() { [ -f "$OUT/.$1.done" ]; }

( while true; do
    U=$(free -g | awk '/Mem:/{print $3}')
    if [ "${U:-0}" -ge "$CEIL_GB" ]; then
      echo "[$(date '+%F %T')] ABORT: RAM ${U}G >= ${CEIL_GB}G, killing inference" >> "$LOG"
      pkill -f "inference.py" 2>/dev/null
      sleep 30
    fi
    sleep 15
  done ) & GUARD=$!
trap 'kill $GUARD 2>/dev/null' EXIT

if ! done_p weights; then
  say "stage 1: pre-trained weights (Zenodo 14193621)"
  mkdir -p "$EXT/train_models/CGTensorProductEquivariantModel"
  curl -L --retry 3 -o "$EXT/train_models/CGTensorProductEquivariantModel/rapidock_local.pt" \
    "https://zenodo.org/records/14193621/files/rapidock_local.pt?download=1" >>"$LOG" 2>&1
  if [ -s "$EXT/train_models/CGTensorProductEquivariantModel/rapidock_local.pt" ]; then
    say "  weights: $(du -h "$EXT/train_models/CGTensorProductEquivariantModel/rapidock_local.pt" | cut -f1)"
    stamp weights
  else
    say "  FAILED to download weights"; exit 1
  fi
fi

if ! done_p env; then
  say "stage 2: python venv at $VENV"
  rm -rf "$VENV"
  "$BASEPY" -m venv "$VENV" >>"$LOG" 2>&1
  "$VENV/bin/pip" install --quiet --upgrade pip >>"$LOG" 2>&1
  "$VENV/bin/pip" install --quiet --upgrade pip >>"$LOG" 2>&1
  "$VENV/bin/pip" install --quiet torch --index-url https://download.pytorch.org/whl/cu128 >>"$LOG" 2>&1
  "$VENV/bin/pip" install --quiet torch_geometric e3nn "MDAnalysis==2.6.1" rdkit \
      "fair-esm==2.0.0" biopython pyyaml scipy pandas networkx >>"$LOG" 2>&1
  TV=$("$VENV/bin/python" -c "import torch;print(torch.__version__.split('+')[0])" 2>/dev/null)
  "$VENV/bin/pip" install --quiet torch-cluster torch-scatter torch-sparse \
      -f "https://data.pyg.org/whl/torch-${TV}+cu128.html" >>"$LOG" 2>&1 || \
      say "  note: no prebuilt pyg wheels for torch ${TV}; RAPiDock may still run without them"
  if "$VENV/bin/python" -c "import torch, e3nn, MDAnalysis; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())" >>"$LOG" 2>&1; then
    say "  env OK: $("$VENV/bin/python" -c 'import torch;print(torch.__version__, torch.cuda.is_available())')"
    stamp env
  else
    say "  ENV BUILD FAILED - see $LOG"; exit 1
  fi
fi

if ! done_p inputs; then
  say "stage 3: pockets + sequences for the expressible complexes"
  conda run -n pepgym-curate python "$REPO/scripts/run/track_b/rapidock_prep.py" >>"$LOG" 2>&1
  say "  prepared: $(wc -l < "$OUT/tasks.csv" 2>/dev/null || echo 0) rows (incl. header)"
  stamp inputs
fi

if ! done_p infer; then
  say "stage 4: RAPiDock inference (local docking, N=10 samples)"
  cd "$EXT" || exit 1
  "$VENV/bin/python" inference.py \
      --protein_peptide_csv "$OUT/tasks.csv" \
      --output_dir "$OUT/out" \
      --model_dir train_models/CGTensorProductEquivariantModel \
      --ckpt rapidock_local.pt \
      --N 10 --batch_size 4 --inference_steps 16 --actual_steps 16 \
      --no_final_step_noise --conformation_partial 1:1:1 --cpu 8 >>"$LOG" 2>&1
  rc=$?
  n=$(ls "$OUT/out" 2>/dev/null | wc -l)
  say "  inference exit=$rc; complexes with output: $n"
  if [ "$n" -gt 0 ]; then stamp infer; else say "  NO OUTPUT - not stamping, see $LOG"; exit 1; fi
fi

say "stage 5: DockQ"
conda run -n pepgym-curate --no-capture-output python "$REPO/scripts/run/track_b/score_fresh.py" \
    --tool rapidock --workers 8 >>"$LOG" 2>&1
say "DONE. leaderboard: $OUT/../rapidock_leaderboard.json"
