#!/bin/bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
FRESH=${FRESH:-$REPO/runs/track_b/fresh}
OUT=$FRESH/haddock3
mkdir -p "$OUT"
HADDOCK="conda run -n haddock3 haddock3"
WORKLIST=${WORKLIST:-$FRESH/worklist.csv}
SAMPLING=${SAMPLING:-1000}
SELECT=${SELECT:-200}

pick_ncores() {
  if [ "$(pgrep -xc gmx_mpi 2>/dev/null || echo 0)" -gt 0 ]; then echo "${NCORES_BUSY:-6}"
  else echo "${NCORES_FREE:-24}"; fi
}

tail -n +2 "$WORKLIST" | while IFS=, read -r cid rest; do
  indir="$FRESH/inputs/$cid"; rundir="$OUT/$cid"
  [ -d "$rundir/run" ] && { echo "skip $cid (done)"; continue; }
  [ -f "$indir/receptor.pdb" ] || { echo "skip $cid (no inputs)"; continue; }

  if grep -q '^HET ' "$indir/peptide.pdb" 2>/dev/null; then
    echo "SKIP $cid (non-canonical residues — topoaa would silently truncate the peptide)"; continue
  fi
  nch=$(awk '/^ATOM/{print substr($0,22,1)}' "$indir/receptor.pdb" | sort -u | wc -l)
  if [ "$nch" -gt 1 ]; then
    echo "SKIP $cid ($nch receptor chains — needs segid prep before CNS will accept it)"; continue
  fi

  NC=$(pick_ncores)
  mkdir -p "$rundir"
  cat > "$rundir/dock.cfg" <<CFG
run_dir = "$rundir/run"
molecules = ["$indir/receptor.pdb", "$indir/peptide.pdb"]
ncores = $NC
mode = "local"
[topoaa]
[rigidbody]
sampling = $SAMPLING
cmrest = true
[seletop]
select = $SELECT
[flexref]
[emref]
[clustfcc]
CFG
  echo "=== docking $cid (recCA chains=$nch, ncores=$NC) $(date '+%F %T') ==="
  $HADDOCK "$rundir/dock.cfg" || echo "FAIL $cid"
done
echo "DONE. Score with: conda run -n pepgym-curate python score_fresh.py --tool haddock3"
