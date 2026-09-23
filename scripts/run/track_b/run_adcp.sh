#!/bin/bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
FRESH=${FRESH:-$REPO/runs/track_b/fresh}
OUT=$FRESH/adcp
mkdir -p "$OUT"
ADFR=${ADFR_HOME:-$HOME/ADFRsuite/install/bin}   # install.sh -d $HOME/ADFRsuite puts binaries here
export PATH="$ADFR:$PATH"

WORKLIST=${WORKLIST:-$FRESH/worklist.csv}
NBRUNS=${NBRUNS:-50}
NSTEPS=${NSTEPS:-2500000}
MAXCORES=${MAXCORES:-24}
PAD=${PAD:-8.0}

box_for() {  # native peptide bbox + PAD, floor 16 A per side
  python3 -c "
import sys
xs=[];ys=[];zs=[]
for l in open(sys.argv[1]):
    if l.startswith(('ATOM','HETATM')):
        xs.append(float(l[30:38])); ys.append(float(l[38:46])); zs.append(float(l[46:54]))
P=float(sys.argv[2])
c=[(min(v)+max(v))/2 for v in (xs,ys,zs)]
s=[max(max(v)-min(v)+2*P, 16.0) for v in (xs,ys,zs)]
print(' '.join(f'{v:.2f}' for v in c+s))" "$1" "$PAD"
}

conda run -n pepgym-curate python -c "
import csv, sys
w = csv.writer(sys.stdout, delimiter='\t')
for r in csv.DictReader(open('$WORKLIST')):
    w.writerow([r['complex_id'], r['peptide_seq'], r['peptide_len'],
                r['is_cyclic'], r['cyclization_type']])
" | while IFS=$'\t' read -r cid pepseq peplen iscyclic cyctype; do
  indir="$FRESH/inputs/$cid"; rundir="$OUT/$cid"
  [ -f "$rundir/${cid}_out_ranked_1.pdb" ] && { echo "skip $cid (done)"; continue; }
  [ -f "$indir/receptor.pdb" ] || { echo "skip $cid (no inputs)"; continue; }

  if grep -q '^HET ' "$indir/peptide.pdb" 2>/dev/null; then
    echo "SKIP $cid (non-canonical residues — ADCP has no rotamers for them)"; continue
  fi
  if [ "$peplen" -lt 5 ]; then
    echo "SKIP $cid (peptide_len=$peplen — below ADCP's supported range)"; continue
  fi

  mkdir -p "$rundir"; cd "$rundir"
  prepare_receptor -r "$indir/receptor.pdb" -o receptor.pdbqt -A checkhydrogens 2>/dev/null
  BOX=$(box_for "$indir/peptide.pdb")
  agfr -r receptor.pdbqt -b user $BOX -o target >agfr.log 2>&1
  [ -f target.trg ] || { echo "FAIL $cid (agfr produced no target; see $rundir/agfr.log)"; cd "$REPO"; continue; }

  CYC=""; [ "$iscyclic" = "True" ] && case "$cyctype" in
    head_to_tail|bicyclic) CYC="-cyc";; disulfide) CYC="-cys";; esac
  echo "=== ADCP $cid len=$peplen ($pepseq) ${CYC:-linear} box=[$BOX] $(date '+%F %T') ==="
  adcp -t target.trg -s "$pepseq" -N "$NBRUNS" -n "$NSTEPS" -c "$MAXCORES" \
       -o "${cid}_out" $CYC -O || echo "FAIL $cid"
  rm -f rigidReceptor.*.map target.trg          # ~10 MB/complex of regenerable grids
  cd "$REPO"
done
echo "DONE. Score with: conda run -n pepgym-curate python score_fresh.py --tool adcp"
