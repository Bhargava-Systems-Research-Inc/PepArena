#!/bin/bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
FRESH=${FRESH:-$REPO/runs/track_b/fresh}
OUT=$FRESH/unidock
mkdir -p "$OUT"
export PATH="${ADFR_HOME:-$HOME/ADFRsuite/install/bin}:$PATH"
UNIDOCK="conda run -n unidock_env --no-capture-output unidock"
PAD=${PAD:-8.0}
NMODES=${NMODES:-9}          # Vina default
EXHAUST=${EXHAUST:-8}        # Vina default
MAX_TORSDOF=${MAX_TORSDOF:-15}
TIMEOUT=${TIMEOUT:-1200}

if [ -n "${ACCEPTED:-}" ]; then
  cp "$ACCEPTED" "$OUT/_accepted_$(basename "$ACCEPTED")"
  ACC="$OUT/_accepted_$(basename "$ACCEPTED")"
else
  ACC="$OUT/_accepted.txt"
  "$HOME/anaconda3/envs/pepgym-curate/bin/python" -c "
import pandas as pd
d = pd.read_csv('$FRESH/coverage.csv')
d = d.loc[(d.unidock == 'ok') & (d.torsdof <= $MAX_TORSDOF)]
print('\n'.join(d.complex_id))
" > "$ACC"
fi
echo "unidock set: $(wc -l < "$ACC") complexes"

while read -r cid; do
  [ -n "$cid" ] || continue
  indir="$FRESH/inputs/$cid"; rundir="$OUT/$cid"
  [ -f "$rundir/${cid}_out_ranked_1.pdb" ] && { echo "skip $cid (done)"; continue; }
  mkdir -p "$rundir"; cd "$rundir"

  cp "$indir/receptor.pdb" rec_in.pdb; cp "$indir/peptide.pdb" lig_in.pdb
  prepare_receptor -r rec_in.pdb -o receptor.pdbqt -A checkhydrogens >prep.log 2>&1
  prepare_ligand   -l lig_in.pdb -o ligand.pdbqt                     >>prep.log 2>&1
  if [ ! -s receptor.pdbqt ] || [ ! -s ligand.pdbqt ]; then
    echo "FAIL $cid (prep produced no pdbqt)"; cd "$REPO"; continue
  fi
  A_IN=$(grep -c '^ATOM\|^HETATM' lig_in.pdb); A_OUT=$(grep -c '^ATOM\|^HETATM' ligand.pdbqt)
  if [ "$A_IN" -ne "$A_OUT" ]; then
    echo "FAIL $cid (ligand prep dropped $((A_IN-A_OUT))/$A_IN atoms)"; cd "$REPO"; continue
  fi

  read -r CX CY CZ SX SY SZ <<<"$(python3 -c "
xs=[];ys=[];zs=[]
for l in open('lig_in.pdb'):
    if l.startswith(('ATOM','HETATM')):
        xs.append(float(l[30:38])); ys.append(float(l[38:46])); zs.append(float(l[46:54]))
P=$PAD
c=[(min(v)+max(v))/2 for v in (xs,ys,zs)]
s=[max(max(v)-min(v)+2*P, 16.0) for v in (xs,ys,zs)]
print(' '.join(f'{v:.2f}' for v in c+s))")"

  echo "=== Uni-Dock $cid box=($CX,$CY,$CZ) size=($SX,$SY,$SZ) $(date '+%F %T') ==="
  timeout "$TIMEOUT" $UNIDOCK --receptor receptor.pdbqt --gpu_batch ligand.pdbqt \
    --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
    --size_x "$SX" --size_y "$SY" --size_z "$SZ" \
    --num_modes "$NMODES" --exhaustiveness "$EXHAUST" --scoring vina \
    --dir . >dock.log 2>&1
  rc=$?; [ $rc -eq 124 ] && echo "TIMEOUT $cid (>${TIMEOUT}s)"
  [ $rc -ne 0 ] && [ $rc -ne 124 ] && echo "FAIL $cid (unidock rc=$rc)"

  python3 - "$cid" <<'PYEOF'
import sys, pathlib
cid = sys.argv[1]
src = next((p for p in (pathlib.Path("ligand_out.pdbqt"), pathlib.Path("ligand.pdbqt_out.pdbqt"))
            if p.exists()), None)
if src is None:
    cands = [p for p in pathlib.Path(".").glob("*out*.pdbqt")]
    src = cands[0] if cands else None
if src is None:
    print("no unidock output pdbqt"); raise SystemExit
n = 0
buf = []
for ln in src.read_text().splitlines():
    if ln.startswith("MODEL"):
        buf = []
    elif ln.startswith("ENDMDL"):
        n += 1
        pathlib.Path(f"{cid}_out_ranked_{n}.pdb").write_text("\n".join(buf) + "\nEND\n")
    elif ln.startswith(("ATOM", "HETATM")):
        buf.append(ln[:66])          # drop the pdbqt charge/type columns
print(f"wrote {n} ranked poses")
PYEOF
  rm -f receptor.pdbqt ligand.pdbqt rec_in.pdb lig_in.pdb
  cd "$REPO"
done < "$ACC"
echo "DONE. Score with: conda run -n pepgym-curate python score_fresh.py --tool unidock"
