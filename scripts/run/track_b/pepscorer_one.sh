#!/usr/bin/env bash
set -uo pipefail
PDB="$1"
REPO=/12TBDrive1/mega_pep_bench
EXT=$REPO/methods/external
PS=$EXT/PepScorer/PepScorerRMSD/PepScorerRMSD
VENV=$EXT/pepscorer_venv
DECOY=$REPO/data/external/MM_PBGBSA-CP/decoy
DSET=$REPO/data/external/MM_PBGBSA-CP/dataset2
OUT=$REPO/runs/track_b/pepscorer
W="$OUT/$PDB"
export VEGADIR=$EXT/Vega
export PATH=$EXT/Vega/Bin/Linux_x64:$PATH
export LD_LIBRARY_PATH=$EXT/Vega/Bin/Linux_x64:${LD_LIBRARY_PATH:-}

cap() { nice -n 19 systemd-run --user --scope -q -p MemoryMax=4G "$@"; }

SCORES="CHARMM22 ELECT ELECTDD MLPINS MLPINS2 MLPINS3 MLPINSF CHEMPLP RPSCORE"

[ -f "$W/predicted.csv" ] && exit 0
mkdir -p "$W/mol2"
REC="$DSET/$PDB/${PDB}_protein.pdb"; [ -f "$REC" ] || exit 1

RECP="$W/receptor_prep.iff"
if [ ! -s "$RECP" ]; then
  cap vega "$REC" -l PROT -c GASTEIGER -p CHARMM22_PROT -f IFF -o "$RECP" >>"$W/prep.log" 2>&1
fi
[ -s "$RECP" ] || exit 1

if [ ! -f "$W/.prepared_vega" ]; then
  rm -f "$W/mol2"/*.mol2 "$W/all.mol2"
  for i in $(seq 1 100); do
    D="$DECOY/$PDB/$i.pdb"; [ -f "$D" ] || continue
    N="${PDB}_$(printf '%03d' $i)"
    cap vega "$D" -l PROT -c GASTEIGER -f MOL2 -o "$W/mol2/$N.mol2" >>"$W/prep.log" 2>&1
  done
  "$VENV/bin/python" "$REPO/scripts/run/track_b/repair_vega_types.py" "$W/mol2" >>"$W/prep.log" 2>&1
  : > "$W/all.mol2"
  for f in "$W/mol2"/*.mol2; do cat "$f" >> "$W/all.mol2"; done
  [ -s "$W/all.mol2" ] && touch "$W/.prepared_vega"
fi
[ -s "$W/all.mol2" ] || exit 1

has_rows() { [ -f "$1" ] && [ "$(wc -l < "$1")" -gt 1 ]; }
if ! has_rows "$W/rescore_raw.csv"; then
  rm -f "$W/rescore_raw.csv"
  cap rescore+ -p "$SCORES" -d "$W/all.mol2" -DDB -r "$RECP" -t 1 \
      "$W/rescore_raw.csv" >>"$W/rescore.log" 2>&1
fi
has_rows "$W/rescore_raw.csv" || exit 1

"$VENV/bin/python" - "$W" <<'PY'
import sys, os, pandas as pd
W = sys.argv[1]
ref = list(pd.read_excel('/12TBDrive1/mega_pep_bench/methods/external/PepScorer/PepScorerRMSD/'
                         'PepScorerRMSD/test/rescore.xlsx').columns)
d = pd.read_csv(f'{W}/rescore_raw.csv', sep=';')
d = d.rename(columns={'CHARMM22': 'CHARMM'})
ff = ['CHARMM', 'Elect', 'ElectDD']
failed = (d[ff] == 0).all(axis=1)
d.loc[failed, ff] = float('nan')
if failed.any():
    print(f'{int(failed.sum())}/{len(d)} poses: force-field block failed, set to NaN')
for c in [c for c in ref if c.startswith('XS_') or c.startswith('APBS_')]:
    d[c] = float('nan')
missing = [c for c in ref if c not in d.columns]
if missing:
    raise SystemExit(f'missing Rescore+ columns: {missing}')
order = [f.rsplit('.', 1)[0] for f in os.listdir(f'{W}/mol2') if f.endswith('.mol2')]
d = d.set_index('Name').reindex(order).reset_index()
if d[ref[1]].isna().all() and d['CHARMM'].isna().all():
    raise SystemExit('row realignment produced no matches -- name mismatch')
if len(d) != len(order):
    raise SystemExit(f'row count {len(d)} != {len(order)} poses')
d[ref].to_csv(f'{W}/rescore.csv', index=False)
print(f'{len(d)} poses, {len(ref)} columns, aligned to listdir order')
PY
[ -s "$W/rescore.csv" ] || exit 1

RUND="$W/psrun"
mkdir -p "$RUND"
ln -sfn "$PS/objects" "$RUND/objects"
cd "$RUND" || exit 1
cap "$VENV/bin/python" "$REPO/scripts/run/track_b/pepscorer_predict.py" -l "$W/mol2/" -r "$W/rescore.csv" >>"$W/predict.log" 2>&1 \
  && mv "$RUND/PepScorerRMSD_output.csv" "$W/predicted.csv"
