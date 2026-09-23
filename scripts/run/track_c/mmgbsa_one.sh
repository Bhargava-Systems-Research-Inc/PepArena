#!/usr/bin/env bash
set -uo pipefail
source ~/anaconda3/etc/profile.d/conda.sh; conda activate gmxmmpbsa
PDB=$1
REPO=/12TBDrive1/mega_pep_bench
SRC=$REPO/data/external/MM_PBGBSA-CP/dataset1/$PDB
D=$REPO/runs/track_c/mmgbsa/$PDB
RES=$REPO/runs/track_c/mmgbsa/results.tsv
[ -f "$D/out.dat" ] && { echo "skip $PDB (done)"; exit 0; }
rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1

python $REPO/scripts/run/track_c/_prep_pdb.py "$SRC/${PDB}_protein.pdb" rec.pdb > prep.log 2>&1 || exit 1
python $REPO/scripts/run/track_c/_prep_pdb.py "$SRC/${PDB}_CP.pdb" lig.pdb >> prep.log 2>&1 || exit 1
grep -v '^TER' rec.pdb > com.pdb; echo TER >> com.pdb; cat lig.pdb >> com.pdb; echo END >> com.pdb

for x in com rec lig; do
  printf 'source leaprc.protein.ff14SB\nm = loadpdb %s.pdb\nsaveamberparm m %s.prmtop %s.inpcrd\nquit\n' $x $x $x > ${x}.leap
  tleap -f ${x}.leap > ${x}.leap.log 2>&1
  [ -s ${x}.prmtop ] || { echo -e "$PDB\tTLEAP_FAIL_$x\tNA\tNA" >> $RES; exit 1; }
done

printf 'minimization\n &cntrl\n  imin=1, maxcyc=1000, ncyc=500, igb=2, cut=999.0, ntb=0, ntpr=200,\n /\n' > min.in
S=$(date +%s)
timeout 3600 sander -O -i min.in -o min.out -p com.prmtop -c com.inpcrd -r min.rst > /dev/null 2>&1
MINRC=$?; M=$(date +%s)
[ $MINRC -ne 0 ] && { echo -e "$PDB\tMIN_FAIL\tNA\tNA" >> $RES; exit 1; }

printf 'parm com.prmtop\ntrajin min.rst\ntrajout snap.mdcrd\nrun\nquit\n' > traj.cpptraj
cpptraj -i traj.cpptraj > cpptraj.log 2>&1
printf 'Single-snapshot MM-GBSA on the minimized complex\n&general\n startframe=1, endframe=1, interval=1,\n/\n&gb\n igb=2, saltcon=0.150,\n/\n' > mmpbsa.in
timeout 1800 MMPBSA.py -O -i mmpbsa.in -o out.dat -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y snap.mdcrd > mmpbsa.log 2>&1
RC=$?; E=$(date +%s)
DG=$(awk '/DELTA TOTAL/{print $3; exit}' out.dat 2>/dev/null)
echo -e "$PDB\trc=$RC\t${DG:-NA}\tmin_s=$((M-S))\ttotal_s=$((E-S))" >> $RES
echo "$PDB dG=${DG:-NA} min=$((M-S))s total=$((E-S))s"
