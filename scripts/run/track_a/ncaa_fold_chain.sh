#!/bin/bash
set -u
cd /12TBDrive1/mega_pep_bench
PY=/home/yash/anaconda3/envs/protenix/bin/python
CHAIN=runs/track_a/logs/ncaa_chain.log
MDRUN=2722655

stamp(){ echo "[$(date +%H:%M:%S)] $*" >> "$CHAIN"; }

while pgrep -f 'ncaa_fold.py --model protenix' >/dev/null; do sleep 30; done
stamp "protenix done: $(find runs/track_a/predictions_ncaa/protenix -path '*seed_*/predictions/*.cif' | wc -l) cifs"

stamp "starting boltz2"
$PY scripts/run/track_a/ncaa_fold.py --model boltz2 >> runs/track_a/logs/ncaa_fold_boltz.log 2>&1
stamp "boltz2 done: $(find runs/track_a/predictions_ncaa/boltz2 -name '*.cif' 2>/dev/null | wc -l) cifs"

stamp "starting af3"
$PY scripts/run/track_a/ncaa_fold.py --model af3 >> runs/track_a/logs/ncaa_fold_af3.log 2>&1
stamp "af3 done: $(find runs/track_a/predictions_ncaa/af3 -name '*_model.cif' 2>/dev/null | wc -l) cifs"

stamp "starting highfold3"
$PY scripts/run/track_a/ncaa_fold.py --model highfold3 >> runs/track_a/logs/ncaa_fold_hf3.log 2>&1
stamp "highfold3 done: $(find runs/track_a/predictions_ncaa/highfold3 -name '*_model.cif' 2>/dev/null | wc -l) cifs"

kill -CONT "$MDRUN" 2>/dev/null
stamp "ALL DONE, resumed MD ($MDRUN state $(ps -o stat= -p $MDRUN))"
