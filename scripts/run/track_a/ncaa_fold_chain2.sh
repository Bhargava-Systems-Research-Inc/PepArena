#!/bin/bash
set -u
cd /12TBDrive1/mega_pep_bench
PY=/home/yash/anaconda3/envs/protenix/bin/python
stamp(){ echo "[$(date +%H:%M:%S)] $*" >> runs/track_a/logs/ncaa_chain.log; }

stamp "af3 (fixed: bare ptmType + empty-MSA modified chains) starting"
$PY scripts/run/track_a/ncaa_fold.py --model af3 >> runs/track_a/logs/ncaa_fold_af3.log 2>&1
stamp "af3 done: $(find runs/track_a/predictions_ncaa/af3 -name '*_model.cif' 2>/dev/null | wc -l) cifs"

stamp "highfold3 (AF3+CycPOEM, 91 cyclics) starting"
$PY scripts/run/track_a/ncaa_fold.py --model highfold3 >> runs/track_a/logs/ncaa_fold_hf3.log 2>&1
stamp "highfold3 done: $(find runs/track_a/predictions_ncaa/highfold3 -name '*_model.cif' 2>/dev/null | wc -l) cifs"

kill -CONT 2722655 2>/dev/null
stamp "ALL DONE, resumed MD (2722655 state $(ps -o stat= -p 2722655))"
