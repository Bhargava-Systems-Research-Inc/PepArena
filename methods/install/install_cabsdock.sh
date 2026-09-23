#!/bin/bash
# CABS-dock: coarse-grained flexible peptide docking with a cyclic protocol (Track B).
# Python package (CABS) installable via pip into a dedicated env. Needs DSSP + modeller for
# all-atom reconstruction (modeller is license-keyed but free for academics).
set -euo pipefail
ENV=cabsdock
if ! conda env list | grep -qE "^\s*${ENV}\b|/${ENV}$"; then
  conda create -y -n "$ENV" -c salilab -c conda-forge python=3.9 dssp modeller
fi
conda run -n "$ENV" pip install pycabs || conda run -n "$ENV" pip install CABS
echo "OK: CABS-dock in env '$ENV'."
echo "NOTE: set your MODELLER license key (KEY_MODELLER) in \$CONDA_PREFIX/lib/modeller-*/modlib/modeller/config.py"
echo "Test: conda run -n $ENV CABSdock -h"
