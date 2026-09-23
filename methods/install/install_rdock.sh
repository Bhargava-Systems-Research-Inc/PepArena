#!/bin/bash
# rDock: open-source docking (Track B). Installable straight from conda-forge. Fast, reliable.
set -euo pipefail
ENV=rdock
if conda env list | grep -qE "^\s*${ENV}\b|/${ENV}$"; then
  echo "env ${ENV} already exists"; exit 0
fi
conda create -y -n "$ENV" -c bioconda -c conda-forge rdock
echo "OK: rDock in env '$ENV'. Test: conda run -n $ENV rbdock -h"
