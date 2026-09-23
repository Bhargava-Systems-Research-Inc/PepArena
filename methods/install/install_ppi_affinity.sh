#!/bin/bash
# PPI-Affinity: SVM protein/peptide binding-affinity predictor (Track C).
# NOTE: PPI-Affinity (Romero-Molina & Ruiz-Blanco, J. Proteome Res. 2022) is distributed as a
# WEB TOOL (https://ppi-affinity.mutabind.org / the group's server), not a maintained public repo.
# Two ways to use it in the harness:
#   (a) script submissions to the web server for the affinity set (flag as web-derived), or
#   (b) use the authors' open sequence-based VS code, SVSBI (github.com/WeilabMSU/SVSBI),
#       which implements the same ESM/transformer + ML affinity recipe locally.
# This sets up an env for option (b) and clones SVSBI.
set -euo pipefail
REPO=/12TBDrive1/mega_pep_bench
DEST=$REPO/data/external/SVSBI
ENV=ppi-affinity
if [ ! -d "$DEST" ]; then
  git clone https://github.com/WeilabMSU/SVSBI "$DEST" || echo "clone failed — verify SVSBI URL"
fi
if ! conda env list | grep -qE "^\s*${ENV}\b|/${ENV}$"; then
  conda create -y -n "$ENV" -c conda-forge python=3.9 scikit-learn pandas numpy biopython
fi
[ -f "$DEST/requirements.txt" ] && conda run -n "$ENV" pip install -r "$DEST/requirements.txt" || true
echo "OK: env '$ENV' ready; SVSBI at $DEST. For the web tool, script submissions and tag as web-derived."
