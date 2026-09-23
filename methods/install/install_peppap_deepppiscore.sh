#!/bin/bash
# PepPAP (peptide binding-affinity/screening) + DeepPpIScore (Track C). Clone + torch(CPU) env.
# Repo URLs may need updating to the authors' current locations.
set -euo pipefail
REPO=/12TBDrive1/mega_pep_bench
ENV=pepgym-trackc-dl
if ! conda env list | grep -qE "^\s*${ENV}\b|/${ENV}$"; then
  conda create -y -n "$ENV" -c conda-forge python=3.10 pytorch cpuonly numpy pandas scikit-learn rdkit
fi
clone() { [ -d "$2" ] || git clone "$1" "$2" || echo "clone failed: $1 (update URL)"; }
clone https://github.com/ChunhuaLab/PepPAP          "$REPO/data/external/PepPAP"        # T100/PDZ/CXCR4 sets
clone https://github.com/zjujdj/DeepPpIScore        "$REPO/data/external/DeepPpIScore"
for d in PepPAP DeepPpIScore; do
  [ -f "$REPO/data/external/$d/requirements.txt" ] && \
    conda run -n "$ENV" pip install -r "$REPO/data/external/$d/requirements.txt" || true
done
echo "OK: PepPAP + DeepPpIScore cloned to data/external/, env '$ENV'."
