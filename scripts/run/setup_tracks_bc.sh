#!/bin/bash
set -uo pipefail
REPO=/12TBDrive1/mega_pep_bench
I=$REPO/methods/install

echo "== quick/reliable (conda-forge + clones) =="
bash "$I/install_rdock.sh"                 || echo "  rdock: see log"
bash "$I/install_ppi_affinity.sh"          || echo "  ppi-affinity: check repo URL"

echo ""
echo "== heavier / interactive (run individually when ready) =="
echo "  bash $I/install_adcp_adfrsuite.sh     # ADCP: ~500 MB CCSB download"
echo "  bash $I/install_cabsdock.sh           # CABS-dock: needs MODELLER key"
echo "  bash $I/install_peppap_deepppiscore.sh# PepPAP + DeepPpIScore: verify repo URLs"
echo "  bash $I/install_rosetta.sh            # Rosetta FlexPepDock: academic license + build"

echo ""
echo "== datasets (deterministic, no GPU) =="
echo "  Track B: python $REPO/scripts/run/track_b/prep_scoring_set.py; prep_fresh_docking.py; prep_complex_pdbs.py"
echo "  Track C: python $REPO/scripts/run/track_c/prep_affinity_set.py; fetch_receptor_seqs.py; materialize_affinity_structures.py"
echo "DONE."
