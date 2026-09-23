#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."                 # repo root
PY="${PEPGYM_PY:-$HOME/anaconda3/envs/pepgym-curate/bin/python}"
C="scripts/curate"

echo "== external source repos =="
mkdir -p data/external
[ -d data/external/PepPCBench ] || git clone --depth 1 https://github.com/zhaisilong/PepPCBench.git data/external/PepPCBench
[ -d data/external/TPepPro ]    || git clone --depth 1 https://github.com/wanglabhku/TPepPro.git    data/external/TPepPro

echo "== affinity source =="
mkdir -p data/raw/affinity
[ -f data/raw/affinity/PPB-Affinity.xlsx ] || \
  curl -sL "https://zenodo.org/api/records/14271435/files/PPB-Affinity.xlsx/content" -o data/raw/affinity/PPB-Affinity.xlsx

echo "== PepBench / LNR structure sources (Zenodo 13373108) =="
mkdir -p data/raw/pepbench
[ -d data/raw/pepbench/LNR ] || { \
  [ -f data/raw/pepbench/LNR.tar.gz ] || curl -sL "https://zenodo.org/api/records/13373108/files/LNR.tar.gz/content" -o data/raw/pepbench/LNR.tar.gz; \
  tar xzf data/raw/pepbench/LNR.tar.gz -C data/raw/pepbench; }
[ -d data/raw/pepbench/train_valid ] || { \
  [ -f data/raw/pepbench/train_valid.tar.gz ] || curl -sL "https://zenodo.org/api/records/13373108/files/train_valid.tar.gz/content" -o data/raw/pepbench/train_valid.tar.gz; \
  tar xzf data/raw/pepbench/train_valid.tar.gz -C data/raw/pepbench; }

echo "== ingest structure (Tracks A/B) =="
"$PY" $C/ingest_peppcbench.py
"$PY" $C/ingest_cyclic_complexes.py
"$PY" $C/ingest_cyclic_monomers.py
"$PY" $C/ingest_pdb_ncaa.py        # non-canonical peptide-protein complexes (PDB mine)
"$PY" $C/ingest_pdb_list.py --tsv data/raw/pepbench/LNR/test.txt --source lnr --out data/curated/structure_lnr.parquet
"$PY" $C/ingest_pepbench.py        # large PepBench/PepGLAD structure pool (~6105)

echo "== ingest affinity (Track C) =="
"$PY" $C/ingest_ppb_affinity.py
"$PY" $C/enrich_affinity.py         # structure-based cyclization + ncAA labels
[ -d data/external/MM_PBGBSA-CP ] || git clone --depth 1 https://github.com/huifengzhao/MM_PBGBSA-CP.git data/external/MM_PBGBSA-CP
[ -d data/external/MM_PBGBSA-CP/dataset1 ] || tar xjf data/external/MM_PBGBSA-CP/dataset1.tar.bz2 -C data/external/MM_PBGBSA-CP
"$PY" $C/ingest_mmpbgbsa_cp.py

echo "== ingest interaction (Track D) =="
"$PY" $C/ingest_tpeppro.py

echo "== combine + splits + native structures =="
"$PY" $C/combine_manifests.py
"$PY" $C/make_splits.py
"$PY" $C/materialize_structures.py

echo "== enrich features (secondary structure + interface) + target class =="
"$PY" $C/enrich_features.py
"$PY" $C/annotate_target_class.py
"$PY" $C/combine_manifests.py        # re-merge now that features + classes exist

echo "== QC (fails build on error) =="
"$PY" $C/validate.py

echo "== DONE =="
