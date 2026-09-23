# Track C — binding affinity + screening

Eval set: **2,045 measurements** (629 complexes, 533 receptors, 470 peptides; log-Kd −12.9..−2.5;
123 cyclic, 66 ncAA), from `affinity_ppb` (1,995) + `affinity_mmpbgbsa_cp` (50 cyclic Kd).
Two leakage-aware split families: `split_random` and `split_receptor` (novel targets).

Metrics: Spearman, Pearson, RMSE/MAE (log-Kd). Screening: NDCG@k, EF@1% (add per target set).

## Prepare (ready)
```
conda run -n pepgym-curate python prep_affinity_set.py           # manifest + splits (done)
conda run -n pepgym-curate python fetch_receptor_seqs.py         # receptor seqs (RCSB API, light)
conda run -n pepgym-curate python materialize_affinity_structures.py   # native PDBs (for scorers)
```

## Methods (cheapest first)
| Method | Status | Compute | How |
|---|---|---|---|
| **ESM-2 + regressor** (this recipe) | **ready, run** | CPU, cheap | `run_affinity_regressor.py` — reuses Track-D embedding cache |
| **Boltz-2 affinity head** | installed (`boltz`) | GPU, moderate | `run_boltz2_affinity.py` (verify YAML schema for boltz 2.2.1) |
| HERMES | installed (`hermes`) | CPU/GPU | scores native complexes; wire per HERMES CLI |
| PPI-Affinity | install | CPU | `methods/install/install_ppi_affinity.sh` |
| PepPAP / DeepPpIScore | install | CPU | clone + env |
| MM-PB(GB)SA | installed (`gmxmmpbsa`) | **expensive (MD)** | single-snapshot MM-GBSA on the 50 cyclic-Kd set only |

## First result (already produced — `leaderboard_regressor.csv`)
ESM-2 (t12_35M) + GBM, peptide-only features:
- `split_random`: Spearman **0.727** — but this is leaky (peptides of the same target on both sides)
- `split_receptor` (novel targets): Spearman **0.336** — the honest number

The random→receptor collapse (0.73 → 0.34) reproduces the same leakage-inflation story as
Tracks A and D. Adding receptor features (`fetch_receptor_seqs.py`) should lift the receptor split.

## Notes
- **No overlap with Track A predictions** by construction: affinity structures are pre-2021,
  the Track-A cofolding tier is post-2023-06-01. Track C runs on its own ~629-complex structure set.
- `split_receptor` groups by receptor PDB (also keeps SKEMPI mutational variants together).
  Upgrade path: mmseqs2 receptor-cluster grouping (as Track A) once receptor seqs are materialized.
