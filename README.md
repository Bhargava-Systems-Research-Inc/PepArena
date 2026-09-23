# PepArena

A leakage-controlled, cyclization-aware benchmark for peptide–protein prediction: **6,933
structural complexes, 356 curated affinity measurements, 19,187 interaction pairs and 80 cyclic
monomers** on one manifest with a shared, geometry-aware cyclization spine.

Companion to *PepArena: A Leakage-Aware Benchmark for Peptide–Protein Structure, Affinity, and
Interaction Prediction* (Bhargava).

## What is here

This repository is the **inputs and the code**: the curated manifests, the leakage-aware splits,
the per-method runners and the evaluation harness. Results, predicted structures and the
per-complex score rows are in the Zenodo deposit, which is far too large for git.

| Path | Contents |
|---|---|
| `data/curated/` | the manifests — structures, affinity, interaction, per-complex features |
| `data/splits/` | leakage-aware splits: random, homology-clustered, novel-receptor, novel-both |
| `data/*.md` | `DATASET.md`, `STRATIFICATION.md`, `SOURCES.md` (per-source licence and redistribution terms), `QC_REPORT.md` |
| `scripts/curate/` | regenerates the dataset from its sources — `run_all.sh`, deterministic |
| `scripts/run/track_{a,b,c,d}/` | the per-method runners and scorers for the four tracks |
| `scripts/release/` | builds the submission bundle and the Zenodo payload |
| `methods/` | environment specs and installers for the external tools |
| `prediction_inputs/` | the input manifest each runner folds from |

## The four tracks

- **A — cofolding.** 8 models on a 388-complex canonical tier and 4 CCD-capable cofolders on a
  225-complex non-canonical tier, both restricted to depositions after **2023-06-01**, the latest
  training cutoff among the models compared — the only boundary that is clean for all of them.
- **B — classic docking.** 4 programs on one 183-complex worklist. The headline is coverage:
  34% of the worklist is accepted by none of them under their default workflows. HADDOCK3 docks
  blind; the other three receive a native-derived box, so their rows are never pooled.
- **C — affinity.** Deliberately the smallest track, 356 measurements, because most source
  records could not be shown to describe the binding reaction they were attributed to.
- **D — interaction.** Frozen language-model embeddings plus a learned head. AUROC falls from
  0.903–0.918 on a random split to 0.663–0.734 when receptor separation is enforced pairwise.

## Two results about benchmark construction

Independent of any model, and the ones we would ask others to adopt:

1. **Define the scoring unit.** Reading DockQ's global average over all chain pairs admits
   interfaces internal to the receptor, which compressed the apparent differences between models.
2. **Report the negative-sampling design with its pool size.** Sharing one decoy pool between
   train and test raises AUROC by 0.239 at a 400-peptide pool and 0.062 at 2,000.

## Reproducing

```bash
bash scripts/curate/run_all.sh     # rebuild the dataset from its sources (no GPU)
```

Every released table regenerates from the per-complex rows in the Zenodo deposit, and the
manuscript build refuses to run unless each one still hashes to the scoring pass that produced it.

## Data and licence

Zenodo deposit (data, predictions, per-complex scores, analysis harness): DOI in the paper's Data
Availability Statement. `PepArena_scores.csv` in that deposit is the flat table — one row per
complex per model, with the stratification axes beside the metrics.

Code is MIT (`LICENSE`). Redistributed source data keeps the licence of its origin; per-source
terms are in `data/LICENSES.md` and provenance in `data/SOURCES.md`.
