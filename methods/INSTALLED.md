# Method inventory (host `yash`)

Snapshot of what is already installed on this box vs. what still needs setup.
Legend: yes = installed, web = web-server only (no local install), no = missing (to install),
🔎 installed but verify weights before first run.

> This session is **setup + curation only** — no methods are run. Heavy/finicky docking
> installs (Rosetta, ADFRsuite/ADCP) are deferred to when Track B starts; they are queued
> below, not done now.

## Track A — Cofolding

| Method | Status | Env / path | Version |
|---|---|---|---|
| AlphaFold3 | 🔎 | env `af3`, `af3x` | alphafold3 3.0.1 (jax 0.9.1) — verify model params |
| Boltz-1/2 | yes | env `boltz`, cache `~/.boltz` | boltz 2.2.1 (torch 2.11) |
| Chai-1 | yes | env `chai`, code `/6TBDrive1/ImmunoBench/Chai_Bench/software/chai-lab` | chai_lab 0.6.1 |
| Protenix | yes | env `protenix` | protenix 2.0.0 |
| Protenix (cyclic) | yes | env `cycproteinix`, code `/12TBDrive1/cycproteinix/Protenix` | 2.0.0 |
| OpenFold | yes | env `openfold` | openfold 2.2.0 |
| RoseTTAFold-AA | 🔎 | env `rfaa` | dgl 1.1.2 / torch 2.0.1 — verify weights |
| HighFold / HighFold3 | yes | repos `~/HighFold3`, `~/HighFold3_onerai`, `/12TBDrive1/cycesmfold2/HighFold3`; env `highfold3` (AF3-based) | — |
| ESMFold / ESM3 | yes | env `cycesmfold2`, `esmfold2_api` | esm 3.3.0 |
| AF2-Multimer / ColabFold | no | — | install localcolabfold |
| HelixFold3 | no | — | PaddlePaddle-based; install later |
| CyclicBoltz1 | no | — | Boltz extension; clone + patch later |

## Track B — Classic docking

> Track B/C harnesses built (`scripts/run/track_b`, `scripts/run/track_c`). Install scripts for
> the missing engines live in `methods/install/`; master: `scripts/run/setup_tracks_bc.sh`.

| Method | Status | Env / path | Notes |
|---|---|---|---|
| HADDOCK3 | yes | env `haddock3` | haddock3 2026.5.0 (incl. cyclic protocol); `run_haddock3.sh` |
| Uni-Dock | yes | env `unidock_env`, repo `~/Uni-Dock` | GPU small-molecule docking; peptide-as-ligand option |
| ADCP (decoy scoring) | yes ready | data `MM_PBGBSA-CP` + env `dockq` | 81 complexes × 100 decoys; `score_docking.py` (no docking compute) |
| AutoDock CrankPep (ADCP, fresh) | no script | — | `methods/install/install_adcp_adfrsuite.sh` (~500 MB); native `-cyc`/`-cys` |
| rDock | yes | env `rdock` (bioconda 24.04) | `rbdock`/`rbcavity` verified; `methods/install/install_rdock.sh` |
| CABS-dock | no script | — | `install_cabsdock.sh` (needs MODELLER key) |
| Rosetta FlexPepDock / PIPER | no script | — | `install_rosetta.sh` (academic license + build) |
| HPEPDOCK2 / MDockPeP2 / GalaxyPepDock / pepATTRACT | web | — | web; reproduce HPEPDOCK2 61% vs ADCP 39% (18-disulfide) as harness check |

## Track C — Affinity + screening

Eval set built: `runs/track_c/affinity_manifest.parquet` (2,045 measurements) + splits.
First result run: ESM-2 + GBM regressor, random Spearman 0.727 → novel-receptor 0.336.

| Method | Status | Env / path | Notes |
|---|---|---|---|
| ESM-2 + regressor | yes ran | env `pepgym-trackd` | `run_affinity_regressor.py`; reuses Track-D embedding cache |
| Boltz-2 affinity head | yes (template) | env `boltz` | `run_boltz2_affinity.py`; verify YAML schema for boltz 2.2.1 |
| gmx_MMPBSA (MM-PBSA/GBSA) | yes | env `gmxmmpbsa` | 1.6.5; single-snapshot MM-GBSA on the 50 cyclic-Kd set |
| HERMES | yes | env `hermes`, `/6TBDrive2/MoggingHERMES` | scores native complexes (materialize first) |
| PPI-Affinity | web/script | env `ppi-affinity` | web tool; open alternative SVSBI via `install_ppi_affinity.sh` |
| PepPAP | no script | — | `install_peppap_deepppiscore.sh` → `ChunhuaLab/PepPAP` (T100/PDZ/CXCR4 sets) |
| DeepPpIScore | no script | — | same script → `zjujdj/DeepPpIScore` |

## Track D — Sequence interaction prediction

| Method | Status | Notes |
|---|---|---|
| CAMP | no | GitHub; lightweight |
| TPepPro | no | GitHub |
| PepBAN | no | GitHub |
| PepInter | no | GitHub |
| PepLand | no | GitHub (canonical + non-canonical reps) |

## Evaluation tools

| Tool | Status | Env | Notes |
|---|---|---|---|
| DockQ | yes | env `dockq` | dockq 2.1.3 |
| OpenStructure (lDDT / lDDT-PLI) | no | — | install via conda/container (curation env) |
| PoseBusters | no | — | pip install (curation env) |
| RDKit | yes | base | 2026.03.1 |

## Install queue (later sessions, per track when we reach it)
- Track A gaps: localColabFold (AF2-Multimer), HelixFold3, CyclicBoltz1.
- Track B: ADFRsuite/ADCP, Rosetta+FlexPepDock+PIPER, CABS-dock, rDock.
- Track C: PPI-Affinity, PepPAP, DeepPpIScore.
- Track D: CAMP, TPepPro, PepBAN, PepInter, PepLand.

## ipSAE (added 2026-09-23)

`methods/external/ipsae/ipsae.py` — Dunbrack lab reference implementation, v4 (Jan 2026, the
version that added Boltz support), MIT. Fetched from
`https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py`. `methods/external/**` is
gitignored, so **refetch it before rerunning** `scripts/run/track_a/run_ipsae.py`. Used rather
than a reimplementation because ipSAE's normalisation is the whole point of the metric and is
easy to get subtly wrong. Also emits pDockQ, pDockQ2 and LIS, which is how we found that the
first two are constant on peptide interfaces.
