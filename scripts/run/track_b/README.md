# Track B — classic docking (poses)

Two sub-tracks. **(1) Scoring/sampling power** on pre-generated ADCP decoys — *zero docking
compute*, ready now. **(2) Fresh docking** (HADDOCK3/ADCP/rDock) on the same leakage-clean
complexes Track A cofolded — heavy, scripts provided.

## Metrics (CASF/CPSet convention)
- **sampling power** — fraction of complexes with any near-native pose (best-of-N DockQ ≥ 0.23 / 0.49)
- **docking power** — fraction whose top-ranked pose is near-native (top-1 DockQ ≥ 0.23)
- **scoring power** — Spearman(pose rank, DockQ): does the tool's own ranking track quality?

## Sub-track 1 — decoy scoring power (READY, no docking needed)
Source: `data/external/MM_PBGBSA-CP` (cloned + extracted) — 81 cyclic complexes, 100 ADCP decoy
poses each (8,100 total; filename N = ADCP rank), 27 with Kd.

```
python prep_scoring_set.py                                   # assemble 81 native complexes + manifest
conda run -n pepgym-curate python score_docking.py           # DockQ every decoy (resumable, ~4 h CPU)
# -> runs/track_b/scoring/{leaderboard.json, power_per_complex.csv, dockq_all.parquet}
```
Smoke-tested (2 complexes): DockQ scores are sane; harness works end-to-end.

## Sub-track 2 — fresh docking (scripts; heavy)
Docks the same complexes Track A cofolded → direct classic-vs-AI contrast on identical,
leakage-controlled targets.

```
python prep_fresh_docking.py            # worklist (default 183: 143 cyclic + 40 linear; --all for 605)
python prep_complex_pdbs.py             # split natives -> inputs/<cid>/{receptor,peptide}.pdb + fasta
nohup bash run_haddock3.sh &            # HADDOCK3 (installed, env haddock3)
nohup bash run_adcp.sh &                # ADCP (install first: methods/install/install_adcp_adfrsuite.sh)
conda run -n pepgym-curate python score_fresh.py --tool haddock3   # DockQ the outputs
```

## Methods & status
| Method | Status | Notes |
|---|---|---|
| ADCP (decoy scoring) | **ready** | reuse MM_PBGBSA-CP decoys; `score_docking.py` |
| HADDOCK3 | **installed** (`haddock3`) | `run_haddock3.sh` (CoM protocol; add cyclic restraints as needed) |
| ADCP (fresh) | install | `methods/install/install_adcp_adfrsuite.sh`; `run_adcp.sh` (native `-cyc`/`-cys`) |
| Uni-Dock | installed (`unidock_env`) | peptide-as-ligand option |
| rDock / CABS-dock | install | `methods/install/install_{rdock,cabsdock}.sh` |
| Rosetta FlexPepDock | install | `methods/install/install_rosetta.sh` (license-gated) |
| HPEPDOCK2 / MDockPeP2 / GalaxyPepDock | web | cite or script submission; reproduce HPEPDOCK2 61% vs ADCP 39% (18-disulfide) as harness check |
