# Dataset provenance ledger

Every source dataset that feeds the unified manifest. Raw bytes live in `data/raw/` or
`data/external/` (gitignored); only this ledger + the curated manifests are tracked.

For each source: what it provides, scale, peptide-type coverage, license, format, and the
manifest type(s) it feeds. **Leakage note** flags whether the set is usable as a clean test
(time-split / dedup status).

Status column: no = not pulled, wip = pulling/inspecting, yes = ingested into manifest.

| # | Source | Type | Provides | Scale | Cyclic? | License | Status |
|---|---|---|---|---|---|---|---|
| 1 | PepPCBench (zhaisilong) | structure | PDB complexes + DockQ/Fnat/L-RMSD, time-split 2023–24 | 261 | linear | open (GitHub) | yes |
| 2 | FoldBench (BEAM-Labs) | structure | all-atom benchmark incl. peptide interfaces, low-homology | ~51 pep | mixed | open (GitHub) | no |
| 3 | PepSet | structure | nonredundant pep-prot complexes, IL_RMSD/Fnat | 185 | linear | open | no |
| 4 | PepPro | structure | nonredundant + unbound receptors, SS-stratified | 89 | linear | open | no |
| 5 | LNR / LNR_clean (Tsaban) | structure | curated nonredundant test set | 96 (93 ingested) | linear | open (Zenodo 13373108) | yes |
| 6 | LEADS-PEP | structure | flexible-peptide docking benchmark | 53 | linear | open | no |
| 7 | CPSet (10-program cyclic, JCIM 2024) | structure | diverse protein–cyclic-peptide complexes | ~ | cyclic | open | no |
| 8 | Cyclic disulfide set (HPEPDOCK2, 2022) | structure | disulfide-cyclized pep complexes | 18 | cyclic(SS) | open | no |
| 9 | Cyclic set (ADCP / CABS-dock) | structure | backbone/disulfide cyclic complexes | 38 | cyclic | open | no |
| 10 | HighFold / AfCycDesign cyclic | structure | cyclic peptide NMR/PDB | 63/80 | cyclic | open | no |
| 11 | HighFold3 ncAA set | structure | cyclic peptides w/ unnatural AAs | ~ | cyclic+ncAA | open | no |
| 12 | CyclicBoltz1 set | structure | cyclic ncAA complexes | ~17 | cyclic+ncAA | open | no |
| 13 | NCPepFold set | structure | non-canonical cyclic peptides | ~ | cyclic+ncAA | open | no |
| 14 | CycPeptMPDB | affinity/struct | cyclic peptide SMILES + permeability | 7,334 | cyclic | web | no |
| 15 | CREMP / CREMP-CycPeptMPDB | struct(conf) | macrocycle conformers (xTB) | 36k | cyclic | CC-BY-4.0 (Zenodo) | no |
| 16 | PPB-Affinity | affinity | pep/prot binding affinities + structures | 12,062 (2,001 peptide) | mixed | CC-BY-4.0 (Zenodo) | yes |
| 17 | PEPBI | affinity | pred+exp complexes w/ ΔG/ΔH/ΔS | 329 | linear | Dryad (auth-gated dl) | wip deferred |
| 18 | PpI[S/A]DS / PpI[S/A]BM | affinity | structure-based pep-prot affinity benchmark | ~ | linear | open | no |
| 19 | MM_PBGBSA-CP datasetI (Zhao 2025) | affinity | protein–cyclic-peptide Kd + receptor/peptide structures | 50 | cyclic | CC-BY (GitHub) | yes |
| 20 | PepPAP sets (T100/PDZ/CXCR4) | affinity/screen | sequence-level affinity + screening | — | linear | open | no |
| 21 | CAMP benchmark | interaction | binary PepPI + binding residues | — | linear | open (reproduce-only) | no |
| 22 | TPepPro pair set | interaction | pep-prot interaction pairs | 19,187 | linear | open | yes |
| 23 | PepBAN / PepInter sets | interaction | PepPI benchmarks | — | linear | open | no |
| 24 | PepLand | interaction/repr | canonical + non-canonical peptide reps | — | mixed+ncAA | open | no |
| 25 | cycproteinix complex34 (local curated) | structure | cyclic pep–protein complexes + release dates + native refs | 34 | cyclic | local | yes |
| 26 | AfCycDesign/HighFold monomers (local) | structure(monomer) | cyclic peptide monomers w/ explicit modes+disulfides | 80 | cyclic | local | yes |
| 27 | 2025+ holdout (local) | structure | post-2025 leakage-clean cyclic+linear complexes | 6+18 | mixed | local | no |
| 28 | PDB ncAA mine | structure | non-canonical peptide–protein complexes (CCD-code search) | 885 | mixed+ncAA | PDB/CC0 | yes |
| 29 | PepBench / PepGLAD | structure | large peptide–protein complex pool (NeurIPS'24) | 6105 | mixed | open (Zenodo 13373108) | yes |

## Detailed records

(Filled in as each source is pulled and inspected — URL, commit/version, files, columns,
counts, and any caveats.)

### 1. PepPCBench  yes
- **Repo:** https://github.com/zhaisilong/PepPCBench (cloned to `data/external/PepPCBench`).
- **License:** see repo LICENSE (open).
- **What we use:** `job_list.csv` = 261 protein–peptide complexes (pdb_id + protein/peptide
  chains). Structures fetched from RCSB (`data/raw/structures/cif/`, gitignored).
  Repo also ships `clusters/` (mmseqs2 + foldseek test-vs-train clustering) and `results/`
  (AF3 DockQ/confidence/ΔG) — reusable later for homology splits and as a prior baseline.
- **Ingested by:** `scripts/curate/ingest_peppcbench.py` → `data/curated/structure_peppcbench.parquet`.
- **Profile (261):** length 5–29 (median 10); bins 1-5:14, 6-10:150, 11-15:63, 16-25:30,
  26-50:4. Methods: xray 242, em 14, nmr 5. Median resolution 2.06 Å. Deposition 2023–2024
  (true post-AF3-cutoff novelty set).
- **Cyclization:** 260 linear + 1 thioether-stapled (8tor, CH3–SG side-chain bond → "other").
  Essentially a linear test set — cyclic coverage comes from sources 7–13.
- **Caveat fixed during ingest:** naïve C–N covale detection false-flagged 14 linear peptides
  (explicit peptide-bond `covale` records). Head-to-tail now requires the bond to span the
  chain's first↔last residue. Disulfide = peptide SG–SG. Validated by manual bond inspection.

### 25. cycproteinix complex34 (local curated cyclic complexes)  yes
- **Source:** `/12TBDrive1/cycproteinix/benchmark/complex34/` (the user's prior cyclic-peptide
  benchmark). `metadata.csv` = 34 cyclic peptide–protein complexes with RCSB release dates,
  peptide/receptor chains + lengths, and native reference PDBs (`refs/*.pdb`).
- **Ingested by:** `scripts/curate/ingest_cyclic_complexes.py` →
  `data/curated/structure_cyclic_complexes.parquet`. Refs copied to
  `data/raw/structures/cyclic_complex34/` (gitignored).
- **Cyclization (re-derived with the shared detector, not trusting source labels):**
  30 head_to_tail + 4 bicyclic (head-to-tail + disulfide). **Validation:** 1SFI → bicyclic,
  correct (SFTI-1 is the canonical head-to-tail + disulfide bicyclic peptide).
- **Profile:** length 6–34; release 1999–2025, **including 5 post-2024 entries**
  (9CDU/9CDZ/9HGC/9HGD/9PHQ) usable as a leakage-clean cyclic novelty slice.
- **Companion local sets (queued):** the post-2025 holdout (6 cyclic + 18 linear).

### 26. AfCycDesign/HighFold cyclic monomers (local)  yes
- **Source:** `/12TBDrive1/cycesmfold2/data/benchmark/monomer_manifest.csv` (80 cyclic peptide
  monomers — the AfCycDesign/HighFold benchmark) + native ref PDBs.
- **Ingested by:** `scripts/curate/ingest_cyclic_monomers.py` →
  `data/curated/structure_monomer_cyclic.parquet` (separate **monomer** sub-track, no receptor).
- **Cyclization (our detector):** 63 bicyclic + 17 head_to_tail. **Cross-check vs the source's
  hand-curated mode labels: 80/80 agree (0 disagreements)** — independent validation of the
  geometric+conn detector.
- **Leakage:** classic pre-2020 PDB structures → training-contaminated for all current methods;
  flagged legacy/not-clean. Use for capability comparison, not novelty claims.
- **ncAA:** these 80 are canonical (disulfide/head-to-tail). Genuine ncAA still needs external
  HighFold3/NCPepFold/CyclicBoltz1/CycPeptMPDB sets.

### 16. PPB-Affinity (peptide subset)  yes
- **Source:** Zenodo 14271435 (CC-BY-4.0). Downloaded `PPB-Affinity.xlsx` (1.2 MB) to
  `data/raw/affinity/` (skipped the 3 GB PDB.zip — not needed for the values).
- **Ingested by:** `scripts/curate/ingest_ppb_affinity.py` → `data/curated/affinity_ppb.parquet`.
  Peptide subset = rows whose ligand is a single polypeptide 4–50 aa bound to a protein
  receptor; ligand lengths + sequences fetched from RCSB GraphQL (the xlsx lacks them).
- **Profile:** 2,001 measurements, 427 unique peptides, 485 PDBs, length 8–50 (median 13),
  log10 Kd −12.9 .. −2.5 (pM→mM, ~10 orders). Sources: SKEMPI (mutational), PDBbind, SAbDab,
  ATLAS, Affinity Benchmark — recorded per-row in `notes` (`src=`, `mut=`).
- **Caveats:** many rows are SKEMPI mutational variants of the same pair (use for ΔΔG /
  ranking; filter `mut=nan` for wildtype-only). **Cyclization = unknown** (no structures
  downloaded yet) — flagged for a later structure-based labeling pass.

### 22. TPepPro interaction pairs  yes
- **Repo:** github.com/wanglabhku/TPepPro (cloned to `data/external/TPepPro`).
- **Files:** `receptor-peptide.actions.tsv` (19,187 pairs, **already balanced** 9594 pos /
  9593 neg) + `receptor(14374)-peptide(9594)_dictionary.xlsx` (id→sequence, Sheet1).
- **Ingested by:** `scripts/curate/ingest_tpeppro.py` → `data/curated/interaction_tpeppro.parquet`.
- **Profile:** 19,187 pairs, 7,892 unique receptors, 5,643 unique peptides, peptide length
  2–50 (median 13). Linear/canonical; cyclization unknown. 0 dropped (all ids resolved).
- **Note:** CAMP (github.com/twopin/CAMP) ships only a 100-row example; the full CAMP
  benchmark must be reproduced via its `data_prepare/` pipeline (RCSB/UniProt/DrugBank) —
  deferred (marked reproduce-only).

### 28. PDB ncAA mine  yes
- **Method:** `scripts/curate/ingest_pdb_ncaa.py`. RCSB count attributes aren't search-enabled,
  but `rcsb_polymer_entity_container_identifiers.chem_comp_monomers` is — so we search for short
  (4–50) polypeptide entities containing any of a curated set of modified-residue CCD codes
  (N-methyl, D-amino acids, Aib, Nle, Orn, Hyp, …), keep entries that also have a protein
  receptor, download the cif, and auto-label with the shared ncAA + cyclization detectors.
- **Yield:** 885 non-canonical peptide–protein complexes (from 1962 candidate entries; 964
  lacked a peptide+receptor pair). Cyclization 421 linear / 197 head_to_tail / 198 side-chain
  macrocycle / 69 disulfide. Top chemistries: DAL, HYP, MVA, SAR, NLE, MLE, AIB, BMT
  (cyclosporin). Length 4–48, deposition 1991–2026.
- **License:** PDB / CC0 (we ship the index + ingester, not bulk structures).
- **Caveat:** broad/heterogeneous (all dates → many pre-cutoff; use leakage time tiers).
  Side-chain macrocycles bucketed as `other`; SMILES-level subtyping (lactam/thioether/staple)
  is the next refinement.

### 5 + 29. PepBench / PepGLAD + LNR (Zenodo 13373108)  yes
- **Source:** PepBench (NeurIPS'24 PepGLAD) on Zenodo: `LNR.tar.gz` (curated 96-complex test,
  93 ingested) + `train_valid.tar.gz` (6105 complexes, pre-cleaned receptor+peptide PDBs;
  chains in filename `{rec}_{pep}_pdb{id}`).
- **Ingested by:** `ingest_pdb_list.py` (LNR; fetches cif) and `ingest_pepbench.py` (PepBench;
  uses local PDBs + batched RCSB dates). LNR: 90 linear/2 disulfide/1 thioether. PepBench:
  5943 linear / 128 disulfide / 31 head_to_tail / 3 bicyclic, **1137 ncAA**, dated 6103/6105.
- **Role:** PepBench is the benchmark's large structure pool (supervised-training scale +
  breadth); the curated test sets (LNR, PepPCBench cyclic) stay distinct. After dedup the
  combined structure set is **6933 complexes** (445 cross-source overlaps removed).
- **Note (PEPBI, source 17):** Dryad file download is auth-gated (bearer token); the 329-row
  thermodynamic affinity set (ΔG/ΔH/ΔS) is deferred until obtained via the Dryad UI.

### 19. MM_PBGBSA-CP — cyclic-peptide affinity benchmark (Zhao 2025)  yes
- **Repo:** github.com/huifengzhao/MM_PBGBSA-CP (CC-BY). `dataset1.csv` + `dataset1.tar.bz2`
  = 50 protein–cyclic-peptide complexes with measured Kd/pKd + cyclization (`Chain type`:
  SS*n / BB / BB+SS), each as separate `{pdb}_protein.pdb` + `{pdb}_CP.pdb`; repo also ships
  minimized + ADCP decoy structures (for physics scoring + docking scoring-power).
- **Ingested by:** `ingest_mmpbgbsa_cp.py` → `data/curated/affinity_mmpbgbsa_cp.parquet`
  (a dedicated **cyclic** affinity manifest for the physics-affinity comparison MM-PBSA/FEP
  vs ML). 50 cyclic: 35 disulfide / 8 bicyclic / 7 head_to_tail; pKd 3.55–10.29.
- **Detector cross-check:** 48/50 agree with the paper's Chain-type (2 head_to_tail edge cases
  use the paper label). Surfaced + fixed an **AMBER residue-naming bug** in `struct_utils`
  (CYX/HID/… were dropped, corrupting seq/length on MD-processed structures) — now CYX→CYS,
  caps ACE/NME/NHE skipped; existing RCSB-based manifests byte-identical (no regression).
- **Companion (same group):** `dataset2` (81 complexes + ADCP decoys, no Kd) is the docking
  scoring-power set for Track B; `huifengzhao/CPSet` is the 10-program cyclic complex set.

### ncAA cofolding-benchmark gap (still tracked)
The dedicated ncAA cyclic *cofolding* benchmarks — CyclicBoltz1 (bioRxiv 2025.02.11.637752),
NCPepFold (2024.12.05.626948), HighFold3 ncAA — have **no public code/data repo** as of
June 2026. The PDB ncAA mine (source 28) is our open substitute for the ncAA axis; CycPeptMPDB/
CREMP remain available for a peptide-level (no-receptor) chemistry sub-track if needed.
