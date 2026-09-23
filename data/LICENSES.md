# Data licensing & redistribution

The PepArena **harness code** is MIT (`/LICENSE`). The **data** is derived from external
sources, each under its own license. We track curated *manifests* (which we author) plus the
*provenance* to regenerate everything; large raw bytes (structures, the affinity workbook) are
not redistributed in-repo — `scripts/curate/run_all.sh` pulls them from the original providers.

| Source | License | In-repo | Redistribution |
|---|---|---|---|
| RCSB PDB structures (PepPCBench ids, materialized natives) | PDB / CC0 | natives gitignored (regenerable) | CC0 — free; we ship the index, not bulk PDBs |
| PepPCBench (ids, clusters) | repo LICENSE (open) | ids only | pull from GitHub |
| cycproteinix complex34 / AfCycDesign monomers | derived from RCSB (CC0) | refs gitignored | regenerate from PDB |
| PPB-Affinity | CC-BY-4.0 (Zenodo 14271435) | values gitignored (xlsx) | CC-BY; attribute |
| TPepPro pairs + dictionary | repo (open) | manifest only | pull from GitHub |
| TPepPro / CAMP code | their repo licenses | not vendored | — |

**Manifests** (`data/curated/*.parquet`) are our derived annotations over public data and are
released with the harness under MIT, with attribution to the underlying sources (see
`SOURCES.md`). When redistributing, cite the original datasets and respect CC-BY attribution
for PPB-Affinity.

> The harness MIT license is a sensible default for research code; change it here if the
> project needs a different license (e.g. Apache-2.0, CC-BY-4.0 for the manifests).
