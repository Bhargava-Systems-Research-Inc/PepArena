from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CUR = ROOT / "data/curated"
OUT = ROOT / "runs/track_a/foldlist.parquet"

CUTOFFS = {
    "af2_multimer_v1": "2018-04-30",
    "rfaa": "2020-05-01",
    "esmfold2": "2020-05-01",
    "chai1": "2021-01-12",
    "af2_multimer_v23": "2021-09-30",
    "alphafold3": "2021-09-30",
    "boltz1": "2021-09-30",
    "protenix": "2021-09-30",
    "helixfold3": "2021-09-30",
    "highfold3": "2021-09-30",
    "boltz2": "2023-06-01",
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-date", default="2023-06-01")
    ap.add_argument("--max-receptor", type=int, default=1800)
    args = ap.parse_args()

    sa = pd.read_parquet(CUR / "structure_all.parquet").drop(columns=["receptor_fasta"], errors="ignore")
    ni = pd.read_parquet(CUR / "structure_native_index.parquet").set_index("complex_id")
    sa = sa.join(ni[["native_path", "receptor_fasta", "n_receptor_res", "n_peptide_res"]], on="complex_id")
    sa["dep"] = pd.to_datetime(sa.deposition_date, errors="coerce")

    fl = sa[sa.dep > pd.Timestamp(args.min_date)].copy()
    fl["receptor_seq"] = fl.receptor_fasta
    fl["receptor_len"] = fl.n_receptor_res
    fl["foldable"] = fl.receptor_len.fillna(9999).le(args.max_receptor) & \
        fl.peptide_seq.fillna("").str.len().gt(0)
    for m, c in CUTOFFS.items():
        fl[f"ok_{m}"] = fl.dep > pd.Timestamp(c)

    cols = ["complex_id", "receptor_pdb_id", "receptor_chains", "peptide_chain", "receptor_seq",
            "peptide_seq", "peptide_len", "receptor_len", "is_cyclic", "cyclization_type",
            "has_ncaa", "deposition_date", "native_path", "foldable"] + [f"ok_{m}" for m in CUTOFFS]
    out = fl[cols].reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)
    print(f"fold-list (deposited > {args.min_date}): {len(out)} complexes -> {OUT}")
    print(f"  foldable (receptor <= {args.max_receptor} aa): {int(out.foldable.sum())} "
          f"(cyclic {int(out[out.foldable].is_cyclic.sum())})")
    print(f"  oversized receptors skipped: {int((~out.foldable).sum())}")
    print("  cyclization (foldable):", out[out.foldable].cyclization_type.value_counts().to_dict())
    print("  per-method eligible (within this fold-list):",
          {m: int(out[f'ok_{m}'].sum()) for m in CUTOFFS})

if __name__ == "__main__":
    main()
