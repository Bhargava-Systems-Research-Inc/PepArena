#!/usr/bin/env python3
import os
import argparse
import pandas as pd

REPO = "/12TBDrive1/mega_pep_bench"
OUT = f"{REPO}/runs/track_b/fresh"
os.makedirs(OUT, exist_ok=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="emit all 605 (default: small tier)")
    ap.add_argument("--n_linear", type=int, default=40, help="linear complexes in default tier")
    a = ap.parse_args()

    fl = pd.read_parquet(f"{REPO}/runs/track_a/foldlist.parquet")
    cols = ["complex_id", "receptor_pdb_id", "receptor_chains", "peptide_chain",
            "peptide_seq", "peptide_len", "is_cyclic", "cyclization_type", "has_ncaa",
            "native_path", "deposition_date"]
    fl = fl[cols]

    if a.all:
        sel = fl
    else:
        cyc = fl[fl.is_cyclic]
        lin = (fl[~fl.is_cyclic].sort_values("peptide_len")
                 .groupby(pd.cut(fl[~fl.is_cyclic].peptide_len, [0, 5, 10, 15, 25, 60]),
                          observed=True)
                 .head(max(1, a.n_linear // 5)))
        sel = pd.concat([cyc, lin]).drop_duplicates("complex_id")

    sel.to_csv(f"{OUT}/worklist.csv", index=False)
    print(f"fresh-docking worklist: {len(sel)} complexes "
          f"({sel.is_cyclic.sum()} cyclic / {(~sel.is_cyclic).sum()} linear) -> {OUT}/worklist.csv")
    print("length bins:", pd.cut(sel.peptide_len, [0, 5, 10, 15, 25, 60]).value_counts().sort_index().to_dict())
    print("next: bash run_haddock3.sh  |  bash run_adcp.sh  (per-complex, reads worklist.csv)")

if __name__ == "__main__":
    main()
