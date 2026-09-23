#!/usr/bin/env python3
import os
import hashlib
import pandas as pd

REPO = "/12TBDrive1/mega_pep_bench"
OUT = f"{REPO}/runs/track_c"
os.makedirs(OUT, exist_ok=True)

KEEP = ["complex_id", "source_dataset", "receptor_pdb_id", "receptor_chains", "peptide_chain",
        "peptide_seq", "peptide_len", "is_cyclic", "cyclization_type", "has_ncaa",
        "target_class", "deposition_date", "measured_value", "value_type", "log_affinity"]

def h(s):
    return int(hashlib.md5(str(s).encode()).hexdigest(), 16)

def main():
    ppb = pd.read_parquet(f"{REPO}/data/curated/affinity_ppb.parquet")
    cp = pd.read_parquet(f"{REPO}/data/curated/affinity_mmpbgbsa_cp.parquet")
    df = pd.concat([ppb[KEEP], cp[KEEP]], ignore_index=True).drop_duplicates("complex_id")
    df = df.dropna(subset=["log_affinity"])

    df["split_random"] = df["complex_id"].map(lambda c: "test" if h(c) % 5 == 0 else "train")
    df["split_receptor"] = df["receptor_pdb_id"].map(
        lambda p: "test" if h(p) % 5 == 0 else "train")

    df.to_parquet(f"{OUT}/affinity_manifest.parquet")
    df[["complex_id", "receptor_pdb_id", "split_random", "split_receptor"]].to_csv(
        f"{OUT}/affinity_splits.csv", index=False)

    print(f"affinity eval set: {len(df)} measurements | "
          f"{df.receptor_pdb_id.nunique()} receptors | {df.peptide_seq.nunique()} peptides")
    print(f"  log_affinity range {df.log_affinity.min():.2f}..{df.log_affinity.max():.2f}")
    print(f"  cyclic {int(df.is_cyclic.sum())} | ncAA {int(df.has_ncaa.sum())}")
    for s in ("split_random", "split_receptor"):
        vc = df[s].value_counts().to_dict()
        print(f"  {s}: train {vc.get('train',0)} / test {vc.get('test',0)}")

if __name__ == "__main__":
    main()
