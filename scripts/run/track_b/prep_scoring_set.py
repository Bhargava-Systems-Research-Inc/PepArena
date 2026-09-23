#!/usr/bin/env python3
import os
import glob
import csv
import pandas as pd

REPO = "/12TBDrive1/mega_pep_bench"
SRC = f"{REPO}/data/external/MM_PBGBSA-CP"
OUT = f"{REPO}/runs/track_b/scoring"
NAT = f"{OUT}/native"
os.makedirs(NAT, exist_ok=True)

def relabel_chain(pdb_path, chain):
    out = []
    with open(pdb_path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                out.append(ln[:21] + chain + ln[22:])
    return out

def main():
    pkd = {}
    with open(f"{SRC}/dataset1.csv") as fh:
        for row in csv.DictReader(fh):
            pkd[row["PDB code"].upper()] = float(row["pKd"])

    rows = []
    for d in sorted(glob.glob(f"{SRC}/dataset2/*")):
        pdb = os.path.basename(d)
        prot = f"{d}/{pdb}_protein.pdb"
        pep = f"{d}/{pdb}_CP.pdb"
        decoys = sorted(glob.glob(f"{SRC}/decoy/{pdb}/*.pdb"))
        if not (os.path.exists(prot) and os.path.exists(pep) and decoys):
            print("skip (missing)", pdb)
            continue
        lines = relabel_chain(prot, "R") + ["TER\n"] + relabel_chain(pep, "P") + ["TER\n", "END\n"]
        with open(f"{NAT}/{pdb}.pdb", "w") as fh:
            fh.writelines(lines)
        pep_len = len({ln[22:26] for ln in relabel_chain(pep, "P")})
        rows.append({"pdb": pdb, "n_decoys": len(decoys), "peptide_len": pep_len,
                     "pKd": pkd.get(pdb.upper())})

    m = pd.DataFrame(rows).sort_values("pdb")
    m.to_parquet(f"{OUT}/manifest.parquet")
    m.to_csv(f"{OUT}/manifest.csv", index=False)
    print(f"prepared {len(m)} native complexes -> {NAT}")
    print(f"  decoys/complex: min {m.n_decoys.min()} max {m.n_decoys.max()} "
          f"(total {int(m.n_decoys.sum())})")
    print(f"  with pKd: {m.pKd.notna().sum()}  | peptide_len {m.peptide_len.min()}-{m.peptide_len.max()}")

if __name__ == "__main__":
    main()
