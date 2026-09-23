from __future__ import annotations
import sys, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import cyclization, struct_utils

AFF = ROOT / "data/curated/affinity_ppb.parquet"
CIF_DIR = ROOT / "data/raw/structures/cif"

def fetch(pdb):
    dst = CIF_DIR / f"{pdb.lower()}.cif.gz"
    if dst.exists() and dst.stat().st_size > 0:
        return dst
    CIF_DIR.mkdir(parents=True, exist_ok=True)
    for _ in range(3):
        try:
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb.lower()}.cif.gz", dst)
            if dst.stat().st_size > 0:
                return dst
        except Exception:
            pass
    return None

def main():
    df = pd.read_parquet(AFF)
    pairs = df[["receptor_pdb_id", "peptide_chain"]].drop_duplicates()
    print(f"labeling {len(pairs)} unique (pdb, peptide_chain) for {len(df)} affinity rows...")
    label = {}
    for i, (_, p) in enumerate(pairs.iterrows()):
        pdb, ch = str(p.receptor_pdb_id), str(p.peptide_chain)
        cif = fetch(pdb)
        if cif is None:
            continue
        try:
            st = struct_utils.read_any(cif)
            ctype, _, _ = cyclization.classify(st, ch)
            _, ncaa, _ = struct_utils.chain_seq_and_ncaa(st, ch)
            label[(pdb, ch)] = (ctype, bool(ncaa), ",".join(sorted(set(ncaa))))
        except Exception:
            pass
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(pairs)}", flush=True)

    def lk(row, j, default):
        return label.get((str(row.receptor_pdb_id), str(row.peptide_chain)), (None, None, None))[j] \
            if (str(row.receptor_pdb_id), str(row.peptide_chain)) in label else default

    df["cyclization_type"] = df.apply(lambda r: lk(r, 0, "unknown"), axis=1)
    df["is_cyclic"] = df.apply(lambda r: (lk(r, 0, None) not in (None, "linear")) if lk(r, 0, None) else None, axis=1)
    df["has_ncaa"] = df.apply(lambda r: lk(r, 1, False), axis=1)
    df["ncaa_list"] = df.apply(lambda r: lk(r, 2, None), axis=1)
    df.to_parquet(AFF, index=False)

    labeled = df[df.cyclization_type != "unknown"]
    print(f"\nlabeled {len(labeled)}/{len(df)} rows ({len(label)} structures resolved)")
    print("cyclization:", df.cyclization_type.value_counts().to_dict())
    print("cyclic rows:", int((df.is_cyclic == True).sum()), "| ncAA rows:", int(df.has_ncaa.sum()))

if __name__ == "__main__":
    main()
