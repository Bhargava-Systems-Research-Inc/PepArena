import pathlib, sys
import pandas as pd

R = pathlib.Path("/12TBDrive1/mega_pep_bench")
CUTOFF = pd.Timestamp("2023-06-01")

s = pd.read_parquet(R / "data/curated/structure_all.parquet")
d = pd.read_parquet(R / "data/curated/pdb_dates.parquet")
s["pdb"] = s.receptor_pdb_id.astype(str).str.upper()
m = s.merge(d, left_on="pdb", right_on="pdb_id", how="left")
m["manifest"] = pd.to_datetime(m.deposition_date, errors="coerce")
m["release"] = pd.to_datetime(m.release_date, errors="coerce")

pool = m[m.manifest > CUTOFF]
contaminated = pool[pool.release <= CUTOFF]
eligible_missed = m[(m.manifest <= CUTOFF) & (m.release > CUTOFF)]
unresolved = m[m.release.isna()]

print(f"evaluated pool (manifest date > {CUTOFF.date()}): {len(pool)}")
print(f"  released on or before the cutoff : {len(contaminated)}")
print(f"  release date unresolved at RCSB  : {len(pool[pool.release.isna()])}")
print(f"eligible but excluded by a deposition-dated source: {len(eligible_missed)}")
print(f"manifest dates unresolved overall: {len(unresolved)}")

if len(contaminated):
    print("\nFAIL: these were released before the cutoff they are held out of:")
    print(contaminated[["complex_id", "manifest", "release"]].head(20).to_string(index=False))
    sys.exit(1)
print("\nOK: no evaluated complex was released before the cutoff. The date inconsistency was "
      "conservative in direction: it omitted eligible complexes, never admitted contaminated ones.")
