import pathlib
import pandas as pd
R = pathlib.Path("/12TBDrive1/mega_pep_bench")
s = pd.read_parquet(R / "data/curated/structure_all.parquet")
d = pd.read_parquet(R / "data/curated/pdb_dates.parquet")
s["pdb"] = s.receptor_pdb_id.astype(str).str.upper()
m = s.merge(d, left_on="pdb", right_on="pdb_id", how="left")
cut = pd.Timestamp("2023-06-01")
newly = m[(pd.to_datetime(m.deposition_date, errors="coerce") <= cut) &
          (pd.to_datetime(m.release_date, errors="coerce") > cut)]
cols = ["complex_id", "receptor_pdb_id", "peptide_chain", "peptide_len", "is_cyclic",
        "has_ncaa", "deposition_date", "release_date"]
out = R / "data/splits/newly_eligible_by_release_date.csv"
newly[cols].sort_values("complex_id").to_csv(out, index=False)
print(f"wrote {out} ({len(newly)} complexes; "
      f"{int((~newly.has_ncaa).sum())} canonical, {int(newly.is_cyclic.sum())} cyclic)")
