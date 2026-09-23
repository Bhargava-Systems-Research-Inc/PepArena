import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
p = R / "runs/track_b/control/worklist.csv"
w = pd.read_csv(p)
m = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
for c in ("peptide_len", "peptide_seq", "has_ncaa"):
    w[f"{c}_orig"] = w[c]
    w[c] = w.complex_id.map(m[c])
w["peptide_len"] = w["peptide_len"].astype(int)
w["has_ncaa"] = w["has_ncaa"].astype(bool)
w["ncaa_codes"] = w.complex_id.map(m["ncaa_codes"]).fillna("")
w.to_csv(p, index=False)
changed = int((w.peptide_len != w.peptide_len_orig).sum())
print(f"control worklist: {len(w)} complexes, {changed} lengths corrected, "
      f"{int(w.has_ncaa.sum())} now carry a non-canonical residue")
print(w.loc[w.has_ncaa, ["complex_id", "peptide_len_orig", "peptide_len", "ncaa_codes"]]
       .to_string(index=False))
