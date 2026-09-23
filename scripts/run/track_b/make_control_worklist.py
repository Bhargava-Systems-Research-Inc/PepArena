import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "runs/track_b/control"
OUT.mkdir(parents=True, exist_ok=True)

fl = pd.read_parquet(R / "runs/track_a/foldlist.parquet")
main = set(pd.read_csv(R / "runs/track_b/fresh/worklist.csv").complex_id)

can = fl[(~fl.is_cyclic) & (~fl.has_ncaa) & (~fl.complex_id.isin(main))].copy()

can = can.sort_values(["peptide_len", "complex_id"])
bins = pd.cut(can.peptide_len, [0, 5, 10, 15, 25, 60])
sel = can.groupby(bins, observed=True).head(8)

cols = ["complex_id", "receptor_pdb_id", "receptor_chains", "peptide_chain", "peptide_seq",
        "peptide_len", "is_cyclic", "cyclization_type", "has_ncaa", "native_path",
        "deposition_date"]
sel[cols].to_csv(OUT / "worklist.csv", index=False)
print(f"control worklist: {len(sel)} canonical linear complexes -> {OUT/'worklist.csv'}")
print(sel.groupby(pd.cut(sel.peptide_len, [0, 5, 10, 15, 25, 60]), observed=True)
      .size().rename("n").to_string())
assert not sel.is_cyclic.any() and not sel.has_ncaa.any(), "control must be canonical and linear"
assert not (set(sel.complex_id) & main), "control must not overlap the main worklist"
