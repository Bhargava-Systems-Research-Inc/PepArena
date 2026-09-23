from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "runs/track_b/fresh"

fl = pd.read_parquet(R / "runs/track_a/foldlist.parquet")
cols = ["complex_id", "receptor_pdb_id", "receptor_chains", "peptide_chain", "peptide_seq",
        "peptide_len", "is_cyclic", "cyclization_type", "has_ncaa", "native_path",
        "deposition_date"]

cyc = fl[fl.is_cyclic]
lin = fl[~fl.is_cyclic].sort_values(["peptide_len", "complex_id"])
lin = lin.groupby(pd.cut(lin.peptide_len, [0, 5, 10, 15, 25, 60]), observed=True).head(8)
sel = pd.concat([cyc, lin]).drop_duplicates("complex_id")[cols]

prev = pd.read_csv(OUT / "worklist.csv")
old, new = set(prev.complex_id), set(sel.complex_id)
sel.to_csv(OUT / "worklist.csv", index=False)
prev.to_csv(OUT / "worklist_preaudit.csv", index=False)

print(f"re-derived worklist: {len(sel)} complexes ({int(sel.is_cyclic.sum())} cyclic, "
      f"{int(sel.has_ncaa.sum())} non-canonical)")
print(f"  carried over from the pre-audit worklist: {len(old & new)}")
print(f"  added: {len(new - old)}   dropped: {len(old - new)}")
print(f"  previous worklist kept as worklist_preaudit.csv")
