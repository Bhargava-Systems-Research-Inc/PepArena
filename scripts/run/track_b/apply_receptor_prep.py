import pandas as pd
from pathlib import Path

F = Path("/12TBDrive1/mega_pep_bench/runs/track_b/fresh")
cov = pd.read_csv(F / "coverage.csv")
rp = pd.read_csv(F / "unidock_receptor_prep.csv").set_index("complex_id")
bad = rp[~rp.receptor_prepares]
before = int((cov.unidock == "ok").sum())
for cid, r in bad.iterrows():
    m = cov.complex_id.astype(str) == str(cid)
    cov.loc[m, "unidock"] = (f"refused: receptor prep failed "
                             f"({int(r.rec_atoms)} atoms exceeds the 99,999 PDB serial limit)")
cov.to_csv(F / "coverage.csv", index=False)

acc = {p: int((cov[p] == "ok").sum()) for p in ("haddock3", "adcp", "unidock", "rapidock")}
anyok = ((cov.haddock3 == "ok") | (cov.adcp == "ok") | (cov.unidock == "ok")
         | (cov.rapidock == "ok"))
print(f"Uni-Dock accepted {before} -> {acc['unidock']}")
print("acceptance:", {k: f"{v} ({100*v/len(cov):.0f}%)" for k, v in acc.items()})
print(f"accepted by none of four: {int((~anyok).sum())} ({100*(~anyok).mean():.1f}%)")
print(f"accepted by any: {int(anyok.sum())} ({100*anyok.mean():.1f}%)")
u = cov[cov.unidock == "ok"]
print("Uni-Dock under caps:", {c: int((u.torsdof <= c).sum()) for c in (15, 20, 25, 32)})
