import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix (v1.0.0)", "af3": "AlphaFold3", "boltz2": "Boltz-2",
       "chai1": "Chai-1", "highfold3": "HighFold3", "helixfold3": "HelixFold3",
       "rfaa": "RoseTTAFold-AA", "esmfold2": "ESMFold2 (ESMC-6B)"}
KEYS = list(LAB)

sc = {k: pd.read_parquet(R / f"runs/track_a/scores/{k}.parquet").set_index("complex_id")
      for k in KEYS}
common = sorted(set.intersection(*(set(v.index) for v in sc.values())))
lk = pd.read_parquet(R / "data/splits/leakage.parquet").set_index("complex_id")
clu = lk.reindex(common)["receptor_cluster"].fillna("none")
dq = pd.DataFrame({k: sc[k].loc[common, "DockQ"].values for k in KEYS}, index=common)
groups = [np.flatnonzero((clu == c).values) for c in clu.unique()]

ref = dq.mean().idxmax()
print(f"reference (highest mean): {LAB[ref]}   common set {len(common)} complexes, "
      f"{len(groups)} receptor clusters")

rng = np.random.default_rng(0)
rows = []
for k in KEYS:
    if k == ref:
        continue
    d = (dq[ref] - dq[k]).values
    obs = d.mean()
    null = []
    for _ in range(10000):
        signs = rng.choice([-1.0, 1.0], size=len(groups))
        v = d.copy()
        for g, sgn in zip(groups, signs):
            v[g] *= sgn
        null.append(v.mean())
    null = np.asarray(null)
    p = float((np.abs(null) >= abs(obs)).mean())
    rows.append({"model": LAB[k], "mean_delta": round(float(obs), 3), "p_raw": p})

out = pd.DataFrame(rows).sort_values("p_raw").reset_index(drop=True)
m = len(out)
adj, running = [], 0.0
for i, p in enumerate(out.p_raw):
    running = max(running, (m - i) * p)
    adj.append(min(1.0, running))
out["p_holm"] = [round(x, 4) for x in adj]
out["survives_holm_0.05"] = out.p_holm < 0.05
out.to_csv(R / "runs/track_a/paired_holm_canonical.csv", index=False)
print(out.to_string(index=False))
print(f"\nsurviving Holm correction at 0.05: {int(out['survives_holm_0.05'].sum())} of {m}")
