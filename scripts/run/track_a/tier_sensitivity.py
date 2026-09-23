import functools
import json
from pathlib import Path

import numpy as np
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix", "af3": "AlphaFold3", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3", "rfaa": "RoseTTAFold-AA",
       "esmfold2": "ESMFold2"}

ids = sorted(set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv").iloc[:, 0].astype(str)))
m = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
moved = [c for c in ids if bool(m.has_ncaa.get(c, False)) and not bool(m.has_ncaa_orig.get(c, False))]
keep = [c for c in ids if c not in moved]
print(f"tier {len(ids)}; reclassified {len(moved)}; strictly canonical {len(keep)}\n")

sc = {k: pd.read_parquet(R / f"runs/track_a/scores/{k}.parquet").set_index("complex_id")
      for k in LAB}
rows = []
for k, lab in LAB.items():
    d = sc[k]
    full = d.DockQ.mean()
    sub = d.loc[[c for c in keep if c in d.index]].DockQ.mean()
    rows.append({"model": lab, "n_full": len(d), "mean_full": round(full, 3),
                 "n_strict": len(set(keep) & set(d.index)), "mean_strict": round(sub, 3),
                 "delta": round(sub - full, 4)})
t = pd.DataFrame(rows).sort_values("mean_full", ascending=False)
print(t.to_string(index=False))
print(f"\nlargest shift in any model mean: {t.delta.abs().max():.4f}")
print(f"ordering unchanged: "
      f"{list(t.sort_values('mean_full', ascending=False).model) == list(t.sort_values('mean_strict', ascending=False).model)}")

common = functools.reduce(set.intersection, (set(v.index) for v in sc.values()))
ck = common & set(keep)
print(f"\ncommon set: {len(common)} -> {len(ck)} on the strictly canonical subset")
lead = t.model.iloc[0]
key = {v: k for k, v in LAB.items()}[lead]
ref = sc[key].loc[sorted(ck)].DockQ
n_excl = 0
for k, lab in LAB.items():
    if lab == lead:
        continue
    diff = (ref - sc[k].loc[sorted(ck)].DockQ).values
    boot = np.array([np.random.default_rng(i).choice(diff, len(diff)).mean() for i in range(2000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    n_excl += (lo > 0 or hi < 0)
print(f"paired intervals excluding zero on the strict common set: {n_excl} of 7")

t.to_csv(R / "runs/track_a/tier_sensitivity.csv", index=False)
order_same = (list(t.sort_values("mean_full", ascending=False).model)
              == list(t.sort_values("mean_strict", ascending=False).model))
json.dump({"tier": len(ids), "reclassified": len(moved), "strictly_canonical": len(keep),
           "largest_mean_shift": round(float(t.delta.abs().max()), 4),
           "ordering_unchanged": bool(order_same),
           "leader": t.model.iloc[0], "trailing": t.model.iloc[-1],
           "common_set": len(common), "common_set_strict": len(ck),
           "paired_intervals_excluding_zero": int(n_excl)},
          open(R / "runs/track_a/tier_sensitivity.json", "w"), indent=1)
print(f"\nwrote {R}/runs/track_a/tier_sensitivity.json")
