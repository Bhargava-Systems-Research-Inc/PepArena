import json
import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
new = pd.read_parquet(R / "runs/track_a/scores_ncaa/highfold3.parquet").set_index("complex_id")
old = pd.read_parquet(R / "runs/track_a/scores_ncaa/highfold3_nocyc.parquet").set_index("complex_id")
af3 = pd.read_parquet(R / "runs/track_a/scores_ncaa/af3.parquet").set_index("complex_id")
sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")

cyc = [c for c in set(new.index) & set(old.index) if bool(sp.is_cyclic.get(c, False))]
cyc = sorted(cyc)
d = (new.loc[cyc, "DockQ"] - old.loc[cyc, "DockQ"]).to_numpy()
rng = np.random.default_rng(0)
bs = rng.choice(d, size=(10000, len(d)), replace=True).mean(axis=1)
lo, hi = np.percentile(bs, [2.5, 97.5])

out = {
    "n_paired_cyclic": len(cyc),
    "mean_with_cycpoem": round(float(new.loc[cyc, "DockQ"].mean()), 3),
    "mean_without_cycpoem": round(float(old.loc[cyc, "DockQ"].mean()), 3),
    "mean_paired_delta": round(float(d.mean()), 3),
    "ci_low": round(float(lo), 3),
    "ci_high": round(float(hi), 3),
    "excludes_zero": bool(lo > 0 or hi < 0),
    "improved_on": int((d > 0).sum()),
    "worsened_on": int((d < 0).sum()),
    "unchanged_on": int((d == 0).sum()),
}
both = [c for c in cyc if c in af3.index]
out["n_vs_af3"] = len(both)
out["af3_mean_same_complexes"] = round(float(af3.loc[both, "DockQ"].mean()), 3)
out["highfold3_mean_same_complexes"] = round(float(new.loc[both, "DockQ"].mean()), 3)
dd = (new.loc[both, "DockQ"] - af3.loc[both, "DockQ"]).to_numpy()
bs2 = rng.choice(dd, size=(10000, len(dd)), replace=True).mean(axis=1)
l2, h2 = np.percentile(bs2, [2.5, 97.5])
out["delta_vs_af3"] = round(float(dd.mean()), 3)
out["delta_vs_af3_ci"] = [round(float(l2), 3), round(float(h2), 3)]

(R / "runs/track_a/cycpoem_effect.json").write_text(json.dumps(out, indent=1))
for k, v in out.items():
    print(f"  {k}: {v}")
