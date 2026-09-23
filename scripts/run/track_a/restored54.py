import json
import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix (v1.0.0)", "af3": "AlphaFold3", "boltz2": "Boltz-2",
       "chai1": "Chai-1", "highfold3": "HighFold3", "helixfold3": "HelixFold3",
       "rfaa": "RoseTTAFold-AA", "esmfold2": "ESMFold2 (ESMC-6B)"}

ids = sorted(set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv").iloc[:, 0].astype(str)))
old = sorted(set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv.bak334").iloc[:, 0].astype(str)))
restored = sorted(set(ids) - set(old))
lk = pd.read_parquet(R / "data/splits/leakage.parquet").set_index("complex_id")
sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
rng = np.random.default_rng(0)

def boot_diff(a, b, ca, cb):
    def draw(v, c):
        groups = [np.flatnonzero(c == g) for g in pd.unique(c)]
        pick = rng.integers(0, len(groups), size=(2000, len(groups)))
        return np.array([v[np.concatenate([groups[j] for j in row])].mean() for row in pick])
    return draw(a, ca) - draw(b, cb)

rows, out = [], {}
for k, lab in LAB.items():
    f = R / f"runs/track_a/scores/{k}.parquet"
    if not f.exists():
        continue
    d = pd.read_parquet(f).set_index("complex_id")
    d["cluster"] = lk.reindex(d.index)["receptor_cluster"].fillna("none")
    r = d[d.index.isin(restored)]
    o = d[d.index.isin(old)]
    if len(r) < 5 or len(o) < 5:
        continue
    diff = boot_diff(r.DockQ.to_numpy(), o.DockQ.to_numpy(),
                     r.cluster.to_numpy(), o.cluster.to_numpy())
    rows.append({"model": lab, "n_restored": len(r), "mean_restored": round(r.DockQ.mean(), 3),
                 "n_original": len(o), "mean_original": round(o.DockQ.mean(), 3),
                 "difference": round(r.DockQ.mean() - o.DockQ.mean(), 3),
                 "ci_low": round(float(np.percentile(diff, 2.5)), 3),
                 "ci_high": round(float(np.percentile(diff, 97.5)), 3)})
    rows[-1]["excludes_zero"] = bool(rows[-1]["ci_low"] > 0 or rows[-1]["ci_high"] < 0)

t = pd.DataFrame(rows)
t.to_csv(R / "runs/track_a/restored_54_comparison.csv", index=False)
lr = sp.reindex(restored).peptide_len.dropna()
lo = sp.reindex(old).peptide_len.dropna()
from scipy.stats import mannwhitneyu
mwu = float(mannwhitneyu(lr, lo).pvalue)
out = {"n_restored": len(restored), "n_original": len(old), "n_tier": len(ids),
       "per_model": rows, "peptide_len_MWU_P": round(mwu, 3),
       "restored_median_len": float(lr.median()), "original_median_len": float(lo.median()),
       "models_with_interval_excluding_zero": int(t.excludes_zero.sum()) if len(t) else 0}
json.dump(out, open(R / "runs/track_a/restored_54_comparison.json", "w"), indent=1)
print(t.to_string(index=False))
print(f"\npeptide length restored vs original: median {lr.median():.0f} vs {lo.median():.0f}, "
      f"Mann-Whitney P = {mwu:.3f}")
print(f"models whose restored-vs-original interval excludes zero: "
      f"{int(t.excludes_zero.sum()) if len(t) else 0} of {len(t)}")
