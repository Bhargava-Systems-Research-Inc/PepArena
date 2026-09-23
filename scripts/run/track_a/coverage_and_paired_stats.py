import pathlib
import numpy as np, pandas as pd

R = pathlib.Path("/12TBDrive1/mega_pep_bench")
KEYS = ["protenix", "af3", "boltz2", "chai1", "highfold3", "helixfold3", "rfaa", "esmfold2"]
LAB = {"protenix": "Protenix (v1.0.0)", "af3": "AlphaFold3", "boltz2": "Boltz-2",
       "chai1": "Chai-1", "highfold3": "HighFold3", "helixfold3": "HelixFold3",
       "rfaa": "RoseTTAFold-AA", "esmfold2": "ESMFold2 (ESMC-6B)"}
TIER = 334

sc = {k: pd.read_parquet(R / f"runs/track_a/scores/{k}.parquet") for k in KEYS}
idx = {k: v.set_index("complex_id") for k, v in sc.items()}
common = sorted(set.intersection(*(set(v.index) for v in idx.values())))

pd.DataFrame([{"model": LAB[k], "scored": len(sc[k]), "tier": TIER,
               "coverage": round(len(sc[k]) / TIER, 3), "missing": TIER - len(sc[k])}
              for k in KEYS]).sort_values("coverage", ascending=False
              ).to_csv(R / "runs/track_a/coverage.csv", index=False)

lk = pd.read_parquet(R / "data/splits/leakage.parquet").set_index("complex_id")
clu = lk.loc[common, "receptor_cluster"]
dq = pd.DataFrame({k: idx[k].loc[common, "DockQ"].values for k in KEYS}, index=common)
groups = {c: np.flatnonzero((clu == c).values) for c in clu.unique()}
gkeys = list(groups)
print(f"common set {len(common)} complexes in {len(gkeys)} receptor clusters "
      f"(largest {max(len(v) for v in groups.values())})")

rng = np.random.default_rng(0)
rows = []
for k in KEYS:
    if k == "protenix":
        continue
    d = (dq["protenix"] - dq[k]).values
    b1 = rng.choice(d, size=(10000, len(d)), replace=True).mean(axis=1)
    l1, h1 = np.percentile(b1, [2.5, 97.5])
    b2 = []
    for _ in range(10000):
        pick = rng.choice(len(gkeys), size=len(gkeys), replace=True)
        sel = np.concatenate([groups[gkeys[i]] for i in pick])
        b2.append(d[sel].mean())
    l2, h2 = np.percentile(b2, [2.5, 97.5])
    rows.append({"model": LAB[k], "n": len(d), "n_clusters": len(gkeys),
                 "mean_delta": round(float(d.mean()), 3),
                 "ci_low": round(float(l2), 3), "ci_high": round(float(h2), 3),
                 "excludes_zero": bool(l2 > 0 or h2 < 0),
                 "ci_low_complex": round(float(l1), 3), "ci_high_complex": round(float(h1), 3)})
out = pd.DataFrame(rows)
out.to_csv(R / "runs/track_a/paired_delta_vs_protenix.csv", index=False)
print("\n=== paired Delta vs Protenix (receptor-cluster bootstrap) ===")
print(out.to_string(index=False))
print(f"\nCI excludes zero (cluster): {int(out.excludes_zero.sum())} of {len(out)}")

s = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
key = s.loc[common].apply(lambda r: (r.receptor_pdb_id, r.peptide_seq), axis=1)
rep = key.reset_index().groupby(key.values)["complex_id"].first().tolist()
print(f"\n=== near-duplicate collapse: {len(common)} -> {len(rep)} representatives ===")
full = dq.mean().sort_values(ascending=False)
coll = dq.loc[rep].mean().sort_values(ascending=False)
print(f"  ranking identical: {list(full.index) == list(coll.index)}")
print(f"  max |mean shift| : {float((full - coll).abs().max()):.4f}")
for k in full.index:
    print(f"    {k:<12} {full[k]:.3f} -> {coll[k]:.3f}")
