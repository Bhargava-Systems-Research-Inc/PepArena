import argparse
import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix (v1.0.0)", "af3": "AlphaFold3", "boltz2": "Boltz-2",
       "chai1": "Chai-1", "highfold3": "HighFold3", "helixfold3": "HelixFold3",
       "rfaa": "RoseTTAFold-AA", "esmfold2": "ESMFold2 (ESMC-6B)"}
CANON = ["protenix", "af3", "boltz2", "chai1", "highfold3", "helixfold3", "rfaa", "esmfold2"]
NCAA = ["protenix", "af3", "boltz2", "highfold3"]
NOT_SUBMITTED = {"canonical": {"af3": "cyclic"}}

NOT_ATTEMPTED_RESTORED = {"canonical": {"rfaa", "esmfold2"}}

def tier_ids(tier):
    if tier == "canonical":
        f = R / "runs/track_a/esmfold2_inputs_canonical.csv"
        ids = pd.read_csv(f).iloc[:, 0].astype(str).tolist()
    else:
        ids = pd.read_csv(R / "ncaa_ccd/manifest.csv").complex_id.astype(str).tolist()
    return sorted(set(ids))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="scores_pep")
    ap.add_argument("--tier", default="canonical", choices=["canonical", "ncaa"])
    a = ap.parse_args()
    models = CANON if a.tier == "canonical" else NCAA
    sdir = R / "runs/track_a" / (a.scores if a.tier == "canonical" else a.scores + "_ncaa")
    ids = tier_ids(a.tier)
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    cyclic = {c for c in ids if bool(sp.is_cyclic.get(c, False))}
    print(f"frozen {a.tier} tier: {len(ids)} complexes ({len(cyclic)} cyclic) from {sdir.name}")

    sc = {}
    for m in models:
        f = sdir / f"{m}.parquet"
        if f.exists():
            d = pd.read_parquet(f)
            sc[m] = d[d.complex_id.isin(ids)].set_index("complex_id")
    if not sc:
        raise SystemExit(f"no score files in {sdir}")

    _bak = R / "runs/track_a/esmfold2_inputs_canonical.csv.bak334"
    restored = (set(ids) - set(pd.read_csv(_bak).iloc[:, 0].astype(str))) if _bak.exists() else set()

    rows = []
    for m, d in sc.items():
        got = set(d.index)
        missing = set(ids) - got
        ns = NOT_SUBMITTED.get(a.tier, {}).get(m)
        not_sub = (missing & cyclic) if ns == "cyclic" else set()
        no_attempt = (missing & restored) if m in NOT_ATTEMPTED_RESTORED.get(a.tier, set()) else set()
        rows.append({"model": LAB.get(m, m), "tier": len(ids), "scored": len(got),
                     "coverage": round(len(got) / len(ids), 3),
                     "missing": len(missing), "not_submitted": len(not_sub),
                     "no_usable_attempt": len(no_attempt),
                     "submitted_but_unreturned": len(missing) - len(not_sub) - len(no_attempt),
                     "zero_scored_no_interface": int((d.get("n_interfaces", pd.Series(dtype=int)) == 0).sum()),
                     "route_merged": int((d.get("route", pd.Series(dtype=object)) == "merged").sum()),
                     "route_per_interface": int((d.get("route", pd.Series(dtype=object)) == "per_interface").sum())})
    cov = pd.DataFrame(rows).sort_values("coverage", ascending=False)
    cov.to_csv(R / f"runs/track_a/coverage_v2_{a.tier}.csv", index=False)
    print("\n=== coverage against the frozen tier")
    print(cov.to_string(index=False))

    common = sorted(set.intersection(*(set(d.index) for d in sc.values())))
    print(f"\ncommon set: {len(common)} complexes returned by all {len(sc)} models")
    lb = []
    for m, d in sc.items():
        q = d.DockQ
        cy = d.index.isin(list(cyclic))
        lb.append({"model": LAB.get(m, m), "n": len(q), "mean_dockq": round(q.mean(), 3),
                   "median": round(q.median(), 3),
                   "acc_0p23": round((q >= 0.23).mean(), 3),
                   "med_0p49": round((q >= 0.49).mean(), 3),
                   "high_0p80": round((q >= 0.80).mean(), 3),
                   "linear_mean": round(q[~cy].mean(), 3) if (~cy).any() else None,
                   "cyclic_mean": round(q[cy].mean(), 3) if cy.any() else None,
                   "cyclic_n": int(cy.sum()),
                   "common_mean": round(d.loc[common, "DockQ"].mean(), 3),
                   "mean_dockq_imputed": round(q.sum() / len(ids), 3),
                   "acc_0p23_imputed": round((q >= 0.23).sum() / len(ids), 3)})
    lb = pd.DataFrame(lb).sort_values("mean_dockq", ascending=False)
    lb.to_csv(R / f"runs/track_a/leaderboard_v2_{a.tier}.csv", index=False)
    print("\n=== leaderboard (conditional on returning, and unconditional over the frozen tier)")
    print(lb.to_string(index=False))

    ref = lb.iloc[0].model
    refkey = [k for k, v in LAB.items() if v == ref][0]
    lk = pd.read_parquet(R / "data/splits/leakage.parquet").set_index("complex_id")
    clu = lk.reindex(common)["receptor_cluster"]
    dq = pd.DataFrame({m: sc[m].loc[common, "DockQ"].values for m in sc}, index=common)
    groups = {c: np.flatnonzero((clu == c).values) for c in clu.dropna().unique()}
    gk = list(groups)
    rng = np.random.default_rng(0)
    rows = []
    for m in sc:
        if m == refkey:
            continue
        d = (dq[refkey] - dq[m]).values
        b = [d[np.concatenate([groups[gk[i]] for i in rng.choice(len(gk), len(gk), True)])].mean()
             for _ in range(2000)]
        lo, hi = np.percentile(b, [2.5, 97.5])
        pair = sorted(set(sc[refkey].index) & set(sc[m].index))
        dp = (sc[refkey].loc[pair, "DockQ"] - sc[m].loc[pair, "DockQ"]).values
        full = np.array([sc[refkey].DockQ.get(c, 0.0) - sc[m].DockQ.get(c, 0.0) for c in ids])
        rows.append({"model": LAB.get(m, m), "n_common": len(d), "n_clusters": len(gk),
                     "delta_common": round(float(d.mean()), 3),
                     "ci_low": round(float(lo), 3), "ci_high": round(float(hi), 3),
                     "excludes_zero": bool(lo > 0 or hi < 0),
                     "n_pairwise": len(dp), "delta_pairwise": round(float(dp.mean()), 3),
                     "n_imputed": len(ids), "delta_imputed": round(float(full.mean()), 3)})
    pv = pd.DataFrame(rows)
    pv.to_csv(R / f"runs/track_a/paired_delta_v2_{a.tier}.csv", index=False)
    print(f"\n=== paired DockQ difference against {ref} (receptor-cluster bootstrap)")
    print(pv.to_string(index=False))

    rest = sorted(set(ids) - set(common))
    comp = pd.DataFrame([
        {"subset": "common set", "n": len(common),
         "cyclic_pct": round(100 * np.mean([c in cyclic for c in common]), 1),
         "median_pep_len": float(sp.reindex(common).peptide_len.median())},
        {"subset": "returned by some or none", "n": len(rest),
         "cyclic_pct": round(100 * np.mean([c in cyclic for c in rest]), 1) if rest else None,
         "median_pep_len": float(sp.reindex(rest).peptide_len.median()) if rest else None}])
    comp.to_csv(R / f"runs/track_a/common_set_composition_v2_{a.tier}.csv", index=False)
    print("\n=== common-set composition (every tier complex is in exactly one row)")
    print(comp.to_string(index=False))

if __name__ == "__main__":
    main()
