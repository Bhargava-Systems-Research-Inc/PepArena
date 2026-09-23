import json, csv, sys
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr, pearsonr

REPO = Path("/12TBDrive1/mega_pep_bench")
OUT = REPO / "runs/track_c/boltz2_affinity"

rows = []
globs = ["out_final/**/affinity_*.json", "out/*/boltz_results_*/predictions/*/affinity*.json"]
files = [f for g in globs for f in OUT.glob(g)]
if any("out_final" in str(f) for f in files):
    files = [f for f in files if "out_final" in str(f)]
for f in files:
    cid = f.parent.name
    try:
        d = json.loads(f.read_text())
    except Exception:
        continue
    v = d.get("affinity_pred_value")
    rows.append({"complex_id": cid, "shard": f.parts[-5],
                 "affinity_pred_value": v,
                 "affinity_probability_binary": d.get("affinity_probability_binary")})
p = pd.DataFrame(rows).drop_duplicates("complex_id")
print(f"collected {len(p)} affinity predictions")
if p.empty:
    sys.exit(0)
n_nan = p.affinity_pred_value.isna().sum() + (p.affinity_pred_value != p.affinity_pred_value).sum()
print(f"  NaN/missing values: {int(p.affinity_pred_value.isna().sum())}")

t = pd.read_csv(REPO / "runs/track_c/affinity_targets.csv")
j = p.merge(t, on="complex_id", how="inner").dropna(subset=["affinity_pred_value", "log_affinity"])
print(f"joined with measured affinity: {len(j)} / {len(t)} targets")
j.to_csv(OUT / "boltz2_affinity_vs_measured.csv", index=False)

if len(j) > 5:
    rho, pr = spearmanr(j.affinity_pred_value, j.log_affinity)
    r, pp = pearsonr(j.affinity_pred_value, j.log_affinity)
    summary = {"n": int(len(j)), "spearman": round(float(rho), 3), "spearman_p": float(pr),
               "pearson": round(float(r), 3), "pearson_p": float(pp),
               "n_over_56_heavy_atoms": int((j.heavy > 56).sum()),
               "note": ("Boltz-2's affinity head was trained on ligands up to 56 heavy atoms; "
                        "peptides here have a median of 77, so nearly every target is out of "
                        "domain and boltz emits that warning per target.")}
    print(f"\nBoltz-2 affinity vs log10 Kd:  Spearman {rho:+.3f} (p={pr:.2g})   Pearson {r:+.3f}")
    print("\nby source dataset:")
    for src, g in j.groupby("source_dataset"):
        if len(g) >= 5:
            rr, _ = spearmanr(g.affinity_pred_value, g.log_affinity)
            print(f"  {src:18s} n={len(g):5d}  rho={rr:+.3f}")
    print("\nby cyclic status:")
    for cy, g in j.groupby("is_cyclic"):
        if len(g) >= 5:
            rr, _ = spearmanr(g.affinity_pred_value, g.log_affinity)
            print(f"  cyclic={str(cy):5s} n={len(g):5d}  rho={rr:+.3f}")
    wr = []
    for rec, g in j.groupby("receptor_pdb_id"):
        if len(g) >= 5 and g.log_affinity.nunique() > 2:
            rr, _ = spearmanr(g.affinity_pred_value, g.log_affinity)
            if rr == rr:
                wr.append((rec, len(g), rr))
    if wr:
        import statistics
        med = statistics.median(r for _, _, r in wr)
        summary["within_receptor_median_spearman"] = round(float(med), 3)
        summary["within_receptor_n_series"] = len(wr)
        print(f"\nwithin-receptor ranking ({len(wr)} series with n>=5): "
              f"median rho {med:+.3f}")
        for rec, n, rr in sorted(wr, key=lambda x: -x[1])[:8]:
            print(f"  {rec:8s} n={n:4d}  rho={rr:+.3f}")

    print("\nby peptide size (heavy atoms, 56 = Boltz-2 affinity training ceiling):")
    j["bin"] = pd.cut(j.heavy, [0, 56, 80, 100, 128])
    for b, g in j.groupby("bin", observed=True):
        if len(g) >= 5:
            rr, _ = spearmanr(g.affinity_pred_value, g.log_affinity)
            print(f"  {str(b):12s} n={len(g):5d}  rho={rr:+.3f}")
    json.dump(summary, open(OUT / "boltz2_affinity_leaderboard.json", "w"), indent=2)
    print(f"\nwrote {OUT/'boltz2_affinity_leaderboard.json'}")
