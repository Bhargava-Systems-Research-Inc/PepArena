import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "runs/track_c/boltz2_affinity_v2_batch"
tg = pd.read_csv(R / "runs/track_c/affinity_targets_v2.csv").set_index("complex_id")

rows = []
for f in sorted(OUT.rglob("affinity*.json")):
    cid = f.parent.name
    if cid not in tg.index:
        continue
    try:
        j = json.loads(f.read_text())
    except Exception:
        continue
    val = j.get("affinity_pred_value")
    if val is None:
        continue
    rows.append({"complex_id": cid, "pred": float(val),
                 "affinity_probability_binary": j.get("affinity_probability_binary"),
                 "log_affinity": float(tg.loc[cid, "log_affinity"]),
                 "heavy_atoms": int(tg.loc[cid, "heavy_atoms"]),
                 "set": tg.loc[cid, "set"]})

d = pd.DataFrame(rows)
d.to_csv(R / "runs/track_c/boltz2_affinity_v2_scored.csv", index=False)
n_all, n_128 = len(tg), int((tg.heavy_atoms <= 128).sum())
print(f"targets in the rebuilt set: {n_all}")
print(f"  within Boltz-2's 128-atom hard limit : {n_128}")
print(f"  within its 56-atom training domain   : {int((tg.heavy_atoms <= 56).sum())}")
print(f"  predictions returned                 : {len(d)}")
if len(d) >= 5:
    rho = spearmanr(d.log_affinity, d.pred)
    print(f"\nSpearman over the {len(d)} scored: {rho.statistic:.3f} (P = {rho.pvalue:.3g})")
    ind = d[d.heavy_atoms <= 56]
    if len(ind) >= 5:
        r2 = spearmanr(ind.log_affinity, ind.pred)
        print(f"  in-domain subset (<=56 atoms, n={len(ind)}): {r2.statistic:.3f}")
    for s, g in d.groupby("set"):
        if len(g) >= 5:
            print(f"  {s} subset (n={len(g)}): {spearmanr(g.log_affinity, g.pred).statistic:.3f}")
    summ = {"n_targets": n_all, "n_within_hard_limit": n_128,
            "n_within_training_domain": int((tg.heavy_atoms <= 56).sum()),
            "n_scored": len(d), "spearman": round(float(rho.statistic), 3),
            "coverage_of_rebuilt_set": round(len(d) / n_all, 3)}
    (R / "runs/track_c/boltz2_affinity_v2_leaderboard.json").write_text(json.dumps(summ, indent=1))
    print("\nwrote boltz2_affinity_v2_leaderboard.json")
