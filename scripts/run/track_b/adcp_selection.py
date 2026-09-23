import json
from math import comb
import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench/runs/track_b/scoring")
THRESH = {"acceptable": 0.23, "medium": 0.49, "high": 0.80}

d = pd.read_parquet(R / "dockq_all.parquet")
print("columns:", d.columns.tolist(), " rows:", len(d))
idcol = "pdb" if "pdb" in d.columns else d.columns[0]
qcol = "DockQ" if "DockQ" in d.columns else [c for c in d.columns if "dockq" in c.lower()][0]
rcol = next((c for c in d.columns if "rank" in c.lower() or c == "pose"), None)
print(f"using id={idcol} quality={qcol} rank={rcol}")

rows = []
for pdb, g in d.groupby(idcol):
    g = g.sort_values(rcol) if rcol else g
    q = g[qcol].to_numpy()
    n = len(q)
    rec = {"pdb": pdb, "n_poses": n, "best": q.max(), "top1": q[0],
           "mean_all": q.mean(), "top5_best": q[:5].max(), "top10_best": q[:10].max()}
    for name, t in THRESH.items():
        rec[f"hit_{name}_any"] = float((q >= t).any())
        rec[f"hit_{name}_top1"] = float(q[0] >= t)
        rec[f"hit_{name}_top5"] = float((q[:5] >= t).any())
        rec[f"hit_{name}_top10"] = float((q[:10] >= t).any())
        rec[f"hit_{name}_rand"] = float((q >= t).mean())
        miss = int((q < t).sum())
        rec[f"hit_{name}_rand5"] = 1.0 - (comb(miss, 5) / comb(n, 5) if miss >= 5 and n >= 5 else 0.0)
    rows.append(rec)

s = pd.DataFrame(rows)
out = R / "selection_quality.csv"
s.to_csv(out, index=False)

print(f"\n=== ADCP pose selection over {len(s)} cyclic complexes, {int(s.n_poses.median())} poses each")
print(f"{'threshold':12s} {'oracle':>7s} {'top-1':>7s} {'top-5':>7s} {'top-10':>7s} {'random-1':>9s} {'random-5':>9s}")
for name in THRESH:
    print(f"{name:12s} {s[f'hit_{name}_any'].mean():7.3f} {s[f'hit_{name}_top1'].mean():7.3f} "
          f"{s[f'hit_{name}_top5'].mean():7.3f} {s[f'hit_{name}_top10'].mean():7.3f} "
          f"{s[f'hit_{name}_rand'].mean():9.3f} {s[f'hit_{name}_rand5'].mean():9.3f}")

print(f"\nmean best available DockQ : {s.best.mean():.3f}")
print(f"mean top-1 DockQ          : {s.top1.mean():.3f}")
print(f"mean quality lost to ranking (best - top1): {(s.best - s.top1).mean():.3f}")
print(f"mean DockQ of a random pose: {s.mean_all.mean():.3f}")
print(f"headroom at acceptable: top-1 {s['hit_acceptable_top1'].mean():.3f} vs oracle "
      f"{s['hit_acceptable_any'].mean():.3f} = {(s['hit_acceptable_any'] - s['hit_acceptable_top1']).mean()*100:.1f} points")
print(f"\nwrote {out}")
