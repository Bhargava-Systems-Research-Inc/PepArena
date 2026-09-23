import json
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
F = R / "runs/track_b/fresh"
C = R / "runs/track_b/control"

cov = pd.read_csv(F / "coverage.csv")
ccov = pd.read_csv(C / "coverage.csv")
print("=== acceptance, enriched worklist (183)")
for p in ("haddock3", "adcp", "rapidock", "unidock"):
    n = int((cov[p] == "ok").sum())
    print(f"  {p:10s} {n:3d}  {n/len(cov):.2f}")
anyok = ((cov.haddock3 == "ok") | (cov.adcp == "ok") | (cov.unidock == "ok")
         | (cov.rapidock == "ok"))
print(f"  refused by all four: {int((~anyok).sum())} ({100*(~anyok).mean():.0f}%)")
u = cov[cov.unidock == "ok"]
print(f"  Uni-Dock caps: " + ", ".join(f"<={c}: {int((u.torsdof <= c).sum())}" for c in (15, 20, 25, 32)))
print(f"  median torsions  ncAA {cov.loc[cov.has_ncaa, 'torsdof'].median():.0f}  "
      f"canonical {cov.loc[~cov.has_ncaa, 'torsdof'].median():.0f}")

print("\n=== acceptance, canonical control (33)")
for p in ("haddock3", "adcp", "unidock"):
    n = int((ccov[p] == "ok").sum())
    s = ccov[~ccov.has_ncaa]
    print(f"  {p:10s} {n:3d}  {n/len(ccov):.2f}   strictly canonical {int((s[p]=='ok').sum())}/"
          f"{len(s)} = {(s[p]=='ok').mean():.2f}")
rj = C / "rapidock_coverage.json"
if rj.exists():
    d = json.loads(rj.read_text())
    print(f"  rapidock   {d['expressible']:3d}  {d['expressible']/d['n']:.2f}   "
          f"caps stripped {d['expressible_ignoring_caps']} = "
          f"{d['expressible_ignoring_caps']/d['n']:.2f}")

print("\n=== accuracy, enriched worklist")
for t in ("haddock3", "adcp", "unidock", "rapidock"):
    f = F / f"{t}_leaderboard.json"
    if not f.exists():
        print(f"  {t}: no leaderboard yet"); continue
    d = json.loads(f.read_text())
    print(f"  {t:10s} n={d['n_complexes']:3d}  best {d['mean_best_dockq']:.3f}  "
          f"top1 {d['mean_top1_dockq']:.3f}  sampling {d['sampling_power_acceptable']:.3f}  "
          f"docking {d['docking_power_top1']:.3f}")
ca = C / "adcp_leaderboard.json"
if ca.exists():
    d = json.loads(ca.read_text())
    print(f"  adcp control n={d['n_complexes']}  best {d['mean_best_dockq']:.3f}  "
          f"top1 {d['mean_top1_dockq']:.3f}  sampling {d['sampling_power_acceptable']:.3f}  "
          f"docking {d['docking_power_top1']:.3f}")

rt = F / "runtime_per_complex.csv"
if rt.exists():
    m = pd.read_csv(rt).groupby("program").seconds.median().round(1)
    print("\n=== median wall-clock per complex\n" + m.to_string())

ts = F / "unidock_torsion_sensitivity.csv"
if ts.exists():
    print("\n=== Uni-Dock torsion-cap sensitivity")
    print(pd.read_csv(ts).to_string(index=False))
