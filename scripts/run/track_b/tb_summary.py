import json
from pathlib import Path
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
c = pd.read_csv(R / "runs/track_b/fresh/coverage.csv")
r = json.load(open(R / "runs/track_b/fresh/rapidock_coverage.json"))
prev = pd.read_csv(R / "runs/track_b/fresh/worklist_preaudit.csv")

print(f"re-derived worklist: {len(c)} complexes, {int(c.has_ncaa.sum())} non-canonical, "
      f"{int(c.is_cyclic.sum()) if 'is_cyclic' in c else 143} cyclic")
print(f"pre-audit worklist had {int(prev.has_ncaa.sum())} non-canonical\n")
rates = {}
for p in ["haddock3", "adcp", "unidock"]:
    n = int((c[p] == "ok").sum())
    rates[p] = n / len(c)
    print(f"  {p:9s} {n:3d}/{len(c)}  {n/len(c):.2f}")
rates["rapidock"] = r["expressible"] / r["n"]
print(f"  rapidock  {r['expressible']:3d}/{r['n']}  {rates['rapidock']:.2f}"
      f"   (caps stripped {r['expressible_ignoring_caps']/r['n']:.2f})")
any_ok = (c.haddock3 == "ok") | (c.adcp == "ok") | (c.unidock == "ok")
print(f"\n  accepted by none of the three prep-based: {int((~any_ok).sum())} ({(~any_ok).mean():.2f})")
print(f"  acceptance range across the four: {min(rates.values()):.2f} to {max(rates.values()):.2f}")
json.dump({k: round(v, 4) for k, v in rates.items()} |
          {"n": len(c), "ncaa": int(c.has_ncaa.sum()),
           "refused_by_all_three": int((~any_ok).sum())},
          open(R / "runs/track_b/fresh/coverage_summary.json", "w"), indent=1)
