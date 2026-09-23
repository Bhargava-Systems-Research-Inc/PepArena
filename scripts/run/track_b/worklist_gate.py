import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
F = R / "runs/track_b/fresh"
cov = pd.read_csv(F / "coverage.csv")

def run_set(tool):
    if tool == "rapidock":
        rd = pd.read_csv(F / "rapidock_coverage.csv")
        return set(rd.loc[rd.ok_ignoring_caps, "complex_id"].astype(str))
    ok = cov[tool] == "ok"
    if tool == "unidock":
        ok &= cov.torsdof <= 32
    return set(cov.loc[ok, "complex_id"].astype(str))

bad = []
for tool in ("haddock3", "adcp", "unidock", "rapidock"):
    f = F / f"{tool}_scored.csv"
    if not f.exists():
        bad.append(f"{tool}: no scored file"); continue
    scored = set(pd.read_csv(f).complex_id.astype(str))
    accepted = run_set(tool)
    stale = scored - accepted
    root = (F / tool / "out") if (F / tool / "out").is_dir() else (F / tool)
    attempted = {d.name for d in root.glob("*") if d.is_dir()} & accepted
    if stale:
        bad.append(f"{tool}: {len(stale)} scored complexes are not on the current worklist "
                   f"(e.g. {sorted(stale)[0]}) -- the run predates the re-derived worklist")
    elif attempted != accepted:
        bad.append(f"{tool}: {len(accepted) - len(attempted)} of {len(accepted)} accepted "
                   f"complexes were never attempted -- the run is unfinished")
    else:
        print(f"  OK    {tool}: {len(accepted)} attempted, {len(scored)} scored, all on the "
              f"current worklist")
if bad:
    print("\nTrack B accuracy does not describe the released worklist:")
    for b in bad:
        print("  " + b)
    sys.exit(1)
print("\nTrack B coverage and accuracy describe the same worklist.")
