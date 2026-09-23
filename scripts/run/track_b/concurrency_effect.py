import re
from pathlib import Path

import pandas as pd

F = Path("/12TBDrive1/mega_pep_bench/runs/track_b/fresh")
new = {}
for log in sorted(Path("/tmp").glob("adcp_shard*.log")):
    txt = log.read_text(errors="ignore")
    ev = sorted([(m.start(), "h", m.group(1)) for m in re.finditer(r"=== ADCP (\S+) ", txt)] +
                [(m.start(), "d", m.group(1))
                 for m in re.finditer(r"Docking performed in ([\d.]+) seconds", txt)])
    cur = None
    for _, k, v in ev:
        if k == "h":
            cur = v
        elif cur:
            new[cur] = float(v)
            cur = None
print(f"complexes timed in the sharded run: {len(new)}")

_base = F / "runtime_per_complex.preaudit.csv"
if not _base.exists():
    print("no preserved unsharded baseline; leaving the existing measurement in place")
    raise SystemExit(0)
old = pd.read_csv(_base)
old = old[old.program == "ADCP"].set_index("complex_id").seconds
shared = [c for c in new if c in old.index]
print(f"also timed in the unsharded run: {len(shared)}")
if shared:
    r = pd.DataFrame({"unsharded": [old[c] for c in shared],
                      "sharded": [new[c] for c in shared]}, index=shared)
    r["ratio"] = (r.sharded / r.unsharded).round(2)
    print(r.sort_values("ratio").to_string())
    print(f"\nmedian slowdown on identical complexes: {r.ratio.median():.2f}x")
import json
out = {"n_paired": len(shared),
       "median_slowdown": round(float(r.ratio.median()), 2) if shared else None,
       "min_slowdown": round(float(r.ratio.min()), 2) if shared else None,
       "max_slowdown": round(float(r.ratio.max()), 2) if shared else None,
       "sharded_median_seconds": round(float(pd.Series(new).median()), 1) if new else None,
       "unsharded_median_seconds": round(float(old.median()), 1),
       "note": ("wall-clock for the re-run includes contention: the worklists were sharded across "
                "concurrent workers so the work would finish in a day. Measured on complexes "
                "docked in both configurations, where the input is identical.")}
json.dump(out, open(F / "adcp_concurrency_effect.json", "w"), indent=1)
if shared:
    r.to_csv(F / "adcp_concurrency_effect.csv")
print(json.dumps(out, indent=1))
