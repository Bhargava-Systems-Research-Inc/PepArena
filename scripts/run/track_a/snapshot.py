import json
import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench/runs/track_a")
FILES = ["leaderboard_v2_canonical.csv", "leaderboard_v2_ncaa.csv", "coverage_v2_canonical.csv",
         "paired_delta_v2_canonical.csv", "paired_holm_canonical.csv", "length_bins_v2_canonical.csv",
         "linear_cyclic_v2_canonical.csv", "length_multivariable.csv", "leaderboard.csv",
         "leaderboard_ncaa.csv", "common_set_composition_v2_canonical.csv"]

def snap():
    out = {}
    for f in FILES:
        p = R / f
        if not p.exists():
            continue
        d = pd.read_csv(p)
        key = d.columns[0]
        for _, row in d.iterrows():
            for c in d.columns[1:]:
                v = row[c]
                if isinstance(v, (int, float)) and pd.notna(v):
                    out[f"{f}:{row[key]}:{c}"] = round(float(v), 4)
    return out

if __name__ == "__main__":
    dest = Path(sys.argv[1])
    if len(sys.argv) > 2 and sys.argv[2] == "--diff":
        before = json.load(open(dest))
        after = snap()
        moved = {k: (before.get(k), v) for k, v in after.items()
                 if k in before and before[k] != v}
        added = sorted(set(after) - set(before))
        print(f"values that moved: {len(moved)} of {len(after)}")
        for k, (a, b) in sorted(moved.items())[:60]:
            print(f"  {k}: {a} -> {b}")
        if added:
            print(f"\nnew values: {len(added)}")
    else:
        json.dump(snap(), open(dest, "w"), indent=0)
        print(f"snapshotted {len(snap())} artifact values -> {dest}")
