from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
man = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
COLS = ["peptide_len", "is_cyclic", "cyclization_type"]
total = 0
for d in ("scores", "scores_ncaa", "scores_pep", "scores_pep_ncaa"):
    for f in sorted((R / "runs/track_a" / d).glob("*.parquet")):
        s = pd.read_parquet(f)
        if "complex_id" not in s.columns or not set(COLS) & set(s.columns):
            continue
        before = s.set_index("complex_id")
        idx = before.index
        changed = 0
        for c in COLS:
            if c in s.columns:
                new = man[c].reindex(idx)
                changed += int((before[c].astype(str) != new.astype(str)).sum())
                s[c] = new.values
        if changed:
            s.to_parquet(f, index=False)
            total += changed
            print(f"  {d}/{f.name}: {changed} cell(s) refreshed")
print(f"\n{total} stale metadata cells refreshed; DockQ values untouched")
