import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
MAN = R / "data/curated/structure_all.parquet"
BAK = R / "runs/qc/structure_all.before_recuration.parquet"
STRAT = R / "data/STRATIFICATION.md"
STRAT_BAK = R / "runs/qc/STRATIFICATION.before_recuration.md"

shutil.copy(MAN, BAK)
shutil.copy(STRAT, STRAT_BAK)
before = pd.read_parquet(BAK)
try:
    r = subprocess.run([sys.executable, str(R / "scripts/curate/combine_manifests.py")],
                       cwd=R, capture_output=True, text=True)
    print(r.stdout[-700:])
    if r.returncode != 0:
        print("curation failed:", r.stderr[-600:])
        raise SystemExit(1)
    after = pd.read_parquet(MAN)
    keys = ["peptide_len", "peptide_seq", "has_ncaa", "length_bin"]
    missing = [c for c in keys + ["peptide_len_orig", "terminal_caps", "ncaa_codes"]
               if c not in after.columns]
    print(f"\nrows {len(before)} -> {len(after)}")
    print(f"non-canonical {int(before.has_ncaa.sum())} -> {int(after.has_ncaa.sum())}")
    print(f"columns missing after re-curation: {missing or 'none'}")
    a = before.set_index("complex_id")[keys].sort_index()
    b = after.set_index("complex_id").reindex(a.index)[keys]
    diff = {c: int((a[c].astype(str) != b[c].astype(str)).sum()) for c in keys}
    print(f"per-column disagreement with the pre-curation manifest: {diff}")
    ok = not missing and all(v == 0 for v in diff.values()) and len(before) == len(after)
    print("\nGUARD HOLDS: a re-curation preserves the repair" if ok
          else "\nGUARD FAILS: a re-curation would change the released chemistry")
finally:
    shutil.copy(BAK, MAN)
    shutil.copy(STRAT_BAK, STRAT)
    print("released manifest and stratification restored from backup")
