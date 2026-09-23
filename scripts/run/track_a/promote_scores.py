import shutil
import sys
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench/runs/track_a")
SPINE = pd.read_parquet("/12TBDrive1/mega_pep_bench/data/curated/structure_all.parquet")
SPINE = SPINE.set_index("complex_id")[["is_cyclic", "cyclization_type", "peptide_len"]]

for suffix in ("", "_ncaa"):
    src, dst = R / f"scores_pep{suffix}", R / f"scores{suffix}"
    if not src.is_dir():
        print(f"missing {src}"); sys.exit(1)
    arch = R / f"_superseded_globaldockq{suffix}"
    if not arch.exists():
        shutil.copytree(dst, arch)
    n = 0
    for f in sorted(src.glob("*.parquet")):
        d = pd.read_parquet(f)
        if d.empty or "complex_id" not in d.columns:
            print(f"  {f.name}: empty, skipped (nothing to promote)")
            continue
        d = d.join(SPINE, on="complex_id")
        missing = d.peptide_len.isna().sum()
        d.to_parquet(dst / f.name, index=False)
        n += 1
        print(f"  {f.name}: {len(d)} rows, mean DockQ {d.DockQ.mean():.3f}"
              + (f"  ({missing} without stratification metadata)" if missing else ""))
    print(f"promoted {n} files -> {dst} (previous archived in {arch.name})")
