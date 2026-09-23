import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/12TBDrive1/mega_pep_bench")
sys.path.insert(0, str(ROOT / "scripts/curate"))
import schema

p = ROOT / "data/curated/structure_all.parquet"
df = pd.read_parquet(p)
old = df["length_bin"].copy()
df["length_bin"] = df["peptide_len"].fillna(0).astype(int).map(schema.length_bin)
moved = int((old != df["length_bin"]).sum())
df.to_parquet(p, index=False)
print(f"length_bin re-derived; {moved} complexes change bin")
print(pd.crosstab(old, df.length_bin).to_string())
