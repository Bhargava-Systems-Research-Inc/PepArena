import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
IDENT, COV = 0.3, 0.8

sp = pd.read_parquet(R / "data/splits/interaction_splits.parquet")
man = pd.read_parquet(R / "data/curated/interaction_tpeppro.parquet")[
    ["complex_id", "receptor_seq"]]
d = man.merge(sp, on="complex_id")

sep = pd.read_csv(R / "runs/track_d/receptor_separation_split_receptor.csv")
sep["mincov"] = sep[["qcov", "tcov"]].min(axis=1)
leaky = set(sep.loc[(sep.fident > IDENT) & (sep.mincov > COV), "test_seq"])

col = "split_receptor"
test_recs = set(d.loc[d[col] == "test", "receptor_seq"])
violate = test_recs & leaky
print(f"receptor-cluster split: {len(test_recs)} test receptors, "
      f"{len(violate)} violate the criterion pairwise ({100*len(violate)/len(test_recs):.1f}%)")

strict = d[col].copy()
strict[d.receptor_seq.isin(violate) & (d[col] == "test")] = "train"
d["split_receptor_strict"] = strict

n_tr = int((strict == "train").sum())
n_te = int((strict == "test").sum())
kept = len(test_recs - violate)
print(f"strict split: {n_tr} train / {n_te} test pairs over {kept} test receptors")
print(f"  pairs moved out of test: {int(((d[col]=='test') & (strict=='train')).sum())}")

out = sp.merge(d[["complex_id", "split_receptor_strict"]], on="complex_id", how="left")
out["split_receptor_strict"] = out["split_receptor_strict"].fillna("none")
out.to_parquet(R / "data/splits/interaction_splits.parquet", index=False)
print(f"\nwrote split_receptor_strict into {R/'data/splits/interaction_splits.parquet'}")

still = set(out.loc[out.split_receptor_strict == "test", "complex_id"])
rem = d[d.complex_id.isin(still)].receptor_seq.unique()
print(f"verification: {sum(1 for r in rem if r in leaky)} of {len(rem)} remaining test receptors "
      f"violate the criterion (must be 0)")
