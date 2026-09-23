from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
MAN = R / "data/curated/structure_all.parquet"

man = pd.read_parquet(MAN)
d = pd.read_csv(R / "runs/qc/manifest_rederived.csv")
man = man.drop(columns=[c for c in ("ncaa_codes", "cap_codes", "terminal_caps",
                                   "n_intra_peptide_bonds") if c in man],
                errors="ignore")
m = man.merge(d, on="complex_id", how="left")

was = lambda c: m[f"{c}_orig"] if f"{c}_orig" in m else m[c]
before = {"ncaa": int(was("has_ncaa").sum()),
          "len_median": float(was("peptide_len").median())}

for c in ("peptide_len", "peptide_seq", "has_ncaa"):
    if f"{c}_orig" not in m:
        m[f"{c}_orig"] = m[c]

got = m.peptide_len_new.notna()
m.loc[got, "peptide_len"] = m.loc[got, "peptide_len_new"].astype(int)
m.loc[got, "peptide_seq"] = m.loc[got, "peptide_seq_new"]
m.loc[got, "has_ncaa"] = m.loc[got, "has_ncaa_new"].astype(bool)
m["ncaa_codes"] = m["ncaa_codes"].fillna("")
m["terminal_caps"] = m["cap_codes"].fillna("")
m["n_intra_peptide_bonds"] = m["n_intra_bonds"].fillna(0).astype(int)
m = m.drop(columns=["peptide_len_new", "peptide_seq_new", "has_ncaa_new", "n_intra_bonds",
                    "cap_codes",
                    "native_path"], errors="ignore")

after = {"ncaa": int(m.has_ncaa.sum()), "len_median": float(m.peptide_len.median())}
m.to_parquet(MAN, index=False)

print(f"re-derived for {int(got.sum())} of {len(m)} complexes")
print(f"  non-canonical complexes: {before['ncaa']} -> {after['ncaa']}")
print(f"  median peptide length:   {before['len_median']:.0f} -> {after['len_median']:.0f}")
print(f"  sequence length now equals peptide_len for "
      f"{int((m.peptide_seq.str.len() == m.peptide_len).sum())} of {len(m)}")
print(f"\noriginals retained as peptide_len_orig / peptide_seq_orig / has_ncaa_orig")
