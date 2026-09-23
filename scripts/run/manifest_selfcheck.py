import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
m = pd.read_parquet(R / "data/curated/structure_all.parquet")
fails = []

def check(name, bad_mask):
    n = int(bad_mask.sum())
    print(f"  {'OK  ' if n == 0 else 'FAIL'}  {name}: {n}")
    if n:
        fails.append(f"{name} ({n})")
        print("        e.g. " + ", ".join(m.loc[bad_mask, "complex_id"].head(3)))

check("len(peptide_seq) != peptide_len", m.peptide_seq.fillna("").str.len() != m.peptide_len)
check("peptide_len <= 0", m.peptide_len <= 0)
check("has_ncaa but no ncaa_codes", m.has_ncaa & (m.ncaa_codes.fillna("") == ""))
check("ncaa_codes but not has_ncaa", (~m.has_ncaa) & (m.ncaa_codes.fillna("") != ""))
check("length_bin disagrees with peptide_len",
      pd.Series([lb != ("1-5" if n <= 5 else "6-10" if n <= 10 else "11-15" if n <= 15
                        else "16-25" if n <= 25 else "26-50" if n <= 50 else "51+")
                 for lb, n in zip(m.length_bin, m.peptide_len)]))
check("duplicate complex_id", m.complex_id.duplicated())
check("is_cyclic but cyclization_type linear", m.is_cyclic & (m.cyclization_type == "linear"))
check("cyclization_type set but not is_cyclic",
      (~m.is_cyclic) & (~m.cyclization_type.isin(["linear"])) & m.cyclization_type.notna())
check("negative intra-peptide bond count", m.n_intra_peptide_bonds < 0)
orig = m.peptide_len_orig.fillna(m.peptide_len)
print(f"\n  {len(m)} complexes; {int((m.peptide_len != orig).sum())} lengths differ from the "
      f"pre-repair value, {int(m.has_ncaa.sum())} non-canonical, "
      f"{int((m.terminal_caps.fillna('') != '').sum())} carry a terminal cap")
if fails:
    print("\nmanifest self-consistency FAILED: " + "; ".join(fails))
    sys.exit(1)
print("\nmanifest is internally consistent")
