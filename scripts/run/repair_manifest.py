import argparse
import json
from pathlib import Path

import gemmi
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
STD3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
CAPS = {"ACE", "NH2", "NME", "FOR", "TFA"}

def parent_letter(name):
    info = gemmi.find_tabulated_residue(name)
    if info and info.is_amino_acid():
        c = info.one_letter_code.upper()
        return c if c.isalpha() else None
    return None

def derive(row):
    p = Path(row.native_path)
    if not p.exists():
        return None
    st = gemmi.read_structure(str(p))
    st.setup_entities()
    ch = next((c for c in st[0] if c.name == row.peptide_chain), None)
    if ch is None:
        return None
    poly = {r.seqid.num for r in ch.get_polymer()}
    seq, ncaa, caps, n_poly = [], [], [], 0
    for res in ch:
        if res.name in CAPS:
            caps.append(res.name)
            continue
        if res.name == "HOH":
            continue
        letter = parent_letter(res.name)
        is_poly = res.seqid.num in poly or letter is not None
        if not is_poly:
            continue
        n_poly += 1
        seq.append(letter or "X")
        if res.name not in STD3:
            ncaa.append(res.name)
    doc = gemmi.cif.read(str(p)).sole_block()
    tab = doc.find("_struct_conn.", ["conn_type_id", "ptnr1_auth_asym_id", "ptnr2_auth_asym_id"])
    intra = sum(1 for r in tab if r[0] in ("covale", "disulf")
                and r[1] == row.peptide_chain and r[2] == row.peptide_chain)
    return {"complex_id": row.complex_id, "peptide_len_new": n_poly,
            "peptide_seq_new": "".join(seq), "has_ncaa_new": bool(ncaa),
            "ncaa_codes": ",".join(sorted(set(ncaa))),
            "cap_codes": ",".join(sorted(set(caps))), "n_intra_bonds": intra}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="change nothing; print the impact")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    man = pd.read_parquet(R / "data/curated/structure_all.parquet")
    idx = pd.read_parquet(R / "data/curated/structure_native_index.parquet")[
        ["complex_id", "native_path"]]
    man = man.merge(idx, on="complex_id", how="left").dropna(subset=["native_path"])
    if a.limit:
        man = man.head(a.limit)

    out, bad = [], 0
    for k, row in enumerate(man.itertuples(), 1):
        d = derive(row)
        if d is None:
            bad += 1
            continue
        out.append(d)
        if k % 1000 == 0:
            print(f"  {k}/{len(man)}", flush=True)
    d = pd.DataFrame(out)
    d.to_csv(R / "runs/qc/manifest_rederived.csv", index=False)

    m = man.merge(d, on="complex_id")
    for c in ("peptide_len", "peptide_seq", "has_ncaa"):
        if f"{c}_orig" in m:
            m[c] = m[f"{c}_orig"]
    len_moved = (m.peptide_len != m.peptide_len_new)
    ncaa_gain = (~m.has_ncaa) & m.has_ncaa_new
    ncaa_lose = m.has_ncaa & (~m.has_ncaa_new)
    bond_linear = (~m.is_cyclic) & (m.n_intra_bonds > 0)

    print(f"\nre-derived {len(d)} complexes ({bad} unreadable)\n")
    print(f"  peptide_len changes:            {int(len_moved.sum())} "
          f"({100*len_moved.mean():.1f}%), median shortfall "
          f"{(m.peptide_len_new - m.peptide_len)[len_moved].median():.0f}")
    print(f"  labelled canonical, carry ncAA: {int(ncaa_gain.sum())} "
          f"({100*ncaa_gain.mean():.1f}%)")
    print(f"  labelled ncAA, none found:      {int(ncaa_lose.sum())}")
    print(f"  flagged linear, intra bond:     {int(bond_linear.sum())}")

    tier = set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv")
               .iloc[:, 0].astype(str))
    t = m[m.complex_id.isin(tier)]
    print(f"\n  canonical Track A tier: {len(t)} complexes")
    print(f"    would become non-canonical: {int(((~t.has_ncaa) & t.has_ncaa_new).sum())}")
    print(f"    peptide_len would change:   {int((t.peptide_len != t.peptide_len_new).sum())}")
    json.dump({"n": len(d), "len_changed": int(len_moved.sum()),
               "canonical_gain_ncaa": int(ncaa_gain.sum()),
               "tier_would_move": int(((~t.has_ncaa) & t.has_ncaa_new).sum())},
              open(R / "runs/qc/manifest_repair_impact.json", "w"), indent=1)
    if not a.report:
        print("\n(--report not given: this run still only wrote the re-derived columns)")

if __name__ == "__main__":
    main()
