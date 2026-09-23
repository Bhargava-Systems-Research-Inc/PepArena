import argparse
import json
import gemmi
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
CAPS = {"ACE", "NH2", "NME", "FOR", "TFA"}

def one_letter(name):
    info = gemmi.find_tabulated_residue(name)
    return info.one_letter_code.upper() if info and info.is_amino_acid() else None

def audit(row):
    d = {"complex_id": row.complex_id}
    p = Path(row.native_path)
    if not p.exists():
        return {**d, "missing_file": True}
    st = gemmi.read_structure(str(p))
    st.setup_entities()
    doc = gemmi.cif.read(str(p))
    model = st[0]
    chains = {c.name: c for c in model}

    pep = chains.get(row.peptide_chain)
    if pep is None:
        return {**d, "peptide_chain_absent": True}

    poly = {r.seqid.num for r in pep.get_polymer()}
    aa, caps, ncaa, hets = [], [], [], []
    for res in pep:
        code = one_letter(res.name)
        if code and res.name not in CAPS:
            aa.append(code)
        elif res.name in CAPS:
            caps.append(res.name)
        elif res.seqid.num in poly:
            ncaa.append(res.name)
        elif res.name != "HOH":
            hets.append(res.name)
    d["ncaa_in_peptide"] = ",".join(sorted(set(ncaa)))
    d["ncaa_unrecorded"] = bool(ncaa) and not bool(row.has_ncaa)
    d["n_polymer"] = len(aa) + len(ncaa)
    d["len_undercount"] = d["n_polymer"] - int(row.peptide_len)

    d["n_aa"] = len(aa)
    d["len_mismatch"] = (len(aa) + len(ncaa)) != int(row.peptide_len)
    d["seq_mismatch"] = len(aa) == int(row.peptide_len) and "".join(aa) != row.peptide_seq
    recorded = getattr(row, "terminal_caps", None)
    d["cap_unrecorded"] = bool(caps) and not bool(
        recorded if recorded is not None else row.has_ncaa)
    d["caps"] = ",".join(sorted(set(caps)))
    d["nonpolymer_in_peptide_chain"] = ",".join(sorted(set(hets)))
    want = [c for c in str(row.receptor_chains).split(",") if c]
    d["receptor_chain_absent"] = ",".join(c for c in want if c not in chains)

    block = doc.sole_block()
    tab = block.find("_struct_conn.", ["conn_type_id", "ptnr1_label_asym_id", "ptnr2_label_asym_id",
                                       "ptnr1_auth_asym_id", "ptnr2_auth_asym_id",
                                       "ptnr1_label_comp_id", "ptnr2_label_comp_id"])
    intra = 0
    for r_ in tab:
        if r_[0] not in ("covale", "disulf"):
            continue
        if r_[3] == row.peptide_chain and r_[4] == row.peptide_chain:
            intra += 1
    d["n_intrapeptide_bonds"] = intra
    d["intra_bond_unrecorded"] = intra > 0 and int(getattr(row, "n_intra_peptide_bonds", 0)) != intra
    d["cyclic_flag_vs_bonds"] = ("bond_but_linear" if intra and not row.is_cyclic else
                                 "cyclic_but_no_bond" if row.is_cyclic and not intra else "")
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="complexes per stratum ceiling")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--orig", action="store_true",
                    help="audit the pre-repair annotation retained as *_orig, so the "
                         "repair's effect is a regenerable number, not a remembered one")
    a = ap.parse_args()

    man = pd.read_parquet(R / "data/curated/structure_all.parquet")
    sfx = "_prerepair" if a.orig else ""
    if a.orig:
        for c in ("peptide_len", "peptide_seq", "has_ncaa"):
            man[c] = man[c + "_orig"]
        man = man.drop(columns=["n_intra_peptide_bonds", "terminal_caps"],
                       errors="ignore")
    man = man.merge(pd.read_parquet(R / "data/curated/structure_native_index.parquet")
                    [["complex_id", "native_path"]], on="complex_id", how="left")
    need = ["complex_id", "peptide_chain", "receptor_chains", "peptide_seq", "peptide_len",
            "is_cyclic", "has_ncaa", "native_path", "cyclization_type",
            "n_intra_peptide_bonds", "terminal_caps"]
    man = man[[c for c in need if c in man.columns]].dropna(subset=["native_path"])

    man["stratum"] = (man.cyclization_type.fillna("linear").astype(str)
                      + "|" + man.has_ncaa.astype(str))
    sample = (man.groupby("stratum", group_keys=False)
                 .apply(lambda g: g.sample(min(len(g), a.n), random_state=a.seed),
                        include_groups=True))
    print(f"auditing {len(sample)} complexes over {sample.stratum.nunique()} strata")

    rows = [audit(r) for r in sample.itertuples()]
    out = pd.DataFrame(rows).merge(sample[["complex_id", "stratum"]], on="complex_id", how="left")
    share = man.stratum.value_counts(normalize=True)
    drawn = out.stratum.value_counts()
    out["w"] = out.stratum.map(share / drawn)
    out.to_csv(R / f"runs/qc/release_audit{sfx}.csv", index=False)

    n = len(out)
    def col(name, default=False):
        return (out[name] if name in out else pd.Series(default, index=out.index)).fillna(default)
    masks = {
        "file missing": col("missing_file"),
        "peptide chain absent": col("peptide_chain_absent"),
        "manifest length != polymer residues present": col("len_mismatch"),
        "intra-peptide bond count not recorded": col("intra_bond_unrecorded"),
        "unmappable residue dropped from length and sequence": col("len_undercount", 0) > 0,
        "manifest sequence != structure": col("seq_mismatch"),
        "terminal cap present, not recorded": col("cap_unrecorded"),
        "free ion or ligand in peptide chain": col("nonpolymer_in_peptide_chain", "") != "",
        "non-canonical residue present, not recorded": col("ncaa_unrecorded"),
        "named receptor chain absent": col("receptor_chain_absent", "") != "",
        "intra-peptide bond present, flagged linear (adjudication, not error)":
            col("cyclic_flag_vs_bonds", "") == "bond_but_linear",
        "flagged cyclic, no deposited bond": col("cyclic_flag_vs_bonds", "") == "cyclic_but_no_bond",
    }
    flags = {k: (int(m.sum()), m) for k, m in masks.items()}
    def rate(mask):
        return 100 * float((out.w * mask.astype(float)).sum() / out.w.sum())
    print(f"\n{'in sample':>17s}  {'release-weighted':>16s}   check")
    for k, (v, mask) in flags.items():
        print(f"  {int(v):4d} / {n} ({100*int(v)/n:4.1f}%)  {rate(mask):15.1f}%   {k}")
    clean = out[[c for c in ["len_mismatch", "seq_mismatch", "cap_unrecorded", "ncaa_unrecorded"]
                 if c in out]].fillna(False).any(axis=1)
    clean = clean | (out.len_undercount.fillna(0) > 0)
    clean = clean | (out.nonpolymer_in_peptide_chain.fillna("") != "") | (out.receptor_chain_absent.fillna("") != "")
    u = out[out.len_undercount.fillna(0) > 0].len_undercount
    if len(u):
        print(f"\n  where a residue was dropped, the shortfall is {u.median():.0f} residues "
              f"(median), up to {u.max():.0f}")
    print(f"\n  complexes with no discrepancy on any check: {int((~clean).sum())} / {n} "
          f"({100*int((~clean).sum())/n:.1f}%) in sample, {100 - rate(clean):.1f}% release-weighted")
    json.dump({"n_audited": n, "n_clean": int((~clean).sum()),
               "in_sample": {k: v for k, (v, _) in flags.items()},
               "release_weighted_pct": {k: round(rate(m), 2) for k, (_, m) in flags.items()},
               "clean_release_weighted_pct": round(100 - rate(clean), 2)},
              open(R / f"runs/qc/release_audit{sfx}.json", "w"), indent=2)
    print(f"\nwrote {R}/runs/qc/release_audit{sfx}.csv")

if __name__ == "__main__":
    main()
