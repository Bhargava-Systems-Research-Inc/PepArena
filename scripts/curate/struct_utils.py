from __future__ import annotations
import gzip
import gemmi

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

AMBER_ALIASES = {
    "CYX": "CYS", "CYM": "CYS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "ARN": "ARG", "TYM": "TYR",
}
CAP_RESIDUES = {"ACE", "NME", "NMA", "NHE", "NH2", "FOR"}

def aa_info(res):
    name = res.name
    if name in CAP_RESIDUES:
        return ("", False, True, name, False)
    canon = AMBER_ALIASES.get(name, name)
    info = gemmi.find_tabulated_residue(canon)
    is_aa = info.is_amino_acid() if info else False
    olc = info.one_letter_code.upper() if info else "X"
    if not olc.isalpha():
        olc = "X"
    is_std = (info.is_standard() if info else False) or (name in AMBER_ALIASES)
    return (olc, is_aa, False, canon, is_std)

def read_any(path):
    p = str(path)
    if p.endswith(".gz"):
        with gzip.open(p, "rt") as fh:
            text = fh.read()
        if ".cif" in p:
            doc = gemmi.cif.read_string(text)
            st = gemmi.make_structure_from_block(doc.sole_block())
        else:
            st = gemmi.read_pdb_string(text)
    elif p.endswith(".cif"):
        doc = gemmi.cif.read(p)
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_structure(p)
    st.setup_entities()
    return st

def chain_seq_and_ncaa(st, chain_id):
    model = st[0]
    chain = model.find_chain(chain_id)
    if chain is None:
        return "", [], []
    seq, ncaa, resnames = [], [], []
    for res in chain:
        olc, is_aa, is_cap, canon, is_std = aa_info(res)
        if is_cap or res.name in ("HOH", "WAT"):
            continue
        if not is_aa:
            continue
        if olc not in STANDARD_AA or not is_std:
            ncaa.append(res.name)
        seq.append(olc if olc in STANDARD_AA else "X")
        resnames.append(res.name)
    return "".join(seq), ncaa, resnames
