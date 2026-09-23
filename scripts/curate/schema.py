from __future__ import annotations
import pandas as pd

CYCLIZATION_TYPES = [
    "linear", "head_to_tail", "disulfide", "lactam", "ester",
    "thioether", "staple", "bicyclic", "nmethyl", "other", "unknown",
]
EXP_METHODS = ["xray", "nmr", "em", "predicted", "unknown"]
LENGTH_BINS = [(1, 5), (6, 10), (11, 15), (16, 25), (26, 50), (51, 10**6)]

SPINE = [
    "complex_id",
    "source_dataset",
    "receptor_pdb_id",
    "receptor_chains",
    "peptide_chain",
    "peptide_seq",
    "peptide_len",
    "is_cyclic",
    "cyclization_type",
    "cyclization_bonds",
    "has_ncaa",
    "ncaa_list",
    "peptide_smiles",
    "peptide_ccd_codes",
    "target_class",
    "sec_struct",
    "n_interface_res",
    "msa_neff_over_l",
    "deposition_date",
    "resolution",
    "exp_method",
    "notes",
]

STRUCTURE_COLS = SPINE + [
    "receptor_path",
    "peptide_path",
    "native_complex_path",
    "receptor_fasta",
]

AFFINITY_COLS = SPINE + [
    "measured_value",
    "value_type",
    "value_units",
    "log_affinity",
    "assay",
    "native_complex_path",
]

INTERACTION_COLS = SPINE + [
    "receptor_seq",
    "label",
    "binding_residues",
    "split_origin",
]

def empty_frame(kind: str) -> pd.DataFrame:
    cols = {"structure": STRUCTURE_COLS, "affinity": AFFINITY_COLS,
            "interaction": INTERACTION_COLS}[kind]
    return pd.DataFrame(columns=cols)

def length_bin(n: int) -> str:
    for lo, hi in LENGTH_BINS:
        if lo <= n <= hi:
            return f"{lo}-{hi}" if hi < 10**6 else f"{lo}+"
    return "na"

if __name__ == "__main__":
    for k in ("structure", "affinity", "interaction"):
        assert list(empty_frame(k).columns)[:3] == ["complex_id", "source_dataset", "receptor_pdb_id"]
    assert length_bin(3) == "1-5" and length_bin(12) == "11-15" and length_bin(99) == "51+"
    assert "head_to_tail" in CYCLIZATION_TYPES
    print("schema self-check ok")
