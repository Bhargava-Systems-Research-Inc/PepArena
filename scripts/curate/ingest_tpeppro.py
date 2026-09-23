from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema

SRC = ROOT / "data/external/TPepPro/data"
ACTIONS = SRC / "receptor-peptide.actions.tsv"
DICT = SRC / "receptor(14374)-peptide(9594)_dictionary.xlsx"
OUT = ROOT / "data/curated/interaction_tpeppro.parquet"

def load_dict():
    d = pd.read_excel(DICT, sheet_name="Sheet1", header=None, usecols=[0, 1],
                      names=["id", "seq"])
    return dict(zip(d.id.astype(str), d.seq.astype(str)))

def main():
    seqs = load_dict()
    print(f"dictionary: {len(seqs)} id->seq entries")
    act = pd.read_csv(ACTIONS, sep="\t", header=None,
                      names=["receptor_id", "peptide_id", "label"])
    rows, miss = [], 0
    for _, a in act.iterrows():
        rid, pid = str(a.receptor_id), str(a.peptide_id)
        rseq, pseq = seqs.get(rid), seqs.get(pid)
        if not rseq or not pseq:
            miss += 1; continue
        row = {c: None for c in schema.INTERACTION_COLS}
        row.update(
            complex_id=f"tpeppro__{rid}__{pid}__{int(a.label)}",
            source_dataset="tpeppro",
            receptor_pdb_id=rid.split("_")[0],
            receptor_chains=rid.split("_")[1] if "_" in rid else None,
            peptide_chain=pid.split("_")[1] if "_" in pid else None,
            peptide_seq=pseq,
            peptide_len=len(pseq),
            is_cyclic=None,
            cyclization_type="unknown",
            has_ncaa=False,
            receptor_seq=rseq,
            label=int(a.label),
            notes="TPepPro benchmark; structure-derived pair",
        )
        rows.append(row)
    df = pd.DataFrame(rows, columns=schema.INTERACTION_COLS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"wrote {len(df)} interaction pairs -> {OUT} ({miss} dropped for missing seq)")
    print("labels:", df.label.value_counts().to_dict())
    print(f"unique receptors: {df.receptor_pdb_id.nunique()}, unique peptides: {df.peptide_seq.nunique()}")
    print(f"peptide length: {df.peptide_len.min()}-{df.peptide_len.max()} (median {int(df.peptide_len.median())})")

if __name__ == "__main__":
    main()
