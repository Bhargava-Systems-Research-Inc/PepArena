from __future__ import annotations
import json, shutil, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema, cyclization, struct_utils

SRC = Path("/12TBDrive1/cycesmfold2/data/benchmark")
MANIFEST = SRC / "monomer_manifest.csv"
RAW = ROOT / "data/raw/structures/cyclic_monomers"
OUT = ROOT / "data/curated/structure_monomer_cyclic.parquet"

MODE_MAP = {"head_to_tail": "head_to_tail", "head_to_tail_and_disulfide": "bicyclic",
            "disulfide": "disulfide"}

def main():
    man = pd.read_csv(MANIFEST)
    RAW.mkdir(parents=True, exist_ok=True)
    rows, disagree = [], 0
    for _, m in man.iterrows():
        tid = str(m["target_id"]).lower()
        ref = SRC / "refs" / f"{tid}.pdb"
        if not ref.exists():
            ref = Path(str(m["reference_path"]))
            if not ref.is_absolute():
                ref = Path("/12TBDrive1/cycesmfold2") / m["reference_path"]
        if not ref.exists():
            print(f"  {tid}: ref missing, skip"); continue
        local = RAW / f"{tid}.pdb"
        if not local.exists():
            shutil.copy(ref, local)
        chain = str(m["chain_id"])
        st = struct_utils.read_any(local)
        seq, ncaa, resnames = struct_utils.chain_seq_and_ncaa(st, chain)
        ctype, bonds, evidence = cyclization.classify(st, chain)
        src_mode = MODE_MAP.get(str(m["mode"]), str(m["mode"]))
        ok = (ctype == src_mode)
        if not ok:
            disagree += 1
        row = {c: None for c in schema.STRUCTURE_COLS}
        row.update(
            complex_id=f"cycmono__{tid}_{chain}",
            source_dataset="afcycdesign_monomer80",
            receptor_pdb_id=tid,
            peptide_chain=chain,
            peptide_seq=seq or str(m["sequence"]),
            peptide_len=len(resnames) or len(str(m["sequence"])),
            is_cyclic=(ctype != "linear"),
            cyclization_type=ctype,
            cyclization_bonds=json.dumps(bonds),
            has_ncaa=bool(ncaa),
            ncaa_list=",".join(sorted(set(ncaa))),
            exp_method="unknown",
            native_complex_path=str(local),
            notes=f"monomer; src_mode={src_mode}; agree={ok}; {evidence}",
        )
        rows.append(row)
    df = pd.DataFrame(rows, columns=schema.STRUCTURE_COLS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"wrote {len(df)} cyclic monomers -> {OUT}")
    print("cyclization (ours):", df.cyclization_type.value_counts().to_dict())
    print(f"source-mode disagreements: {disagree}/{len(df)}")
    if disagree:
        bad = df[df.notes.str.contains("agree=False")]
        for _, r in bad.iterrows():
            print("  DISAGREE", r.receptor_pdb_id, "->", r.cyclization_type, "|", r.notes)

if __name__ == "__main__":
    main()
