from __future__ import annotations
import json, shutil, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema, cyclization, struct_utils

SRC34 = Path("/12TBDrive1/cycproteinix/benchmark/complex34")
RAW = ROOT / "data/raw/structures/cyclic_complex34"
OUT = ROOT / "data/curated/structure_cyclic_complexes.parquet"

def ingest_complex34():
    meta = pd.read_csv(SRC34 / "metadata.csv")
    RAW.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, m in meta.iterrows():
        target = str(m["target"])
        ref = SRC34 / "refs" / f"{target}.pdb"
        if not ref.exists():
            print(f"  {target}: ref missing, skip"); continue
        local = RAW / f"{target}.pdb"
        if not local.exists():
            shutil.copy(ref, local)
        pep_chain = str(m["pep_chain"])
        st = struct_utils.read_any(local)
        seq, ncaa, resnames = struct_utils.chain_seq_and_ncaa(st, pep_chain)
        ctype, bonds, evidence = cyclization.classify(st, pep_chain)
        row = {c: None for c in schema.STRUCTURE_COLS}
        row.update(
            complex_id=f"cyc34__{target}_{pep_chain}",
            source_dataset="cycproteinix_complex34",
            receptor_pdb_id=target.lower(),
            receptor_chains=str(m["rec_chain"]),
            peptide_chain=pep_chain,
            peptide_seq=seq,
            peptide_len=len(resnames) or int(m["pep"]),
            is_cyclic=(ctype != "linear"),
            cyclization_type=ctype,
            cyclization_bonds=json.dumps(bonds),
            has_ncaa=bool(ncaa),
            ncaa_list=",".join(sorted(set(ncaa))),
            deposition_date=str(m["release"]) if pd.notna(m.get("release")) else None,
            exp_method="unknown",
            native_complex_path=str(local),
            notes=evidence,
        )
        rows.append(row)
        print(f"  {target} {pep_chain} len={row['peptide_len']} {ctype}"
              f"{' ncAA' if ncaa else ''} rel={row['deposition_date']}", flush=True)
    return rows

def main():
    rows = ingest_complex34()
    df = pd.DataFrame(rows, columns=schema.STRUCTURE_COLS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"\nwrote {len(df)} cyclic complexes -> {OUT}")
    print("cyclization:", df.cyclization_type.value_counts().to_dict())

if __name__ == "__main__":
    main()
