from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema, cyclization, struct_utils
from ingest_peppcbench import fetch_cif, read_structure

def ingest_pdb_list(rows, source, out_path):
    out, fail = [], 0
    for i, r in enumerate(rows):
        pdb = str(r["pdb_id"]).lower()
        pep = str(r["peptide_chain"])
        rec = str(r["receptor_chains"]).replace(":", ",")
        try:
            st, dep, res, method = read_structure(fetch_cif(pdb))
            seq, ncaa, resn = struct_utils.chain_seq_and_ncaa(st, pep)
            if not resn:
                fail += 1; continue
            ctype, bonds, ev = cyclization.classify(st, pep)
            row = {c: None for c in schema.STRUCTURE_COLS}
            row.update(
                complex_id=f"{source}__{pdb}_{pep}", source_dataset=source,
                receptor_pdb_id=pdb, receptor_chains=rec, peptide_chain=pep,
                peptide_seq=seq, peptide_len=len(resn),
                is_cyclic=(ctype != "linear"), cyclization_type=ctype,
                cyclization_bonds=json.dumps(bonds), has_ncaa=bool(ncaa),
                ncaa_list=",".join(sorted(set(ncaa))), deposition_date=dep or None,
                resolution=res, exp_method=method, notes=ev)
            out.append(row)
        except Exception as e:
            fail += 1
            print(f"  {pdb} {pep}: FAIL {e}")
        if (i + 1) % 25 == 0:
            print(f"  {i+1} done", flush=True)
    df = pd.DataFrame(out, columns=schema.STRUCTURE_COLS)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    print(f"wrote {len(df)} -> {out_path} ({fail} failed)")
    if len(df):
        print("cyclization:", df.cyclization_type.value_counts().to_dict(),
              "| ncAA:", int(df.has_ncaa.sum()))
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    d = pd.read_csv(args.tsv, sep="\t", header=None,
                    names=["pdb_id", "receptor_chains", "peptide_chain", "extra"],
                    usecols=[0, 1, 2])
    ingest_pdb_list(d.to_dict("records"), args.source, args.out)

if __name__ == "__main__":
    main()
