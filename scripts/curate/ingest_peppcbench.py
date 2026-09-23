from __future__ import annotations
import argparse, gzip, json, os, sys, urllib.request
from pathlib import Path
import pandas as pd
import gemmi

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema
import cyclization
import struct_utils

JOB_LIST = ROOT / "data/external/PepPCBench/job_list.csv"
CIF_DIR = ROOT / "data/raw/structures/cif"
OUT = ROOT / "data/curated/structure_peppcbench.parquet"

def fetch_cif(pdb_id: str) -> Path:
    pdb_id = pdb_id.lower()
    dst = CIF_DIR / f"{pdb_id}.cif.gz"
    if dst.exists() and dst.stat().st_size > 0:
        return dst
    CIF_DIR.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    urllib.request.urlretrieve(url, dst)
    return dst

def read_structure(cif_gz: Path):
    with gzip.open(cif_gz, "rt") as fh:
        doc = gemmi.cif.read_string(fh.read())
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    dep = block.find_value("_pdbx_database_status.recvd_initial_deposition_date") \
          or block.find_value("_database_PDB_rev.date_original") or ""
    dep = dep.strip().strip("'\"") if dep else ""
    res = block.find_value("_refine.ls_d_res_high") or block.find_value("_em_3d_reconstruction.resolution") or ""
    try:
        res = float(str(res).strip().strip("'\""))
    except ValueError:
        res = None
    method = (block.find_value("_exptl.method") or "").strip().strip("'\"").lower()
    if "x-ray" in method:
        method = "xray"
    elif "nmr" in method:
        method = "nmr"
    elif "microscop" in method or "em" == method:
        method = "em"
    else:
        method = method[:10] or "unknown"
    return st, dep, res, method

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    jobs = pd.read_csv(JOB_LIST)
    if args.limit:
        jobs = jobs.head(args.limit)

    rows = []
    n_fail = 0
    for i, job in jobs.iterrows():
        pdb_id = str(job["pdb_id"]).lower()
        pep_chain = str(job["peptide_chains"]).split(":")[0]
        rec_chains = str(job["protein_chains"]).replace(":", ",")
        try:
            cif = fetch_cif(pdb_id)
            st, dep, res, method = read_structure(cif)
            seq, ncaa, resnames = struct_utils.chain_seq_and_ncaa(st, pep_chain)
            plen = len(resnames)
            ctype, bonds, evidence = cyclization.classify(st, pep_chain)
            row = {c: None for c in schema.STRUCTURE_COLS}
            row.update(
                complex_id=f"peppcbench__{pdb_id}_{pep_chain}",
                source_dataset="peppcbench",
                receptor_pdb_id=pdb_id,
                receptor_chains=rec_chains,
                peptide_chain=pep_chain,
                peptide_seq=seq,
                peptide_len=plen,
                is_cyclic=(ctype != "linear"),
                cyclization_type=ctype,
                cyclization_bonds=json.dumps(bonds),
                has_ncaa=bool(ncaa),
                ncaa_list=",".join(sorted(set(ncaa))),
                deposition_date=dep or None,
                resolution=res,
                exp_method=method,
                notes=evidence,
            )
            rows.append(row)
            print(f"[{i+1}/{len(jobs)}] {pdb_id} {pep_chain} len={plen} "
                  f"{ctype}{' ncAA' if ncaa else ''} {dep}", flush=True)
        except Exception as e:
            n_fail += 1
            print(f"[{i+1}/{len(jobs)}] {pdb_id} FAILED: {e}", flush=True)

    df = pd.DataFrame(rows, columns=schema.STRUCTURE_COLS)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    print(f"\nwrote {len(df)} rows -> {args.out}  ({n_fail} failed)")

if __name__ == "__main__":
    main()
