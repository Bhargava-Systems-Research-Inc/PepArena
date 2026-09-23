from __future__ import annotations
import json, sys, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema, cyclization, struct_utils

SRC = ROOT / "data/raw/pepbench/train_valid"
OUT = ROOT / "data/curated/structure_pepbench.parquet"
GRAPHQL = "https://data.rcsb.org/graphql"

def graphql_dates(pdb_ids, batch=80):
    Q = """{ entries(entry_ids:[%s]) { rcsb_id
      rcsb_accession_info { initial_release_date }
      refine { ls_d_res_high } exptl { method } } }"""
    info = {}
    ids = sorted(pdb_ids)
    for i in range(0, len(ids), batch):
        q = Q % ",".join(f'"{p}"' for p in ids[i:i + batch])
        req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                entries = json.load(r)["data"]["entries"]
        except Exception:
            entries = []
        for e in entries or []:
            method = ((e.get("exptl") or [{}])[0] or {}).get("method", "") or ""
            m = "xray" if "X-RAY" in method else "nmr" if "NMR" in method else \
                "em" if "MICROSCOPY" in method else "unknown"
            res = ((e.get("refine") or [{}])[0] or {}).get("ls_d_res_high")
            info[e["rcsb_id"].upper()] = {
                "date": (e.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
                "res": float(res) if res else None, "method": m}
        print(f"  dates {min(i+batch,len(ids))}/{len(ids)}", flush=True)
    return info

def main():
    entries = []
    for line in (SRC / "all.txt").read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name, rec, pep = parts[0], parts[1], parts[2]
        pdb = name.split("pdb")[-1].lower()
        entries.append((name, pdb, rec, pep))
    pdbs = {e[1].upper() for e in entries}
    print(f"{len(entries)} complexes, {len(pdbs)} unique PDBs; fetching dates...")
    meta = graphql_dates(pdbs)

    rows, fail = [], 0
    for n, (name, pdb, rec, pep) in enumerate(entries):
        pdbfile = SRC / "pdbs" / f"{name}.pdb"
        if not pdbfile.exists():
            fail += 1; continue
        try:
            st = struct_utils.read_any(pdbfile)
            seq, ncaa, resn = struct_utils.chain_seq_and_ncaa(st, pep)
            if not resn:
                fail += 1; continue
            ctype, bonds, ev = cyclization.classify(st, pep)
            md = meta.get(pdb.upper(), {})
            row = {c: None for c in schema.STRUCTURE_COLS}
            row.update(
                complex_id=f"pepbench__{pdb}_{rec}_{pep}", source_dataset="pepbench",
                receptor_pdb_id=pdb, receptor_chains=rec, peptide_chain=pep,
                peptide_seq=seq, peptide_len=len(resn),
                is_cyclic=(ctype != "linear"), cyclization_type=ctype,
                cyclization_bonds=json.dumps(bonds), has_ncaa=bool(ncaa),
                ncaa_list=",".join(sorted(set(ncaa))), deposition_date=md.get("date") or None,
                resolution=md.get("res"), exp_method=md.get("method", "unknown"),
                native_complex_path=str(pdbfile), notes=f"pepbench; {ev}")
            rows.append(row)
        except Exception:
            fail += 1
        if (n + 1) % 1000 == 0:
            print(f"  parsed {n+1}/{len(entries)}", flush=True)

    df = pd.DataFrame(rows, columns=schema.STRUCTURE_COLS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"\nwrote {len(df)} -> {OUT} ({fail} skipped)")
    print("cyclization:", df.cyclization_type.value_counts().to_dict())
    print("ncAA:", int(df.has_ncaa.sum()), "| length median", int(df.peptide_len.median()),
          "| dated:", int(df.deposition_date.notna().sum()))

if __name__ == "__main__":
    main()
