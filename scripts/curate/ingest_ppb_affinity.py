from __future__ import annotations
import argparse, json, math, sys, time, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema

XLSX = ROOT / "data/raw/affinity/PPB-Affinity.xlsx"
OUT = ROOT / "data/curated/affinity_ppb.parquet"
GRAPHQL = "https://data.rcsb.org/graphql"

Q = """{ entries(entry_ids:[%s]) { rcsb_id
  polymer_entities { entity_poly { type rcsb_sample_sequence_length pdbx_seq_one_letter_code_can }
    rcsb_polymer_entity_container_identifiers { auth_asym_ids } } } }"""

def fetch_chain_info(pdb_ids, batch=70):
    info = {}
    for i in range(0, len(pdb_ids), batch):
        chunk = pdb_ids[i:i + batch]
        q = Q % ",".join(f'"{p}"' for p in chunk)
        req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json"})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = json.load(r)["data"]["entries"]
                break
            except Exception as e:
                if attempt == 2:
                    print(f"  batch {i} failed: {e}"); data = []
                time.sleep(2)
        for entry in data or []:
            pdb = entry["rcsb_id"].upper()
            chains = {}
            for pe in entry.get("polymer_entities") or []:
                ep = pe.get("entity_poly") or {}
                ids = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []
                for c in ids:
                    chains[c] = {"type": ep.get("type", ""),
                                 "len": ep.get("rcsb_sample_sequence_length"),
                                 "seq": (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")}
            info[pdb] = chains
        print(f"  fetched {min(i+batch,len(pdb_ids))}/{len(pdb_ids)} PDBs", flush=True)
    return info

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-pep-len", type=int, default=50)
    ap.add_argument("--min-pep-len", type=int, default=4)
    args = ap.parse_args()

    df = pd.read_excel(XLSX)
    pdbs = sorted(df.PDB.dropna().astype(str).str.upper().unique())
    print(f"{len(df)} affinity rows over {len(pdbs)} PDBs; fetching chain info from RCSB...")
    info = fetch_chain_info(pdbs)

    def is_pep(pdb, ch):
        c = info.get(pdb, {}).get(ch)
        return c and "polypeptide" in c["type"].lower() and c["len"] \
            and args.min_pep_len <= c["len"] <= args.max_pep_len

    def is_prot(pdb, ch):
        c = info.get(pdb, {}).get(ch)
        return c and "polypeptide" in c["type"].lower() and c["len"] and c["len"] > args.max_pep_len

    rows, kept = [], 0
    for _, r in df.iterrows():
        pdb = str(r.PDB).upper()
        lig = [c.strip() for c in str(r["Ligand Chains"]).split(",") if c.strip()]
        rec = [c.strip() for c in str(r["Receptor Chains"]).split(",") if c.strip()]
        pep_ch = [c for c in lig if is_pep(pdb, c)]
        if len(pep_ch) != 1 or not any(is_prot(pdb, c) for c in rec):
            continue
        kd = r["KD(M)"]
        if not (isinstance(kd, (int, float)) and kd > 0):
            continue
        ch = pep_ch[0]
        seq = info[pdb][ch]["seq"]
        row = {c: None for c in schema.AFFINITY_COLS}
        row.update(
            complex_id=f"ppb__{pdb}_{ch}",
            source_dataset="ppb_affinity",
            receptor_pdb_id=pdb.lower(),
            receptor_chains=",".join(rec),
            peptide_chain=ch,
            peptide_seq=seq,
            peptide_len=len(seq),
            is_cyclic=None,
            cyclization_type="unknown",
            has_ncaa=False,
            measured_value=float(kd),
            value_type="Kd",
            value_units="M",
            log_affinity=math.log10(float(kd)),
            assay=str(r.get("Affinity Method") or ""),
            deposition_date=str(r.get("PDB Release Date") or "")[:10] or None,
            resolution=r.get("Resolution(Å)"),
            exp_method="xray",
            notes=f"src={r.get('Source Data Set')}; mut={r.get('Mutations')}",
        )
        rows.append(row); kept += 1

    out = pd.DataFrame(rows, columns=schema.AFFINITY_COLS)
    before = len(out)
    out = out.drop_duplicates(
        subset=["receptor_pdb_id", "peptide_chain", "measured_value", "value_type", "assay", "notes"]
    ).reset_index(drop=True)
    out["complex_id"] = out["complex_id"] + "__m" + \
        out.groupby(["receptor_pdb_id", "peptide_chain"]).cumcount().astype(str)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)
    print(f"\nkept {len(out)} affinity measurements ({before - len(out)} exact-dup rows dropped) -> {OUT}")
    if len(out):
        print(f"  unique complex_id: {out.complex_id.nunique()} (== rows), "
              f"distinct complexes: {out.groupby(['receptor_pdb_id','peptide_chain']).ngroups}")
        print(f"  unique peptides: {out.peptide_seq.nunique()}, unique PDBs: {out.receptor_pdb_id.nunique()}")
        print(f"  peptide length: {out.peptide_len.min()}-{out.peptide_len.max()} (median {int(out.peptide_len.median())})")
        print(f"  log10 Kd range: {out.log_affinity.min():.1f} .. {out.log_affinity.max():.1f}")

if __name__ == "__main__":
    main()
