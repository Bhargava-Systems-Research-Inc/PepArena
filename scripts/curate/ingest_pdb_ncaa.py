from __future__ import annotations
import argparse, gzip, json, sys, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd
import gemmi

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema, cyclization, struct_utils

CIF_DIR = ROOT / "data/raw/structures/cif"
OUT = ROOT / "data/curated/structure_pdb_ncaa.parquet"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query?json="
GRAPHQL = "https://data.rcsb.org/graphql"

NCAA_CODES = [
    "MLE", "MVA", "MAA", "SAR", "AIB", "MK8", "MP8",
    "DAL", "DAR", "DAS", "DCY", "DGL", "DGN", "DHI", "DIL", "DLE", "DLY",
    "DPN", "DPR", "DSN", "DTH", "DTR", "DTY", "DVA", "MED", "DNE", "DSG",
    "ORN", "HYP", "NLE", "ABA", "DAB", "PCA", "BMT", "PHI",
]

def search_entities(codes, lo=4, hi=50):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "entity_poly.rcsb_sample_sequence_length",
            "operator": "range", "value": {"from": lo, "to": hi}}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "entity_poly.rcsb_entity_polymer_type",
            "operator": "exact_match", "value": "Protein"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.chem_comp_monomers",
            "operator": "in", "value": codes}},
    ]}, "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True, "results_content_type": ["experimental"]}}
    url = SEARCH + urllib.parse.quote(json.dumps(q))
    with urllib.request.urlopen(url, timeout=60) as r:
        d = json.load(r)
    return [x["identifier"] for x in d.get("result_set", [])]

def graphql_chains(pdb_ids, batch=60):
    Q = """{ entries(entry_ids:[%s]) { rcsb_id
      rcsb_accession_info { initial_release_date }
      polymer_entities { entity_poly { type rcsb_sample_sequence_length }
        rcsb_polymer_entity_container_identifiers { auth_asym_ids } } } }"""
    info = {}
    for i in range(0, len(pdb_ids), batch):
        q = Q % ",".join(f'"{p}"' for p in pdb_ids[i:i + batch])
        req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            entries = json.load(r)["data"]["entries"]
        for e in entries or []:
            chains = {}
            for pe in e.get("polymer_entities") or []:
                ep = pe.get("entity_poly") or {}
                for c in (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []:
                    chains[c] = {"type": ep.get("type", ""), "len": ep.get("rcsb_sample_sequence_length")}
            info[e["rcsb_id"].upper()] = {
                "chains": chains,
                "date": (e.get("rcsb_accession_info") or {}).get("initial_release_date", "")}
        print(f"  graphql {min(i+batch,len(pdb_ids))}/{len(pdb_ids)}", flush=True)
    return info

def fetch_cif(pdb_id):
    dst = CIF_DIR / f"{pdb_id.lower()}.cif.gz"
    if dst.exists() and dst.stat().st_size > 0:
        return dst
    CIF_DIR.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.cif.gz"
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(url, dst)
            if dst.stat().st_size > 0:
                return dst
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)
    return dst

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-rec-len", type=int, default=50)
    ap.add_argument("--max-pep-len", type=int, default=50)
    args = ap.parse_args()

    ents = search_entities(NCAA_CODES, hi=args.max_pep_len)
    pdbs = sorted({e.split("_")[0].upper() for e in ents})
    print(f"search: {len(ents)} entities across {len(pdbs)} PDB entries")
    info = graphql_chains(pdbs)
    if args.limit:
        pdbs = pdbs[:args.limit]

    rows, kept, no_rec = [], 0, 0
    for n, pdb in enumerate(pdbs):
        meta = info.get(pdb.upper())
        if not meta:
            continue
        chains = meta["chains"]
        peps = sorted(c for c, v in chains.items()
                      if "polypeptide" in v["type"].lower() and v["len"] and 4 <= v["len"] <= args.max_pep_len)
        recs = sorted(c for c, v in chains.items()
                      if "polypeptide" in v["type"].lower() and v["len"] and v["len"] >= args.min_rec_len)
        if not peps or not recs:
            no_rec += 1; continue
        try:
            st = struct_utils.read_any(fetch_cif(pdb))
        except Exception as e:
            print(f"  {pdb}: cif fail {e}"); continue
        chosen = None
        for pc in peps:
            seq, ncaa, resn = struct_utils.chain_seq_and_ncaa(st, pc)
            if ncaa and len(resn) >= 4:
                chosen = (pc, seq, ncaa, resn); break
        if not chosen:
            continue
        pc, seq, ncaa, resn = chosen
        rec = sorted(c for c in recs if c != pc)
        if not rec:
            no_rec += 1; continue
        ctype, bonds, evidence = cyclization.classify(st, pc)
        row = {c: None for c in schema.STRUCTURE_COLS}
        row.update(
            complex_id=f"pdbncaa__{pdb.lower()}_{pc}",
            source_dataset="pdb_ncaa_mine",
            receptor_pdb_id=pdb.lower(),
            receptor_chains=",".join(rec),
            peptide_chain=pc,
            peptide_seq=seq,
            peptide_len=len(resn),
            is_cyclic=(ctype != "linear"),
            cyclization_type=ctype,
            cyclization_bonds=json.dumps(bonds),
            has_ncaa=True,
            ncaa_list=",".join(sorted(set(ncaa))),
            deposition_date=(meta["date"] or "")[:10] or None,
            exp_method="unknown",
            notes=f"ncAA mine; {evidence}",
        )
        rows.append(row); kept += 1
        if kept % 25 == 0:
            print(f"  [{n+1}/{len(pdbs)}] kept {kept}", flush=True)

    df = pd.DataFrame(rows, columns=schema.STRUCTURE_COLS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"\nkept {kept} ncAA peptide-protein complexes -> {OUT} ({no_rec} lacked peptide+receptor)")
    if kept:
        print("cyclization:", df.cyclization_type.value_counts().to_dict())
        print(f"length {df.peptide_len.min()}-{df.peptide_len.max()} (median {int(df.peptide_len.median())})")
        from collections import Counter
        cc = Counter(c for lst in df.ncaa_list for c in lst.split(","))
        print("top ncAA codes:", cc.most_common(12))

if __name__ == "__main__":
    main()
