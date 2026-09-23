import json
import time
import urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/track_c"
GQL = "https://data.rcsb.org/graphql"

Q = """query($id:String!){entry(entry_id:$id){polymer_entities{
 entity_poly{pdbx_seq_one_letter_code_can}
 rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}"""

def fetch(pdb):
    req = urllib.request.Request(GQL, data=json.dumps({"query": Q, "variables": {"id": pdb.upper()}}).encode(),
                                 headers={"Content-Type": "application/json"})
    d = json.loads(urllib.request.urlopen(req, timeout=30).read())
    ents = (d.get("data") or {}).get("entry", {}).get("polymer_entities") or []
    out = {}
    for e in ents:
        seq = (e["entity_poly"] or {}).get("pdbx_seq_one_letter_code_can", "")
        for ch in (e["rcsb_polymer_entity_container_identifiers"] or {}).get("auth_asym_ids", []) or []:
            out[ch] = seq.replace("\n", "")
    return out

def main():
    df = pd.read_parquet(OUT / "affinity_manifest.parquet")
    want = df[["receptor_pdb_id", "receptor_chains"]].drop_duplicates()
    done = {}
    csv = OUT / "receptor_seqs.csv"
    if csv.exists():
        done = dict(zip(pd.read_csv(csv).receptor_pdb_id, pd.read_csv(csv).receptor_seq))
    rows = [{"receptor_pdb_id": k, "receptor_seq": v} for k, v in done.items()]
    for i, r in enumerate(want.itertuples(index=False)):
        if r.receptor_pdb_id in done:
            continue
        try:
            chains = fetch(r.receptor_pdb_id)
            first = str(r.receptor_chains).split(",")[0]
            seq = chains.get(first) or (max(chains.values(), key=len) if chains else "")
            rows.append({"receptor_pdb_id": r.receptor_pdb_id, "receptor_seq": seq})
        except Exception as e:
            print("FAIL", r.receptor_pdb_id, e)
        if (i + 1) % 50 == 0:
            pd.DataFrame(rows).to_csv(csv, index=False)
            print(f"  {i+1}/{len(want)}", flush=True)
        time.sleep(0.05)
    pd.DataFrame(rows).to_csv(csv, index=False)
    print(f"receptor seqs -> {csv} ({len(rows)} PDBs)")

if __name__ == "__main__":
    main()
