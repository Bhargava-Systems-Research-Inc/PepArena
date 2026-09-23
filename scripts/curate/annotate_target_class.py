from __future__ import annotations
import json, sys, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CUR = ROOT / "data/curated"
OUT = CUR / "complex_target_class.parquet"
GRAPHQL = "https://data.rcsb.org/graphql"

KEYWORDS = [
    ("mhc", "MHC"), ("hla", "MHC"), ("histocompat", "MHC"), ("beta-2 microglob", "MHC"),
    ("pdz", "PDZ"),
    ("sh3", "SH3"), ("sh2", "SH2"), ("ww domain", "WW"),
    ("bromodomain", "bromodomain"), ("14-3-3", "14-3-3"),
    ("protease", "protease"), ("peptidase", "protease"), ("trypsin", "protease"),
    ("proteinase", "protease"), ("caspase", "protease"), ("elastase", "protease"),
    ("kinase", "kinase"), ("phosphatase", "phosphatase"),
    ("ubiquitin ligase", "ubiquitin-system"), ("e3 ", "ubiquitin-system"),
    ("ubiquitin", "ubiquitin-system"), ("sumo", "ubiquitin-system"),
    ("bcl-2", "apoptosis"), ("bcl2", "apoptosis"), ("mdm2", "apoptosis"),
    ("integrin", "integrin"), ("g protein-coupled", "GPCR"), ("gpcr", "GPCR"),
    ("receptor", "receptor"), ("antibody", "antibody"), ("immunoglobulin", "antibody"),
    ("fab ", "antibody"), ("nanobody", "antibody"),
    ("chaperone", "chaperone"), ("heat shock", "chaperone"), ("hsp", "chaperone"),
    ("transcription factor", "transcription-factor"), ("bromo", "bromodomain"),
    ("calmodulin", "calmodulin"), ("importin", "transport"), ("karyopherin", "transport"),
]

def classify(text):
    t = (text or "").lower()
    for kw, cls in KEYWORDS:
        if kw in t:
            return cls
    return "other"

def fetch(pdb_ids, batch=60):
    Q = """{ entries(entry_ids:[%s]) { rcsb_id
      struct { title }
      polymer_entities { rcsb_polymer_entity { pdbx_description }
        rcsb_polymer_entity_container_identifiers { auth_asym_ids } } } }"""
    info = {}
    for i in range(0, len(pdb_ids), batch):
        q = Q % ",".join(f'"{p}"' for p in pdb_ids[i:i + batch])
        req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            entries = json.load(r)["data"]["entries"]
        for e in entries or []:
            ents = {}
            for pe in e.get("polymer_entities") or []:
                desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
                for c in (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []:
                    ents[c] = desc
            info[e["rcsb_id"].upper()] = {"title": (e.get("struct") or {}).get("title", ""),
                                          "entities": ents}
        print(f"  {min(i+batch,len(pdb_ids))}/{len(pdb_ids)}", flush=True)
    return info

def main():
    df = pd.read_parquet(CUR / "structure_all.parquet")
    pdbs = sorted({str(p).upper() for p in df.receptor_pdb_id.dropna().unique()})
    print(f"fetching receptor descriptions for {len(pdbs)} PDBs...")
    info = fetch(pdbs)
    rows = []
    for _, r in df.iterrows():
        pdb = str(r.receptor_pdb_id).upper()
        meta = info.get(pdb, {})
        rec_chains = [c.strip() for c in str(r.receptor_chains or "").replace(":", ",").split(",") if c.strip()]
        desc = next((meta.get("entities", {}).get(c) for c in rec_chains if meta.get("entities", {}).get(c)), "")
        cls = classify(desc) if desc else classify(meta.get("title"))
        rows.append({"complex_id": r.complex_id, "target_class": cls,
                     "receptor_name": (desc or meta.get("title") or "")[:80]})
    out = pd.DataFrame(rows)
    out.to_parquet(OUT, index=False)
    print(f"\nannotated {len(out)} complexes -> {OUT}")
    print(out.target_class.value_counts().to_string())

if __name__ == "__main__":
    main()
