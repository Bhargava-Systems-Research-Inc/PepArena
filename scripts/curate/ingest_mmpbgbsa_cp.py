from __future__ import annotations
import json, math, sys, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import schema, cyclization, struct_utils

SRC = ROOT / "data/external/MM_PBGBSA-CP"
CSV = SRC / "dataset1.csv"
OUT = ROOT / "data/curated/affinity_mmpbgbsa_cp.parquet"
GRAPHQL = "https://data.rcsb.org/graphql"

def chaintype_to_cyc(ct):
    ct = str(ct).upper()
    if "BB" in ct and "SS" in ct:
        return "bicyclic"
    if ct.startswith("BB"):
        return "head_to_tail"
    if ct.startswith("SS"):
        return "disulfide"
    return "unknown"

def fetch_dates(pdb_ids):
    Q = '{ entries(entry_ids:[%s]) { rcsb_id rcsb_accession_info { initial_release_date } } }'
    q = Q % ",".join(f'"{p}"' for p in pdb_ids)
    req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            es = json.load(r)["data"]["entries"]
        return {e["rcsb_id"].upper(): (e.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10]
                for e in es or []}
    except Exception:
        return {}

def peptide_chain(st):
    for c in st[0]:
        if any((i := struct_utils.gemmi.find_tabulated_residue(r.name)) and i.is_amino_acid() for r in c):
            return c.name
    return None

def main():
    df = pd.read_csv(CSV)
    dates = fetch_dates([str(p).upper() for p in df["PDB code"]])
    rows, disagree = [], 0
    for _, r in df.iterrows():
        pdb = str(r["PDB code"])
        rec, pep = (str(r["Chain id"]).split(":") + [""])[:2]
        cp = SRC / "dataset1" / pdb / f"{pdb}_CP.pdb"
        prot = SRC / "dataset1" / pdb / f"{pdb}_protein.pdb"
        if not cp.exists():
            print(f"  {pdb}: CP structure missing"); continue
        st = struct_utils.read_any(cp)
        ch = peptide_chain(st)
        seq, ncaa, resn = struct_utils.chain_seq_and_ncaa(st, ch)
        ctype, _, ev = cyclization.classify(st, ch)
        expected = chaintype_to_cyc(r["Chain type"])
        ok = (ctype == expected) or (expected == "disulfide" and ctype == "disulfide")
        if not ok:
            disagree += 1
        kd = float(r["Kd (M)"])
        row = {c: None for c in schema.AFFINITY_COLS}
        row.update(
            complex_id=f"mmpbgbsa_cp__{pdb.lower()}_{pep or ch}",
            source_dataset="mmpbgbsa_cp", receptor_pdb_id=pdb.lower(),
            receptor_chains=rec, peptide_chain=pep or ch, peptide_seq=seq,
            peptide_len=len(resn), is_cyclic=True,
            cyclization_type=ctype if ok else expected,
            has_ncaa=bool(ncaa), ncaa_list=",".join(sorted(set(ncaa))),
            measured_value=kd, value_type="Kd", value_units="M",
            log_affinity=math.log10(kd), assay="literature (MM_PBGBSA-CP datasetI)",
            deposition_date=dates.get(pdb.upper()) or None, exp_method="xray",
            native_complex_path=str(prot),
            notes=f"chain_type={r['Chain type']}; det={ctype}; agree={ok}; pep={cp}; {ev}")
        rows.append(row)
    out = pd.DataFrame(rows, columns=schema.AFFINITY_COLS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)
    print(f"wrote {len(out)} cyclic affinity measurements -> {OUT}")
    print("cyclization:", out.cyclization_type.value_counts().to_dict())
    print(f"detector vs paper Chain-type disagreements: {disagree}/{len(out)}")
    print(f"pKd(=-log10 Kd) range: {(-out.log_affinity).min():.2f} .. {(-out.log_affinity).max():.2f}; ncAA: {int(out.has_ncaa.sum())}")

if __name__ == "__main__":
    main()
