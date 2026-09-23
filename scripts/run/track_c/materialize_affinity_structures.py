import os
from pathlib import Path
import pandas as pd
import gemmi

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/track_c/structures"
CP = ROOT / "data/external/MM_PBGBSA-CP/dataset1"
OUT.mkdir(parents=True, exist_ok=True)

def keep_chains(cif_path, chains, out_pdb):
    st = gemmi.read_structure(str(cif_path))
    sel = gemmi.Structure()
    m = gemmi.Model("1")
    for ch in st[0]:
        if ch.name in chains:
            m.add_chain(ch.clone())
    sel.add_model(m); sel.setup_entities(); sel.write_pdb(str(out_pdb))

def main():
    df = pd.read_parquet(ROOT / "runs/track_c/affinity_manifest.parquet")
    ok = miss = 0
    for r in df.itertuples(index=False):
        out = OUT / f"{r.complex_id}.pdb"
        if out.exists():
            ok += 1; continue
        local = CP / r.receptor_pdb_id.upper()
        if r.source_dataset == "mmpbgbsa_cp" and (local / f"{r.receptor_pdb_id.upper()}_protein.pdb").exists():
            prot = (local / f"{r.receptor_pdb_id.upper()}_protein.pdb").read_text()
            pep = (local / f"{r.receptor_pdb_id.upper()}_CP.pdb").read_text()
            out.write_text(prot + "TER\n" + pep + "END\n"); ok += 1; continue
        try:
            path = _fetch(r.receptor_pdb_id)
            chains = set(str(r.receptor_chains).split(",")) | {str(r.peptide_chain)}
            keep_chains(path, chains, out); ok += 1
        except Exception as e:
            miss += 1
            if miss <= 10:
                print("MISS", r.complex_id, e)
    print(f"materialized {ok} | missing {miss} -> {OUT}")

def _fetch(pdb):
    import urllib.request
    dst = OUT / f"_{pdb}.cif"
    if not dst.exists():
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb.upper()}.cif", dst)
    return dst

if __name__ == "__main__":
    main()
