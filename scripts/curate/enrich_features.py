from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
CUR = ROOT / "data/curated"
OUT = CUR / "complex_features.parquet"
SSE_MAP = {"a": "H", "b": "E", "c": "C"}

def features(native_path, pep_chain, rec_chains):
    arr = pdbx.get_structure(pdbx.CIFFile.read(native_path), model=1)
    aa = arr[struc.filter_amino_acids(arr)]
    pep = aa[aa.chain_id == str(pep_chain)]
    rec = aa[np.isin(aa.chain_id, [c for c in rec_chains])] if rec_chains else aa[aa.chain_id != str(pep_chain)]
    if pep.array_length() == 0 or rec.array_length() == 0:
        return None
    fr = {"H": 0.0, "E": 0.0, "C": 1.0}
    dom = "C"
    try:
        sse = struc.annotate_sse(pep)
        if len(sse):
            labs = [SSE_MAP.get(s, "C") for s in sse]
            n = len(labs)
            fr = {k: labs.count(k) / n for k in ("H", "E", "C")}
            dom = max(fr, key=fr.get)
    except Exception:
        pass
    rtree = cKDTree(rec.coord)
    dp, _ = rtree.query(pep.coord, distance_upper_bound=5.0)
    n_pep_if = len(set(pep.res_id[dp < 5.0]))
    ptree = cKDTree(pep.coord)
    dr, _ = ptree.query(rec.coord, distance_upper_bound=5.0)
    n_rec_if = len(set(rec.res_id[dr < 5.0]))
    return {"sec_struct": dom, "frac_helix": round(fr["H"], 3),
            "frac_sheet": round(fr["E"], 3), "frac_coil": round(fr["C"], 3),
            "n_pep_interface": n_pep_if, "n_rec_interface": n_rec_if,
            "n_interface_res": n_pep_if + n_rec_if}

def main():
    idx = pd.read_parquet(CUR / "structure_native_index.parquet").set_index("complex_id")
    meta = pd.read_parquet(CUR / "structure_all.parquet").set_index("complex_id")
    rows, fail = [], 0
    for cid, ni in idx.iterrows():
        m = meta.loc[cid] if cid in meta.index else None
        if m is None:
            continue
        recc = [c.strip() for c in str(m.get("receptor_chains") or "").replace(":", ",").split(",") if c.strip()]
        try:
            f = features(ni["native_path"], m["peptide_chain"], recc)
        except Exception:
            f = None
        if f is None:
            fail += 1; continue
        f["complex_id"] = cid
        rows.append(f)
    df = pd.DataFrame(rows)
    df.to_parquet(OUT, index=False)
    print(f"enriched {len(df)} complexes -> {OUT} ({fail} skipped)")
    print("sec_struct:", df.sec_struct.value_counts().to_dict())
    print(f"interface residues: median {int(df.n_interface_res.median())}, "
          f"range {df.n_interface_res.min()}-{df.n_interface_res.max()}")

if __name__ == "__main__":
    main()
