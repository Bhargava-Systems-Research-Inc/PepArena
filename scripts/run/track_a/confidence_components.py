import json
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
PRED = R / "runs/track_a/predictions"

def pep_index(cif):
    st = gemmi.read_structure(str(cif)); st.setup_entities(); st.remove_ligands_and_waters()
    lens = [len(c) for c in st[0]]
    return int(np.argmin(lens)), len(lens)

def pair_mean(mat, i, n):
    others = [j for j in range(n) if j != i]
    return float(np.mean([mat[i][j] for j in others])) if others else np.nan

def protenix(cid):
    js = sorted(PRED.glob(f"protenix/{cid}/seed_*/predictions/*summary_confidence_sample_*.json"))
    if not js:
        return None
    best = max(js, key=lambda p: json.load(open(p))["ranking_score"])
    d = json.load(open(best))
    cif = best.parent / best.name.replace("summary_confidence_sample_", "sample_").replace(".json", ".cif")
    if not cif.exists():
        return None
    i, n = pep_index(cif)
    return {"score": d["ranking_score"], "iptm": d["iptm"], "ptm": d["ptm"],
            "iptm_pep": pair_mean(np.array(d["chain_pair_iptm"]), i, n)}

def boltz2(cid):
    js = sorted(PRED.glob(f"boltz2/boltz_results_*/predictions/{cid}/confidence_{cid}_model_0.json"))
    cifs = sorted(PRED.glob(f"boltz2/boltz_results_*/predictions/{cid}/{cid}_model_0.cif"))
    if not js or not cifs:
        return None
    d = json.load(open(js[0]))
    i, n = pep_index(cifs[0])
    pc = d["pair_chains_iptm"]
    mat = [[pc[str(a)][str(b)] for b in range(n)] for a in range(n)] if len(pc) >= n else None
    return {"score": d["confidence_score"], "iptm": d["iptm"], "ptm": d["ptm"],
            "iptm_pep": pair_mean(mat, i, n) if mat else np.nan,
            "ipde": d.get("complex_ipde", np.nan)}

def chai1(cid):
    npzs = sorted(PRED.glob(f"chai1/{cid}/scores.model_idx_*.npz"))
    cifs = sorted(PRED.glob(f"chai1/{cid}/pred.model_idx_*.cif"))
    if not npzs or not cifs:
        return None
    best = max(range(len(npzs)), key=lambda k: float(np.load(npzs[k])["aggregate_score"][0]))
    d = np.load(npzs[best])
    i, n = pep_index(cifs[best])
    return {"score": float(d["aggregate_score"][0]), "iptm": float(d["iptm"][0]),
            "ptm": float(d["ptm"][0]),
            "iptm_pep": pair_mean(d["per_chain_pair_iptm"][0], i, n)}

EXTRACT = {"protenix": protenix, "boltz2": boltz2, "chai1": chai1}

def main():
    out = []
    for m, fn in EXTRACT.items():
        sc = pd.read_parquet(R / f"runs/track_a/scores_pep/{m}.parquet").set_index("complex_id")
        for cid, dq in sc.DockQ.items():
            try:
                c = fn(cid)
            except Exception:
                c = None
            if c:
                out.append({"model": m, "complex_id": cid, "DockQ": float(dq), **c})
        print(f"  {m}: {sum(1 for r in out if r['model']==m)}/{len(sc)}", flush=True)
    d = pd.DataFrame(out)
    dest = R / "runs/track_a/confidence_components.csv"
    d.to_csv(dest, index=False)

    rows = []
    for m, g in d.groupby("model"):
        r = {"model": m, "n": len(g)}
        for c in ("score", "iptm", "iptm_pep", "ptm", "ipde"):
            if c in g and g[c].notna().sum() > 10:
                r[f"rho_{c}"] = round(g[c].corr(g.DockQ, method="spearman"), 3)
        rows.append(r)
    S = pd.DataFrame(rows)
    S.to_csv(R / "runs/track_a/confidence_rho.csv", index=False)
    print(f"\nwrote {dest}")
    print(S.to_string(index=False))

if __name__ == "__main__":
    main()
