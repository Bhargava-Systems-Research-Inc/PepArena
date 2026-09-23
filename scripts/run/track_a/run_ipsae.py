import shutil
import subprocess
import tempfile
from pathlib import Path

import gemmi
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
IPSAE = R / "methods/external/ipsae/ipsae.py"
PY = Path.home() / "anaconda3/envs/pepgym-curate/bin/python"
PRED = R / "runs/track_a/predictions/boltz2"

def peptide_chain(cif):
    st = gemmi.read_structure(str(cif)); st.setup_entities(); st.remove_ligands_and_waters()
    return min(st[0], key=lambda c: len(c)).name

def one(cid, tmp):
    pae = next(PRED.glob(f"*/predictions/{cid}/pae_{cid}_model_0.npz"), None)
    cif = next(PRED.glob(f"*/predictions/{cid}/{cid}_model_0.cif"), None)
    conf = next(PRED.glob(f"*/predictions/{cid}/confidence_{cid}_model_0.json"), None)
    if not (pae and cif):
        return None
    d = Path(tempfile.mkdtemp(dir=tmp))
    for f in (pae, cif, conf):
        if f:
            shutil.copy(f, d)
    r = subprocess.run([str(PY), str(IPSAE), str(d / pae.name), str(d / cif.name), "10", "10"],
                       capture_output=True, text=True, timeout=600)
    out = next(d.glob("*_10_10.txt"), None)
    if out is None:
        shutil.rmtree(d, ignore_errors=True)
        return None
    pep = peptide_chain(cif)
    rows = []
    for line in out.read_text().splitlines()[1:]:
        p = line.split()
        if len(p) < 13 or p[4] != "max":
            continue
        if pep not in (p[0], p[1]):
            continue
        rows.append({"ipSAE": float(p[5]), "ipSAE_d0chn": float(p[6]),
                     "ipTM_af": float(p[8]), "pDockQ": float(p[10]),
                     "pDockQ2": float(p[11]), "LIS": float(p[12])})
    shutil.rmtree(d, ignore_errors=True)
    if not rows:
        return None
    best = max(rows, key=lambda r: r["ipSAE"])
    return {"complex_id": cid, **best}

def main():
    sc = pd.read_parquet(R / "runs/track_a/scores_pep/boltz2.parquet").set_index("complex_id")
    dest = R / "runs/track_a/ipsae_boltz2.csv"
    done = set(pd.read_csv(dest).complex_id) if dest.exists() else set()
    rows = pd.read_csv(dest).to_dict("records") if dest.exists() else []
    with tempfile.TemporaryDirectory() as tmp:
        ids = [c for c in sc.index if c not in done]
        print(f"ipsae: {len(ids)} to run ({len(done)} done)", flush=True)
        for k, cid in enumerate(ids, 1):
            try:
                r = one(cid, tmp)
            except Exception as e:
                print(f"  ERR {cid} {type(e).__name__}", flush=True); r = None
            if r:
                r["DockQ"] = float(sc.DockQ[cid]); rows.append(r)
            if k % 25 == 0:
                pd.DataFrame(rows).to_csv(dest, index=False)
                print(f"  {k}/{len(ids)}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(dest, index=False)
    print(f"\nwrote {dest}  n={len(d)}")
    rows = []
    for c in ("ipSAE", "ipSAE_d0chn", "ipTM_af", "pDockQ", "pDockQ2", "LIS"):
        if c in d:
            rows.append({"metric": c, "n": len(d), "n_distinct": int(d[c].nunique()),
                         "rho": round(float(d[c].corr(d.DockQ, method="spearman")), 3),
                         "median": round(float(d[c].median()), 4)})
    S = pd.DataFrame(rows)
    S.to_csv(R / "runs/track_a/confidence_rho_boltz2.csv", index=False)
    print(S.to_string(index=False))

if __name__ == "__main__":
    main()
