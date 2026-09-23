#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge

ROOT = Path("/12TBDrive1/mega_pep_bench")
OUT = ROOT / "runs/track_c"
MODEL = "esm2_t12_35M_UR50D"
CACHE = ROOT / f"runs/track_d/cache/emb_{MODEL}"
sha1 = lambda s: hashlib.sha1(s.encode()).hexdigest()

def feats(rows):
    X = []
    for r, p in zip(rows.receptor_seq, rows.peptide_seq):
        a = np.load(CACHE / f"{sha1(r)}.npy")
        b = np.load(CACHE / f"{sha1(p)}.npy")
        X.append(np.concatenate([a, b, np.abs(a - b), a * b]))
    return np.vstack(X)

def main():
    t = pd.read_csv(OUT / "affinity_targets.csv").dropna(subset=["log_affinity"])
    have = lambda s: (CACHE / f"{sha1(s)}.npy").exists()
    t = t[[have(r) and have(p) for r, p in zip(t.receptor_seq, t.peptide_seq)]]
    train = t[t.source_dataset == "ppb_affinity"]
    test = t[t.source_dataset == "mmpbgbsa_cp"]
    print(f"train (PPB-Affinity, out-of-domain for the target set): {len(train)}")
    print(f"test  (cyclic MM_PBGBSA-CP complexes):                  {len(test)}\n")

    Xtr, ytr = feats(train), train.log_affinity.to_numpy()
    Xte, yte = feats(test), test.log_affinity.to_numpy()

    rows = []
    for name, mk in (("ESM-2 + Ridge (learned)", lambda: Ridge(alpha=1.0)),
                     ("ESM-2 + GBM (learned)", HistGradientBoostingRegressor)):
        m = mk().fit(Xtr, ytr)
        rho = float(spearmanr(m.predict(Xte), yte).statistic)
        rows.append({"method": name, "family": "ML", "trained_on": "PPB-Affinity (held out)",
                     "n": len(test), "spearman": round(rho, 3)})

    rng = np.random.default_rng(0)
    idx = rng.permutation(len(train))
    cut = int(0.8 * len(train))
    m = HistGradientBoostingRegressor().fit(Xtr[idx[:cut]], ytr[idx[:cut]])
    rho_in = float(spearmanr(m.predict(Xtr[idx[cut:]]), ytr[idx[cut:]]).statistic)
    rows.append({"method": "ESM-2 + GBM (learned)", "family": "ML",
                 "trained_on": "PPB-Affinity (within-source random split)",
                 "n": len(train) - cut, "spearman": round(rho_in, 3)})

    mm = pd.read_csv(OUT / "mmgbsa_scored.csv")
    mm["cid"] = mm.complex_id.str.replace("mmpbgbsa_cp__", "", regex=False).str.lower()
    tst = test.assign(cid=test.complex_id.str.replace("mmpbgbsa_cp__", "", regex=False).str.lower())
    j = tst.merge(mm[["cid", "dG_mmgbsa"]], on="cid", how="inner")
    if len(j) >= 5:
        rows.append({"method": "MM-GBSA (physics)", "family": "physics", "trained_on": "not applicable (physics-based score)",
                     "n": len(j),
                     "spearman": round(float(spearmanr(j.dG_mmgbsa, j.log_affinity).statistic), 3)})
    bz = OUT / "boltz2_affinity/boltz2_affinity_vs_measured.csv"
    if bz.exists():
        b = pd.read_csv(bz)
        b = b[b.source_dataset == "mmpbgbsa_cp"]
        if len(b) >= 5:
            rows.append({"method": "Boltz-2 affinity head (learned)", "family": "ML",
                         "trained_on": "its own (undisclosed) corpus", "notes": "applied out of domain: >56 heavy atoms",
                         "n": len(b),
                         "spearman": round(float(spearmanr(b.affinity_pred_value,
                                                           b.log_affinity).statistic), 3)})

    R = pd.DataFrame(rows)
    R.to_csv(OUT / "e1_ml_vs_physics.csv", index=False)
    print(R.to_string(index=False))
    json.dump(rows, open(OUT / "e1_ml_vs_physics.json", "w"), indent=2)
    print(f"\nwrote {OUT/'e1_ml_vs_physics.csv'}")

if __name__ == "__main__":
    main()
