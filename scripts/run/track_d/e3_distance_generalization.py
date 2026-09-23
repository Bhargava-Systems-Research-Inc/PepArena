#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.neural_network import MLPClassifier

ROOT = Path("/12TBDrive1/mega_pep_bench")
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
SPLITS = ROOT / "data/splits/interaction_splits.parquet"
OUT = ROOT / "runs/track_d/e3_distance"
SPLIT = "split_peptide"
BINS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.01]
FRACS = [0.0, 0.2, 0.4, 0.6, 0.8]

def sha1(s):
    return hashlib.sha1(s.encode()).hexdigest()

def pair_features(emb, rseq, pseq):
    r, p = emb[sha1(rseq)], emb[sha1(pseq)]
    return np.concatenate([r, p, np.abs(r - p), r * p])

def nearest_distance(test_peps, train_peps):
    tr = list(train_peps)
    out = {}
    for t in test_peps:
        best = min(Levenshtein.distance(t, s) / max(len(t), len(s)) for s in tr)
        out[t] = best
    return out

def fit_auroc(X, y, tr, te, head):
    if tr.sum() < 50 or len(np.unique(y[te])) < 2:
        return float("nan")
    clf = (MLPClassifier(hidden_layer_sizes=(256,), max_iter=300, random_state=0)
           if head == "mlp" else
           LogisticRegression(max_iter=2000, n_jobs=-1))
    clf.fit(X[tr], y[tr])
    return float(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="composition,esm2_t12_35M_UR50D,esm2_t33_650M_UR50D")
    ap.add_argument("--head", default="mlp")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(MAN)[["complex_id", "receptor_seq", "peptide_seq", "label"]]
    sp = pd.read_parquet(SPLITS).set_index("complex_id")
    df = df.join(sp[[SPLIT]], on="complex_id")
    tr_mask = (df[SPLIT] == "train").to_numpy()
    te_mask = (df[SPLIT] == "test").to_numpy()

    test_peps = sorted(set(df.loc[te_mask, "peptide_seq"]))
    train_peps = sorted(set(df.loc[tr_mask, "peptide_seq"]))
    print(f"{len(test_peps)} unique test peptides vs {len(train_peps)} train peptides", flush=True)
    dist = nearest_distance(test_peps, train_peps)
    pd.Series(dist).to_csv(OUT / "test_peptide_nearest_distance.csv", header=["distance"])
    ds = pd.Series(dist)
    print("nearest-distance distribution:\n" + ds.describe().to_string(), flush=True)

    train_to_test = nearest_distance(train_peps, test_peps)
    tr_rank = sorted(train_peps, key=lambda s: train_to_test[s])

    rows = []
    for model in a.models.split(","):
        cache = ROOT / f"runs/track_d/cache/emb_{model}"
        if not cache.exists():
            print(f"[skip] no cache for {model}"); continue
        emb = {p.stem: np.load(p) for p in cache.glob("*.npy")}
        X = np.vstack([pair_features(emb, r, p)
                       for r, p in zip(df.receptor_seq, df.peptide_seq)]).astype(np.float32)
        y = df.label.to_numpy()
        print(f"\n=== {model}  X{X.shape} ===", flush=True)

        d_col = df.peptide_seq.map(dist)
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            sub = te_mask & (d_col >= lo).to_numpy() & (d_col < hi).to_numpy()
            if sub.sum() < 30:
                continue
            au = fit_auroc(X, y, tr_mask, sub, a.head)
            rows.append({"model": model, "panel": "A_test_stratified", "bin": f"[{lo},{hi})",
                         "removed_frac": 0.0, "n_train": int(tr_mask.sum()),
                         "n_test": int(sub.sum()), "auroc": au})
            print(f"  A dist [{lo},{hi}): n={int(sub.sum()):5d} AUROC {au:.3f}", flush=True)

        for panel, order in (("B_remove_closest", tr_rank),
                             ("C_remove_furthest", tr_rank[::-1])):
            for f in FRACS:
                k = int(len(order) * f)
                drop = set(order[:k])
                keep = tr_mask & ~df.peptide_seq.isin(drop).to_numpy()
                au = fit_auroc(X, y, keep, te_mask, a.head)
                rows.append({"model": model, "panel": panel, "bin": "",
                             "removed_frac": f, "n_train": int(keep.sum()),
                             "n_test": int(te_mask.sum()), "auroc": au})
                print(f"  {panel[0]} removed {f:.0%} of train peptides "
                      f"(n_train={int(keep.sum()):5d}): AUROC {au:.3f}", flush=True)
        del X, emb

    R = pd.DataFrame(rows)
    R.to_csv(OUT / "e3_results.csv", index=False)
    json.dump({"split": SPLIT, "head": a.head, "n_rows": len(R)},
              open(OUT / "e3_meta.json", "w"), indent=2)
    print(f"\nwrote {OUT/'e3_results.csv'} ({len(R)} rows)")

if __name__ == "__main__":
    main()
