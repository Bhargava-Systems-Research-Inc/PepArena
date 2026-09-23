#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.neural_network import MLPClassifier

ROOT = Path("/12TBDrive1/mega_pep_bench")
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
SPLITS = ROOT / "data/splits/interaction_splits.parquet"
OUT = ROOT / "runs/track_d/e4_negatives"
SPLIT = "split_receptor"
RNG = np.random.default_rng(0)

def sha1(s):
    return hashlib.sha1(s.encode()).hexdigest()

def pair_features(emb, r, p):
    a, b = emb[sha1(r)], emb[sha1(p)]
    return np.concatenate([a, b, np.abs(a - b), a * b])

def build(pos, design, pool):
    rows = [{"receptor_seq": r, "peptide_seq": p, "label": 1}
            for r, p in zip(pos.receptor_seq, pos.peptide_seq)]
    true_pairs = set(zip(pos.receptor_seq, pos.peptide_seq))
    if design == "shared_pool":
        cand = pool
    else:
        cand = sorted(set(pos.peptide_seq))
    cand = list(cand)
    by_len = {}
    for s in cand:
        by_len.setdefault(len(s), []).append(s)
    for r, p in zip(pos.receptor_seq, pos.peptide_seq):
        choices = by_len.get(len(p), cand) if design == "length_matched" else cand
        for _ in range(12):
            d = choices[RNG.integers(len(choices))]
            if (r, d) not in true_pairs and d != p:
                rows.append({"receptor_seq": r, "peptide_seq": d, "label": 0})
                break
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--head", default="mlp")
    ap.add_argument("--seed", type=int, default=0,
                    help="seeds the negative draw and the head; sweep it to bound the variance")
    a = ap.parse_args()
    global RNG
    RNG = np.random.default_rng(a.seed)
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(MAN)[["complex_id", "receptor_seq", "peptide_seq", "label"]]
    sp = pd.read_parquet(SPLITS).set_index("complex_id")
    df = df.join(sp[[SPLIT]], on="complex_id")
    pos = df[df.label == 1]
    tr_pos = pos[pos[SPLIT] == "train"]
    te_pos = pos[pos[SPLIT] == "test"]
    print(f"positives: {len(tr_pos)} train / {len(te_pos)} test", flush=True)

    pool = sorted(set(tr_pos.peptide_seq))
    pool = [pool[i] for i in RNG.choice(len(pool), size=min(400, len(pool)), replace=False)]

    cache = ROOT / f"runs/track_d/cache/emb_{a.model}"
    emb = {p.stem: np.load(p) for p in cache.glob("*.npy")}

    rows = []
    for design in ("shared_pool", "within_split", "length_matched"):
        tr = build(tr_pos, design, pool)
        te = build(te_pos, design, pool)
        Xtr = np.vstack([pair_features(emb, r, p)
                         for r, p in zip(tr.receptor_seq, tr.peptide_seq)]).astype(np.float32)
        Xte = np.vstack([pair_features(emb, r, p)
                         for r, p in zip(te.receptor_seq, te.peptide_seq)]).astype(np.float32)
        clf = (MLPClassifier(hidden_layer_sizes=(256,), max_iter=300, random_state=a.seed)
               if a.head == "mlp" else LogisticRegression(max_iter=2000, n_jobs=-1))
        clf.fit(Xtr, tr.label.to_numpy())
        s = clf.predict_proba(Xte)[:, 1]
        au = float(roc_auc_score(te.label, s))
        apr = float(average_precision_score(te.label, s))
        rows.append({"design": design, "n_train": len(tr), "n_test": len(te),
                     "auroc": round(au, 4), "aupr": round(apr, 4)})
        print(f"  {design:15s} n_train={len(tr):6d} n_test={len(te):5d}  "
              f"AUROC {au:.3f}  AUPR {apr:.3f}", flush=True)
        del Xtr, Xte

    R = pd.DataFrame(rows)
    suffix = "" if a.seed == 0 else f"_seed{a.seed}"
    R.to_csv(OUT / f"e4_results{suffix}.csv", index=False)
    base = R.loc[R.design == "within_split", "auroc"].iloc[0]
    infl = {r.design: round(r.auroc - base, 4) for r in R.itertuples()}
    json.dump({"model": a.model, "head": a.head, "split": SPLIT,
               "auroc_by_design": {r.design: r.auroc for r in R.itertuples()},
               "seed": a.seed,
               "inflation_vs_within_split": infl}, open(OUT / f"e4_summary{suffix}.json", "w"), indent=2)
    print("\ninflation vs the honest within-split design:", json.dumps(infl))

if __name__ == "__main__":
    main()
