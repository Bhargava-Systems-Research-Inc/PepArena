import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef

ROOT = Path("/12TBDrive1/mega_pep_bench")
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
SPLITS = ROOT / "data/splits/interaction_splits.parquet"
SPLIT_COLS = ["split_random", "split_peptide", "split_receptor", "split_both",
              "split_receptor_strict"]

def head(seed):
    return make_pipeline(
        StandardScaler(),
        MLPClassifier(hidden_layer_sizes=(256,), activation="relu", solver="adam",
                      learning_rate_init=1e-3, early_stopping=True, validation_fraction=0.1,
                      max_iter=300, random_state=seed))

def sha1(s):
    return hashlib.sha1(s.encode()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--splits", default="")
    a = ap.parse_args()

    df = pd.read_parquet(MAN)[["complex_id", "receptor_seq", "peptide_seq", "label"]]
    sp = pd.read_parquet(SPLITS).set_index("complex_id")
    df = df.join(sp[SPLIT_COLS], on="complex_id")

    cache_root = ROOT / "runs/track_d/cache"
    models = [d.name.replace("emb_", "") for d in sorted(cache_root.glob("emb_*"))]
    if a.models:
        models = [m for m in models if m in a.models.split(",")]

    out = ROOT / "runs/track_d/leaderboard_v2.csv"
    rows = []
    if out.exists():
        rows = pd.read_csv(out).to_dict("records")
    done = {(r["method"], r["split"], r["seed"]) for r in rows}

    for m in models:
        emb = {p.stem: np.load(p) for p in (cache_root / f"emb_{m}").glob("*.npy")}
        if not emb:
            continue
        keys = df.receptor_seq.map(sha1), df.peptide_seq.map(sha1)
        ok = keys[0].isin(emb) & keys[1].isin(emb)
        d = df[ok].copy()
        X = np.vstack([np.concatenate([emb[r], emb[p], np.abs(emb[r] - emb[p]), emb[r] * emb[p]])
                       for r, p in zip(keys[0][ok], keys[1][ok])]).astype(np.float32)
        for split in ([x for x in SPLIT_COLS if x in a.splits.split(",")] if a.splits else SPLIT_COLS):
            if (m, split, a.seed) in done:
                continue
            tr, te = (d[split] == "train").to_numpy(), (d[split] == "test").to_numpy()
            if tr.sum() == 0 or te.sum() == 0:
                continue
            clf = head(a.seed).fit(X[tr], d.label.to_numpy()[tr])
            s = clf.predict_proba(X[te])[:, 1]
            y = d.label.to_numpy()[te]
            rows.append({"method": m, "split": split, "head": "mlp_locked", "seed": a.seed,
                         "n_train": int(tr.sum()), "n_test": int(te.sum()), "n": int(te.sum()),
                         "auroc": roc_auc_score(y, s), "aupr": average_precision_score(y, s),
                         "mcc": matthews_corrcoef(y, (s > 0.5).astype(int))})
            pd.DataFrame(rows).to_csv(out, index=False)
            print(f"  {m:26s} {split:16s} seed={a.seed} AUROC={rows[-1]['auroc']:.3f}", flush=True)
        del emb, X
    print(f"\nwrote {out} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
