from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef,
                             f1_score, accuracy_score, precision_score, recall_score)

ROOT = Path(__file__).resolve().parents[3]
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
SPLITS = ROOT / "data/splits/interaction_splits.parquet"
import hashlib
sha1 = lambda s: hashlib.sha1(s.encode()).hexdigest()

HEADS = {
    "logreg": lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
    "mlp": lambda: make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(256, 64),
                                                                 max_iter=60, early_stopping=True)),
    "rf": lambda: RandomForestClassifier(n_estimators=300, n_jobs=-1, max_depth=None),
}
SPLIT_COLS = ["split_random", "split_receptor", "split_peptide", "split_both"]

def pair_features(emb, rseq, pseq):
    r, p = emb[sha1(rseq)], emb[sha1(pseq)]
    return np.concatenate([r, p, np.abs(r - p), r * p])

def metrics(y, score):
    pred = (score >= 0.5).astype(int)
    return {"auroc": float(roc_auc_score(y, score)), "aupr": float(average_precision_score(y, score)),
            "mcc": float(matthews_corrcoef(y, pred)), "f1": float(f1_score(y, pred)),
            "accuracy": float(accuracy_score(y, pred)), "precision": float(precision_score(y, pred)),
            "recall": float(recall_score(y, pred)), "n": int(len(y)), "pos_rate": float(y.mean())}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t12_35M_UR50D")
    args = ap.parse_args()
    cache = ROOT / f"runs/track_d/cache/emb_{args.model}"

    df = pd.read_parquet(MAN)[["complex_id", "receptor_seq", "peptide_seq", "label"]]
    sp = pd.read_parquet(SPLITS).set_index("complex_id")
    df = df.join(sp[SPLIT_COLS], on="complex_id")

    print("loading cached embeddings...")
    emb = {p.stem: np.load(p) for p in cache.glob("*.npy")}
    print(f"  {len(emb)} embeddings; building features for {len(df)} pairs...")
    X = np.vstack([pair_features(emb, r, p) for r, p in zip(df.receptor_seq, df.peptide_seq)])
    y = df.label.to_numpy()
    print(f"  feature matrix {X.shape}")

    rows = []
    for split in SPLIT_COLS:
        tr = (df[split] == "train").to_numpy(); te = (df[split] == "test").to_numpy()
        for head_name, make in HEADS.items():
            run = ROOT / f"runs/track_d/{args.model}__{head_name}__{split}"
            mfile = run / "metrics.json"
            if mfile.exists():
                rows.append({"split": split, "head": head_name, **json.loads(mfile.read_text())["test"]})
                print(f"  [skip] {run.name} (done)"); continue
            run.mkdir(parents=True, exist_ok=True)
            clf = make().fit(X[tr], y[tr])
            score = clf.predict_proba(X[te])[:, 1]
            m = metrics(y[te], score)
            pred = df[te][["complex_id", "label"]].copy(); pred["score"] = score
            pred.to_parquet(run / "predictions.parquet", index=False)
            meta = {"model": args.model, "head": head_name, "split": split,
                    "n_train": int(tr.sum()), "test": m}
            mfile.write_text(json.dumps(meta, indent=2))
            rows.append({"split": split, "head": head_name, **m})
            print(f"  {run.name}: AUROC {m['auroc']:.3f} AUPR {m['aupr']:.3f} "
                  f"MCC {m['mcc']:.3f} F1 {m['f1']:.3f}")

    lb = pd.DataFrame(rows).sort_values(["split", "auroc"], ascending=[True, False])
    lb.to_csv(ROOT / "runs/track_d/leaderboard.csv", index=False)
    print("\n=== leaderboard ===")
    print(lb[["split", "head", "auroc", "aupr", "mcc", "f1", "accuracy"]].to_string(index=False))

if __name__ == "__main__":
    main()
