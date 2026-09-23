import argparse
import hashlib
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

ROOT = Path("/12TBDrive1/mega_pep_bench")
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
SPLITS = ROOT / "data/splits/interaction_splits.parquet"
IDENT, COV = 0.3, 0.8

def sha1(s):
    return hashlib.sha1(s.encode()).hexdigest()

def head(seed):
    return make_pipeline(StandardScaler(),
                         MLPClassifier(hidden_layer_sizes=(256,), activation="relu", solver="adam",
                                       learning_rate_init=1e-3, early_stopping=True,
                                       validation_fraction=0.1, max_iter=300, random_state=seed))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="esm2_t33_650M_UR50D,esm2_t36_3B_UR50D,esmc_300m,prot_t5_xl,composition")
    ap.add_argument("--split", default="split_receptor")
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()

    sep = pd.read_csv(ROOT / f"runs/track_d/receptor_separation_{a.split}.csv")
    sep["mincov"] = sep[["qcov", "tcov"]].min(axis=1)
    leaky = set(sep.loc[(sep.fident > IDENT) & (sep["mincov"] > COV), "test_seq"])
    print(f"test receptors with a training neighbour over {IDENT:.0%} identity and {COV:.0%} "
          f"coverage: {len(leaky)}")

    df = pd.read_parquet(MAN)[["complex_id", "receptor_seq", "peptide_seq", "label"]]
    sp = pd.read_parquet(SPLITS).set_index("complex_id")
    df = df.join(sp[[a.split]], on="complex_id")

    rows = []
    for m in a.models.split(","):
        cache = ROOT / f"runs/track_d/cache/emb_{m}"
        emb = {p.stem: np.load(p) for p in cache.glob("*.npy")}
        if not emb:
            print(f"  {m}: no cached embeddings"); continue
        rk, pk = df.receptor_seq.map(sha1), df.peptide_seq.map(sha1)
        ok = rk.isin(emb) & pk.isin(emb)
        d = df[ok].copy()
        X = np.vstack([np.concatenate([emb[r], emb[p], np.abs(emb[r] - emb[p]), emb[r] * emb[p]])
                       for r, p in zip(rk[ok], pk[ok])]).astype(np.float32)
        tr = (d[a.split] == "train").to_numpy()
        te = (d[a.split] == "test").to_numpy()
        strict = te & ~d.receptor_seq.isin(leaky).to_numpy()
        for seed in range(a.seeds):
            clf = head(seed).fit(X[tr], d.label.to_numpy()[tr])
            s = clf.predict_proba(X)[:, 1]
            y = d.label.to_numpy()
            rows.append({"method": m, "seed": seed,
                         "n_test": int(te.sum()), "auroc_all": roc_auc_score(y[te], s[te]),
                         "n_strict": int(strict.sum()),
                         "auroc_strict": roc_auc_score(y[strict], s[strict]),
                         "n_leaky": int((te & ~strict).sum()),
                         "auroc_leaky": roc_auc_score(y[te & ~strict], s[te & ~strict])
                         if len(set(y[te & ~strict])) > 1 else float("nan")})
            print(f"  {m:24s} seed={seed} all={rows[-1]['auroc_all']:.3f} "
                  f"strict={rows[-1]['auroc_strict']:.3f} "
                  f"leaky={rows[-1]['auroc_leaky']:.3f}", flush=True)
        del emb, X

    out = pd.DataFrame(rows)
    out.to_csv(ROOT / f"runs/track_d/strict_receptor_eval_{a.split}.csv", index=False)
    if len(out):
        g = out.groupby("method")[["auroc_all", "auroc_strict", "auroc_leaky"]].mean().round(3)
        print("\n=== mean over seeds")
        print(g.to_string())
        print(f"\nmean drop from all test pairs to strictly separated ones: "
              f"{(out.auroc_all - out.auroc_strict).mean():+.3f}")

if __name__ == "__main__":
    main()
