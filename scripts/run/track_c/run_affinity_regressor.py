from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import esm
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/track_c"

def sha1(s): return hashlib.sha1(s.encode()).hexdigest()

def embed_all(seqs, model_name, max_len=1022, threads=24):
    cache = ROOT / f"runs/track_d/cache/emb_{model_name}"
    cache.mkdir(parents=True, exist_ok=True)
    todo = sorted({s for s in seqs if not (cache / f"{sha1(s)}.npy").exists()})
    if todo:
        torch.set_num_threads(threads)
        model, alphabet = getattr(esm.pretrained, model_name)()
        model = model.eval()
        bc = alphabet.get_batch_converter()
        layer = model.num_layers
        print(f"embedding {len(todo)} new seqs (model {model_name})", flush=True)
        with torch.no_grad():
            for i, s in enumerate(todo):
                _, _, toks = bc([("x", s[:max_len])])
                rep = model(toks, repr_layers=[layer])["representations"][layer]
                v = rep[0, 1:len(s[:max_len]) + 1].mean(0).numpy().astype("float32")
                np.save(cache / f"{sha1(s)}.npy", v)
                if (i + 1) % 200 == 0:
                    print(f"  {i+1}/{len(todo)}", flush=True)
    return lambda s: np.load(cache / f"{sha1(s)}.npy")

def metrics(y, p):
    return {"spearman": round(spearmanr(y, p).correlation, 3),
            "pearson": round(pearsonr(y, p)[0], 3),
            "rmse": round(mean_squared_error(y, p) ** 0.5, 3),
            "mae": round(mean_absolute_error(y, p), 3)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t12_35M_UR50D")
    a = ap.parse_args()
    df = pd.read_parquet(OUT / "affinity_manifest.parquet")

    recseq = {}
    rs = OUT / "receptor_seqs.csv"
    if rs.exists():
        r = pd.read_csv(rs)
        recseq = dict(zip(r.receptor_pdb_id, r.receptor_seq))
        print(f"using receptor seqs for {len(recseq)} PDBs (peptide+receptor features)")
    else:
        print("no receptor_seqs.csv -> peptide-only features (run fetch_receptor_seqs.py to upgrade)")

    seqs = list(df.peptide_seq) + [recseq[p] for p in df.receptor_pdb_id if p in recseq]
    get = embed_all([s for s in seqs if isinstance(s, str) and s], a.model)

    X, keep = [], []
    for row in df.itertuples(index=False):
        p = get(row.peptide_seq)
        if recseq:
            r = get(recseq[row.receptor_pdb_id]) if row.receptor_pdb_id in recseq else np.zeros_like(p)
            feat = np.concatenate([p, r, np.abs(p - r), p * r])
        else:
            feat = p
        X.append(feat); keep.append(True)
    X = np.array(X); y = df.log_affinity.values

    rows = []
    for split in ("split_random", "split_receptor"):
        tr = (df[split] == "train").values; te = (df[split] == "test").values
        for name, reg in [("ridge", Ridge(alpha=10.0)),
                          ("gbm", HistGradientBoostingRegressor(max_iter=400))]:
            reg.fit(X[tr], y[tr])
            m = metrics(y[te], reg.predict(X[te]))
            rows.append({"method": f"esm2+{name}", "features": "pep+rec" if recseq else "pep",
                         "split": split, "n_test": int(te.sum()), **m})
            print(f"{split:15s} {name:5s}: Spearman {m['spearman']} | Pearson {m['pearson']} | "
                  f"RMSE {m['rmse']} (n_test={int(te.sum())})")

    lb = pd.DataFrame(rows)
    lb.to_csv(OUT / "leaderboard_regressor.csv", index=False)
    print("\nwrote", OUT / "leaderboard_regressor.csv")

if __name__ == "__main__":
    main()
