import argparse
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge
from sklearn.ensemble import HistGradientBoostingRegressor

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "runs/track_c"
MODEL = "esm2_t12_35M_UR50D"

def embed(seqs, batch=4, max_len=1022):
    import esm
    model, alphabet = getattr(esm.pretrained, MODEL)()
    model.eval()
    bc = alphabet.get_batch_converter()
    out = {}
    seqs = list(seqs)
    for i in range(0, len(seqs), batch):
        chunk = [(f"s{j}", s[:max_len]) for j, s in enumerate(seqs[i:i + batch])]
        _, _, toks = bc(chunk)
        with torch.no_grad():
            rep = model(toks, repr_layers=[12])["representations"][12]
        for k, (_, s) in enumerate(chunk):
            n = (toks[k] != alphabet.padding_idx).sum().item()
            out[seqs[i + k]] = rep[k, 1:n - 1].mean(0).numpy()
        if (i // batch) % 20 == 0:
            print(f"  embedded {min(i + batch, len(seqs))}/{len(seqs)}", flush=True)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default=str(OUT / "affinity_targets_v2.csv"))
    a = ap.parse_args()
    d = pd.read_csv(a.targets)
    d = d[d.receptor_seq.notna() & d.peptide_seq.notna()].copy()
    print(f"rows with both sequences: {len(d)} of {len(pd.read_csv(a.targets))} "
          f"(cyclic rows carry SMILES, not a peptide sequence, and are excluded here)")

    uniq = sorted(set(d.receptor_seq) | set(d.peptide_seq))
    print(f"unique sequences to embed: {len(uniq)}", flush=True)
    emb = embed(uniq)

    X = np.vstack([np.concatenate([emb[p], emb[r], np.abs(emb[p] - emb[r]), emb[p] * emb[r]])
                   for p, r in zip(d.peptide_seq, d.receptor_seq)]).astype(np.float32)
    y = d.log_affinity.to_numpy()

    rows = []
    for split in ("split_random", "split_receptor_cluster"):
        tr = (d[split] == "train").to_numpy()
        te = (d[split] == "test").to_numpy()
        if te.sum() < 10:
            continue
        for name, reg in [("ridge", Ridge(alpha=10.0)),
                          ("gbm", HistGradientBoostingRegressor(max_iter=400, random_state=0))]:
            reg.fit(X[tr], y[tr])
            p = reg.predict(X[te])
            rows.append({"method": f"esm2+{name}", "split": split, "n_train": int(tr.sum()),
                         "n_test": int(te.sum()), "spearman": float(spearmanr(y[te], p).statistic),
                         "pearson": float(pearsonr(y[te], p)[0]),
                         "rmse": float(np.sqrt(((y[te] - p) ** 2).mean()))})
            print(f"  {split:24s} {name:5s} n={int(te.sum()):3d} "
                  f"Spearman {rows[-1]['spearman']:.3f} RMSE {rows[-1]['rmse']:.3f}", flush=True)
    out = OUT / "leaderboard_regressor_v2.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")

if __name__ == "__main__":
    main()
