import argparse
import hashlib
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path("/12TBDrive1/mega_pep_bench")
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
SPLITS = ROOT / "data/splits/interaction_splits.parquet"
SPLIT = "split_receptor"
OUT = ROOT / "runs/track_d/e4_negatives"

def sha1(s):
    return hashlib.sha1(s.encode()).hexdigest()

def pair_features(emb, r, p):
    a, b = emb[sha1(r)], emb[sha1(p)]
    return np.concatenate([a, b, np.abs(a - b), a * b])

def negatives(pos, pool, rng, true_pairs):
    rows = [{"receptor_seq": r, "peptide_seq": p, "label": 1}
            for r, p in zip(pos.receptor_seq, pos.peptide_seq)]
    pool = list(pool)
    for r, p in zip(pos.receptor_seq, pos.peptide_seq):
        for _ in range(12):
            d = pool[rng.integers(len(pool))]
            if (r, d) not in true_pairs and d != p:
                rows.append({"receptor_seq": r, "peptide_seq": d, "label": 0})
                break
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--sizes", default="200,400,1000,2000")
    a = ap.parse_args()
    sizes = [int(s) for s in a.sizes.split(",")]

    df = pd.read_parquet(MAN)[["complex_id", "receptor_seq", "peptide_seq", "label"]]
    sp = pd.read_parquet(SPLITS).set_index("complex_id")
    df = df.join(sp[[SPLIT]], on="complex_id")
    pos = df[df.label == 1]
    tr_pos, te_pos = pos[pos[SPLIT] == "train"], pos[pos[SPLIT] == "test"]
    true_pairs = set(zip(pos.receptor_seq, pos.peptide_seq))
    tr_peps = sorted(set(tr_pos.peptide_seq))
    te_peps = sorted(set(te_pos.peptide_seq))
    print(f"positives: {len(tr_pos)} train / {len(te_pos)} test; "
          f"unique peptides {len(tr_peps)} train / {len(te_peps)} test", flush=True)

    cache = ROOT / f"runs/track_d/cache/emb_{a.model}"
    emb = {p.stem: np.load(p) for p in cache.glob("*.npy")}
    print(f"embeddings cached: {len(emb)}", flush=True)

    rows = []
    for N in sizes:
        for seed in range(a.seeds):
            rng = np.random.default_rng(seed)
            shared = [tr_peps[i] for i in rng.choice(len(tr_peps), min(N, len(tr_peps)), replace=False)]
            wtr = [tr_peps[i] for i in rng.choice(len(tr_peps), min(N, len(tr_peps)), replace=False)]
            wte = [te_peps[i] for i in rng.choice(len(te_peps), min(N, len(te_peps)), replace=False)]
            for design, ptr, pte in (("shared", shared, shared), ("within", wtr, wte)):
                tr = negatives(tr_pos, ptr, rng, true_pairs)
                te = negatives(te_pos, pte, rng, true_pairs)
                Xtr = np.vstack([pair_features(emb, r, p)
                                 for r, p in zip(tr.receptor_seq, tr.peptide_seq)]).astype(np.float32)
                Xte = np.vstack([pair_features(emb, r, p)
                                 for r, p in zip(te.receptor_seq, te.peptide_seq)]).astype(np.float32)
                clf = MLPClassifier(hidden_layer_sizes=(256,), max_iter=300, random_state=seed)
                clf.fit(Xtr, tr.label.to_numpy())
                s = clf.predict_proba(Xte)[:, 1]
                rows.append({"pool_size": N, "seed": seed, "design": design,
                             "n_train": len(tr), "n_test": len(te),
                             "auroc": roc_auc_score(te.label, s),
                             "aupr": average_precision_score(te.label, s)})
                print(f"  N={N:5d} seed={seed} {design:7s} AUROC={rows[-1]['auroc']:.3f}", flush=True)

    d = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT / "e4b_controlled.csv", index=False)
    piv = d.groupby(["pool_size", "design"]).auroc.agg(["mean", "std"]).round(3)
    print("\n=== AUROC by pool size and sharing (pool size now held fixed within a row pair)")
    print(piv.to_string())
    gap = (d.pivot_table(index=["pool_size", "seed"], columns="design", values="auroc")
             .assign(inflation=lambda x: x["shared"] - x["within"]))
    print("\n=== sharing effect at matched pool size")
    print(gap.groupby("pool_size").inflation.agg(["mean", "std"]).round(3).to_string())
    print(f"\nwrote {OUT / 'e4b_controlled.csv'}")

if __name__ == "__main__":
    main()
