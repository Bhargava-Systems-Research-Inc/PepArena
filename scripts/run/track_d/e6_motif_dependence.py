#!/usr/bin/env python3
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path("/12TBDrive1/mega_pep_bench")
OUT = ROOT / "runs/track_d/e6_motif"
MIN_PAIRS = 10
MIN_BINDERS = 3
AA = "ACDEFGHIKLMNPQRSTVWY"

BEST = {"split_random": "prot_t5_xl__rf__split_random",
        "split_receptor": "esm2_t36_3B_UR50D__mlp__split_receptor"}

def self_similarity(seqs):
    pairs = list(combinations(seqs, 2))
    if not pairs:
        return np.nan
    return float(np.mean([1 - Levenshtein.distance(a, b) / max(len(a), len(b))
                          for a, b in pairs]))

def composition_ic(seqs, background):
    joined = "".join(seqs)
    counts = np.array([joined.count(a) for a in AA], dtype=float)
    if counts.sum() == 0:
        return np.nan
    p = (counts + 1) / (counts.sum() + len(AA))
    return float(np.sum(p * np.log2(p / background)))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    man = pd.read_parquet(ROOT / "data/curated/interaction_tpeppro.parquet")[
        ["complex_id", "peptide_seq"]]
    splits = pd.read_parquet(ROOT / "data/splits/interaction_splits.parquet")[
        ["complex_id", "receptor_cluster"]]

    allb = "".join(man.peptide_seq.astype(str))
    bg_counts = np.array([allb.count(a) for a in AA], dtype=float)
    background = (bg_counts + 1) / (bg_counts.sum() + len(AA))

    summary, frames = {}, []
    for split, run in BEST.items():
        pred = pd.read_parquet(ROOT / "runs/track_d" / run / "predictions.parquet")
        d = pred.merge(splits, on="complex_id").merge(man, on="complex_id")
        rows = []
        for cl, g in d.groupby("receptor_cluster"):
            if len(g) < MIN_PAIRS or g.label.nunique() < 2:
                continue
            binders = [s for s in g.loc[g.label == 1, "peptide_seq"].astype(str) if s]
            if len(binders) < MIN_BINDERS:
                continue
            rows.append({
                "split": split, "receptor_cluster": cl, "n_pairs": len(g),
                "n_binders": len(binders),
                "auroc": float(roc_auc_score(g.label, g.score)),
                "binder_self_similarity": self_similarity(binders),
                "composition_ic": composition_ic(binders, background),
            })
        c = pd.DataFrame(rows)
        if c.empty:
            summary[split] = {"run": run, "n_clusters": 0}
            continue
        frames.append(c)

        res = {"run": run, "n_clusters": int(len(c)),
               "n_pairs_covered": int(c.n_pairs.sum()),
               "median_cluster_auroc": round(float(c.auroc.median()), 3)}
        for metric in ["binder_self_similarity", "composition_ic"]:
            v = c[[metric, "auroc"]].dropna()
            rho, p = spearmanr(v[metric], v.auroc)
            res[f"rho_auroc_vs_{metric}"] = round(float(rho), 3)
            res[f"p_{metric}"] = round(float(p), 4)
            hi = c[c[metric] >= c[metric].median()]
            lo = c[c[metric] < c[metric].median()]
            res[f"auroc_high_{metric}"] = round(float(hi.auroc.mean()), 3)
            res[f"auroc_low_{metric}"] = round(float(lo.auroc.mean()), 3)
        v = c[["n_binders", "auroc"]].dropna()
        res["rho_auroc_vs_n_binders"] = round(float(spearmanr(v.n_binders, v.auroc)[0]), 3)
        summary[split] = res
        print(f"\n{split} ({run})")
        print(json.dumps(res, indent=2))

    if frames:
        pd.concat(frames).to_csv(OUT / "motif_per_cluster.csv", index=False)
    (OUT / "e6_motif.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {OUT}/motif_per_cluster.csv and e6_motif.json")

if __name__ == "__main__":
    main()
