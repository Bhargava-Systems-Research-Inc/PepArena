import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
MMSEQS = "mmseqs"

def search(train, test, tmp):
    tdir = Path(tmp)
    fa_tr, fa_te = tdir / "train.fasta", tdir / "test.fasta"
    for path, seqs in ((fa_tr, train), (fa_te, test)):
        with open(path, "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">s{i}\n{s}\n")
    out = tdir / "hits.tsv"
    cmd = [MMSEQS, "easy-search", str(fa_te), str(fa_tr), str(out), str(tdir / "tmp"),
           "-s", "7.5", "--max-seqs", "2000", "-e", "10000",
           "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"mmseqs easy-search failed:\n{r.stderr[-800:]}")
    cols = ["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"]
    if not out.exists() or out.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    return pd.read_csv(out, sep="\t", names=cols)

def report(label, train, test):
    print(f"\n=== {label}: {len(test)} test receptors against {len(train)} training receptors")
    with tempfile.TemporaryDirectory() as tmp:
        hits = search(train, test, tmp)
    if not len(hits):
        print("  no hits at all — every test receptor is unmatched in the training set")
        return pd.DataFrame()
    hits["score"] = hits.fident * hits[["qcov", "tcov"]].min(axis=1)
    best = hits.sort_values("score", ascending=False).groupby("query", as_index=False).first()
    idx = {f"s{i}": s for i, s in enumerate(test)}
    best["test_seq"] = best["query"].map(idx)
    unmatched = len(test) - len(best)
    ident = best.fident
    print(f"  test receptors with no hit found: {unmatched}")
    print(f"  nearest-training identity: median {ident.median():.3f}, "
          f"90th pct {ident.quantile(0.9):.3f}, max {ident.max():.3f}")
    for t in (0.3, 0.4, 0.5, 0.7, 0.9):
        n = int((ident > t).sum())
        print(f"    above {int(t*100)}% identity: {n} of {len(test)} ({100*n/len(test):.1f}%)")
    hi = best[best.fident > 0.3]
    if len(hi):
        print(f"  among those above 30%: median coverage "
              f"{hi[['qcov','tcov']].min(axis=1).median():.2f}, median alignment length "
              f"{hi.alnlen.median():.0f}")
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", default="d", choices=["d", "c"])
    a = ap.parse_args()

    if a.track == "d":
        man = pd.read_parquet(R / "data/curated/interaction_tpeppro.parquet")
        sp = pd.read_parquet(R / "data/splits/interaction_splits.parquet").set_index("complex_id")
        man = man.join(sp[["split_receptor", "split_both"]], on="complex_id")
        for col in ("split_receptor", "split_both"):
            tr = sorted(set(man[man[col] == "train"].receptor_seq))
            te = sorted(set(man[man[col] == "test"].receptor_seq))
            b = report(f"Track D, {col}", tr, te)
            if len(b):
                b.to_csv(R / f"runs/track_d/receptor_separation_{col}.csv", index=False)
    else:
        t = pd.read_csv(R / "runs/track_c/affinity_targets_v2.csv")
        tr = sorted(set(t[t.split_receptor_cluster == "train"].receptor_seq))
        te = sorted(set(t[t.split_receptor_cluster == "test"].receptor_seq))
        b = report("Track C, receptor-cluster holdout", tr, te)
        if len(b):
            b.to_csv(R / "runs/track_c/receptor_separation.csv", index=False)

if __name__ == "__main__":
    sys.exit(main())
