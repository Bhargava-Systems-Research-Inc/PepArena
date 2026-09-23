import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
IDENT, COV = 0.3, 0.8
COLS = ["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"]
MMSEQS = shutil.which("mmseqs") or os.path.expanduser(
    "~/anaconda3/envs/pepgym-curate/bin/mmseqs")

def search(train, test):
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        for path, seqs in ((t / "train.fasta", train), (t / "test.fasta", test)):
            path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
        out = t / "hits.tsv"
        r = subprocess.run([MMSEQS, "easy-search", str(t / "test.fasta"), str(t / "train.fasta"),
                            str(out), str(t / "tmp"), "-s", "7.5", "--max-seqs", "2000",
                            "-e", "10000", "--format-output", ",".join(COLS)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(f"mmseqs easy-search failed:\n{r.stderr[-800:]}")
        if not out.exists() or out.stat().st_size == 0:
            return pd.DataFrame(columns=COLS + ["test_seq"])
        hits = pd.read_csv(out, sep="\t", names=COLS)
    hits["mincov"] = hits[["qcov", "tcov"]].min(axis=1)
    hits["test_seq"] = hits["query"].map({f"s{i}": s for i, s in enumerate(test)})
    return hits

def main():
    sp = pd.read_parquet(R / "data/splits/interaction_splits.parquet")
    man = pd.read_parquet(R / "data/curated/interaction_tpeppro.parquet")[
        ["complex_id", "receptor_seq"]]
    d = man.merge(sp, on="complex_id")
    col = "split_receptor_strict"
    rounds, moved_total = [], 0
    for it in range(1, 6):
        train = sorted(set(d.loc[d[col] == "train", "receptor_seq"]))
        test = sorted(set(d.loc[d[col] == "test", "receptor_seq"]))
        hits = search(train, test)
        bad = set(hits.loc[(hits.fident > IDENT) & (hits.mincov > COV), "test_seq"])
        n_pairs = int((d[col] == "test").sum())
        print(f"round {it}: {len(test)} test receptors, {n_pairs} test pairs, "
              f"{len(bad)} receptors violate against the current training set", flush=True)
        rounds.append({"round": it, "test_receptors": len(test), "test_pairs": n_pairs,
                       "violating_receptors": len(bad)})
        if not bad:
            break
        mv = d.receptor_seq.isin(bad) & (d[col] == "test")
        moved_total += int(mv.sum())
        d.loc[mv, col] = "train"
    else:
        raise SystemExit("strict split did not converge in 5 rounds")

    if moved_total:
        out = sp.drop(columns=[col]).merge(d[["complex_id", col]], on="complex_id", how="left")
        out[col] = out[col].fillna("none")
        out.to_parquet(R / "data/splits/interaction_splits.parquet", index=False)
        print(f"\nmoved a further {moved_total} pairs to training; split rewritten")

    best = hits.sort_values("fident", ascending=False).groupby("query", as_index=False).first()
    rep = {"criterion": {"identity": IDENT, "mutual_coverage": COV},
           "rounds": rounds, "extra_pairs_moved": moved_total,
           "final_test_pairs": int((d[col] == "test").sum()),
           "final_train_pairs": int((d[col] == "train").sum()),
           "final_test_receptors": int(d.loc[d[col] == "test", "receptor_seq"].nunique()),
           "max_identity_in_test": round(float(best.fident.max()), 3) if len(best) else None,
           "median_identity_in_test": round(float(best.fident.median()), 3) if len(best) else None,
           "violating_receptors_final": 0}
    json.dump(rep, open(R / "runs/track_d/strict_split_report.json", "w"), indent=1)
    print(json.dumps(rep, indent=1))

if __name__ == "__main__":
    main()
