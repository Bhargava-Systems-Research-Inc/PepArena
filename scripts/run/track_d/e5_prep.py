#!/usr/bin/env python3
import argparse
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

REPO = Path("/12TBDrive1/mega_pep_bench")
OUT = REPO / "runs/track_d/e5_zeroshot"
MSA, YML = OUT / "msa", OUT / "yaml"
SEED = 16
MAX_RLEN = 400

def sample(n_receptors, per_class):
    man = pd.read_parquet(REPO / "data/curated/interaction_tpeppro.parquet")
    spl = pd.read_parquet(REPO / "data/splits/interaction_splits.parquet")[
        ["complex_id", "receptor_cluster"]]
    d = man.merge(spl, on="complex_id")
    d["rlen"] = d.receptor_seq.astype(str).str.len()
    d = d[(d.rlen <= MAX_RLEN) & d.peptide_seq.astype(str).str.len().gt(0)]
    g = d.groupby("receptor_pdb_id").agg(pos=("label", "sum"), n=("label", "size"),
                                         cluster=("receptor_cluster", "first"))
    g["neg"] = g.n - g.pos
    ok = g[(g.pos >= per_class) & (g.neg >= per_class)].reset_index()
    ok = ok.sample(frac=1, random_state=SEED).drop_duplicates("cluster")
    keep = ok.sample(n=min(n_receptors, len(ok)), random_state=SEED).receptor_pdb_id
    rows = []
    for r in keep:
        sub = d[d.receptor_pdb_id == r]
        for lab in (1, 0):
            rows.append(sub[sub.label == lab].sample(n=per_class, random_state=SEED))
    return pd.concat(rows).reset_index(drop=True)

def fetch_msas(uniq):
    from boltz.data.msa.mmseqs2 import run_mmseqs2
    todo = [(p, s) for p, s in uniq.items() if not (MSA / f"{p}.a3m").exists()]
    print(f"unique receptors={len(uniq)} todo={len(todo)}", flush=True)
    chunk_n, lock, done = int(os.environ.get("MSA_CHUNK", 8)), threading.Lock(), [0]
    chunks = [todo[i:i + chunk_n] for i in range(0, len(todo), chunk_n)]

    def fetch(arg):
        idx, chunk = arg
        t0 = time.time()
        try:
            a3ms = run_mmseqs2([s for _, s in chunk], str(MSA / f"_tmp_{idx}"),
                               use_env=True, use_pairing=False)
        except Exception as e:
            return f"chunk {idx}: FAILED {type(e).__name__}: {e}"
        for (pdb, _), a3m in zip(chunk, a3ms):
            (MSA / f"{pdb}.a3m").write_text(a3m)
        with lock:
            done[0] += len(chunk)
            return f"chunk {idx}: +{len(chunk)} ({done[0]}/{len(todo)}) in {time.time()-t0:.0f}s"

    if chunks:
        with ThreadPoolExecutor(max_workers=int(os.environ.get("MSA_PAR", 4))) as ex:
            for f in as_completed([ex.submit(fetch, (i, c)) for i, c in enumerate(chunks)]):
                print(" ", f.result(), flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-receptors", type=int, default=40)
    ap.add_argument("--per-class", type=int, default=2)
    ap.add_argument("--sample-only", action="store_true",
                    help="write targets.csv and stop (for the env that has pyarrow)")
    a = ap.parse_args()
    for d in (MSA, YML):
        d.mkdir(parents=True, exist_ok=True)

    tgt = OUT / "targets.csv"
    if tgt.exists():
        s = pd.read_csv(tgt)
        print(f"reusing {tgt} ({len(s)} pairs)", flush=True)
    else:
        s = sample(a.n_receptors, a.per_class)
        s.to_csv(tgt, index=False)
        print(f"sampled {len(s)} pairs over {s.receptor_pdb_id.nunique()} receptors "
              f"({int(s.label.sum())} binders / {int((s.label == 0).sum())} non-binders)",
              flush=True)
    if a.sample_only:
        return

    fetch_msas(dict(zip(s.receptor_pdb_id, s.receptor_seq)))
    have = {p.stem for p in MSA.glob("*.a3m") if p.stat().st_size > 0}

    written = skipped = 0
    for _, r in s.iterrows():
        if r.receptor_pdb_id not in have:
            skipped += 1
            continue
        (YML / f"{r.complex_id.replace('/', '_')}.yaml").write_text(
            "version: 1\nsequences:\n  - protein:\n      id: R\n"
            f"      sequence: {r.receptor_seq}\n"
            f"      msa: {MSA}/{r.receptor_pdb_id}.a3m\n"
            "  - protein:\n      id: P\n"
            f"      sequence: {r.peptide_seq}\n"
            "      msa: empty\n")
        written += 1
    print(f"YAMLs written: {written}  skipped(no MSA): {skipped} -> {YML}", flush=True)

if __name__ == "__main__":
    main()
