from __future__ import annotations
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
AA = "ACDEFGHIKLMNPQRSTVWY"
AI = {a: i for i, a in enumerate(AA)}
DI = {a + b: i for i, (a, b) in enumerate((x, y) for x in AA for y in AA)}
sha1 = lambda s: hashlib.sha1(s.encode()).hexdigest()

def featurize(seq):
    seq = "".join(c for c in seq if c in AI)
    aac = np.zeros(20, np.float32)
    for c in seq:
        aac[AI[c]] += 1
    dip = np.zeros(400, np.float32)
    for i in range(len(seq) - 1):
        dip[DI[seq[i:i + 2]]] += 1
    if len(seq):
        aac /= len(seq)
    if len(seq) > 1:
        dip /= (len(seq) - 1)
    return np.concatenate([aac, dip])

def main():
    cache = ROOT / "runs/track_d/cache/emb_composition"
    cache.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(MAN)
    seqs = sorted(set(df.receptor_seq) | set(df.peptide_seq))
    n = 0
    for s in seqs:
        f = cache / f"{sha1(s)}.npy"
        if f.exists():
            continue
        np.save(f, featurize(s)); n += 1
    print(f"composition: {n} embedded ({len(seqs)-n} cached) -> {cache}")

if __name__ == "__main__":
    main()
