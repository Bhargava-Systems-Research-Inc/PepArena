from __future__ import annotations
import argparse, glob
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CSV = ROOT / "runs/track_a/esmfold2_inputs_canonical.csv"
OUTROOT = ROOT / "runs/track_a/predictions/esmfold2/esmfold2-2026-05_loops20_steps100"

def has_struct(cid):
    return bool(glob.glob(str(OUTROOT / cid / f"{cid}_model_0.cif")) or
               glob.glob(str(OUTROOT / cid / f"{cid}_model_0.pdb")))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="fold all, not just missing")
    ap.add_argument("--min-len", type=int, default=0)
    args = ap.parse_args()
    df = pd.read_csv(CSV)
    chain_cols = [c for c in df.columns if c.startswith("chain_") and c.endswith("_sequence")]
    rows = []
    for _, r in df.iterrows():
        cid = str(r["pdb_id"])
        if not args.all and has_struct(cid):
            continue
        seqs = [str(r[c]) for c in chain_cols if isinstance(r[c], str) and r[c].strip()]
        if not seqs or int(r["total_length"]) < args.min_len:
            continue
        rows.append((cid, ":".join(seqs)))
    print(f"{len(rows)} complexes to fold locally (full 6.6B ESMFold2)", flush=True)

    import os
    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    import torch
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    ESMFOLD2_REVISION = "1ebf0e3481a5184eb6171d40615c79e384b48796"
    model = (ESMFold2Model
             .from_pretrained("biohub/ESMFold2", revision=ESMFOLD2_REVISION)
             .cuda().eval())
    done = 0
    for cid, seq in rows:
        try:
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
                pdb = model.infer_protein_as_pdb(seq)
            d = OUTROOT / cid
            d.mkdir(parents=True, exist_ok=True)
            (d / f"{cid}_model_0.pdb").write_text(pdb)
            done += 1
            print(f"  {cid} (len {len(seq)}): ok", flush=True)
        except Exception as e:
            print(f"  {cid} (len {len(seq)}): FAIL {type(e).__name__}: {str(e)[:120]}", flush=True)
    print(f"DONE: {done}/{len(rows)} folded locally", flush=True)

if __name__ == "__main__":
    main()
