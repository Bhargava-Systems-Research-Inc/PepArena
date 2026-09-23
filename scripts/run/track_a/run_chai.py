from __future__ import annotations
import argparse, os
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "prediction_inputs/chai1"
PRED = ROOT / "runs/track_a/predictions/chai1"
MSA_DIR = ROOT / "runs/track_a/_chai_msa"

def predicted(cid):
    d = PRED / cid
    return d.exists() and any(d.glob("*.cif"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    import torch
    from chai_lab.chai1 import run_inference

    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    todo = [c for c in man.complex_id if (INPUTS / f"{c}.fasta").exists() and not predicted(c)]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(man)} {args.tier} complexes; {len(todo)} to fold "
          f"({len(man)-len(todo)} already done).", flush=True)

    dev = torch.device("cuda:0")
    done, fail = 0, 0
    for cid in todo:
        out = PRED / cid
        out.mkdir(parents=True, exist_ok=True)
        try:
            run_inference(fasta_file=INPUTS / f"{cid}.fasta", output_dir=out,
                          num_diffn_samples=5, num_diffn_timesteps=200, num_trunk_recycles=10,
                          use_esm_embeddings=True, use_msa_server=False,
                          msa_directory=MSA_DIR, seed=1, device=dev)
            done += 1
        except Exception as e:
            fail += 1
            print(f"  {cid}: FAIL {type(e).__name__}: {str(e)[:200]}", flush=True)
        if (done + fail) % 10 == 0:
            print(f"  {done} folded, {fail} failed ({done+fail}/{len(todo)})", flush=True)
    print(f"DONE: {sum(predicted(c) for c in man.complex_id)}/{len(man)} predicted "
          f"({done} new, {fail} failed)", flush=True)

if __name__ == "__main__":
    main()
