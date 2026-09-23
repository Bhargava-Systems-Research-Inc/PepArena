from __future__ import annotations
import argparse, shutil, subprocess, time
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PROTENIX = str(Path.home() / "anaconda3/envs/protenix/bin/protenix")
INPUTS = ROOT / "prediction_inputs/protenix_msa"
PRED = ROOT / "runs/track_a/predictions/protenix"
STAGE = ROOT / "runs/track_a/_protenix_batch"
MODEL = "protenix_base_default_v1.0.0"
import os
MSA_MODE = os.environ.get("PROTENIX_MSA_MODE", "colabfold")

def predicted(cid):
    d = PRED / cid
    return d.exists() and any(d.glob("**/*.cif"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    todo = [c for c in man.complex_id if (INPUTS / f"{c}.json").exists() and not predicted(c)]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(man)} {args.tier} complexes; {len(todo)} to fold "
          f"({len(man)-len(todo)} already done). batch={args.batch}", flush=True)

    PRED.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    for i in range(0, len(todo), args.batch):
        batch = todo[i:i + args.batch]
        shutil.rmtree(STAGE, ignore_errors=True)
        STAGE.mkdir(parents=True)
        for c in batch:
            (STAGE / f"{c}.json").symlink_to((INPUTS / f"{c}.json").resolve())
        subprocess.run([PROTENIX, "pred", "-i", str(STAGE), "-o", str(PRED),
                        "-n", MODEL, "-s", "101", "--use_msa", "true",
                        "--use_default_params", "true",
                        "--trimul_kernel", "torch", "--triatt_kernel", "torch"],
                       env={**os.environ, "LAYERNORM_TYPE": "torch"},
                       check=False)
        n = min(i + args.batch, len(todo))
        rate = n / (time.time() - t0 + 1e-9)
        print(f"== batch done: {n}/{len(todo)} folded "
              f"({rate*60:.1f}/min, ETA {(len(todo)-n)/rate/60:.0f} min) ==", flush=True)
    shutil.rmtree(STAGE, ignore_errors=True)
    print(f"DONE: {sum(predicted(c) for c in man.complex_id)}/{len(man)} predicted", flush=True)

if __name__ == "__main__":
    main()
