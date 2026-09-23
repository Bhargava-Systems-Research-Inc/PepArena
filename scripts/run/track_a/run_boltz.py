from __future__ import annotations
import argparse, shutil, subprocess, sys, time
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BOLTZ = str(Path.home() / "anaconda3/envs/boltz/bin/boltz")
PRED = ROOT / "runs/track_a/predictions/boltz2"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=12)
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--inputs", default="boltz2", help="prediction_inputs subdir (boltz2 | boltz2_msa)")
    ap.add_argument("--diffusion-samples", type=int, default=1)
    ap.add_argument("--recycling-steps", type=int, default=10)
    ap.add_argument("--sampling-steps", type=int, default=200)
    ap.add_argument("--batch-name", default="_boltz_batch", help="staging dir -> boltz_results_<name>")
    args = ap.parse_args()
    inputs = ROOT / "prediction_inputs" / args.inputs
    batchdir = ROOT / "runs/track_a" / args.batch_name
    done = PRED / f"boltz_results_{args.batch_name}" / "predictions"
    use_server = args.inputs == "boltz2"

    def predicted(cid):
        return (done / cid / f"{cid}_model_0.cif").exists()

    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    todo = [c for c in man.complex_id if (inputs / f"{c}.yaml").exists() and not predicted(c)]
    print(f"{len(man)} {args.tier} complexes; {len(todo)} to fold "
          f"({len(man)-len(todo)} already done). batch={args.batch} samples={args.diffusion_samples} "
          f"msa_server={use_server}", flush=True)

    t0 = time.time()
    for i in range(0, len(todo), args.batch):
        batch = todo[i:i + args.batch]
        shutil.rmtree(batchdir, ignore_errors=True)
        batchdir.mkdir(parents=True)
        for c in batch:
            (batchdir / f"{c}.yaml").symlink_to((inputs / f"{c}.yaml").resolve())
        cmd = [BOLTZ, "predict", str(batchdir), "--out_dir", str(PRED),
               "--accelerator", "gpu", "--devices", "1", "--num_workers", "2",
               "--diffusion_samples", str(args.diffusion_samples),
               "--recycling_steps", str(args.recycling_steps),
               "--sampling_steps", str(args.sampling_steps)]
        if use_server:
            cmd.append("--use_msa_server")
        subprocess.run(cmd, check=False)
        n = min(i + args.batch, len(todo))
        rate = n / (time.time() - t0 + 1e-9)
        print(f"== batch done: {n}/{len(todo)} folded "
              f"({rate*60:.1f}/min, ETA {(len(todo)-n)/rate/60:.0f} min) ==", flush=True)
    shutil.rmtree(batchdir, ignore_errors=True)
    print(f"DONE: {sum(predicted(c) for c in man.complex_id)}/{len(man)} predicted", flush=True)

if __name__ == "__main__":
    main()
