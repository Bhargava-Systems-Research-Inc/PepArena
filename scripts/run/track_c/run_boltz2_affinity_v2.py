import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/12TBDrive1/mega_pep_bench")
OUT = ROOT / "runs/track_c/boltz2_affinity_v2"

def yaml_for(recseq, smiles):
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: R\n"
        f"      sequence: {recseq}\n"
        "      msa: empty\n"
        "  - ligand:\n"
        "      id: L\n"
        f"      smiles: '{smiles}'\n"
        "properties:\n"
        "  - affinity:\n"
        "      binder: L\n"
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--targets", default=str(ROOT / "runs/track_c/affinity_targets_v2.csv"))
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.targets)
    done = 0
    for row in df.itertuples(index=False):
        if a.limit and done >= a.limit:
            break
        d = OUT / row.complex_id
        pred = d / "predictions" / row.complex_id / f"affinity_{row.complex_id}.json"
        if pred.exists() or (d / "DONE").exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / "input.yaml").write_text(yaml_for(row.receptor_seq, row.smiles))
        r = subprocess.run([str(Path.home() / "anaconda3/envs/boltz/bin/boltz"), "predict", str(d / "input.yaml"), "--out_dir", str(d),
                            "--diffusion_samples", "1", "--sampling_steps", "100"],
                           capture_output=True, text=True)
        (d / "run.log").write_text(r.stdout[-4000:] + "\n--- stderr ---\n" + r.stderr[-4000:])
        (d / "DONE").write_text("ok\n")
        done += 1
        if done % 10 == 0:
            print(f"  {done} predicted", flush=True)
    print(f"ran {done} Boltz-2 affinity predictions -> {OUT}")

if __name__ == "__main__":
    main()
