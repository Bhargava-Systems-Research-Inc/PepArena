import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "runs/track_c/boltz2_affinity_v2_batch"
BOLTZ = str(Path.home() / "anaconda3/envs/boltz/bin/boltz")

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
    ap.add_argument("--targets", default=str(R / "runs/track_c/affinity_targets_v2_boltzable.csv"))
    a = ap.parse_args()
    d = pd.read_csv(a.targets)
    stage = OUT / "yaml"
    shutil.rmtree(stage, ignore_errors=True)
    stage.mkdir(parents=True, exist_ok=True)
    for row in d.itertuples(index=False):
        (stage / f"{row.complex_id}.yaml").write_text(yaml_for(row.receptor_seq, row.smiles))
    print(f"staged {len(d)} YAMLs -> {stage}", flush=True)
    r = subprocess.run([BOLTZ, "predict", str(stage), "--out_dir", str(OUT),
                        "--diffusion_samples", "1", "--sampling_steps", "100"],
                       capture_output=True, text=True)
    (OUT / "run.log").write_text(r.stdout[-20000:] + "\n--- stderr ---\n" + r.stderr[-20000:])
    n = len(list(OUT.rglob("affinity*.json")))
    print(f"affinity predictions returned: {n} of {len(d)}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
