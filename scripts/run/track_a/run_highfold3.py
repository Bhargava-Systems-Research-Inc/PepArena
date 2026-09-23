from __future__ import annotations
import argparse, json, os, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
AF3_ENV_PY = str(Path.home() / "anaconda3/envs/af3/bin/python")
AF3_CODE = "/data/alphafold3"
INPUTS = ROOT / "prediction_inputs/highfold3_msa"
PRED = ROOT / "runs/track_a/predictions/highfold3"
CYC_MATS = ROOT / "runs/track_a/_hf3_cyc_mats"
SP = Path.home() / "anaconda3/envs/af3/lib/python3.12/site-packages/nvidia"

def ld_path():
    libs = [str(p) for p in SP.glob("*/lib")]
    return ":".join(libs) + ":" + os.environ.get("LD_LIBRARY_PATH", "")

def outfile(cid):
    for name in (cid, cid.lower()):
        p = PRED / name / f"{name}_model.cif"
        if p.exists():
            return p
    return PRED / cid / f"{cid}_model.cif"

def cyc_env(cid, row):
    try:
        sys.path.insert(0, "/12TBDrive1/cycesmfold2")
        from cyc.cyclic_peptide_matrix import build_offset_matrix
        d = json.loads((INPUTS / f"{cid}.json").read_text())
        seqs = d["sequences"]
        R = sum(len(s["protein"]["sequence"]) for s in seqs[:-1])
        P = len(seqs[-1]["protein"]["sequence"])
        pep_id = seqs[-1]["protein"]["id"]
        pairs = []
        for (a1, a2) in d.get("bondedAtomPairs", []):
            if a1[0] == pep_id and a2[0] == pep_id:
                pairs.append((a1[1] - 1, a2[1] - 1))
        CYC_MATS.mkdir(parents=True, exist_ok=True)
        m = build_offset_matrix(P, True, pairs).astype(np.int32)
        np.save(CYC_MATS / f"{cid}.npy", m)
        return {"HF3_CYC_MATRIX": str(CYC_MATS / f"{cid}.npy"), "HF3_CYC_OFFSET": str(R)}
    except Exception as e:
        print(f"  {cid}: cyclic matrix FAILED ({e}) -> running as linear", flush=True)
        return {}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", choices=["linear", "cyclic"])
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    if args.only == "linear":
        man = man[~man.is_cyclic]
    elif args.only == "cyclic":
        man = man[man.is_cyclic]
    todo = [(c, cyc) for c, cyc in zip(man.complex_id, man.is_cyclic)
            if (INPUTS / f"{c}.json").exists() and not outfile(c).exists()]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(man)} complexes; {len(todo)} to fold. LD set: {'nvidia' in ld_path()}", flush=True)
    PRED.mkdir(parents=True, exist_ok=True)
    base_env = dict(os.environ, LD_LIBRARY_PATH=ld_path())
    done = 0
    for cid, is_cyc in todo:
        row = man[man.complex_id == cid].iloc[0]
        env = dict(base_env)
        if is_cyc:
            env.update(cyc_env(cid, row))
        r = subprocess.run([AF3_ENV_PY, "run_alphafold.py",
                            f"--json_path={(INPUTS / f'{cid}.json').resolve()}",
                            f"--output_dir={PRED.resolve()}",
                            "--model_dir=/data/af3_models", "--norun_data_pipeline", "--run_inference"],
                           cwd=AF3_CODE, env=env)
        ok = outfile(cid).exists()
        done += ok
        print(f"  {cid} ({'cyclic' if is_cyc else 'linear'}): {'ok' if ok else 'FAIL rc=%d' % r.returncode}", flush=True)
    print(f"DONE: {sum(outfile(c).exists() for c in man.complex_id)}/{len(man)} predicted", flush=True)

if __name__ == "__main__":
    main()
