from __future__ import annotations
import argparse, glob, json, os, shutil, subprocess, time
from pathlib import Path
import pandas as pd

ROOT = Path("/12TBDrive1/mega_pep_bench")
NC = ROOT / "ncaa_ccd"
PRED = ROOT / "runs/track_a/predictions_ncaa"
HOME = Path.home()

def cids():
    return list(pd.read_csv(NC / "manifest.csv").complex_id)

def protenix_done(c):
    return bool(glob.glob(str(PRED / "protenix" / c / "seed_*/predictions/*.cif")))

def run_protenix(todo, batch):
    PY = str(HOME / "anaconda3/envs/protenix/bin/protenix")
    STAGE = Path("/12TBDrive1/tmp_ncaa_ptx")
    t0 = time.time()
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        shutil.rmtree(STAGE, ignore_errors=True); STAGE.mkdir(parents=True)
        for c in chunk:
            (STAGE / f"{c}.json").symlink_to((NC / "protenix_msa" / f"{c}.json").resolve())
        subprocess.run([PY, "pred", "-i", str(STAGE), "-o", str(PRED / "protenix"),
                        "-n", "protenix_base_default_v1.0.0", "-s", "101", "--use_msa", "true",
                        "--use_default_params", "true", "--trimul_kernel", "torch",
                        "--triatt_kernel", "torch"],
                       env={**os.environ, "LAYERNORM_TYPE": "torch"}, check=False)
        n = min(i + batch, len(todo)); r = n / (time.time() - t0 + 1e-9)
        print(f"== protenix {n}/{len(todo)} ({r*60:.1f}/min, ETA {(len(todo)-n)/r/60:.0f}m) ==", flush=True)

def boltz_done(c):
    return bool(glob.glob(str(PRED / "boltz2/boltz_results__ncaa/predictions" / c / "*.cif")))

def run_boltz(todo, batch):
    PY = str(HOME / "anaconda3/envs/boltz/bin/boltz")
    STAGE = Path("/12TBDrive1/tmp_ncaa_boltz/_ncaa")
    t0 = time.time()
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        shutil.rmtree(STAGE.parent, ignore_errors=True); STAGE.mkdir(parents=True)
        for c in chunk:
            (STAGE / f"{c}.yaml").symlink_to((NC / "boltz2_msa" / f"{c}.yaml").resolve())
        subprocess.run([PY, "predict", str(STAGE), "--out_dir", str(PRED / "boltz2"),
                        "--accelerator", "gpu", "--devices", "1", "--num_workers", "2",
                        "--diffusion_samples", "5", "--recycling_steps", "10",
                        "--sampling_steps", "200", "--output_format", "mmcif"], check=False)
        n = min(i + batch, len(todo)); r = n / (time.time() - t0 + 1e-9)
        print(f"== boltz {n}/{len(todo)} ({r*60:.1f}/min, ETA {(len(todo)-n)/r/60:.0f}m) ==", flush=True)

def _af3_done(sub, c):
    for name in (c, c.lower()):
        if (PRED / sub / name / f"{name}_model.cif").exists():
            return True
    return False

def af3_done(c):
    return _af3_done("af3", c)

def hf3_done(c):
    return _af3_done("highfold3", c)

def _af3_env():
    PY = str(HOME / "anaconda3/envs/af3/bin/python")
    SP = HOME / "anaconda3/envs/af3/lib/python3.12/site-packages/nvidia"
    libs = ":".join(str(p) for p in SP.glob("*/lib"))
    return PY, dict(os.environ, LD_LIBRARY_PATH=libs + ":" + os.environ.get("LD_LIBRARY_PATH", ""))

def _run_af3_one(c, sub, extra_env, PY, base_env):
    env = dict(base_env); env.update(extra_env)
    subprocess.run([PY, "run_alphafold.py",
                    f"--json_path={(NC / 'af3_msa' / f'{c}.json').resolve()}",
                    f"--output_dir={(PRED / sub).resolve()}",
                    "--model_dir=/data/af3_models", "--norun_data_pipeline", "--run_inference"],
                   cwd="/data/alphafold3", env=env, check=False)

def run_af3(todo, batch):
    import shutil
    PY, base_env = _af3_env()
    STAGE = Path("/12TBDrive1/tmp_ncaa_af3")
    shutil.rmtree(STAGE, ignore_errors=True); STAGE.mkdir(parents=True)
    for c in todo:
        (STAGE / f"{c}.json").symlink_to((NC / "af3_msa" / f"{c}.json").resolve())
    subprocess.run([PY, "run_alphafold.py",
                    f"--input_dir={STAGE.resolve()}",
                    f"--output_dir={(PRED / 'af3').resolve()}",
                    "--model_dir=/data/af3_models", "--norun_data_pipeline", "--run_inference"],
                   cwd="/data/alphafold3", env=base_env, check=False)
    print(f"  af3: {sum(af3_done(c) for c in todo)}/{len(todo)} folded", flush=True)

class Inexpressible(Exception):
    pass

def cyc_env(cid):
    import numpy as np
    import sys
    sys.path.insert(0, "/12TBDrive1/cycesmfold2")
    try:
        from cyc.cyclic_peptide_matrix import build_offset_matrix
    except ImportError as e:
        raise SystemExit(f"ABORT: CycPOEM module unavailable ({e}). HighFold3 without its "
                         f"cyclic encoding is AlphaFold3, not HighFold3.")
    d = json.loads((NC / "af3_msa" / f"{cid}.json").read_text())
    seqs = d["sequences"]
    R = sum(len(s["protein"]["sequence"]) for s in seqs[:-1])
    P = len(seqs[-1]["protein"]["sequence"])
    pep_id = seqs[-1]["protein"]["id"]
    pairs = [(a1[1] - 1, a2[1] - 1) for a1, a2 in d.get("bondedAtomPairs", [])
             if a1[0] == pep_id and a2[0] == pep_id]
    if not pairs:
        raise Inexpressible(f"{cid}: ring not expressible in the AF3 JSON "
                            f"(no intra-peptide bondedAtomPairs)")
    mats = PRED / "_hf3_cyc_mats"; mats.mkdir(parents=True, exist_ok=True)
    m = build_offset_matrix(P, True, pairs).astype(np.int32)
    lin = build_offset_matrix(P, False, []).astype(np.int32)
    if m.shape != lin.shape or (m == lin).all():
        raise SystemExit(f"ABORT: {cid} CycPOEM matrix is identical to the linear encoding "
                         f"(pairs={pairs}).")
    np.save(mats / f"{cid}.npy", m)
    return {"HF3_CYC_MATRIX": str(mats / f"{cid}.npy"), "HF3_CYC_OFFSET": str(R)}

def run_highfold3(todo, batch):
    import shutil
    PY, base_env = _af3_env()
    man = pd.read_csv(NC / "manifest.csv").set_index("complex_id")
    t0 = time.time()
    skipped = []
    for k, c in enumerate(todo, 1):
        try:
            env = cyc_env(c)
        except Inexpressible as e:
            skipped.append({"complex_id": c, "reason": str(e)})
            print(f"  hf3 [{k}/{len(todo)}] SKIP {e}", flush=True)
            continue
        _run_af3_one(c, "highfold3", env, PY, base_env)
        r = k / (time.time() - t0 + 1e-9)
        print(f"  hf3 [{k}/{len(todo)}] {c}: {'ok' if hf3_done(c) else 'FAIL'} "
              f"({r*60:.1f}/min, ETA {(len(todo)-k)/r/60:.0f}m)", flush=True)
    for c in man.index:
        if man.loc[c, "cyclization_type"] == "linear" and af3_done(c) and not hf3_done(c):
            src = next((PRED / "af3" / n for n in (c, c.lower()) if (PRED / "af3" / n).exists()), None)
            if src:
                shutil.copytree(src, PRED / "highfold3" / src.name, dirs_exist_ok=True)
    if skipped:
        out = ROOT / "runs/track_a/hf3_ncaa_inexpressible.csv"
        pd.DataFrame(skipped).to_csv(out, index=False)
        print(f"  hf3: {len(skipped)} cyclic complexes could not be encoded -> {out}", flush=True)
    print(f"  hf3: shared linear from af3 ({sum(hf3_done(c) for c in man.index)}/{len(man)} total)", flush=True)

MODELS = {"protenix": (protenix_done, run_protenix, 12),
          "boltz2": (boltz_done, run_boltz, 12),
          "af3": (af3_done, run_af3, 1),
          "highfold3": (hf3_done, run_highfold3, 1)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=list(MODELS))
    ap.add_argument("--batch", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    done_fn, run_fn, defbatch = MODELS[args.model]
    (PRED / args.model).mkdir(parents=True, exist_ok=True)
    indir = {"protenix": ("protenix_msa", "json"), "boltz2": ("boltz2_msa", "yaml"),
             "af3": ("af3_msa", "json"), "highfold3": ("af3_msa", "json")}[args.model]
    pool = cids()
    if args.model == "highfold3":
        man = pd.read_csv(NC / "manifest.csv")
        pool = list(man[man.cyclization_type != "linear"].complex_id)
    todo = [c for c in pool if (NC / indir[0] / f"{c}.{indir[1]}").exists() and not done_fn(c)]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{args.model}: {len(cids())} ncAA complexes, {len(todo)} to fold", flush=True)
    if todo:
        run_fn(todo, args.batch or defbatch)
    print(f"DONE {args.model}: {sum(done_fn(c) for c in cids())}/{len(cids())} predicted", flush=True)

if __name__ == "__main__":
    main()
