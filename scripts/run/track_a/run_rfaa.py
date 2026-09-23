from __future__ import annotations
import argparse, os, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RFAA_REPO = "/data/RoseTTAFold-All-Atom"
PY = str(Path.home() / "anaconda3/envs/rfaa_bw/bin/python")
HHSEARCH = str(Path.home() / "anaconda3/envs/rfaa_bw/bin/hhsearch")
WORK = ROOT / "prediction_inputs/rfaa_work"
PRED = ROOT / "runs/track_a/predictions/rfaa"
PDB100 = "/12TBDrive1/rfaa_db/pdb100_2021Mar03/pdb100_2021Mar03"

def outpdb(cid):
    p = PRED / f"{cid}.pdb"
    return p if p.exists() else None

def gen_templates(cid):
    for chdir in sorted((PRED / cid).iterdir()):
        if not chdir.is_dir():
            continue
        a3m, hhr, atab = chdir / "t000_.msa0.a3m", chdir / "t000_.hhr", chdir / "t000_.atab"
        if a3m.exists() and not (hhr.exists() and atab.exists()):
            subprocess.run([HHSEARCH, "-b", "50", "-B", "500", "-z", "50", "-Z", "500",
                            "-mact", "0.05", "-cpu", "4", "-aliw", "100000", "-e", "100",
                            "-p", "5.0", "-d", PDB100, "-i", str(a3m),
                            "-o", str(hhr), "-atab", str(atab), "-v", "0"],
                           capture_output=True, timeout=1200)

def run_one(cid):
    if outpdb(cid):
        return True
    gen_templates(cid)
    env = dict(os.environ)
    r = subprocess.run([PY, "-m", "rf2aa.run_inference",
                        "--config-dir", str(WORK.resolve()), "--config-name", cid],
                       cwd=RFAA_REPO, env=env, capture_output=True, timeout=3600, text=True)
    ok = outpdb(cid) is not None
    if not ok:
        log = WORK / f"{cid}.err.log"
        log.write_text(f"exit={r.returncode}\n\n--- stderr tail ---\n"
                       + "\n".join((r.stderr or "").splitlines()[-40:])
                       + "\n\n--- stdout tail ---\n"
                       + "\n".join((r.stdout or "").splitlines()[-20:]))
    return ok

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    todo = [c for c in man.complex_id if (WORK / f"{c}.yaml").exists() and not outpdb(c)]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(man)} complexes; {len(todo)} to fold. workers={args.workers}", flush=True)
    done, t0 = 0, time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_one, c): c for c in todo}
        for i, f in enumerate(as_completed(futs), 1):
            cid = futs[f]
            try:
                ok = f.result()
            except Exception:
                ok = False
            done += ok
            r = i / (time.time() - t0 + 1e-9)
            print(f"  [{i}/{len(todo)}] {cid}: {'ok' if ok else 'FAIL'} "
                  f"({r*60:.1f}/min, ETA {(len(todo)-i)/r/60:.0f} min)", flush=True)
    print(f"DONE: {sum(bool(outpdb(c)) for c in man.complex_id)}/{len(man)} predicted", flush=True)

if __name__ == "__main__":
    main()
