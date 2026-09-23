#!/usr/bin/env python3
import os
import re
import glob
import gzip
import json
import argparse
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import pandas as pd

REPO = "/12TBDrive1/mega_pep_bench"
FRESH = os.environ.get("FRESH", f"{REPO}/runs/track_b/fresh")
DOCKQ = str(Path.home() / "anaconda3/envs/dockq/bin/DockQ")
CACHE = f"{FRESH}/_scored_cache"
os.makedirs(CACHE, exist_ok=True)

def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def _is_h(ln):
    return (ln[76:78].strip() or ln[12:16].strip()[:1]) == "H"

def lines(pdb, chain, keep=("ATOM", "HETATM")):
    with _open(pdb) as fh:
        return [ln[:21] + chain + ln[22:] for ln in fh
                if ln.startswith(keep) and not _is_h(ln)]

def model_lines(tool, pose, rec):
    if tool in ("adcp", "unidock", "rapidock"):
        pep = lines(pose, "P")
        pep.sort(key=lambda ln: int(ln[22:26]))
        return rec + pep
    with _open(pose) as fh:
        atoms = [ln for ln in fh if ln.startswith(("ATOM", "HETATM")) and not _is_h(ln)]
    res = defaultdict(set)
    for ln in atoms:
        res[ln[21]].add(ln[22:27])
    if len(res) < 2:
        return []
    pep = min(res, key=lambda c: len(res[c]))
    return [ln[:21] + ("P" if ln[21] == pep else "R") + ln[22:] for ln in atoms]

def ranked_poses(tool, rundir):
    if tool == "rapidock":
        hits = [f for f in glob.glob(f"{FRESH}/rapidock/out/{os.path.basename(rundir)}/rank*.pdb")
                if "reverseprocess" not in f]
        return sorted(hits, key=lambda f: int(re.sub(r"\D", "", os.path.basename(f).split("_")[0]) or 0))
    if tool in ("adcp", "unidock"):
        return sorted(glob.glob(f"{rundir}/*_out_ranked_*.pdb"),
                      key=lambda f: int(f.rsplit("_", 1)[1].split(".")[0]))
    tsv = f"{rundir}/run/5_clustfcc/clustfcc.tsv"
    if not os.path.exists(tsv):
        return []
    t = pd.read_csv(tsv, sep="\t").sort_values("rank")
    out = []
    for name in t.model_name:
        for cand in (f"{rundir}/run/4_emref/{name}.gz", f"{rundir}/run/4_emref/{name}"):
            if os.path.exists(cand):
                out.append(cand)
                break
    return out

def dockq(model_lines, native):
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tf:
        tf.writelines(model_lines + ["TER\n", "END\n"]); model = tf.name
    outj = model + ".json"
    try:
        subprocess.run([DOCKQ, model, native, "--capri_peptide", "--short", "--json", outj],
                       capture_output=True, text=True, timeout=120)
        return json.load(open(outj)).get("GlobalDockQ")
    except Exception:
        return None
    finally:
        for p in (model, outj):
            try: os.unlink(p)
            except OSError: pass

def score_one(args):
    tool, cid = args
    cache = f"{CACHE}/{tool}__{cid}.json"
    if os.path.exists(cache):
        return json.load(open(cache))
    indir, rundir = f"{FRESH}/inputs/{cid}", f"{FRESH}/{tool}/{cid}"
    poses = ranked_poses(tool, rundir)
    if not poses or not os.path.exists(f"{indir}/receptor.pdb"):
        return None
    rec = lines(f"{indir}/receptor.pdb", "R")
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as nf:
        nf.writelines(rec + ["TER\n"] + lines(f"{indir}/peptide.pdb", "P") + ["TER\n", "END\n"])
        native = nf.name
    try:
        dqs = [d for d in (dockq(model_lines(tool, p, rec), native) for p in poses) if d is not None]
    finally:
        os.unlink(native)
    if not dqs:
        return None
    rec_out = {"complex_id": cid, "n_poses": len(dqs),
               "best_dockq": max(dqs), "top1_dockq": dqs[0]}
    json.dump(rec_out, open(cache, "w"))
    print(f"{cid}: best {max(dqs):.3f} top1 {dqs[0]:.3f} (n={len(dqs)})", flush=True)
    return rec_out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=["haddock3", "adcp", "unidock", "rapidock"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    wl = pd.read_csv(f"{FRESH}/worklist.csv")
    cids = wl.complex_id.tolist()[: a.limit] if a.limit else wl.complex_id.tolist()

    with Pool(a.workers) as pool:
        per = [r for r in pool.map(score_one, [(a.tool, c) for c in cids]) if r]

    P = pd.DataFrame(per)
    P.to_csv(f"{FRESH}/{a.tool}_scored.csv", index=False)
    if len(P):
        summary = {
            "tool": a.tool, "n_complexes": int(len(P)),
            "sampling_power_acceptable": round(float((P.best_dockq >= 0.23).mean()), 3),
            "sampling_power_medium": round(float((P.best_dockq >= 0.49).mean()), 3),
            "docking_power_top1": round(float((P.top1_dockq >= 0.23).mean()), 3),
            "mean_best_dockq": round(float(P.best_dockq.mean()), 3),
            "mean_top1_dockq": round(float(P.top1_dockq.mean()), 3),
        }
        json.dump(summary, open(f"{FRESH}/{a.tool}_leaderboard.json", "w"), indent=2)
        print("\n" + json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
