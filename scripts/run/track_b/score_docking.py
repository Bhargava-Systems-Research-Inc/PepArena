#!/usr/bin/env python3
import os
import glob
import json
import argparse
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = "/12TBDrive1/mega_pep_bench"
SRC = f"{REPO}/data/external/MM_PBGBSA-CP/decoy"
SCOR = f"{REPO}/runs/track_b/scoring"
NAT = f"{SCOR}/native"
CACHE = f"{SCOR}/dockq_cache"
DOCKQ = str(Path.home() / "anaconda3/envs/dockq/bin/DockQ")
os.makedirs(CACHE, exist_ok=True)

def protein_lines(native_pdb):
    return [ln for ln in open(native_pdb) if ln.startswith(("ATOM", "HETATM")) and ln[21] == "R"]

def decoy_as_P(decoy_pdb):
    return [ln[:21] + "P" + ln[22:] for ln in open(decoy_pdb)
            if ln.startswith(("ATOM", "HETATM"))]

def run_dockq(model_lines, native_pdb):
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tf:
        tf.writelines(model_lines + ["TER\n", "END\n"])
        model = tf.name
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as jf:
            outj = jf.name
        r = subprocess.run([DOCKQ, model, native_pdb, "--capri_peptide", "--short",
                            "--json", outj], capture_output=True, text=True, timeout=120)
        d = json.load(open(outj))
        best = d.get("GlobalDockQ")
        if best is None and "best_result" in d:
            best = next(iter(d["best_result"].values())).get("DockQ")
        return best
    except (subprocess.TimeoutExpired, json.JSONDecodeError, StopIteration, KeyError):
        return None
    finally:
        for p in (model, outj):
            try:
                os.unlink(p)
            except OSError:
                pass

def score_complex(pdb):
    cache = f"{CACHE}/{pdb}.json"
    if os.path.exists(cache):
        return json.load(open(cache))
    native = f"{NAT}/{pdb}.pdb"
    prot = protein_lines(native)
    rows = []
    for decoy in sorted(glob.glob(f"{SRC}/{pdb}/*.pdb"),
                        key=lambda p: int(Path(p).stem) if Path(p).stem.isdigit() else 9999):
        rank = int(Path(decoy).stem) if Path(decoy).stem.isdigit() else None
        dq = run_dockq(prot + decoy_as_P(decoy), native)
        if dq is not None:
            rows.append({"rank": rank, "dockq": dq})
    json.dump(rows, open(cache, "w"))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--thresh", type=float, default=0.23)
    a = ap.parse_args()

    man = pd.read_csv(f"{SCOR}/manifest.csv")
    pdbs = man["pdb"].tolist()[: a.limit] if a.limit else man["pdb"].tolist()

    per, allrows = [], []
    for i, pdb in enumerate(pdbs):
        rows = score_complex(pdb)
        if not rows:
            continue
        dq = np.array([r["dockq"] for r in rows])
        ranks = np.array([r["rank"] for r in rows], dtype=float)
        top1 = next((r["dockq"] for r in rows if r["rank"] == 1), dq[0])
        rho = spearmanr(-np.asarray(ranks), dq).correlation if len(dq) > 2 else np.nan
        per.append({"pdb": pdb, "n": len(dq), "best_dockq": dq.max(),
                    "top1_dockq": top1, "frac_acceptable": (dq >= a.thresh).mean(),
                    "scoring_rho": rho})
        for r in rows:
            allrows.append({"pdb": pdb, **r})
        print(f"[{i+1}/{len(pdbs)}] {pdb}: best {dq.max():.3f} top1 {top1:.3f} "
              f"n={len(dq)} rho={rho:.2f}")

    P = pd.DataFrame(per)
    P.to_csv(f"{SCOR}/power_per_complex.csv", index=False)
    pd.DataFrame(allrows).to_parquet(f"{SCOR}/dockq_all.parquet")

    summary = {
        "method": "ADCP (decoys, MM_PBGBSA-CP dataset2)",
        "n_complexes": len(P),
        "sampling_power_acceptable": round((P.best_dockq >= a.thresh).mean(), 3),
        "sampling_power_medium": round((P.best_dockq >= 0.49).mean(), 3),
        "docking_power_top1": round((P.top1_dockq >= a.thresh).mean(), 3),
        "scoring_power_mean_rho": round(P.scoring_rho.mean(), 3),
        "mean_best_dockq": round(P.best_dockq.mean(), 3),
    }
    json.dump(summary, open(f"{SCOR}/leaderboard.json", "w"), indent=2)
    print("\n=== Track B scoring-power summary ===")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
