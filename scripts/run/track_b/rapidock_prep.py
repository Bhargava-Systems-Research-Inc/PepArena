#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path("/12TBDrive1/mega_pep_bench")
EXT = REPO / "methods/external/RAPiDock"
FRESH = REPO / "runs/track_b/fresh"
OUT = FRESH / "rapidock"
POCK = OUT / "pockets"
CAPS = {"ACE", "NH2", "NME", "NHE"}
THRESH = 10.0

def strip_caps(seq):
    for c in CAPS:
        seq = seq.replace(f"[{c}]", "")
    return seq

def main():
    cov = pd.read_csv(FRESH / "rapidock_coverage.csv")
    use = cov[cov.ok_ignoring_caps & (cov.n_res > 0)].copy()
    POCK.mkdir(parents=True, exist_ok=True)
    rows, skipped = [], []
    for _, r in use.iterrows():
        cid = r.complex_id
        indir = FRESH / "inputs" / cid
        rec, pep = indir / "receptor.pdb", indir / "peptide.pdb"
        if not rec.exists() or not pep.exists():
            skipped.append((cid, "missing prepared inputs")); continue
        pocket = POCK / f"{cid}_pocket.pdb"
        if not pocket.exists():
            cmd = [sys.executable, str(EXT / "pocket_trunction.py"),
                   "--protein_path", str(rec), "--peptide_path", str(pep),
                   "--threshold", str(THRESH), "--level", "Residue",
                   "--save_name", str(pocket)]
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if not pocket.exists() or pocket.stat().st_size == 0:
                skipped.append((cid, f"pocket failed: {(p.stderr or '').strip()[:80]}")); continue
        seq = strip_caps(str(r.seq))
        if not seq:
            skipped.append((cid, "empty sequence after cap stripping")); continue
        rows.append({"complex_name": cid, "protein_description": str(pocket),
                     "peptide_description": seq})
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "tasks.csv", index=False)
    print(f"prepared {len(d)} RAPiDock tasks -> {OUT/'tasks.csv'}")
    if skipped:
        print(f"skipped {len(skipped)}:")
        for cid, why in skipped[:15]:
            print(f"  {cid}: {why}")
    pd.DataFrame(skipped, columns=["complex_id", "reason"]).to_csv(OUT / "skipped.csv", index=False)

if __name__ == "__main__":
    main()
