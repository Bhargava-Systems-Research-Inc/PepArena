import os
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

F = Path("/12TBDrive1/mega_pep_bench/runs/track_b/fresh")
ADFR = Path(os.environ.get("ADFR_HOME", Path.home() / "ADFRsuite/install/bin"))
cov = pd.read_csv(F / "coverage.csv")
todo = cov.loc[cov.unidock == "ok", "complex_id"].astype(str).tolist()
print(f"testing receptor preparation for {len(todo)} complexes Uni-Dock accepts")

rows = []
for i, cid in enumerate(todo, 1):
    rec = F / "inputs" / cid / "receptor.pdb"
    n_atoms = sum(1 for l in open(rec) if l.startswith(("ATOM", "HETATM"))) if rec.exists() else -1
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        (t / "rec_in.pdb").write_bytes(rec.read_bytes()) if rec.exists() else None
        r = subprocess.run([str(ADFR / "prepare_receptor"), "-r", "rec_in.pdb",
                            "-o", "receptor.pdbqt", "-A", "checkhydrogens"],
                           cwd=t, capture_output=True, text=True)
        out = t / "receptor.pdbqt"
        ok = out.exists() and out.stat().st_size > 0
        err = "" if ok else (r.stderr.strip().splitlines() or ["no output"])[-1][:120]
    rows.append({"complex_id": cid, "rec_atoms": n_atoms, "receptor_prepares": ok, "error": err})
    if i % 20 == 0:
        print(f"  {i}/{len(todo)}", flush=True)

d = pd.DataFrame(rows)
d.to_csv(F / "unidock_receptor_prep.csv", index=False)
bad = d[~d.receptor_prepares]
print(f"\nreceptor preparation fails for {len(bad)} of {len(d)}")
if len(bad):
    print(bad[["complex_id", "rec_atoms", "error"]].to_string(index=False))
print(f"\nUni-Dock acceptance: {len(d)} -> {int(d.receptor_prepares.sum())} of {len(cov)}")
