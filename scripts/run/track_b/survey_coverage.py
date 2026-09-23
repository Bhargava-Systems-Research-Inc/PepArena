#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path("/12TBDrive1/mega_pep_bench")
FRESH = Path(os.environ.get("FRESH", str(REPO / "runs/track_b/fresh")))
ADFR = Path.home() / "ADFRsuite/install/bin"
VINA_MAX_TORSION = 32

def atom_count(p: Path):
    return sum(1 for l in p.read_text(errors="ignore").splitlines()
               if l.startswith(("ATOM", "HETATM")))

def prepare_receptor(rec: Path):
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        (t / "rec.pdb").write_bytes(rec.read_bytes())
        r = subprocess.run([str(ADFR / "prepare_receptor"), "-r", "rec.pdb",
                            "-o", "rec.pdbqt", "-A", "checkhydrogens"],
                           cwd=t, capture_output=True, text=True)
        out = t / "rec.pdbqt"
        if out.exists() and out.stat().st_size > 0:
            return ""
        last = (r.stderr.strip().splitlines() or ["no output"])[-1]
        return last[:90]

def prepare_ligand(pep: Path):
    n_in = atom_count(pep)
    with tempfile.TemporaryDirectory() as td:
        shutil.copy(pep, Path(td) / "lig.pdb")
        env = dict(os.environ, PATH=f"{ADFR}:{os.environ['PATH']}")
        r = subprocess.run([str(ADFR / "prepare_ligand"), "-l", "lig.pdb", "-o", "lig.pdbqt"],
                           cwd=td, capture_output=True, text=True, timeout=300, env=env)
        out = Path(td) / "lig.pdbqt"
        if not out.exists():
            return n_in, 0, None, (r.stderr.strip().splitlines() or ["no output"])[-1][:120]
        txt = out.read_text()
        n_out = sum(1 for l in txt.splitlines() if l.startswith(("ATOM", "HETATM")))
        tors = next((int(l.split()[1]) for l in txt.splitlines() if l.startswith("TORSDOF")), None)
        return n_in, n_out, tors, None

def main():
    wl = pd.read_csv(FRESH / "worklist.csv")
    rows = []
    for _, w in wl.iterrows():
        cid = w.complex_id
        pep = FRESH / "inputs" / cid / "peptide.pdb"
        rec = FRESH / "inputs" / cid / "receptor.pdb"
        if not pep.exists():
            continue
        text = pep.read_text(errors="ignore")
        has_het = any(l.startswith("HET ") for l in text.splitlines())
        nrec_chains = len({l[21] for l in rec.read_text(errors="ignore").splitlines()
                           if l.startswith("ATOM")})

        if has_het:
            h3 = "refused: non-canonical residues (topoaa would truncate)"
        elif nrec_chains > 1:
            h3 = f"refused: {nrec_chains}-chain receptor (CNS mol_fix_origin)"
        else:
            h3 = "ok"

        if has_het:
            ad = "refused: non-canonical residues (no rotamers)"
        elif int(w.peptide_len) < 5:
            ad = f"refused: peptide_len={w.peptide_len} below supported range"
        else:
            ad = "ok"

        n_in, n_out, tors, err = prepare_ligand(pep)
        rec_err = prepare_receptor(rec) if not err else ""
        if err:
            ud = f"refused: ligand prep failed ({err})"
        elif rec_err:
            ud = f"refused: receptor prep failed ({rec_err})"
        elif n_out != n_in:
            ud = f"refused: prep silently dropped {n_in - n_out}/{n_in} atoms"
        elif tors is not None and tors > VINA_MAX_TORSION:
            ud = f"refused: {tors} torsions > Vina ceiling {VINA_MAX_TORSION}"
        else:
            ud = "ok"

        rows.append({"complex_id": cid, "peptide_len": w.peptide_len,
                     "cyclization_type": w.cyclization_type, "has_ncaa": bool(w.has_ncaa),
                     "n_receptor_chains": nrec_chains,
                     "lig_atoms_in": n_in, "lig_atoms_out": n_out, "torsdof": tors,
                     "haddock3": h3, "adcp": ad, "unidock": ud})
        print(f"{cid:28s} h3={h3[:28]:28s} adcp={ad[:28]:28s} ud={ud[:44]}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(FRESH / "coverage.csv", index=False)
    print(f"\nwrote {FRESH/'coverage.csv'} ({len(d)} complexes)\n")
    for tool in ("haddock3", "adcp", "unidock"):
        ok = (d[tool] == "ok").sum()
        print(f"{tool:10s} accepts {ok:3d}/{len(d)}  ({100*ok/len(d):.0f}%)")
        for reason, k in d.loc[d[tool] != "ok", tool].str.split("(").str[0].value_counts().items():
            print(f"    {k:3d}  {reason.strip()}")
    print("\nacceptance by ncAA status:")
    for tool in ("haddock3", "adcp", "unidock"):
        g = d.assign(ok=d[tool] == "ok").groupby("has_ncaa").ok.agg(["sum", "size"])
        parts = [f"ncaa={k}: {int(r['sum'])}/{int(r['size'])}" for k, r in g.iterrows()]
        print(f"  {tool:10s} " + "   ".join(parts))

if __name__ == "__main__":
    main()
