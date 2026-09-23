import importlib.util
import sys
import tempfile
import warnings
from pathlib import Path

import gemmi
import pandas as pd
from rdkit import Chem, RDLogger
from posebusters import PoseBusters

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

R = Path("/12TBDrive1/mega_pep_bench")
NATIVE = R / "data/curated/structures/native"

def prepare(path, d):
    st = gemmi.read_structure(str(path)); st.setup_entities(); st.remove_ligands_and_waters()
    m = st[0]
    if len(m) < 2:
        return None, None
    pep = min(m, key=lambda c: len(c))
    ps = gemmi.Structure(); ps.add_model(gemmi.Model("1"))
    pc = pep.clone(); pc.name = "P"; ps[0].add_chain(pc)
    ps.setup_entities(); ps.write_pdb(str(d / "pep.pdb"))
    rs = gemmi.Structure(); rs.add_model(gemmi.Model("1"))
    for i, c in enumerate(m):
        if c.name != pep.name:
            rc = c.clone(); rc.name = chr(ord("A") + i % 26); rs[0].add_chain(rc)
    rs.setup_entities(); rs.write_pdb(str(d / "rec.pdb"))
    mol = Chem.MolFromPDBFile(str(d / "pep.pdb"), proximityBonding=True, removeHs=False,
                              sanitize=False)
    if mol is None:
        return None, None
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    w = Chem.SDWriter(str(d / "pep.sdf")); w.write(mol); w.close()
    return d / "pep.sdf", d / "rec.pdb"

def buster():
    pb = PoseBusters("dock")
    pb.config = {**pb.config,
                 "modules": [m for m in pb.config["modules"] if m["name"] != "Energy ratio"]}
    return pb

def main():
    which = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    pb = buster()
    ids = sorted(pd.read_parquet(R / "runs/track_a/scores_pep/protenix.parquet").complex_id)
    spec = importlib.util.spec_from_file_location("sd", R / "scripts/run/track_a/score_dockq.py")
    sd = importlib.util.module_from_spec(spec); spec.loader.exec_module(sd)

    out = R / "runs/track_a/pose_validity"; out.mkdir(parents=True, exist_ok=True)
    dest = out / f"{which}.csv"
    done = set(pd.read_csv(dest).complex_id) if dest.exists() else set()
    rows = pd.read_csv(dest).to_dict("records") if dest.exists() else []
    todo = [c for c in ids if c not in done]
    if limit:
        todo = todo[:limit]
    print(f"{which}: {len(todo)} to check ({len(done)} done)", flush=True)

    with tempfile.TemporaryDirectory() as tmp:
        for k, cid in enumerate(todo, 1):
            src = (NATIVE / f"{cid}.cif") if which == "native" else sd.find_pred(which, cid)
            if not src or not Path(src).exists():
                continue
            try:
                d = Path(tempfile.mkdtemp(dir=tmp))
                p, r = prepare(Path(src), d)
                if p is None:
                    rows.append({"complex_id": cid, "prepared": False}); continue
                res = pb.bust([p], None, r).iloc[0]
                rows.append({"complex_id": cid, "prepared": True,
                             **{c: bool(res[c]) for c in res.index
                                if c not in ("file", "molecule") and res[c] in (True, False)}})
            except Exception as e:
                rows.append({"complex_id": cid, "prepared": False,
                             "error": f"{type(e).__name__}: {str(e)[:120]}"})
            if k % 25 == 0:
                pd.DataFrame(rows).to_csv(dest, index=False)
                print(f"  {which}: {k}/{len(todo)}", flush=True)
    pd.DataFrame(rows).to_csv(dest, index=False)
    print(f"{which}: wrote {len(rows)} -> {dest}", flush=True)

if __name__ == "__main__":
    main()
