import pandas as pd
from pathlib import Path
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

R = Path("/12TBDrive1/mega_pep_bench")
PT = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15, "F": 9, "CL": 17, "BR": 35, "I": 53}

def read_prmtop(path):
    blocks, flag = {}, None
    for line in open(path):
        if line.startswith("%FLAG"):
            flag = line.split()[1]; blocks[flag] = []
        elif line.startswith("%"):
            continue
        elif flag:
            blocks[flag].append(line.rstrip("\n"))

    def ints(f):
        return [int(x) for x in " ".join(blocks.get(f, [])).split()]

    names = [s for ln in blocks.get("ATOM_NAME", []) for s in
             [ln[i:i + 4].strip() for i in range(0, len(ln.rstrip()), 4)] if s]
    z = ints("ATOMIC_NUMBER")
    if not z:
        z = [PT.get(n[:2].upper(), PT.get(n[:1].upper(), 6)) for n in names]
    bonds = []
    for f in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        v = ints(f)
        for i in range(0, len(v), 3):
            bonds.append((v[i] // 3, v[i + 1] // 3))
    return names, z, bonds

def build(pdb):
    d = R / "runs/track_c/mmgbsa" / pdb.upper()
    top = d / "lig.prmtop"
    if not top.exists():
        return None, "no ligand topology"
    names, z, bonds = read_prmtop(top)
    if not names or not bonds:
        return None, "empty topology"
    m = Chem.RWMol()
    for num in z:
        m.AddAtom(Chem.Atom(int(num)))
    seen = set()
    for a, b in bonds:
        if a >= len(z) or b >= len(z):
            continue
        k = tuple(sorted((a, b)))
        if k in seen or a == b:
            continue
        seen.add(k)
        m.AddBond(a, b, Chem.BondType.SINGLE)
    mol = m.GetMol()
    mol = Chem.RemoveHs(mol, sanitize=False)
    try:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception as e:
        return None, f"sanitize failed: {str(e)[:50]}"
    return mol, f"{len(seen)} bonds from topology"

def summarise(m):
    ri = m.GetRingInfo()
    macro = max([len(r) for r in ri.AtomRings()] + [0])
    ss = sum(1 for b in m.GetBonds()
             if b.GetBeginAtom().GetSymbol() == "S" and b.GetEndAtom().GetSymbol() == "S")
    return macro, ss

if __name__ == "__main__":
    cyc = pd.read_parquet(R / "data/curated/affinity_mmpbgbsa_cp.parquet")
    rows, refused = [], []
    for r in cyc.itertuples():
        m, note = build(r.receptor_pdb_id)
        if m is None:
            refused.append((r.receptor_pdb_id, note)); continue
        macro, ss = summarise(m)
        expect_ss = r.cyclization_type in ("disulfide", "bicyclic")
        ok = (ss > 0) if expect_ss else (macro > 7)
        if not ok:
            refused.append((r.receptor_pdb_id, f"{r.cyclization_type}: ring {macro}, S-S {ss}"))
            continue
        rows.append({"complex_id": r.complex_id, "receptor_pdb_id": r.receptor_pdb_id,
                     "cyclization_type": r.cyclization_type, "log_affinity": r.log_affinity,
                     "smiles": Chem.MolToSmiles(m), "heavy_atoms": m.GetNumAtoms(),
                     "largest_ring": macro, "n_ss": ss, "note": note})
    d = pd.DataFrame(rows)
    d.to_csv(R / "runs/track_c/cyclic_ligands_rebuilt.csv", index=False)
    print(f"rebuilt {len(d)}/{len(cyc)} cyclic ligands")
    for p, why in refused:
        print(f"  REFUSED {p}: {why}")
    if len(d):
        print(f"\n  with a macrocycle (>7 atoms): {(d.largest_ring > 7).sum()}")
        print(f"  with at least one S-S       : {(d.n_ss > 0).sum()}")
        print(f"  by chemistry: {d.cyclization_type.value_counts().to_dict()}")
        print(f"  heavy atoms: min {d.heavy_atoms.min()} med {int(d.heavy_atoms.median())} max {d.heavy_atoms.max()}")
