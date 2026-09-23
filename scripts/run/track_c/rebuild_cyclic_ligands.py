import gzip
import subprocess
import sys
import pandas as pd
import gemmi
from pathlib import Path
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

R = Path("/12TBDrive1/mega_pep_bench")
CIF = R / "data/raw/structures/cif"

def fetch(pdb):
    f = CIF / f"{pdb.lower()}.cif.gz"
    if f.exists():
        return f
    CIF.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb.upper()}.cif.gz"
    r = subprocess.run(["curl", "-sfL", "--max-time", "60", "-o", str(f), url])
    return f if r.returncode == 0 and f.exists() and f.stat().st_size > 0 else None

def read(pdb):
    f = fetch(pdb)
    if not f:
        return None, None
    with gzip.open(f, "rt") as fh:
        txt = fh.read()
    st = gemmi.read_structure_string(txt)
    st.setup_entities()
    doc = gemmi.cif.read_string(txt)
    return st, doc

def conns(doc, chain):
    out = []
    for block in doc:
        t = block.find("_struct_conn.", ["conn_type_id", "ptnr1_auth_asym_id", "ptnr1_auth_seq_id",
                                         "ptnr1_label_atom_id", "ptnr2_auth_asym_id",
                                         "ptnr2_auth_seq_id", "ptnr2_label_atom_id"])
        for row in t:
            kind = row[0].strip().strip("'\"").lower()
            if kind not in ("disulf", "covale"):
                continue
            c1, r1, a1, c2, r2, a2 = (row[i].strip().strip("'\"") for i in (1, 2, 3, 4, 5, 6))
            if c1 != chain or c2 != chain:
                continue
            try:
                r1i, r2i = int(r1), int(r2)
            except ValueError:
                continue
            if r1i == r2i:
                continue
            out.append((r1i, a1, r2i, a2, kind))
    return out

def peptide_chain(st, want_len):
    best = None
    for ch in st[0]:
        aa = [r for r in ch if gemmi.find_tabulated_residue(r.name)
              and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        if not aa:
            continue
        if best is None or abs(len(aa) - want_len) < abs(len(best[1]) - want_len):
            best = (ch.name, aa)
    return best

def build(pdb, want_len, chem):
    st, doc = read(pdb)
    if st is None:
        return None, "could not fetch deposited structure"
    hit = peptide_chain(st, want_len)
    if hit is None:
        return None, "no amino-acid chain"
    cname, res = hit
    sel = gemmi.Structure()
    sel.add_model(gemmi.Model("1"))
    ch = gemmi.Chain(cname)
    for r in res:
        ch.add_residue(r)
    sel[0].add_chain(ch)
    block = sel.make_pdb_string()
    mol = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False, proximityBonding=True)
    if mol is None:
        return None, "unparsable peptide chain"

    idx = {}
    for i, a in enumerate(mol.GetAtoms()):
        pi = a.GetPDBResidueInfo()
        if pi:
            idx[(pi.GetResidueNumber(), pi.GetName().strip())] = i

    rw = Chem.RWMol(mol)
    added = []
    for r1, a1, r2, a2, kind in conns(doc, cname):
        i, j = idx.get((r1, a1)), idx.get((r2, a2))
        if i is None or j is None or i == j:
            continue
        if rw.GetBondBetweenAtoms(i, j) is None:
            rw.AddBond(i, j, Chem.BondType.SINGLE)
            added.append(f"{kind}:{r1}{a1}-{r2}{a2}")
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception as e:
        return None, f"sanitize failed: {str(e)[:50]}"

    ri = m.GetRingInfo()
    macro = max([len(x) for x in ri.AtomRings()] + [0])
    ss = sum(1 for b in m.GetBonds()
             if b.GetBeginAtom().GetSymbol() == "S" and b.GetEndAtom().GetSymbol() == "S")
    need_ss = chem in ("disulfide", "bicyclic")
    if need_ss and ss == 0:
        return None, f"{chem} but no S-S recovered (added: {added or 'none'})"
    if not need_ss and macro <= 7:
        return None, f"{chem} but largest ring is {macro} (added: {added or 'none'})"
    return {"smiles": Chem.MolToSmiles(m), "heavy_atoms": m.GetNumAtoms(),
            "largest_ring": macro, "n_ss": ss, "closures": ";".join(added) or "none",
            "chain": cname}, "ok"

def main():
    cyc = pd.read_parquet(R / "data/curated/affinity_mmpbgbsa_cp.parquet")
    rows, refused = [], []
    for r in cyc.itertuples():
        want = len(str(r.peptide_seq)) if getattr(r, "peptide_seq", None) else 10
        d, why = build(r.receptor_pdb_id, want, r.cyclization_type)
        if d is None:
            refused.append((r.receptor_pdb_id, r.cyclization_type, why))
            continue
        rows.append({"complex_id": r.complex_id, "receptor_pdb_id": r.receptor_pdb_id,
                     "cyclization_type": r.cyclization_type, "log_affinity": r.log_affinity, **d})
        print(f"  {r.receptor_pdb_id} [{r.cyclization_type:12s}] ring={d['largest_ring']:3d} "
              f"SS={d['n_ss']} closures={d['closures'][:48]}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(R / "runs/track_c/cyclic_ligands_v2.csv", index=False)
    print(f"\nrebuilt {len(d)}/{len(cyc)} cyclic ligands with real ring closure")
    for p, c, why in refused:
        print(f"  REFUSED {p} [{c}]: {why}")
    if len(d):
        print(f"\n  macrocycle (>7): {(d.largest_ring > 7).sum()}   with S-S: {(d.n_ss > 0).sum()}")
        print(f"  heavy atoms: min {d.heavy_atoms.min()} med {int(d.heavy_atoms.median())} "
              f"max {d.heavy_atoms.max()}")

if __name__ == "__main__":
    sys.exit(main())
