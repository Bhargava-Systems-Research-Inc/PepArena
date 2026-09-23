from __future__ import annotations
import argparse, json, string, time, urllib.request
from pathlib import Path
import gemmi
import pandas as pd

ROOT = Path("/12TBDrive1/mega_pep_bench")
OUT = ROOT / "ncaa_ccd"
CCDCACHE = OUT / "_userccd"
BOLTZ_MOLS = Path.home() / ".boltz/mols"
STD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
SOLV = set("HOH EDO GOL SO4 PO4 CL NA MG ZN CA ACT PEG DMS MPD BME TRS IPA NO3 ACE NH2".split())

def is_aa(name):
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())

def parent1(name):
    if name in STD:
        return gemmi.find_tabulated_residue(name).one_letter_code.upper()
    info = gemmi.find_tabulated_residue(name)
    c = info.one_letter_code.upper() if info else ""
    return c if c and c.isalpha() and c != "X" else "G"

def fetch_userccd(code):
    CCDCACHE.mkdir(parents=True, exist_ok=True)
    f = CCDCACHE / f"{code}.cif"
    if not f.exists():
        try:
            urllib.request.urlretrieve(f"https://files.rcsb.org/ligands/download/{code}.cif", f)
            time.sleep(0.1)
        except Exception:
            return None
    txt = f.read_text()
    return txt if txt.strip() else None

def peptide_and_receptors(cid):
    st = gemmi.read_structure(str(ROOT / f"data/curated/structures/native/{cid}.cif"))
    chains = []
    for c in st[0]:
        res = [r for r in c if r.name not in SOLV]
        naa = sum(1 for r in res if is_aa(r.name) or r.name in STD)
        if res and naa >= max(1, 0.3 * len(res)):
            chains.append((c.name, res))
    if not chains:
        return None, None
    pep_cid, pep = min(chains, key=lambda x: len(x[1]))
    rec = [(cn, res) for cn, res in chains if cn != pep_cid]
    order = {str(r.seqid.num): i + 1 for i, r in enumerate(pep)}
    return (pep, order), rec

def build(cid, bondsjson):
    (pep, order), rec = peptide_and_receptors(cid) or (None, None)
    if pep is None:
        return None
    pep_seq = "".join(parent1(r.name) for r in pep)
    mods = [(i + 1, r.name) for i, r in enumerate(pep) if r.name not in STD]
    rec_seqs = []
    for _, aa in rec:
        rec_seqs.append("".join(gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                                if r.name in STD else "X" for r in aa))
    bonds = []
    for b in (json.loads(bondsjson) if isinstance(bondsjson, str) else []):
        i, j = order.get(str(b.get("res1"))), order.get(str(b.get("res2")))
        if not (i and j and i != j):
            continue
        if {b["atom1"], b["atom2"]} == {"C", "N"} and abs(i - j) == 1:
            continue
        bonds.append((i, b["atom1"], j, b["atom2"]))
    bonds = sorted({tuple(x) for x in bonds})
    exotic = sorted({c for _, c in mods if not (BOLTZ_MOLS / f"{c}.pkl").exists()})
    return dict(pep_seq=pep_seq, rec_seqs=rec_seqs, mods=mods, bonds=bonds, exotic=exotic)

def labels(n):
    return list(string.ascii_uppercase)[:n]

def write_af3(d, cid, path):
    labs = labels(len(d["rec_seqs"]) + 1)
    seqs = [{"protein": {"id": labs[i], "sequence": s}} for i, s in enumerate(d["rec_seqs"])]
    plab = labs[-1]
    seqs.append({"protein": {"id": plab, "sequence": d["pep_seq"],
                             "modifications": [{"ptmType": c, "ptmPosition": p} for p, c in d["mods"]]}})
    out = {"name": cid, "modelSeeds": [1], "dialect": "alphafold3", "version": 2, "sequences": seqs}
    if d["bonds"]:
        out["bondedAtomPairs"] = [[[plab, i, a1], [plab, j, a2]] for i, a1, j, a2 in d["bonds"]]
    userccd = [t for t in (fetch_userccd(c) for c in d["exotic"]) if t]
    if userccd:
        out["userCCD"] = "\n".join(userccd)
    path.write_text(json.dumps(out, indent=2))

def write_protenix(d, cid, path):
    seqs = [{"proteinChain": {"sequence": s, "count": 1}} for s in d["rec_seqs"]]
    seqs.append({"proteinChain": {"sequence": d["pep_seq"], "count": 1,
                                  "modifications": [{"ptmType": f"CCD_{c}", "ptmPosition": p} for p, c in d["mods"]]}})
    pep_entity = len(seqs)
    out = [{"name": cid, "sequences": seqs}]
    if d["bonds"]:
        out[0]["covalent_bonds"] = [
            {"entity1": pep_entity, "position1": i, "atom1": a1, "copy1": 1,
             "entity2": pep_entity, "position2": j, "atom2": a2, "copy2": 1}
            for i, a1, j, a2 in d["bonds"]]
    path.write_text(json.dumps(out, indent=2))

def write_boltz(d, cid, path):
    labs = labels(len(d["rec_seqs"]) + 1)
    lines = ["version: 1", "sequences:"]
    for i, s in enumerate(d["rec_seqs"]):
        lines += ["  - protein:", f"      id: {labs[i]}", f"      sequence: {s}"]
    plab = labs[-1]
    lines += ["  - protein:", f"      id: {plab}", f"      sequence: {d['pep_seq']}"]
    if d["mods"]:
        lines.append("      modifications:")
        for p, c in d["mods"]:
            lines += [f"        - position: {p}", f"          ccd: {c}"]
    if d["bonds"]:
        lines.append("constraints:")
        for i, a1, j, a2 in d["bonds"]:
            lines += ["  - bond:", f"      atom1: [{plab}, {i}, {a1}]", f"      atom2: [{plab}, {j}, {a2}]"]
    path.write_text("\n".join(lines) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    ncaa = man[man.has_ncaa]
    bondmap = pd.read_parquet(ROOT / "data/curated/structure_all.parquet").set_index("complex_id")["cyclization_bonds"].to_dict()
    cids = list(ncaa.complex_id)[:args.limit] if args.limit else list(ncaa.complex_id)
    for m in ("af3", "protenix", "boltz2"):
        (OUT / m).mkdir(parents=True, exist_ok=True)
    rows, fail = [], 0
    for cid in cids:
        try:
            d = build(cid, bondmap.get(cid))
            if not d or not d["pep_seq"]:
                fail += 1; continue
            write_af3(d, cid, OUT / "af3" / f"{cid}.json")
            write_protenix(d, cid, OUT / "protenix" / f"{cid}.json")
            write_boltz(d, cid, OUT / "boltz2" / f"{cid}.yaml")
            row = ncaa[ncaa.complex_id == cid].iloc[0]
            rows.append({"complex_id": cid, "cyclization_type": row.cyclization_type,
                         "peptide_len": len(d["pep_seq"]), "n_ncaa": len(d["mods"]),
                         "n_bonds": len(d["bonds"]), "n_exotic": len(d["exotic"]),
                         "exotic_codes": ";".join(d["exotic"]),
                         "tier": "userCCD" if d["exotic"] else "clean"})
        except Exception as e:
            fail += 1; print(f"  {cid}: FAIL {e}")
    idx = pd.DataFrame(rows)
    idx.to_csv(OUT / "manifest.csv", index=False)
    print(f"wrote {len(idx)} ncAA CCD inputs x 3 models -> {OUT} ({fail} failed)")
    print(f"  clean (no userCCD): {(idx.tier=='clean').sum()} | needs userCCD: {(idx.tier=='userCCD').sum()}")
    print(f"  by cyclization: {idx.cyclization_type.value_counts().to_dict()}")

if __name__ == "__main__":
    main()
