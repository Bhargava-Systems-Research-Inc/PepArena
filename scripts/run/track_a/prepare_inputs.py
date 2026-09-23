from __future__ import annotations
import argparse, json, string, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import struct_utils, cyclization

FOLDLIST = ROOT / "runs/track_a/foldlist.parquet"
OUTROOT = ROOT / "prediction_inputs"
MODELS = ["af3", "protenix", "helixfold3", "highfold3", "boltz2", "chai1", "esmfold2", "rfaa"]

def extract_chains(native_path, receptor_chains, peptide_chain, manifest_bonds_json):
    st = struct_utils.read_any(native_path)
    rec_ids = [c.strip() for c in str(receptor_chains).replace(":", ",").split(",") if c.strip()]
    labels = iter(string.ascii_uppercase)
    receptors = []
    for rc in rec_ids:
        seq, _, resn = struct_utils.chain_seq_and_ncaa(st, rc)
        if resn:
            receptors.append((next(labels), seq))
    pep_label = next(labels)
    pseq, _, presn = struct_utils.chain_seq_and_ncaa(st, peptide_chain)
    order = {str(r.seqid): i + 1 for i, r in enumerate(struct_utils_aa(st[0].find_chain(peptide_chain)))}
    bonds = []
    for b in (json.loads(manifest_bonds_json) if isinstance(manifest_bonds_json, str) else []):
        i, j = order.get(b.get("res1")), order.get(b.get("res2"))
        if not (i and j and i != j):
            continue
        if {b["atom1"], b["atom2"]} == {"C", "N"} and abs(i - j) == 1:
            continue
        bonds.append((pep_label, i, b["atom1"], pep_label, j, b["atom2"]))
    bonds = sorted({tuple(x) for x in bonds})
    return receptors, (pep_label, pseq), bonds

def struct_utils_aa(chain):
    return [r for r in chain if struct_utils.aa_info(r)[1]] if chain else []

def write_af3(path, cid, receptors, peptide, bonds):
    seqs = [{"protein": {"id": lab, "sequence": s}} for lab, s in receptors + [peptide]]
    d = {"name": cid, "modelSeeds": [1], "sequences": seqs, "dialect": "alphafold3", "version": 2}
    bp = [[[c1, r1, a1], [c2, r2, a2]] for (c1, r1, a1, c2, r2, a2) in bonds]
    if bp:
        d["bondedAtomPairs"] = bp
    path.write_text(json.dumps(d, indent=2))

def write_protenix(path, cid, receptors, peptide, bonds):
    seqs = [{"proteinChain": {"sequence": s, "count": 1}} for _, s in receptors + [peptide]]
    d = [{"name": cid, "sequences": seqs}]
    if bonds:
        d[0]["covalent_bonds"] = [{"left": [c1, r1, a1], "right": [c2, r2, a2]}
                                  for (c1, r1, a1, c2, r2, a2) in bonds]
    path.write_text(json.dumps(d, indent=2))

def write_boltz(path, cid, receptors, peptide, bonds):
    lines = ["version: 1", "sequences:"]
    for lab, s in receptors + [peptide]:
        lines += [f"  - protein:", f"      id: {lab}", f"      sequence: {s}"]
    if bonds:
        lines.append("constraints:")
        for (c1, r1, a1, c2, r2, a2) in bonds:
            lines += ["  - bond:", f"      atom1: [{c1}, {r1}, {a1}]",
                      f"      atom2: [{c2}, {r2}, {a2}]"]
    path.write_text("\n".join(lines) + "\n")

def write_fasta(path, receptors, peptide):
    blocks = [f">protein|{lab}\n{s}" for lab, s in receptors + [peptide]]
    path.write_text("\n".join(blocks) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap-receptor", type=int, default=1800)
    args = ap.parse_args()
    fl = pd.read_parquet(FOLDLIST)
    fl = fl[fl.foldable & (fl.receptor_len.fillna(9999) <= args.cap_receptor)]
    bondmap = pd.read_parquet(ROOT / "data/curated/structure_all.parquet").set_index("complex_id")["cyclization_bonds"].to_dict()
    for m in MODELS:
        (OUTROOT / m).mkdir(parents=True, exist_ok=True)
    rows, fail = [], 0
    for _, r in fl.iterrows():
        cid = r.complex_id
        try:
            receptors, peptide, bonds = extract_chains(r.native_path, r.receptor_chains,
                                                       r.peptide_chain, bondmap.get(cid))
            if not receptors or not peptide[1]:
                fail += 1; continue
        except Exception as e:
            fail += 1; print(f"  {cid}: FAIL {e}"); continue
        write_af3(OUTROOT / "af3" / f"{cid}.json", cid, receptors, peptide, bonds)
        write_af3(OUTROOT / "highfold3" / f"{cid}.json", cid, receptors, peptide, bonds)
        write_af3(OUTROOT / "helixfold3" / f"{cid}.json", cid, receptors, peptide, bonds)
        write_protenix(OUTROOT / "protenix" / f"{cid}.json", cid, receptors, peptide, bonds)
        write_boltz(OUTROOT / "boltz2" / f"{cid}.yaml", cid, receptors, peptide, bonds)
        for m in ("chai1", "esmfold2", "rfaa"):
            write_fasta(OUTROOT / m / f"{cid}.fasta", receptors, peptide)
        rows.append({"complex_id": cid, "n_receptor_chains": len(receptors),
                     "peptide_len": len(peptide[1]), "n_bonds": len(bonds),
                     "is_cyclic": bool(r.is_cyclic), "cyclization_type": r.cyclization_type,
                     "has_ncaa": bool(r.has_ncaa)})
    idx = pd.DataFrame(rows)
    idx.to_csv(OUTROOT / "inputs_manifest.csv", index=False)
    print(f"wrote inputs for {len(idx)} complexes x {len(MODELS)} models ({fail} failed)")
    print(f"  with cyclization bonds: {int((idx.n_bonds>0).sum())} | ncAA (X positions): {int(idx.has_ncaa.sum())}")
    print(f"  -> {OUTROOT}")

if __name__ == "__main__":
    main()
