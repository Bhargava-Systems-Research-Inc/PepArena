import gzip
import re
import sys
import pandas as pd
import gemmi
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
CIF = R / "data/raw/structures/cif"
CONTACT = 5.0

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

def load(pdb):
    f = CIF / f"{pdb.lower()}.cif.gz"
    if not f.exists():
        return None
    with gzip.open(f, "rt") as fh:
        st = gemmi.read_structure_string(fh.read())
    st.setup_entities()
    return st

def chain_map(st):
    out = {}
    for ch in st[0]:
        res = [r for r in ch if r.name in AA3 or gemmi.find_tabulated_residue(r.name)
               and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        if not res:
            continue
        out[ch.name] = {
            "num2aa": {r.seqid.num: AA3.get(r.name, "X") for r in res},
            "order": [r.seqid.num for r in res],
            "seq": "".join(AA3.get(r.name, "X") for r in res),
        }
    return out

MIN_CONTACT_RES = 3

def contacts(st, pep_chain, cut=CONTACT):
    ns = gemmi.NeighborSearch(st, cut).populate()
    hit = {}
    pep = None
    for ch in st[0]:
        if ch.name == pep_chain:
            pep = ch
    if pep is None:
        return None
    for res in pep:
        for atom in res:
            if atom.element == gemmi.Element("H"):
                continue
            for m in ns.find_atoms(atom.pos, "\0", radius=cut):
                cra = m.to_cra(st[0])
                if cra.chain.name != pep_chain and cra.residue.name not in ("HOH", "DOD"):
                    if gemmi.find_tabulated_residue(cra.residue.name) and \
                       gemmi.find_tabulated_residue(cra.residue.name).is_amino_acid():
                        hit.setdefault(cra.chain.name, set()).add(cra.residue.seqid.num)
    return {c for c, res in hit.items() if len(res) >= MIN_CONTACT_RES}

def parse_mut(s):
    toks = [t.strip() for t in re.split(r"[,;]", str(s)) if t.strip()]
    out = []
    for t in toks:
        g = re.match(r"^([A-Za-z0-9])_([A-Z])(-?\d+)([A-Z])$", t)
        if not g:
            return None
        out.append((g.group(1), g.group(2), int(g.group(3)), g.group(4)))
    return out

def main():
    ppb = pd.read_parquet(R / "data/curated/affinity_ppb.parquet")
    ppb["mut_raw"] = ppb.notes.str.extract(r"mut=(.*)$")[0].fillna("")
    ppb["mut_raw"] = ppb.mut_raw.replace({"nan": "", "None": ""})

    rows, ledger = [], []
    cache = {}
    for r in ppb.itertuples():
        pdb = r.receptor_pdb_id.lower()
        rec_chains = [c.strip() for c in str(r.receptor_chains).split(",") if c.strip()]
        pep_ch = str(r.peptide_chain)
        construct = set(rec_chains) | {pep_ch}

        key = (pdb, pep_ch)
        if key not in cache:
            st = load(pdb)
            cache[key] = ((chain_map(st), contacts(st, pep_ch)) if st else (None, None))
            del st
        cm, touch = cache[key]
        if cm is None:
            ledger.append((r.complex_id, "no deposited structure locally")); continue
        if pep_ch not in cm:
            ledger.append((r.complex_id, f"peptide chain {pep_ch} not in structure")); continue

        if touch is None:
            ledger.append((r.complex_id, "could not compute contacts")); continue
        extra = touch - construct
        if extra:
            ledger.append((r.complex_id,
                           f"incomplete construct: chains {sorted(extra)} contact the peptide "
                           f"but are not in the curated complex")); continue

        muts = parse_mut(r.mut_raw) if r.mut_raw else []
        if muts is None:
            ledger.append((r.complex_id, f"unparseable mutation '{r.mut_raw}'")); continue

        off = [m for m in muts if m[0] != pep_ch]
        if off:
            ledger.append((r.complex_id,
                           f"mutation on non-peptide chain(s) {sorted({m[0] for m in off})}")); continue

        seq = list(cm[pep_ch]["seq"])
        num2idx = {n: i for i, n in enumerate(cm[pep_ch]["order"])}
        bad = None
        for ch, wt, pos, mt in muts:
            i = num2idx.get(pos)
            if i is None:
                bad = f"position {pos} absent from chain {ch}"; break
            if seq[i] != wt:
                bad = f"wild-type mismatch at {pos}: structure has {seq[i]}, mutation says {wt}"; break
            seq[i] = mt
        if bad:
            ledger.append((r.complex_id, bad)); continue

        rec_seq = "".join(cm[c]["seq"] for c in rec_chains if c in cm)
        if not rec_seq:
            ledger.append((r.complex_id, "no receptor sequence")); continue

        rows.append({"complex_id": r.complex_id, "receptor_pdb_id": pdb,
                     "receptor_chains": ",".join(rec_chains), "peptide_chain": pep_ch,
                     "receptor_seq": rec_seq, "peptide_seq": "".join(seq),
                     "peptide_seq_wt": cm[pep_ch]["seq"], "n_mutations": len(muts),
                     "mutations": r.mut_raw, "log_affinity": r.log_affinity,
                     "measured_value": r.measured_value, "assay": r.assay,
                     "source": str(r.notes).split(";")[0].replace("src=", "").strip()})

    d = pd.DataFrame(rows)
    led = pd.DataFrame(ledger, columns=["complex_id", "reason"])
    out = R / "runs/track_c/affinity_manifest_v2.parquet"
    d.to_parquet(out, index=False)
    led.to_csv(R / "runs/track_c/affinity_excluded_v2.csv", index=False)

    print(f"kept {len(d)} of {len(ppb)} PPB affinity rows")
    print("\nexclusion reasons:")
    print(led.reason.str.replace(r"\d+", "N", regex=True).str.slice(0, 70)
             .value_counts().head(12).to_string())
    if len(d):
        combos = d.groupby(["receptor_seq", "peptide_seq"]).ngroups
        print(f"\ndistinct molecular inputs : {combos}  (was 136 for 1,436 rows)")
        print(f"rows                      : {len(d)}")
        print(f"rows carrying a mutation  : {int((d.n_mutations > 0).sum())}")
        dup = d.groupby(["receptor_seq", "peptide_seq"]).log_affinity.nunique()
        print(f"inputs with >1 affinity   : {int((dup > 1).sum())}")
        print(f"unique PDB entries        : {d.receptor_pdb_id.nunique()}")
        print(f"peptide length: min {d.peptide_seq.str.len().min()} "
              f"med {int(d.peptide_seq.str.len().median())} max {d.peptide_seq.str.len().max()}")
    print(f"\nwrote {out}")

if __name__ == "__main__":
    sys.exit(main())
