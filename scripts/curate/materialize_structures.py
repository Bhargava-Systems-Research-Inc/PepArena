from __future__ import annotations
import hashlib, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import struct_utils

CUR = ROOT / "data/curated"
CIF_DIR = ROOT / "data/raw/structures/cif"
NATIVE = CUR / "structures/native"
OUT = CUR / "structure_native_index.parquet"

def split_chains(val):
    return [c.strip() for c in str(val).replace(":", ",").split(",") if c.strip()]

def materialize(row):
    pep = str(row["peptide_chain"])
    rec = split_chains(row.get("receptor_chains"))
    keep = set(rec + [pep])
    src = row.get("native_complex_path")
    if isinstance(src, str) and Path(src).exists():
        st = struct_utils.read_any(src)
    elif pd.notna(row.get("receptor_pdb_id")):
        cif = CIF_DIR / f"{str(row['receptor_pdb_id']).lower()}.cif.gz"
        if not cif.exists():
            return None
        st = struct_utils.read_any(cif)
    else:
        return None

    model = st[0]
    for cn in [c.name for c in model]:
        if cn not in keep:
            model.remove_chain(cn)
    st.remove_waters()
    st.remove_empty_chains()

    n_pep = sum(1 for c in st[0] if c.name == pep
                for r in c if (gi := struct_utils.gemmi.find_tabulated_residue(r.name)) and gi.is_amino_acid())
    n_rec = sum(1 for c in st[0] if c.name in rec
                for r in c if (gi := struct_utils.gemmi.find_tabulated_residue(r.name)) and gi.is_amino_acid())
    rec_fasta = "".join(struct_utils.chain_seq_and_ncaa(st, c)[0] for c in rec)

    NATIVE.mkdir(parents=True, exist_ok=True)
    out = NATIVE / f"{row['complex_id']}.cif"
    st.setup_entities()
    cif_str = st.make_mmcif_document().as_string()
    out.write_text(cif_str)
    sha = hashlib.sha256(cif_str.encode()).hexdigest()
    return str(out), sha, n_rec, n_pep, rec_fasta

def main():
    df = pd.read_parquet(CUR / "structure_all.parquet")
    rows, fail = [], 0
    for _, r in df.iterrows():
        res = materialize(r)
        if res is None:
            fail += 1
            print(f"  {r['complex_id']}: FAILED"); continue
        path, sha, n_rec, n_pep, fasta = res
        rows.append({"complex_id": r["complex_id"], "native_path": path,
                     "sha256": sha, "n_receptor_res": n_rec, "n_peptide_res": n_pep,
                     "receptor_fasta": fasta})
    idx = pd.DataFrame(rows)
    idx.to_parquet(OUT, index=False)
    keep = set(idx.complex_id)
    stale = [p for p in NATIVE.glob("*.cif") if p.stem not in keep]
    for p in stale:
        p.unlink()
    if stale:
        print(f"removed {len(stale)} stale native files")
    print(f"\nmaterialized {len(idx)} native complexes -> {NATIVE} ({fail} failed)")
    print(f"index -> {OUT}")
    bad = idx[(idx.n_peptide_res == 0) | (idx.n_receptor_res == 0)]
    print(f"sanity: {len(bad)} complexes with an empty receptor or peptide chain (should be 0)")
    if len(bad):
        print(bad[["complex_id", "n_receptor_res", "n_peptide_res"]].to_string())

if __name__ == "__main__":
    main()
