import subprocess
import sys
import pandas as pd
from pathlib import Path
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "runs/track_c"
MMSEQS = "mmseqs"

def cluster(seqs, tmp, ident=0.3, cov=0.8):
    tmp.mkdir(parents=True, exist_ok=True)
    fa = tmp / "rec.fasta"
    with open(fa, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">s{i}\n{s}\n")
    pref = tmp / "clu"
    r = subprocess.run([MMSEQS, "easy-cluster", str(fa), str(pref), str(tmp / "tmp"),
                        "--min-seq-id", str(ident), "-c", str(cov)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"ABORT: mmseqs failed, no cluster split can be built:\n{r.stderr[-800:]}")
    rep = {}
    for line in open(f"{pref}_cluster.tsv"):
        a, b = line.split()
        rep[b] = a
    return {seqs[int(k[1:])]: rep[k] for k in rep}

def main():
    man = pd.read_parquet(OUT / "affinity_manifest_v2.parquet")
    cyc = pd.read_csv(OUT / "cyclic_ligands_v2.csv")

    rows, refused = [], []
    for r in man.itertuples():
        pep = str(r.peptide_seq)
        if "X" in pep:
            refused.append((r.complex_id, "peptide carries a residue one-letter chemistry cannot express"))
            continue
        m = Chem.MolFromSequence(pep)
        if m is None:
            refused.append((r.complex_id, "RDKit could not build the peptide")); continue
        rows.append({"complex_id": r.complex_id, "receptor_pdb_id": r.receptor_pdb_id,
                     "receptor_seq": r.receptor_seq, "peptide_seq": pep,
                     "smiles": Chem.MolToSmiles(m), "heavy_atoms": m.GetNumAtoms(),
                     "log_affinity": r.log_affinity, "set": "binary", "is_cyclic": False,
                     "n_mutations": r.n_mutations})
    for r in cyc.itertuples():
        rows.append({"complex_id": r.complex_id, "receptor_pdb_id": r.receptor_pdb_id,
                     "receptor_seq": None, "peptide_seq": None, "smiles": r.smiles,
                     "heavy_atoms": r.heavy_atoms, "log_affinity": r.log_affinity,
                     "set": "cyclic", "is_cyclic": True, "n_mutations": 0})
    d = pd.DataFrame(rows)

    rs = pd.read_csv(OUT / "receptor_seqs.csv")
    seqmap = dict(zip(rs.receptor_pdb_id.astype(str).str.lower(), rs.receptor_seq))
    d.loc[d.receptor_seq.isna(), "receptor_seq"] = \
        d.loc[d.receptor_seq.isna(), "receptor_pdb_id"].str.lower().map(seqmap)
    missing = d.receptor_seq.isna().sum()
    if missing:
        print(f"  dropping {missing} rows with no receptor sequence")
        d = d[d.receptor_seq.notna()].copy()

    uniq = sorted(set(d.receptor_seq))
    cl = cluster(uniq, OUT / "_mmseqs_v2")
    d["receptor_cluster"] = d.receptor_seq.map(cl)
    n_cl = d.receptor_cluster.nunique()
    order = sorted(d.receptor_cluster.unique())
    test = set(order[::5])
    d["split_receptor_cluster"] = ["test" if c in test else "train" for c in d.receptor_cluster]
    d = d.sort_values("complex_id").reset_index(drop=True)
    d["split_random"] = ["test" if i % 5 == 0 else "train" for i in range(len(d))]

    d.to_csv(OUT / "affinity_targets_v2.csv", index=False)
    pd.DataFrame(refused, columns=["complex_id", "reason"]).to_csv(
        OUT / "affinity_targets_v2_refused.csv", index=False)

    print(f"exported {len(d)} affinity targets ({(d.set == 'binary').sum()} binary, "
          f"{(d.set == 'cyclic').sum()} cyclic)")
    print(f"  refused at ligand build: {len(refused)}")
    print(f"  receptor clusters: {n_cl}  (30% identity, 80% coverage)")
    print(f"  cluster split: {(d.split_receptor_cluster == 'train').sum()} train / "
          f"{(d.split_receptor_cluster == 'test').sum()} test")
    tr = set(d[d.split_receptor_cluster == 'train'].receptor_cluster)
    te = set(d[d.split_receptor_cluster == 'test'].receptor_cluster)
    print(f"  clusters spanning the boundary: {len(tr & te)} (must be 0)")
    print(f"  heavy atoms: median {int(d.heavy_atoms.median())}, "
          f"over Boltz-2's 56-atom training domain: {(d.heavy_atoms > 56).sum()}/{len(d)}")
    print(f"\nwrote {OUT / 'affinity_targets_v2.csv'}")

if __name__ == "__main__":
    sys.exit(main())
