import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
REPO = "/12TBDrive1/mega_pep_bench"
m = pd.read_parquet(f"{REPO}/runs/track_c/affinity_manifest.parquet")
rs = pd.read_csv(f"{REPO}/runs/track_c/receptor_seqs.csv")
seqcol = [c for c in rs.columns if "seq" in c.lower()][0]
idcol = [c for c in rs.columns if "pdb" in c.lower() or "id" in c.lower()][0]
seqs = dict(zip(rs[idcol].astype(str).str.lower(), rs[seqcol]))
print("receptor_seqs columns:", rs.columns.tolist(), "->", idcol, seqcol, "n=", len(seqs))

rows = []
drop_nosmiles = drop_toobig = drop_noseq = 0
for r in m.itertuples(index=False):
    mol = Chem.MolFromSequence(str(r.peptide_seq or ""))
    if mol is None:
        drop_nosmiles += 1; continue
    heavy = mol.GetNumAtoms()
    if heavy > 128:
        drop_toobig += 1; continue
    rseq = seqs.get(str(r.receptor_pdb_id).lower())
    if not isinstance(rseq, str) or not rseq:
        drop_noseq += 1; continue
    rows.append({"complex_id": r.complex_id, "receptor_pdb_id": r.receptor_pdb_id,
                 "receptor_seq": rseq, "peptide_seq": r.peptide_seq, "heavy": heavy,
                 "smiles": Chem.MolToSmiles(mol), "log_affinity": r.log_affinity,
                 "source_dataset": r.source_dataset, "is_cyclic": r.is_cyclic})
d = pd.DataFrame(rows)
d.to_csv(f"{REPO}/runs/track_c/affinity_targets.csv", index=False)
print(f"targets={len(d)}  dropped: no_smiles={drop_nosmiles} >128heavy={drop_toobig} no_recseq={drop_noseq}")
print(f"unique receptors needing MSA: {d.receptor_pdb_id.nunique()}")
print(f"receptor length: min={d.receptor_seq.str.len().min()} med={int(d.receptor_seq.str.len().median())} max={d.receptor_seq.str.len().max()}")
