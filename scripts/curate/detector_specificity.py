import pathlib
import gemmi
import pandas as pd

R = pathlib.Path("/12TBDrive1/mega_pep_bench")
s = pd.read_parquet(R / "data/curated/structure_all.parquet")
s = s[s.native_complex_path.notna()]

def probe(path, pep_chain):
    try:
        st = gemmi.read_structure(str(path))
    except Exception:
        return None, None
    any_bond, pep_bond = False, False
    for c in st.connections:
        if c.type not in (gemmi.ConnectionType.Disulf, gemmi.ConnectionType.Covale):
            continue
        a, b = c.partner1, c.partner2
        if a.chain_name == b.chain_name:
            try:
                if abs(a.residue_id.seqid.num - b.residue_id.seqid.num) == 1:
                    continue
            except Exception:
                pass
        any_bond = True
        if str(pep_chain) in {a.chain_name, b.chain_name}:
            pep_bond = True
    return any_bond, pep_bond

cyc = s[s.is_cyclic.astype(bool)]
lin = s[~s.is_cyclic.astype(bool)].sample(n=min(800, (~s.is_cyclic.astype(bool)).sum()), random_state=0)

rows = []
for frame in (cyc, lin):
    for _, r in frame.iterrows():
        any_b, pep_b = probe(r.native_complex_path, r.peptide_chain)
        if any_b is None:
            continue
        rows.append({"complex_id": r.complex_id, "detector_cyclic": bool(r.is_cyclic),
                     "file_has_records": any_b, "deposited_pep_bond": pep_b})
d = pd.DataFrame(rows)
d.to_csv(R / "runs/track_a/cyclization_detector_confusion.csv", index=False)

lab = d[d.file_has_records]
print(f"probed {len(d)} complexes; {len(lab)} in files that preserve connectivity records")
tp = int((lab.deposited_pep_bond & lab.detector_cyclic).sum())
fp = int((~lab.deposited_pep_bond & lab.detector_cyclic).sum())
fn = int((lab.deposited_pep_bond & ~lab.detector_cyclic).sum())
tn = int((~lab.deposited_pep_bond & ~lab.detector_cyclic).sum())
print(f"  TP {tp}  FP {fp}  FN {fn}  TN {tn}")
if tp + fn: print(f"  sensitivity {tp/(tp+fn):.3f}")
if tn + fp: print(f"  specificity {tn/(tn+fp):.3f}")
if tp + fp: print(f"  precision   {tp/(tp+fp):.3f}")
stripped = d[~d.file_has_records]
print(f"\nfiles with connectivity stripped: {len(stripped)} "
      f"({int(stripped.detector_cyclic.sum())} of them detector-cyclic) -- excluded from the counts above")
