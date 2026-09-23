from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
from make_splits import mmseqs_clusters, _bucket

MAN = ROOT / "data/curated/interaction_tpeppro.parquet"
OUT = ROOT / "data/splits/interaction_splits.parquet"

def main():
    df = pd.read_parquet(MAN)
    rec_seq = dict(zip(df.complex_id, df.receptor_seq))
    pep_seq = dict(zip(df.complex_id, df.peptide_seq))
    uniq_rec = {f"r{i}": s for i, s in enumerate(sorted(set(df.receptor_seq)))}
    uniq_pep = {f"p{i}": s for i, s in enumerate(sorted(set(df.peptide_seq)))}
    print(f"clustering {len(uniq_rec)} receptor + {len(uniq_pep)} peptide sequences...")
    rclu = mmseqs_clusters(uniq_rec)
    pclu = mmseqs_clusters(uniq_pep)
    rec_to_clu = {uniq_rec[k]: rclu[k] for k in uniq_rec}
    pep_to_clu = {uniq_pep[k]: pclu[k] for k in uniq_pep}

    out = pd.DataFrame({"complex_id": df.complex_id, "label": df.label})
    out["receptor_cluster"] = df.receptor_seq.map(rec_to_clu)
    out["peptide_cluster"] = df.peptide_seq.map(pep_to_clu)
    out["split_random"] = out.complex_id.map(lambda c: "test" if _bucket(c) == 0 else "train")
    out["split_receptor"] = out.receptor_cluster.map(lambda c: "test" if _bucket(c) == 0 else "train")
    out["split_peptide"] = out.peptide_cluster.map(lambda c: "test" if _bucket(c) == 0 else "train")
    out["split_both"] = [
        "test" if (r == "test" and p == "test")
        else ("train" if (r == "train" and p == "train") else "none")
        for r, p in zip(out.split_receptor, out.split_peptide)
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)
    print(f"wrote {OUT}")
    print("receptor clusters:", out.receptor_cluster.nunique(), "| peptide clusters:", out.peptide_cluster.nunique())
    for s in ("split_random", "split_receptor", "split_peptide", "split_both"):
        vc = out[s].value_counts()
        bal = out.groupby(s).label.mean()
        print(f"  {s}: train {int(vc.get('train',0))} / test {int(vc.get('test',0))} "
              f"| test pos-rate {bal.get('test',float('nan')):.2f}")
    for col, key in (("split_receptor", "receptor_cluster"), ("split_peptide", "peptide_cluster")):
        span = out.groupby(key)[col].nunique()
        print(f"  {key} spanning train+test in {col} (must be 0):", int((span > 1).sum()))
    tb = set(out.loc[out.split_both == "test", "receptor_cluster"])
    trb = set(out.loc[out.split_both == "train", "receptor_cluster"])
    pb = set(out.loc[out.split_both == "test", "peptide_cluster"])
    ptrb = set(out.loc[out.split_both == "train", "peptide_cluster"])
    print("  split_both receptor-cluster overlap (must be 0):", len(tb & trb))
    print("  split_both peptide-cluster overlap  (must be 0):", len(pb & ptrb))
    print("  split_both dropped (mixed cells):", int((out.split_both == "none").sum()))

if __name__ == "__main__":
    main()
