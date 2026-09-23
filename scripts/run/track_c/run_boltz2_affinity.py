import os
import argparse
import subprocess
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/track_c/boltz2_affinity"
OUT.mkdir(parents=True, exist_ok=True)

def yaml_for(row, recseq):
    return f"""version: 1
sequences:
  - protein:
      id: R
      sequence: {recseq}
  - protein:
      id: P
      sequence: {row.peptide_seq}
properties:
  - affinity:
      binder: P
"""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    df = pd.read_parquet(ROOT / "runs/track_c/affinity_manifest.parquet")
    rs = ROOT / "runs/track_c/receptor_seqs.csv"
    if not rs.exists():
        raise SystemExit("run fetch_receptor_seqs.py first (Boltz-2 affinity needs the receptor sequence)")
    recseq = dict(zip(pd.read_csv(rs).receptor_pdb_id, pd.read_csv(rs).receptor_seq))

    rows = df.itertuples(index=False)
    done = 0
    for row in rows:
        if a.limit and done >= a.limit:
            break
        if row.receptor_pdb_id not in recseq:
            continue
        d = OUT / row.complex_id
        pred = d / "predictions" / row.complex_id / "affinity.json"
        if pred.exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / "input.yaml").write_text(yaml_for(row, recseq[row.receptor_pdb_id]))
        subprocess.run(["boltz", "predict", str(d / "input.yaml"), "--out_dir", str(d),
                        "--use_msa_server"], check=False)
        done += 1
    print(f"ran {done} Boltz-2 affinity predictions -> {OUT}")
    print("collect with: python collect_boltz2_affinity.py  (parse affinity.json -> leaderboard)")

if __name__ == "__main__":
    main()
