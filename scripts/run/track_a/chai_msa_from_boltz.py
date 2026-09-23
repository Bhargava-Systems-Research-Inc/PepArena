from __future__ import annotations
import argparse, csv, hashlib, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]

def _boltz_msa_dirs():
    root = ROOT / "runs/track_a/predictions/boltz2"
    ds = [p / "msa" for p in sorted(root.glob("boltz_results__*")) if (p / "msa").is_dir()]
    return sorted(ds, key=lambda p: p.stat().st_mtime, reverse=True)

def _boltz_msa_glob(pattern):
    out = []
    for d in _boltz_msa_dirs():
        out.extend(sorted(d.glob(pattern)))
    return out

OUT = ROOT / "runs/track_a/_chai_msa"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    OUT.mkdir(parents=True, exist_ok=True)
    seen, n_files, n_depth1 = set(), 0, 0
    for cid in man.complex_id:
        for f in sorted(_boltz_msa_glob(f"{cid}_*.csv")):
            if not re.match(rf"^{re.escape(cid)}_\d+\.csv$", f.name):
                continue
            rows = [r["sequence"] for r in csv.DictReader(f.open()) if r.get("sequence")]
            if not rows:
                continue
            q = rows[0]
            h = hashlib.sha256(q.upper().encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h)
            recs = [{"sequence": q, "source_database": "query", "pairing_key": "", "comment": "query"}]
            for s in rows[1:]:
                recs.append({"sequence": s, "source_database": "uniref90",
                             "pairing_key": "", "comment": ""})
            pd.DataFrame.from_records(recs).to_parquet(OUT / f"{h}.aligned.pqt", index=False)
            n_files += 1
            n_depth1 += len(rows) == 1
    print(f"wrote {n_files} unique-chain .aligned.pqt -> {OUT}")
    print(f"  ({n_depth1} query-only / single-seq chains, rest have full Boltz MSA depth)")

if __name__ == "__main__":
    main()
