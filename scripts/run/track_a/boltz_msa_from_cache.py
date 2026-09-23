from __future__ import annotations
import argparse, csv, re
from pathlib import Path
import yaml
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

BOLTZ_IN = ROOT / "prediction_inputs/boltz2"

OUT = ROOT / "prediction_inputs/boltz2_msa"

def chain_csvs(cid):
    out = {}
    for f in sorted(_boltz_msa_glob(f"{cid}_*.csv")):
        if not re.match(rf"^{re.escape(cid)}_\d+\.csv$", f.name):
            continue
        with f.open() as fh:
            r = csv.DictReader(fh)
            first = next(r, None)
        if first and first.get("sequence"):
            out[first["sequence"]] = str(f.resolve())
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    OUT.mkdir(parents=True, exist_ok=True)
    n, n_msa, n_empty = 0, 0, 0
    for cid in man.complex_id:
        src = BOLTZ_IN / f"{cid}.yaml"
        if not src.exists():
            continue
        d = yaml.safe_load(src.read_text())
        csvs = chain_csvs(cid)
        for s in d["sequences"]:
            if "protein" not in s:
                continue
            seq = s["protein"]["sequence"]
            if seq in csvs:
                s["protein"]["msa"] = csvs[seq]; n_msa += 1
            else:
                s["protein"]["msa"] = "empty"; n_empty += 1
        (OUT / f"{cid}.yaml").write_text(yaml.safe_dump(d, sort_keys=False))
        n += 1
    print(f"wrote {n} boltz-with-MSA YAMLs -> {OUT}")
    print(f"  chains with cached MSA: {n_msa} | empty (no cache): {n_empty}")

if __name__ == "__main__":
    main()
