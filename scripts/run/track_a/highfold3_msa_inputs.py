from __future__ import annotations
import argparse, csv, json, re
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

AF3_CLI = ROOT / "prediction_inputs/af3"

OUT = ROOT / "prediction_inputs/highfold3_msa"

def chain_a3m(cid, seq):
    for f in sorted(_boltz_msa_glob(f"{cid}_*.csv")):
        if not re.match(rf"^{re.escape(cid)}_\d+\.csv$", f.name):
            continue
        rows = [r["sequence"] for r in csv.DictReader(f.open()) if r.get("sequence")]
        if rows and rows[0] == seq:
            lines = [">query", rows[0]]
            for i, s in enumerate(rows[1:]):
                lines += [f">seq_{i+1}", s]
            return "\n".join(lines) + "\n"
    return f">query\n{seq}\n"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    OUT.mkdir(parents=True, exist_ok=True)
    n, n_msa, n_single = 0, 0, 0
    for cid in man.complex_id:
        src = AF3_CLI / f"{cid}.json"
        if not src.exists():
            continue
        cli = json.loads(src.read_text())
        seqs = []
        for s in cli["sequences"]:
            seq = s["protein"]["sequence"]
            a3m = chain_a3m(cid, seq)
            depth = a3m.count(">")
            n_msa += depth > 1; n_single += depth == 1
            seqs.append({"protein": {"id": s["protein"]["id"], "sequence": seq,
                                     "modifications": [], "unpairedMsa": a3m,
                                     "pairedMsa": "", "templates": []}})
        d = {"dialect": "alphafold3", "version": 4, "name": cid, "modelSeeds": [1],
             "sequences": seqs}
        if cli.get("bondedAtomPairs"):
            d["bondedAtomPairs"] = cli["bondedAtomPairs"]
        (OUT / f"{cid}.json").write_text(json.dumps(d))
        n += 1
    print(f"wrote {n} HighFold3 data.json (MSA embedded) -> {OUT}")
    print(f"  chains with MSA: {n_msa} | single-seq: {n_single}")

if __name__ == "__main__":
    main()
