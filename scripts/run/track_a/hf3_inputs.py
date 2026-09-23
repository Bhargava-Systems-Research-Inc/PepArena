from __future__ import annotations
import argparse, csv, json, re, string
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

INP = ROOT / "prediction_inputs/helixfold3_hf"
MSAROOT = ROOT / "runs/track_a/_hf3_msa"

def chain_rows(cid, seq):
    for f in sorted(_boltz_msa_glob(f"{cid}_*.csv")):
        if re.match(rf"^{re.escape(cid)}_\d+\.csv$", f.name):
            rows = [r["sequence"] for r in csv.DictReader(f.open()) if r.get("sequence")]
            if rows and rows[0] == seq:
                return rows
    return [seq]

def to_sto(rows):
    lines = ["# STOCKHOLM 1.0"]
    for i, r in enumerate(rows):
        aligned = "".join(c for c in r if not c.islower())
        lines.append(f"seq{i} {aligned}")
    return "\n".join(lines) + "\n//\n"

def write_msa(dest, rows):
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "small_bfd_hits.sto").write_text(to_sto(rows))
    for nm in ("uniref90_hits.sto", "mgnify_hits.sto", "uniprot_hits.sto"):
        (dest / nm).write_text(to_sto(rows[:1]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    INP.mkdir(parents=True, exist_ok=True)
    n = 0
    for cid in man.complex_id:
        src = AF3_CLI / f"{cid}.json"
        if not src.exists():
            continue
        cli = json.loads(src.read_text())
        ents, labels = [], iter(string.ascii_uppercase)
        for s in cli["sequences"]:
            seq = s["protein"]["sequence"]
            ents.append({"type": "protein", "sequence": seq, "count": 1})
            lab = next(labels)
            write_msa(MSAROOT / cid / "msas" / f"protein_{lab}", chain_rows(cid, seq))
        (INP / f"{cid}.json").write_text(json.dumps({"entities": ents}))
        n += 1
    print(f"wrote {n} HelixFold3 entities JSONs -> {INP}")
    print(f"  precomputed MSA (a3m+sto per chain) -> {MSAROOT}/<cid>/msas/protein_<chain>/")

if __name__ == "__main__":
    main()
