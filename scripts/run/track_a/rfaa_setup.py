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

WORK = ROOT / "prediction_inputs/rfaa_work"
PRED = ROOT / "runs/track_a/predictions/rfaa"
WEIGHTS = "/data/RoseTTAFold-All-Atom/weights/RFAA_paper_weights.pt"
PDB100 = "/12TBDrive1/rfaa_db/pdb100_2021Mar03/pdb100_2021Mar03"

def chain_a3m(cid, seq):
    for f in sorted(_boltz_msa_glob(f"{cid}_*.csv")):
        if re.match(rf"^{re.escape(cid)}_\d+\.csv$", f.name):
            rows = [r["sequence"] for r in csv.DictReader(f.open()) if r.get("sequence")]
            if rows and rows[0] == seq:
                return ">query\n" + rows[0] + "\n" + "".join(
                    f">s{i}\n{s}\n" for i, s in enumerate(rows[1:]))
    return f">query\n{seq}\n"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    cids = list(man.complex_id)
    if args.limit:
        cids = cids[:args.limit]
    WORK.mkdir(parents=True, exist_ok=True)
    n = 0
    for cid in cids:
        src = AF3_CLI / f"{cid}.json"
        if not src.exists():
            continue
        chains = json.loads(src.read_text())["sequences"]
        pinputs = {}
        for ch in chains:
            cidlab, seq = ch["protein"]["id"], ch["protein"]["sequence"]
            fa = WORK / cid / f"{cidlab}.fasta"
            fa.parent.mkdir(parents=True, exist_ok=True)
            fa.write_text(f">{cid}_{cidlab}\n{seq}\n")
            msa_dir = PRED / cid / cidlab
            msa_dir.mkdir(parents=True, exist_ok=True)
            (msa_dir / "t000_.msa0.a3m").write_text(chain_a3m(cid, seq))
            pinputs[cidlab] = {"fasta_file": str(fa.resolve())}
        cfg = ["defaults:", "  - base", "", f'job_name: "{cid}"',
               f'output_path: "{PRED.resolve()}"',
               f'checkpoint_path: {WEIGHTS}', "protein_inputs:"]
        for lab, d in pinputs.items():
            cfg += [f"  {lab}:", f"    fasta_file: {d['fasta_file']}"]
        cfg += ["database_params:", f"  hhdb: {PDB100}", "  command: make_msa.sh",
                "  num_cpus: 4", "  mem: 64", '  sequencedb: ""']
        (WORK / f"{cid}.yaml").write_text("\n".join(cfg) + "\n")
        n += 1
    print(f"set up {n} RFAA complexes (reused MSA, empty templates) -> configs in {WORK}, MSA in {PRED}")

if __name__ == "__main__":
    main()
