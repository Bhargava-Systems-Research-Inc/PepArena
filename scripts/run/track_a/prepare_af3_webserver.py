from __future__ import annotations
import argparse, json, math
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
AF3_CLI = ROOT / "prediction_inputs/af3"
OUT = ROOT / "prediction_inputs/af3_webserver"

def server_job(cid, seed):
    cli = json.loads((AF3_CLI / f"{cid}.json").read_text())
    chains = []
    for s in cli["sequences"]:
        seq = s["protein"]["sequence"]
        chains.append({"proteinChain": {"sequence": seq, "count": 1,
                                        "useStructureTemplate": True}})
    return {"name": cid, "modelSeeds": [str(seed)], "sequences": chains,
            "dialect": "alphafoldserver", "version": 2}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=30)
    ap.add_argument("--seed", default="1")
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    lin = man[(~man.has_ncaa) & (~man.is_cyclic)]
    cids = [c for c in lin.complex_id if (AF3_CLI / f"{c}.json").exists()]
    OUT.mkdir(parents=True, exist_ok=True)
    nb = math.ceil(len(cids) / args.batch_size)
    manifest_rows = []
    for b in range(nb):
        batch = cids[b * args.batch_size:(b + 1) * args.batch_size]
        bdir = OUT / f"batch_{b + 1:04d}"
        bdir.mkdir(parents=True, exist_ok=True)
        merged = []
        for cid in batch:
            job = server_job(cid, args.seed)
            (bdir / f"{cid}.json").write_text(json.dumps([job], indent=2) + "\n")
            merged.append(job)
            manifest_rows.append({"batch": b + 1, "complex_id": cid,
                                  "n_chains": len(job["sequences"]),
                                  "pep_len": int(lin.loc[lin.complex_id == cid, "peptide_len"].iloc[0])})
        (bdir / "merged_batch.json").write_text(json.dumps(merged, indent=2) + "\n")
    pd.DataFrame(manifest_rows).to_csv(OUT / "manifest.csv", index=False)
    print(f"wrote {len(cids)} linear canonical jobs in {nb} batches of {args.batch_size} -> {OUT}")
    print(f"  each batch_NNNN/ has single-job <cid>.json files + merged_batch.json (upload either)")

if __name__ == "__main__":
    main()
