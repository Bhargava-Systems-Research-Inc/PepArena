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

PROT_IN = ROOT / "prediction_inputs/protenix"

OUT_JSON = ROOT / "prediction_inputs/protenix_msa"
MSA_ROOT = ROOT / "runs/track_a/_protenix_msa"

def boltz_chain_msas(cid):
    out = {}
    for f in sorted(_boltz_msa_glob(f"{cid}_*.csv")):
        if not re.match(rf"^{re.escape(cid)}_\d+\.csv$", f.name):
            continue
        rows = [r["sequence"] for r in csv.DictReader(f.open()) if r.get("sequence")]
        if rows:
            out[rows[0]] = rows
    return out

def write_non_pairing(dest_dir, rows):
    dest_dir.mkdir(parents=True, exist_ok=True)
    with (dest_dir / "non_pairing.a3m").open("w") as f:
        f.write(f">query\n{rows[0]}\n")
        for i, s in enumerate(rows[1:]):
            f.write(f">hit_{i}\n{s}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    OUT_JSON.mkdir(parents=True, exist_ok=True)
    n_ok, n_chains_msa, n_chains_query_only, n_missing = 0, 0, 0, 0
    for cid in man.complex_id:
        src = PROT_IN / f"{cid}.json"
        if not src.exists():
            n_missing += 1; continue
        d = json.loads(src.read_text())
        msas = boltz_chain_msas(cid)
        for ci, ch in enumerate(d[0]["sequences"]):
            seq = ch["proteinChain"]["sequence"]
            dest = MSA_ROOT / cid / str(ci)
            rows = msas.get(seq, [seq])
            write_non_pairing(dest, rows)
            n_chains_msa += len(rows) > 1
            n_chains_query_only += len(rows) == 1
            ch["proteinChain"]["msa"] = {"precomputed_msa_dir": str(dest.resolve()),
                                         "pairing_db": "uniref100"}
        (OUT_JSON / f"{cid}.json").write_text(json.dumps(d, indent=2))
        n_ok += 1
    print(f"wrote {n_ok} protenix-with-MSA JSONs -> {OUT_JSON}")
    print(f"  chains with Boltz MSA: {n_chains_msa} | query-only (peptides/no-hit): "
          f"{n_chains_query_only} | complexes missing input: {n_missing}")

if __name__ == "__main__":
    main()
