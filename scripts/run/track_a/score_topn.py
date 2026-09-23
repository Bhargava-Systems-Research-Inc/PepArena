import importlib.util
import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
NATIVE = R / "data/curated/structures/native"

spec = importlib.util.spec_from_file_location("spd", R / "scripts/run/track_a/score_peptide_dockq.py")
spd = importlib.util.module_from_spec(spec); spec.loader.exec_module(spd)

PRED = R / "runs/track_a/predictions"

def samples(model, cid):
    b = PRED / model
    if model == "protenix":
        return sorted(b.glob(f"{cid}/seed_*/predictions/{cid}_sample_*.cif"))
    if model == "boltz2":
        return sorted(b.glob(f"boltz_results_*/predictions/{cid}/{cid}_model_*.cif"))
    if model == "chai1":
        return sorted(b.glob(f"{cid}/pred.model_idx_*.cif"))
    if model == "highfold3":
        return sorted(b.glob(f"{cid}/seed-*_sample-*/{cid}_seed-*_sample-*_model.cif"))
    if model == "helixfold3":
        return sorted(b.glob(f"{cid}/{cid}-rank*/predicted_structure.cif"))
    return []

def main():
    model = sys.argv[1]
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    top1 = pd.read_parquet(R / f"runs/track_a/scores_pep/{model}.parquet").set_index("complex_id")
    out = R / "runs/track_a/scores_topn"; out.mkdir(parents=True, exist_ok=True)
    dest = out / f"{model}.parquet"
    done = (set(pd.read_parquet(dest).complex_id) if dest.exists() else set())
    rows = [] if not dest.exists() else pd.read_parquet(dest).to_dict("records")

    ids = [c for c in top1.index if c not in done]
    print(f"{model}: {len(ids)} complexes to score ({len(done)} already done)", flush=True)
    for k, cid in enumerate(ids, 1):
        nat = NATIVE / f"{cid}.cif"
        if not nat.exists():
            continue
        plen = int(sp.peptide_len.get(cid, 10))
        vals = []
        for f in samples(model, cid):
            s = spd.score(f, nat, plen) or spd.score_fallback(f, nat, plen)
            if s:
                vals.append(s["DockQ"])
        if vals:
            rows.append({"complex_id": cid, "n_samples": len(vals),
                         "top1": float(top1.DockQ[cid]), "oracle": max(vals),
                         "worst": min(vals)})
        if k % 20 == 0:
            pd.DataFrame(rows).to_parquet(dest, index=False)
            print(f"  {model}: {k}/{len(ids)}", flush=True)
    pd.DataFrame(rows).to_parquet(dest, index=False)
    print(f"{model}: wrote {len(rows)} rows -> {dest}", flush=True)

if __name__ == "__main__":
    main()
