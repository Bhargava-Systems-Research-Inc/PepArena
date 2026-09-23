import json
import subprocess
import sys
import tempfile
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
DOCKQ = str(Path.home() / "anaconda3/envs/dockq/bin/DockQ")
sys.path.insert(0, str(R / "scripts/run/track_a"))
import importlib.util
spec = importlib.util.spec_from_file_location("sd", R / "scripts/run/track_a/score_dockq.py")
sd = importlib.util.module_from_spec(spec); spec.loader.exec_module(sd)

def run(model, native):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
        out = tf.name
    subprocess.run([DOCKQ, str(model), str(native), "--capri_peptide",
                    "--allowed_mismatches", "10", "--json", out, "--short"],
                   capture_output=True, text=True, timeout=900)
    try:
        return json.loads(Path(out).read_text())
    except Exception:
        return None

def main(model_name="protenix", tier="canonical", limit=40):
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    scores = pd.read_parquet(R / f"runs/track_a/scores/{model_name}.parquet")
    sd.PRED_BASE = R / "runs/track_a/predictions"
    NATIVE = R / "data/curated/structures/native"

    rows = []
    for cid in scores.complex_id:
        if len(rows) >= limit:
            break
        nat = NATIVE / f"{cid}.cif"
        pred = sd.find_pred(model_name, cid)
        if not pred or not nat.exists():
            continue
        d = run(pred, nat)
        if not d:
            continue
        best = d.get("best_result") or {}
        if len(best) < 2:
            continue
        g = d.get("GlobalDockQ")
        per = {k: v.get("DockQ") for k, v in best.items()}
        first = next(iter(best.values()), {})
        rows.append({"complex_id": cid, "n_interfaces": len(best), "GlobalDockQ": g,
                     "first_iface_DockQ": first.get("DockQ"),
                     "min_iface": min(per.values()), "max_iface": max(per.values()),
                     "interfaces": ";".join(f"{k}={v:.3f}" for k, v in per.items())})
        print(f"  {cid}: {len(best)} interfaces, global={g:.3f}, per-interface {per}", flush=True)

    if not rows:
        print("no multi-interface complexes found in the sample")
        return
    d = pd.DataFrame(rows)
    d.to_csv(R / "runs/track_a/dockq_interface_audit.csv", index=False)
    d["gap"] = (d.GlobalDockQ - d.first_iface_DockQ).abs()
    print(f"\nmulti-interface complexes sampled: {len(d)}")
    print(f"  |GlobalDockQ - first interface| : median {d.gap.median():.3f} max {d.gap.max():.3f}")
    print(f"  complexes where they differ >0.05: {(d.gap > 0.05).sum()}")
    print(f"  spread across interfaces (max-min): median {(d.max_iface - d.min_iface).median():.3f}")

if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else "protenix", a[1] if len(a) > 1 else "canonical", int(a[2]) if len(a) > 2 else 40)
