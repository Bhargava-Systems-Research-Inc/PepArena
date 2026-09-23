import importlib.util
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
spec = importlib.util.spec_from_file_location("sp2", R / "scripts/run/track_a/score_peptide_dockq.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

ARCH = R / "runs/track_a/predictions_ncaa/_highfold3_nocycpoem_ARCHIVED"
NATIVE = R / "data/curated/structures/native"
sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")

rows = []
for d in sorted(ARCH.iterdir()):
    if not d.is_dir():
        continue
    cid = d.name
    pred = None
    for name in (cid, cid.lower()):
        for cand in (d / f"{name}_model.cif", d / f"{name}_model_0.cif"):
            if cand.exists():
                pred = cand
                break
        if pred:
            break
    nat = NATIVE / f"{cid}.cif"
    if not pred or not nat.exists():
        continue
    plen = int(sp.peptide_len.get(cid, 10))
    s = m.score(pred, nat, plen) or m.score_fallback(pred, nat, plen)
    if s is None:
        continue
    rows.append({"complex_id": cid, **s})

out = pd.DataFrame(rows)
out.to_parquet(R / "runs/track_a/scores_ncaa/highfold3_nocyc.parquet", index=False)
print(f"scored {len(out)} archived predictions; mean DockQ {out.DockQ.mean():.3f}")
