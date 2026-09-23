#!/usr/bin/env python3
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO = Path("/12TBDrive1/mega_pep_bench")
E5 = REPO / "runs/track_d/e5_zeroshot"
CHANNELS = {"iptm": 1, "ptm": 1, "complex_plddt": 1, "complex_iplddt": 1,
            "complex_pde": -1, "complex_ipde": -1, "ligand_iptm": 1, "protein_iptm": 1}

def main():
    tgt = pd.read_csv(E5 / "targets.csv")
    key = {r.complex_id.replace("/", "_"): r for _, r in tgt.iterrows()}
    rows = []
    for f in sorted((E5 / "out").rglob("confidence_*.json")):
        stem = f.stem[len("confidence_"):].rsplit("_model_", 1)[0]
        if stem not in key:
            continue
        c = json.loads(f.read_text())
        r = key[stem]
        rows.append({"complex_id": r.complex_id, "receptor": r.receptor_pdb_id,
                     "label": int(r.label), "peptide_len": len(str(r.peptide_seq)),
                     **{k: c.get(k) for k in CHANNELS}})
    if not rows:
        print("no confidence files found -- has e5_run_boltz.sh finished?")
        return
    d = pd.DataFrame(rows)
    d.to_csv(E5 / "e5_confidence.csv", index=False)

    summary = {"n_pairs": int(len(d)), "n_receptors": int(d.receptor.nunique()),
               "n_binders": int(d.label.sum()), "n_nonbinders": int((d.label == 0).sum()),
               "channels": {}}
    for ch, sign in CHANNELS.items():
        v = d[["receptor", "label", ch]].dropna()
        if v[ch].nunique() < 2 or v.label.nunique() < 2:
            continue
        per = [roc_auc_score(g.label, sign * g[ch])
               for _, g in v.groupby("receptor") if g.label.nunique() == 2]
        summary["channels"][ch] = {
            "pooled_auroc": round(float(roc_auc_score(v.label, sign * v[ch])), 3),
            "within_receptor_auroc": round(float(np.mean(per)), 3),
            "n_receptors_scored": len(per),
            "mean_binder": round(float(v.loc[v.label == 1, ch].mean()), 3),
            "mean_nonbinder": round(float(v.loc[v.label == 0, ch].mean()), 3),
        }
    (E5 / "e5_zeroshot.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {E5}/e5_confidence.csv and e5_zeroshot.json")

if __name__ == "__main__":
    main()
