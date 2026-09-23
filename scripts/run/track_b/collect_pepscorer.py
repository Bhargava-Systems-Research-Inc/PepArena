#!/usr/bin/env python3
import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

REPO = Path("/12TBDrive1/mega_pep_bench")
OUT = REPO / "runs/track_b/pepscorer"
ACC = 0.23

def main():
    truth = pd.read_parquet(REPO / "runs/track_b/scoring/dockq_all.parquet")
    rows = []
    for pred_f in sorted(OUT.glob("*/predicted.csv")):
        pdb = pred_f.parent.name
        p = pd.read_csv(pred_f)
        p["rank"] = p["Name"].map(lambda s: int(re.sub(r"\D", "", str(s).split("_")[-1])))
        t = truth[truth.pdb == pdb][["rank", "dockq"]]
        j = p.merge(t, on="rank", how="inner").dropna(subset=["Predicted RMSD", "dockq"])
        if len(j) < 5:
            continue
        adcp_top = j.loc[j["rank"].idxmin(), "dockq"]
        ps_top = j.loc[j["Predicted RMSD"].idxmin(), "dockq"]
        rho_ps, _ = spearmanr(-j["Predicted RMSD"], j["dockq"])
        rho_adcp, _ = spearmanr(-j["rank"], j["dockq"])
        rows.append({"pdb": pdb, "n_poses": len(j), "best_dockq": j.dockq.max(),
                     "adcp_top1_dockq": adcp_top, "pepscorer_top1_dockq": ps_top,
                     "pepscorer_rho": rho_ps, "adcp_rho": rho_adcp})
    if not rows:
        print("no PepScorer predictions found"); return
    d = pd.DataFrame(rows)
    d.to_csv(REPO / "runs/track_b/pepscorer_scored.csv", index=False)
    summary = {
        "method": "PepScorer::RMSD re-ranking of ADCP decoys",
        "set": "MM_PBGBSA-CP dataset2, complexes within PepScorer's 3-10 residue domain",
        "n_complexes": int(len(d)),
        "sampling_power_acceptable": round(float((d.best_dockq >= ACC).mean()), 3),
        "docking_power_adcp_ranking": round(float((d.adcp_top1_dockq >= ACC).mean()), 3),
        "docking_power_pepscorer_ranking": round(float((d.pepscorer_top1_dockq >= ACC).mean()), 3),
        "mean_top1_dockq_adcp": round(float(d.adcp_top1_dockq.mean()), 3),
        "mean_top1_dockq_pepscorer": round(float(d.pepscorer_top1_dockq.mean()), 3),
        "scoring_power_pepscorer_rho": round(float(d.pepscorer_rho.mean()), 3),
        "scoring_power_adcp_rho": round(float(d.adcp_rho.mean()), 3),
    }
    json.dump(summary, open(REPO / "runs/track_b/pepscorer_leaderboard.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
