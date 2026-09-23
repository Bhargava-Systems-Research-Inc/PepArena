#!/usr/bin/env python3
import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr, pearsonr

REPO = Path("/12TBDrive1/mega_pep_bench")
MM = REPO / "runs/track_c/mmgbsa"
OUT = REPO / "runs/track_c"

def parse_delta_total(dat: Path):
    val = None
    for ln in dat.read_text(errors="ignore").splitlines():
        if ln.startswith("DELTA TOTAL"):
            m = re.findall(r"[-+]?\d+\.\d+", ln)
            if m:
                val = float(m[0])
    return val

def main():
    truth = pd.read_parquet(REPO / "data/curated/affinity_mmpbgbsa_cp.parquet")
    truth["pdb"] = truth.receptor_pdb_id.str.upper()
    keep = ["pdb", "complex_id", "log_affinity", "cyclization_type", "peptide_len", "has_ncaa"]
    truth = truth[keep]

    rows = []
    for d in sorted(MM.iterdir()):
        dat = d / "out.dat"
        if not dat.is_file():
            continue
        v = parse_delta_total(dat)
        rows.append({"pdb": d.name.upper(), "dG_mmgbsa": v})
    got = pd.DataFrame(rows)
    df = truth.merge(got, on="pdb", how="inner").dropna(subset=["dG_mmgbsa", "log_affinity"])
    df.to_csv(OUT / "mmgbsa_scored.csv", index=False)

    rho, prho = spearmanr(df.dG_mmgbsa, df.log_affinity)
    r, pr = pearsonr(df.dG_mmgbsa, df.log_affinity)
    per_chem = {}
    for chem, g in df.groupby("cyclization_type"):
        if len(g) >= 4:
            cr, cp = spearmanr(g.dG_mmgbsa, g.log_affinity)
            per_chem[chem] = {"n": int(len(g)), "spearman": round(float(cr), 3),
                              "p": round(float(cp), 4)}

    summary = {
        "method": "MM-GBSA (single minimized structure, igb=2)",
        "set": "MM_PBGBSA-CP dataset I (cyclic peptide-protein, measured Kd)",
        "n_complexes": int(len(df)),
        "spearman": round(float(rho), 3), "spearman_p": round(float(prho), 4),
        "pearson": round(float(r), 3), "pearson_p": round(float(pr), 4),
        "spearman_within_chemistry": per_chem,
        "note": ("Pooled rho is inflated by between-chemistry offsets; the within-chemistry "
                 "values are what a ranking user would experience."),
    }
    json.dump(summary, open(OUT / "mmgbsa_leaderboard.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {OUT/'mmgbsa_scored.csv'} ({len(df)} rows)")

if __name__ == "__main__":
    main()
