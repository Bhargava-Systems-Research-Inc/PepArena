from pathlib import Path
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
OUT = R / "dist/PepArena_Proteins_Submission/zenodo"
FRESH = R / "runs/track_b/fresh"
PROG = {"haddock3": "HADDOCK3", "adcp": "ADCP", "unidock": "Uni-Dock", "rapidock": "RAPiDock"}


def docking():
    cov = pd.read_csv(FRESH / "coverage.csv").set_index("complex_id")
    rows = []
    for key, lab in PROG.items():
        f = FRESH / f"{key}_scored.csv"
        sc = pd.read_csv(f).set_index("complex_id") if f.exists() else pd.DataFrame()
        for cid, c in cov.iterrows():
            r = {"complex_id": cid, "program": lab,
                 "site": "blind" if key == "haddock3" else "native_box",
                 "accepted": c.get(key) == "ok", "refusal_reason": None if c.get(key) == "ok" else c.get(key),
                 "peptide_len": c.peptide_len, "cyclization_type": c.cyclization_type,
                 "has_ncaa": c.has_ncaa, "torsdof": c.torsdof}
            if cid in sc.index:
                s = sc.loc[cid]
                r["top1_dockq"] = s.get("top1_dockq")
                r["best_dockq"] = s.get("best_dockq")
                r["n_poses"] = s.get("n_poses")
            rows.append(r)
    d = pd.DataFrame(rows).sort_values(["program", "complex_id"])
    for c in d.select_dtypes("float").columns:
        d[c] = d[c].round(4)
    d.to_csv(OUT / "PepArena_docking.csv", index=False)
    return d


def affinity():
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    t = pd.read_csv(R / "runs/track_c/affinity_targets_v2.csv")
    b = R / "runs/track_c/boltz2_affinity_v2_scored.csv"
    if b.exists():
        sc = pd.read_csv(b)[["complex_id", "pred", "log_affinity", "heavy_atoms", "set"]]
        t = t.merge(sc, on="complex_id", how="left", suffixes=("", "_boltz"))
    m = sp.reindex(t.complex_id)[["peptide_len", "is_cyclic", "cyclization_type", "has_ncaa"]]
    t = pd.concat([t.reset_index(drop=True), m.reset_index(drop=True)], axis=1)
    for c in t.select_dtypes("float").columns:
        t[c] = t[c].round(4)
    t.to_csv(OUT / "PepArena_affinity.csv", index=False)
    return t


def interaction():
    d = pd.read_csv(R / "runs/track_d/leaderboard_v2.csv")
    for c in d.select_dtypes("float").columns:
        d[c] = d[c].round(4)
    d.to_csv(OUT / "PepArena_interaction.csv", index=False)
    return d


if __name__ == "__main__":
    for name, fn in (("docking", docking), ("affinity", affinity), ("interaction", interaction)):
        try:
            t = fn()
            print(f"PepArena_{name}.csv  {len(t)} rows x {len(t.columns)} cols")
        except Exception as e:
            print(f"{name}: {type(e).__name__}: {e}")
