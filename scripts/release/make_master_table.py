from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix", "af3": "AlphaFold3", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3", "rfaa": "RoseTTAFold-AA",
       "esmfold2": "ESMFold2"}
META = ["receptor_pdb_id", "peptide_chain", "peptide_seq", "peptide_len", "length_bin",
        "is_cyclic", "cyclization_type", "has_ncaa", "sec_struct", "target_class",
        "n_interface_res", "resolution", "deposition_date"]

def capri(d):
    return pd.cut(d, [-0.01, 0.23, 0.49, 0.80, 1.01],
                  labels=["Incorrect", "Acceptable", "Medium", "High"])

def main():
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    rows = []
    for tier, sub in (("canonical", "scores_pep"), ("ncAA", "scores_pep_ncaa")):
        for key, lab in LAB.items():
            f = R / f"runs/track_a/{sub}/{key}.parquet"
            if not f.exists():
                continue
            d = pd.read_parquet(f).set_index("complex_id")
            m = sp.reindex(d.index)[META]
            t = d.join(m)
            t["tier"] = tier
            t["model"] = lab
            tn = R / f"runs/track_a/scores_topn/{key}.parquet"
            if tier == "canonical" and tn.exists():
                s = pd.read_parquet(tn).set_index("complex_id")
                t = t.join(s[["oracle", "worst", "n_samples"]].rename(
                    columns={"oracle": "oracle_top5_dockq", "worst": "worst_sample_dockq"}))
            pv = R / f"runs/track_a/pose_validity/{key}.csv"
            if tier == "canonical" and pv.exists():
                v = pd.read_csv(pv)
                v = v[v.prepared.astype(str) == "True"].set_index("complex_id")
                checks = [c for c in v.columns
                          if c not in ("complex_id", "prepared", "error",
                                       "minimum_distance_to_protein")]
                ok = v[checks].apply(lambda s: s.map({True: True, False: False,
                                                      "True": True, "False": False}))
                t = t.join(ok.all(axis=1).rename("pose_passes_validity"))
            rows.append(t.reset_index())

    cc = R / "runs/track_a/confidence_components.csv"
    M = pd.concat(rows, ignore_index=True)
    if cc.exists():
        c = pd.read_csv(cc)
        c["model"] = c.model.map(LAB)
        M = M.merge(c[["complex_id", "model", "score", "iptm", "iptm_pep", "ptm"]],
                    on=["complex_id", "model"], how="left")
        M = M.rename(columns={"score": "model_confidence", "iptm_pep": "iptm_peptide_receptor"})

    M = M.drop(columns=["fnonnat", "F1", "n_interfaces", "route"], errors="ignore")
    M.insert(M.columns.get_loc("DockQ") + 1, "capri_class", capri(M.DockQ).astype(str))
    front = ["complex_id", "tier", "model"] + META
    M = M[front + [c for c in M.columns if c not in front]]
    M = M.sort_values(["tier", "model", "complex_id"])
    for c in M.select_dtypes("float").columns:
        M[c] = M[c].round(4)

    for dest in (R / "dist/zenodo_upload/PepArena_scores.csv",
                 R / "dist/PepArena_Proteins_Submission/zenodo/PepArena_scores.csv"):
        dest.parent.mkdir(parents=True, exist_ok=True)
        M.to_csv(dest, index=False)
    print(f"{len(M)} rows x {len(M.columns)} cols -> PepArena_scores.csv")
    print(f"  tiers: {dict(M.tier.value_counts())}")
    print(f"  models: {M.model.nunique()}   complexes: {M.complex_id.nunique()}")
    print("  columns:", ", ".join(M.columns))

if __name__ == "__main__":
    main()
