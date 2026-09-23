import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
MODELS = {"protenix": "Protenix", "boltz2": "Boltz-2", "af3": "AlphaFold3",
          "highfold3": "HighFold3"}
CHEMS = ["disulfide", "ester", "head_to_tail", "lactam", "staple", "thioether"]

man = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
rows, counts = [], {}
for key, lab in MODELS.items():
    f = R / f"runs/track_a/scores_pep_ncaa/{key}.parquet"
    if not f.exists():
        continue
    sc = pd.read_parquet(f).set_index("complex_id").join(man[["cyclization_type"]])
    sc = sc[sc.cyclization_type.isin(CHEMS)]
    g = sc.groupby("cyclization_type").DockQ.agg(["count", "mean", "median"])
    rows.append({"model": lab, **{c: round(float(g["mean"][c]), 3) if c in g.index else None
                                 for c in CHEMS}})
    if key == "protenix":
        counts = {c: g.loc[c] for c in g.index}

chem = pd.DataFrame(rows)
chem.loc[len(chem)] = ["n (complexes)"] + [int(counts[c]["count"]) if c in counts else 0
                                           for c in CHEMS]
chem.to_csv(R / "runs/track_a/ncaa_chemistry_by_model.csv", index=False)

pd.DataFrame([{"cyclization_type": c, "count": int(v["count"]),
               "mean": round(float(v["mean"]), 3), "median": round(float(v["median"]), 3)}
              for c, v in sorted(counts.items())]).to_csv(
    R / "runs/track_a/ncaa_chemistry_counts.csv", index=False)

print(chem.to_string(index=False))
print(f"\nwrote ncaa_chemistry_by_model.csv and ncaa_chemistry_counts.csv")
