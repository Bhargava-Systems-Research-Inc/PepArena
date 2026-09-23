import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix", "af3": "AlphaFold3", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3", "rfaa": "RoseTTAFold-AA",
       "esmfold2": "ESMFold2"}
BINS, NAMES = [0, 5, 10, 15, 25, 60], ["<=5", "6-10", "11-15", "16-25", "26+"]

def load(tier):
    d = R / "runs/track_a" / ("scores_pep" if tier == "canonical" else "scores_pep_ncaa")
    return {k: pd.read_parquet(d / f"{k}.parquet").set_index("complex_id")
            for k in LAB if (d / f"{k}.parquet").exists()}

def main():
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet")
    sp = sp.set_index("complex_id")
    bond = None
    ba = R / "runs/track_a/ncaa_cyclic_bond_audit.csv"
    if ba.exists():
        bond = pd.read_csv(ba).set_index("complex_id")["af3_pairs"]

    for tier in ("canonical", "ncaa"):
        sc = load(tier)
        if not sc:
            print(f"[{tier}] no corrected scores yet"); continue
        print(f"\n===== {tier} tier ({len(sc)} models) =====")

        rows = []
        for m, d in sc.items():
            f = sp.reindex(d.index)
            b = pd.cut(f.peptide_len, BINS, labels=NAMES)
            r = {"model": LAB.get(m, m)}
            for n in NAMES:
                sel = (b == n).values
                r[n] = round(d.DockQ[sel].mean(), 3) if sel.sum() else None
                r[f"n_{n}"] = int(sel.sum())
            rows.append(r)
        L = pd.DataFrame(rows)
        L.to_csv(R / f"runs/track_a/length_bins_v2_{tier}.csv", index=False)
        print("\n-- mean DockQ by peptide-length bin")
        print(L.to_string(index=False))
        for r in rows:
            vals = [r[n] for n in NAMES if r[n] is not None and r[f"n_{n}"] >= 3]
            mono = all(x >= y for x, y in zip(vals, vals[1:]))
            print(f"   {r['model']:<12} monotonic over bins with n>=3: {mono}")

        rows = []
        for m, d in sc.items():
            f = sp.reindex(d.index)
            cy = f.is_cyclic.fillna(False).values
            rows.append({"model": LAB.get(m, m),
                         "linear_mean": round(d.DockQ[~cy].mean(), 3), "n_linear": int((~cy).sum()),
                         "cyclic_mean": round(d.DockQ[cy].mean(), 3) if cy.sum() else None,
                         "n_cyclic": int(cy.sum())})
        C = pd.DataFrame(rows)
        C.to_csv(R / f"runs/track_a/linear_cyclic_v2_{tier}.csv", index=False)
        print("\n-- linear vs cyclic")
        print(C.to_string(index=False))

        rows = []
        for m, d in sc.items():
            f = sp.reindex(d.index)
            for code, name in (("C", "coil"), ("H", "helix"), ("E", "sheet")):
                sel = (f.sec_struct == code).values
                rows.append({"model": LAB.get(m, m), "sec_struct": name, "n": int(sel.sum()),
                             "mean_dockq": round(d.DockQ[sel].mean(), 3) if sel.sum() else None})
        S = pd.DataFrame(rows)
        S.to_csv(R / f"runs/track_a/sec_struct_v2_{tier}.csv", index=False)
        print("\n-- mean DockQ by peptide secondary structure")
        print(S.pivot_table(index="sec_struct", columns="model", values="mean_dockq").round(3)
              .to_string())

        rows = []
        for m, d in sc.items():
            f = sp.reindex(d.index)
            for cls, g in d.join(f[["target_class"]]).groupby("target_class"):
                rows.append({"model": LAB.get(m, m), "target_class": cls, "n": len(g),
                             "mean_dockq": round(g.DockQ.mean(), 3)})
        T = pd.DataFrame(rows)
        T.to_csv(R / f"runs/track_a/target_class_v2_{tier}.csv", index=False)
        print("\n-- mean DockQ by receptor target class (n>=15 shown)")
        big = T.groupby("target_class").n.min()
        print(T[T.target_class.isin(big[big >= 15].index)]
              .pivot_table(index="target_class", columns="model", values="mean_dockq")
              .round(3).to_string())

        if tier == "ncaa" and bond is not None:
            rows = []
            for m, d in sc.items():
                f = sp.reindex(d.index)
                cy = f.is_cyclic.fillna(False)
                sub = d[cy.values].join(f[["cyclization_type"]]).join(bond.rename("pairs"))
                for chem, g in sub.groupby("cyclization_type"):
                    for has, gg in g.groupby(g.pairs.fillna(0) > 0):
                        rows.append({"model": LAB.get(m, m), "chemistry": chem,
                                     "ring_in_input": bool(has), "n": len(gg),
                                     "mean_dockq": round(gg.DockQ.mean(), 3)})
            B = pd.DataFrame(rows)
            B.to_csv(R / "runs/track_a/chemistry_by_bond_v2.csv", index=False)
            print("\n-- ring chemistry, split by whether the ring bond was in the input")
            piv = B.pivot_table(index=["chemistry", "ring_in_input"], columns="model",
                                values="mean_dockq")
            print(piv.round(3).to_string())
            print("\n   counts:")
            print(B.pivot_table(index=["chemistry", "ring_in_input"], columns="model",
                                values="n").to_string())

if __name__ == "__main__":
    main()
