import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix", "af3": "AlphaFold3", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3", "rfaa": "RoseTTAFold-AA",
       "esmfold2": "ESMFold2"}

feat = pd.read_parquet(R / "data/curated/complex_features.parquet").set_index("complex_id")
lk = pd.read_parquet(R / "data/splits/leakage.parquet").set_index("complex_id")
ncaa = pd.read_parquet(R / "data/curated/structure_pdb_ncaa.parquet")
ncaa_ids = set(ncaa.complex_id) if "complex_id" in ncaa.columns else set()

def ols_cluster(X, y, g):
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ (X.T @ y)
    e = y - X @ b
    meat = np.zeros((X.shape[1], X.shape[1]))
    for c in np.unique(g):
        m = g == c
        s = X[m].T @ e[m]
        meat += np.outer(s, s)
    V = XtX_inv @ meat @ XtX_inv
    return b, np.sqrt(np.diag(V))

rows = []
for k, lab in LAB.items():
    f = R / f"runs/track_a/scores/{k}.parquet"
    if not f.exists():
        continue
    d = pd.read_parquet(f).set_index("complex_id")
    d = d.join(feat[["n_interface_res", "frac_helix"]]).join(lk[["receptor_cluster"]])
    d = d.dropna(subset=["DockQ", "peptide_len", "n_interface_res"])
    d["ncaa"] = d.index.isin(ncaa_ids).astype(float)
    d["cyc"] = d.is_cyclic.astype(float)
    cols = [c for c in ["peptide_len", "cyc", "ncaa", "n_interface_res", "frac_helix"]
            if d[c].astype(float).std() > 0]
    X = np.column_stack([np.ones(len(d))] + [d[c].astype(float).to_numpy() for c in cols])
    for j, c in enumerate(cols, start=1):
        if c in ("peptide_len", "n_interface_res"):
            X[:, j] = (X[:, j] - X[:, j].mean()) / X[:, j].std()
    y = d.DockQ.to_numpy(float)
    g = d.receptor_cluster.fillna("none").to_numpy()
    b, se = ols_cluster(X, y, g)
    uni = np.column_stack([np.ones(len(d)), X[:, 1]])
    bu, seu = ols_cluster(uni, y, g)
    rows.append({"model": lab, "n": len(d), "clusters": len(np.unique(g)),
                 "len_alone": round(bu[1], 4), "len_alone_se": round(seu[1], 4),
                 "len_adj": round(b[1], 4), "len_adj_se": round(se[1], 4),
                 "len_adj_t": round(b[1] / se[1], 2),
                 **{c: round(b[j], 3) for j, c in enumerate(cols, start=1) if c != "peptide_len"},
                 "covariates": "+".join(cols)})

out = pd.DataFrame(rows)
out.to_csv(R / "runs/track_a/length_multivariable.csv", index=False)
print(out.to_string(index=False))
print("\ncoefficients for peptide_len and n_interface_res are per standard deviation of DockQ units")
print(f"models where the adjusted length coefficient stays negative with |t| > 2: "
      f"{int(((out.len_adj < 0) & (out.len_adj_t.abs() > 2)).sum())} of {len(out)}")
