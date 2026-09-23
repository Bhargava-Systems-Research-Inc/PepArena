#!/usr/bin/env python3
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

PS = Path("/12TBDrive1/mega_pep_bench/methods/external/PepScorer/PepScorerRMSD/PepScorerRMSD")
OUT = Path("/12TBDrive1/mega_pep_bench/runs/track_b/pepscorer")
XS = ["XS_HPScore", "XS_HSScore"]
APBS = ["APBS_Ligand"]
GONE = XS + APBS
FF = ["CHARMM", "ElectDD"]
FF_FAIL_RATE = 1348 / 3981

X = pd.read_csv(PS / "Data/features_train_and_evaluation_set.csv")
meta = pd.read_excel(PS / "Data/Poses.xlsx")
assert len(X) == len(meta) == 6854, (len(X), len(meta))
assert all(c in X.columns for c in GONE) and X.shape[1] == 30, X.shape

y = meta.RMSD.to_numpy()
tr, ev = (meta.Split == "Train").to_numpy(), (meta.Split == "Evaluation").to_numpy()
published = joblib.load(PS / "objects/PepScorerRMSD.joblib")

def score(name, pred):
    d = pd.DataFrame({"protein": meta.Protein[ev].to_numpy(), "true": y[ev], "pred": pred})
    top = d.loc[d.groupby("protein").pred.idxmin()]
    per = [spearmanr(g.pred, g.true).statistic for _, g in d.groupby("protein") if len(g) > 2]
    r = dict(model=name, n_eval=int(ev.sum()), n_proteins=int(d.protein.nunique()),
             rmse=float(np.sqrt(np.mean((pred - y[ev]) ** 2))),
             pearson=float(np.corrcoef(pred, y[ev])[0, 1]),
             spearman_pooled=float(spearmanr(pred, y[ev]).statistic),
             spearman_within_protein=float(np.nanmean(per)),
             top1_under_2A=float((top.true < 2.0).mean()))
    print(f"{name:26s} RMSE {r['rmse']:.3f}  rho(pooled) {r['spearman_pooled']:.3f}  "
          f"rho(within) {r['spearman_within_protein']:.3f}  top1<2A {r['top1_under_2A']:.3f}")
    return r

print(f"evaluation split: {ev.sum()} poses over {meta.Protein[ev].nunique()} proteins\n")
rows = [score("published (30 features)", published.predict(X[ev]))]

Xn = X.copy()
Xn[XS] = np.nan
rows.append(score("published + XScore=NaN", published.predict(Xn[ev])))

Xna = X.copy()
Xna[GONE] = np.nan
rows.append(score("published + XScore,APBS=NaN", published.predict(Xna[ev])))

Xff = Xna.copy()
rng = np.random.default_rng(16)
hit = rng.random(len(Xff)) < FF_FAIL_RATE
Xff.loc[hit, FF] = np.nan
rows.append(score(f"+ force field NaN on {FF_FAIL_RATE:.0%}", published.predict(Xff[ev])))

X28 = X.drop(columns=XS)
m28 = HistGradientBoostingRegressor(random_state=16).fit(X28[tr], y[tr])
rows.append(score("retrained (28 features)", m28.predict(X28[ev])))

X27 = X.drop(columns=GONE)
m27 = HistGradientBoostingRegressor(random_state=16).fit(X27[tr], y[tr])
rows.append(score("retrained (27 features)", m27.predict(X27[ev])))

base = rows[0]
for r in rows[1:]:
    r["top1_delta_vs_published"] = round(r["top1_under_2A"] - base["top1_under_2A"], 4)
    r["rho_within_delta_vs_published"] = round(
        r["spearman_within_protein"] - base["spearman_within_protein"], 4)

joblib.dump(m28, OUT / "PepScorerRMSD_28feature.joblib")
joblib.dump(m27, OUT / "PepScorerRMSD_27feature.joblib")
(OUT / "missing_feature_ablation.json").write_text(json.dumps(rows, indent=2) + "\n")
pd.DataFrame(rows).to_csv(OUT / "missing_feature_ablation.csv", index=False)
print(f"\nwrote {OUT}/missing_feature_ablation.{{json,csv}} and the reduced models")
