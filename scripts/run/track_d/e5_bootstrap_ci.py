import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

c = pd.read_csv("/12TBDrive1/mega_pep_bench/runs/track_d/e5_zeroshot/e5_confidence.csv")
ORIENT = {"iptm": 1, "ptm": 1, "complex_plddt": 1, "complex_iplddt": 1,
          "complex_pde": -1, "complex_ipde": -1}
chans = [x for x in c.columns if x in ORIENT]
for ch in chans:
    c[ch] = ORIENT[ch] * c[ch]
recs = c.receptor.unique()
rng = np.random.default_rng(0)

def within(frame, ch):
    vals = [roc_auc_score(g.label, g[ch]) for _, g in frame.groupby("receptor")
            if g.label.nunique() == 2]
    return float(np.mean(vals)) if vals else np.nan

rows = []
for ch in chans:
    pooled, wr = roc_auc_score(c.label, c[ch]), within(c, ch)
    bp, bw = [], []
    for _ in range(4000):
        pick = rng.choice(recs, size=len(recs), replace=True)
        sub = pd.concat([c[c.receptor == r] for r in pick])
        if sub.label.nunique() == 2:
            bp.append(roc_auc_score(sub.label, sub[ch]))
        w = within(sub, ch)
        if not np.isnan(w):
            bw.append(w)
    pl, ph = np.percentile(bp, [2.5, 97.5])
    wl, wh = np.percentile(bw, [2.5, 97.5])
    rows.append({"channel": ch, "pooled_auroc": round(pooled, 3),
                 "pooled_lo": round(pl, 3), "pooled_hi": round(ph, 3),
                 "within_receptor_auroc": round(wr, 3),
                 "within_lo": round(wl, 3), "within_hi": round(wh, 3),
                 "covers_half": bool(wl <= 0.5 <= wh)})
out = pd.DataFrame(rows)
out.to_csv("/12TBDrive1/mega_pep_bench/runs/track_d/e5_zeroshot/e5_bootstrap_ci.csv", index=False)
print(out.to_string(index=False))
print(f"\nchannels whose within-receptor interval covers 0.5: {int(out.covers_half.sum())} of {len(out)}")
