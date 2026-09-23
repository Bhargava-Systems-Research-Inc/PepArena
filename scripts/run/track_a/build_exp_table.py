from pathlib import Path
import pandas as pd, numpy as np
ROOT = Path("/12TBDrive1/mega_pep_bench"); SC = ROOT / "runs/track_a/scores"

def load(m):
    p = SC / f"{m}.parquet"
    return pd.read_parquet(p).dropna(subset=["DockQ"]).set_index("complex_id") if p.exists() else None

def paired(b, f, col="DockQ"):
    j = pd.concat([b[col].rename("base"), f[col].rename("ft")], axis=1).dropna()
    d = j.ft - j.base
    return len(j), j.base.mean(), j.ft.mean(), d.mean(), (d > 0.01).sum(), (d < -0.01).sum()

configs = [("loops20 / 1 sample  (reference)", "esmfold2_local", "esmfold2_ft"),
           ("loops8  / 1 sample  (Exp1)",      "esmfold2_local_l8", "esmfold2_ft_l8"),
           ("loops20 / 5 samp ipTM (Exp2)",    "esmfold2_local_5s", "esmfold2_ft_5s")]
out = ["# Base vs fine-tuned ESMFold2 across recipes (paired, DockQ --capri_peptide)", "",
       "| Config | n | base DockQ | ft DockQ | Δ(ft-base) | ft>base | base>ft | base fnat | ft fnat |",
       "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
for name, bm, fm in configs:
    b, f = load(bm), load(fm)
    if b is None or f is None:
        out.append(f"| {name} | missing ({bm if b is None else fm}) |"); continue
    n, bmn, fmn, dl, w, l = paired(b, f)
    _, bf, ff, _, _, _ = paired(b, f, "fnat")
    out.append(f"| {name} | {n} | {bmn:.3f} | {fmn:.3f} | {dl:+.3f} | {w} | {l} | {bf:.3f} | {ff:.3f} |")

b20, b8, b5 = load("esmfold2_local"), load("esmfold2_local_l8"), load("esmfold2_local_5s")
out += ["", "## Recipe effect on BASE ESMFold2 (paired)"]
if b8 is not None and b20 is not None:
    j = pd.concat([b8.DockQ.rename("l8"), b20.DockQ.rename("l20")], axis=1).dropna()
    out.append(f"- loops 8 -> 20: {j.l8.mean():.3f} -> {j.l20.mean():.3f} (Δ {(j.l20-j.l8).mean():+.3f}, n={len(j)})")
if b5 is not None and b20 is not None:
    j = pd.concat([b20.DockQ.rename("s1"), b5.DockQ.rename("s5")], axis=1).dropna()
    out.append(f"- 1 -> 5 samples (ipTM): {j.s1.mean():.3f} -> {j.s5.mean():.3f} (Δ {(j.s5-j.s1).mean():+.3f}, n={len(j)})")

txt = "\n".join(out) + "\n"
(ROOT / "runs/track_a/EXP_base_vs_ft.md").write_text(txt)
print(txt)
