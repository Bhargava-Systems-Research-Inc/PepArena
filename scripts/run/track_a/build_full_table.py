from pathlib import Path
import pandas as pd, numpy as np

ROOT = Path("/12TBDrive1/mega_pep_bench")
SC = ROOT / "runs/track_a/scores"
LABEL = {"protenix": "Protenix (v1.0.0)", "boltz2": "Boltz-2", "chai1": "Chai-1",
         "esmfold2": "ESMFold2 (ESMC-6B, API)",
         "af3": "AlphaFold3", "helixfold3": "HelixFold3", "highfold3": "HighFold3",
         "rfaa": "RoseTTAFold-AA",
         "esmfold2_local": "ESMFold2 (local exx, this run)",
         "esmfold2_ft": "ESMFold2 fine-tuned (this work)"}
ORDER = ["protenix", "af3", "boltz2", "chai1", "highfold3", "helixfold3", "rfaa",
         "esmfold2", "esmfold2_local", "esmfold2_ft"]
LBINS, LNAMES = [0, 5, 10, 15, 25, 60], ["<=5", "6-10", "11-15", "16-25", "26+"]

loaded = {m: pd.read_parquet(SC / f"{m}.parquet").dropna(subset=["DockQ"])
          for m in ORDER if (SC / f"{m}.parquet").exists()}
FULL = {m: d for m, d in loaded.items() if len(d) >= 250}
common = set.intersection(*[set(d.complex_id) for d in FULL.values()]) if FULL else set()
PARTIAL = {m for m, d in loaded.items() if len(d) < 250}

def stats(dq):
    dq = dq.dropna()
    return dict(n=len(dq), mean=round(dq.mean(), 3), median=round(dq.median(), 3),
                s23=round((dq >= 0.23).mean(), 3), a49=round((dq >= 0.49).mean(), 3),
                h80=round((dq >= 0.80).mean(), 3))

rows = []
for m, d in loaded.items():
    s = stats(d.DockQ)
    cyc = d.groupby("is_cyclic").DockQ.mean().round(3)
    dc = d[d.complex_id.isin(common)]
    nm = LABEL.get(m, m) + (" [PARTIAL, still folding]" if m in PARTIAL else "")
    rows.append({"model": nm, "key": m, **s,
                 "linear": cyc.get(False, np.nan), "cyclic": cyc.get(True, np.nan),
                 "iRMSD": round(d.iRMSD.dropna().mean(), 2), "LRMSD": round(d.LRMSD.dropna().mean(), 2),
                 "fnat": round(d.fnat.dropna().mean(), 3),
                 "common_n": len(dc), "common_mean": round(dc.DockQ.mean(), 3)})
lb = pd.DataFrame(rows).sort_values("common_mean", ascending=False)

L = [f"# Track A cofolding leaderboard — updated with local + fine-tuned ESMFold2", "",
     "Canonical tier, post-2023-06-01, leakage-clean. num_loops=20/steps=100, best sample, "
     "DockQ `--capri_peptide`. Local + fine-tuned ESMFold2 folded on exx (bf16, single-seq, no MSA).", "",
     "| Model | n | mean DockQ | median | succ@0.23 | acc@0.49 | high@0.80 | linear | cyclic | iRMSD | LRMSD | fnat |",
     "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
for _, r in lb.iterrows():
    L.append(f"| {r.model} | {r.n} | {r['mean']} | {r['median']} | {r.s23} | {r.a49} | {r.h80} "
             f"| {r.linear} | {r.cyclic} | {r.iRMSD} | {r.LRMSD} | {r.fnat} |")
L += ["", f"## Common set — {len(common)} complexes ALL listed models predicted", "",
      "| Model | mean DockQ (common) |", "|---|--:|"]
for _, r in lb.iterrows():
    L.append(f"| {r.model} | {r.common_mean} |")
L += ["", "## Mean DockQ by peptide length", "",
      "| Model | " + " | ".join(LNAMES) + " |", "|---|" + "--:|" * len(LNAMES)]
for m in ORDER:
    if m not in loaded: continue
    d = loaded[m].assign(b=pd.cut(loaded[m].peptide_len, LBINS, labels=LNAMES))
    by = d.groupby("b", observed=True).DockQ.mean().round(2)
    L.append(f"| {LABEL.get(m,m)} | " + " | ".join(str(by.get(n, '-')) for n in LNAMES) + " |")

if "esmfold2_local" in loaded and "esmfold2_ft" in loaded:
    b = loaded["esmfold2_local"].set_index("complex_id").DockQ
    f = loaded["esmfold2_ft"].set_index("complex_id").DockQ
    j = pd.concat([b.rename("base"), f.rename("ft")], axis=1).dropna()
    d = j.ft - j.base
    L += ["", f"## Paired local base vs fine-tuned (n={len(j)}, identical exx path)", "",
          f"- mean DockQ: base {j.base.mean():.3f} -> ft {j.ft.mean():.3f} (delta {d.mean():+.3f})",
          f"- ft better: {(d>0.01).sum()}, base better: {(d<-0.01).sum()}, tie: {((d.abs()<=0.01)).sum()}",
          f"- median delta {d.median():+.3f}"]

(ROOT / "runs/track_a/LEADERBOARD_with_ft.md").write_text("\n".join(L) + "\n")
lb.to_csv(ROOT / "runs/track_a/leaderboard_with_ft.csv", index=False)
print(lb[["model", "n", "mean", "median", "s23", "a49", "h80", "linear", "cyclic", "common_mean"]].to_string(index=False))
print(f"\ncommon={len(common)} -> LEADERBOARD_with_ft.md")
