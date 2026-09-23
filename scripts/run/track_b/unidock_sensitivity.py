import json
import re
import numpy as np
import pandas as pd
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
CAPS = (15, 20, 25, 32)

cov = pd.read_csv(R / "runs/track_b/fresh/coverage.csv")
cov = cov[cov.unidock == "ok"][["complex_id", "torsdof"]]
sc = pd.read_csv(R / "runs/track_b/fresh/unidock_scored.csv")
d = cov.merge(sc, on="complex_id", how="left")

RUN = R / "runs/track_b/fresh/unidock"

def outcome(cid):
    dd = RUN / cid
    if not dd.is_dir():
        return "not attempted"
    if list(dd.glob("*_out_ranked_1.pdb")):
        return "scored"
    log = dd / "dock.log"
    if not log.exists():
        return "preparation failed"
    txt = log.read_text(errors="ignore")
    if "Structure parsing error" in txt:
        return "PDBQT parse error"
    if not txt.strip():
        return "no output"
    return "no scorable pose"

d["outcome"] = d.complex_id.map(outcome)
d["timed_out"] = d.outcome == "no output"
print(d.outcome.value_counts().to_string())
d[["complex_id", "torsdof", "outcome"]].to_csv(
    R / "runs/track_b/fresh/unidock_outcomes.csv", index=False)

rt = pd.read_csv(R / "runs/track_b/fresh/runtime_per_complex.csv")
rt = rt[rt.program == "Uni-Dock"][["complex_id", "seconds"]]
d = d.merge(rt, on="complex_id", how="left")

rows = []
for c in CAPS:
    m = d.torsdof <= c
    s = d[m & d.top1_dockq.notna()]
    rows.append({"cap": c, "accepted": int(m.sum()), "scored": len(s),
                 "timed_out": int(d[m].timed_out.sum()),
                 "mean_top1_dockq": round(float(s.top1_dockq.mean()), 3) if len(s) else np.nan,
                 "acceptable_0p23": round(float((s.top1_dockq >= 0.23).mean()), 3) if len(s) else np.nan,
                 "mean_best_dockq": round(float(s.best_dockq.mean()), 3) if len(s) else np.nan,
                 "median_seconds": round(float(s.seconds.median()), 0) if s.seconds.notna().any() else np.nan})

for a, b in zip(CAPS, CAPS[1:]):
    s = d[(d.torsdof > a) & (d.torsdof <= b) & d.top1_dockq.notna()]
    _band = (d.torsdof > a) & (d.torsdof <= b)
    rows.append({"cap": f"{a}<t<={b} (added)", "accepted": int(_band.sum()),
                 "timed_out": int(d[_band].timed_out.sum()),
                 "scored": len(s),
                 "mean_top1_dockq": round(float(s.top1_dockq.mean()), 3) if len(s) else np.nan,
                 "acceptable_0p23": round(float((s.top1_dockq >= 0.23).mean()), 3) if len(s) else np.nan,
                 "mean_best_dockq": round(float(s.best_dockq.mean()), 3) if len(s) else np.nan,
                 "median_seconds": round(float(s.seconds.median()), 0) if s.seconds.notna().any() else np.nan})

t = pd.DataFrame(rows)
t.to_csv(R / "runs/track_b/fresh/unidock_torsion_sensitivity.csv", index=False)
print(t.to_string(index=False))

nested = t[t.cap.isin(CAPS)]
if nested.mean_top1_dockq.notna().all():
    json.dump({"caps": {str(r.cap): {"accepted": int(r.accepted), "scored": int(r.scored),
                                     "mean_top1_dockq": float(r.mean_top1_dockq),
                                     "acceptable_0p23": float(r.acceptable_0p23),
                                     "median_seconds": None if pd.isna(r.median_seconds) else float(r.median_seconds)}
                        for r in nested.itertuples()}},
              open(R / "runs/track_b/fresh/unidock_torsion_sensitivity.json", "w"), indent=1)
    print(f"\ncap 15 -> 32: accepted {int(nested.accepted.iloc[0])} -> {int(nested.accepted.iloc[-1])}, "
          f"mean top-1 DockQ {nested.mean_top1_dockq.iloc[0]:.3f} -> {nested.mean_top1_dockq.iloc[-1]:.3f}")
