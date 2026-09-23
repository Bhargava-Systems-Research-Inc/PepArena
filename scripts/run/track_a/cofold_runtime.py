import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
LOGS = [Path(p) for p in ("/tmp/peparena_run54.log", "/tmp/peparena_run54b.log",
                          "/tmp/peparena_run54c.log")]
LAB = {"protenix": "Protenix (v1.0.0)", "af3": "AlphaFold3", "boltz2": "Boltz-2",
       "chai1": "Chai-1", "highfold3": "HighFold3", "helixfold3": "HelixFold3",
       "rfaa": "RoseTTAFold-AA", "esmfold2": "ESMFold2 (ESMC-6B)"}

present = [p for p in LOGS if p.exists()]
if not present:
    sys.exit(f"no step logs at {', '.join(str(p) for p in LOGS)}")
LOG_DAY = {}
chunks = []
for _p in present:
    d = datetime.fromtimestamp(_p.stat().st_mtime).date()
    for line in _p.read_text().splitlines():
        if line.startswith("=== "):
            LOG_DAY[line.strip()] = d
    chunks.append(_p.read_text())
text = "\n".join(chunks)
day = datetime.fromtimestamp(max(p.stat().st_mtime for p in present)).date()

def stamp(hhmm, on=None):
    h, m = map(int, hhmm.split(":"))
    base = on or day
    return datetime.combine(base, datetime.min.time()) + timedelta(hours=h, minutes=m)

def day_for(model, kind):
    for line, d in LOG_DAY.items():
        if line.startswith(f"=== {model} {kind}"):
            return d
    return day

pairs = re.findall(r"^=== (\w+) starting (\d\d:\d\d)$.*?^=== \1 finished (\d\d:\d\d) rc=",
                   text, re.M | re.S)
windows = {}
for k, t0, t1 in pairs:
    windows.setdefault(k, []).append((t0, t1))
starts = {k: v[0][0] for k, v in windows.items()}
ends = {k: v[0][1] for k, v in windows.items()}

rows = []
for k, wins in windows.items():
    d = R / "runs/track_a/predictions" / k
    if not d.is_dir():
        continue
    tops = [x for x in d.iterdir() if x.is_dir()]
    nested = [x / "predictions" for x in tops if (x / "predictions").is_dir()]
    pool = [y for nd in nested for y in nd.iterdir() if y.is_dir()] if nested else tops

    best = None
    for t0, t1 in wins:
        a = stamp(t0, day_for(k, 'starting'))
        b = stamp(t1, day_for(k, 'finished'))
        if b < a:
            b += timedelta(days=1)
        secs = (b - a).total_seconds()
        if secs <= 0:
            continue
        n = sum(1 for x in pool
                if a.timestamp() <= x.stat().st_mtime <= b.timestamp() + 60)
        if n >= 3 and (best is None or n > best[0]):
            best = (n, secs)
    if best is None:
        continue
    n, secs = best
    rows.append({"model": LAB.get(k, k), "window_minutes": round(secs / 60, 1),
                 "complexes_folded": n, "seconds_per_complex": round(secs / n, 1),
                 "minutes_per_complex": round(secs / n / 60, 2)})

t = pd.DataFrame(rows).sort_values("seconds_per_complex") if rows else pd.DataFrame()
out = R / "runs/track_a/runtime_per_complex.csv"
if len(t):
    t.to_csv(out, index=False)
    print(t.to_string(index=False))
    lo, hi = t.seconds_per_complex.min(), t.seconds_per_complex.max()
    json.dump({"models_timed": int(len(t)), "seconds_per_complex_min": float(lo),
               "seconds_per_complex_max": float(hi),
               "seconds_per_complex_median": float(t.seconds_per_complex.median()),
               "minutes_per_complex_min": round(float(lo) / 60, 2),
               "minutes_per_complex_max": round(float(hi) / 60, 2),
               "source": "completion run step log, one GPU, one manifest"},
              open(out.with_suffix(".json"), "w"), indent=1)
    print(f"\ncofolding on one GPU: {lo/60:.1f}-{hi/60:.1f} min per complex "
          f"({len(t)} models measured this way)")
    b = pd.read_csv(R / "runs/track_b/fresh/runtime_per_complex.csv")
    print("\ndocking programs against the slowest measured cofolder:")
    for p, v in b.groupby("program").seconds.median().sort_values().items():
        print(f"  {p:10s} {v:7.0f} s ({v/60:5.1f} min) = {v/hi:5.1f}x")
else:
    print("no completed steps with enough complexes yet")
