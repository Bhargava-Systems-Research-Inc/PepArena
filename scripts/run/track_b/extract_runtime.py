import csv, datetime as dt, pathlib, re

F = pathlib.Path("/12TBDrive1/mega_pep_bench/runs/track_b/fresh")
rows = []

for log in sorted(F.glob("haddock3/*/run/log")):
    cid = log.parent.parent.name
    m = re.search(r"This HADDOCK3 run took:\s*(?:(\d+) hours?[, ]+)?(?:(\d+) minutes? and )?(\d+) seconds",
                  log.read_text(errors="ignore"))
    if m:
        h, mi, s = (int(x) if x else 0 for x in m.groups())
        rows.append(("HADDOCK3", cid, h*3600 + mi*60 + s, "exact"))

def _logs(stem):
    got = sorted(q for q in F.glob(f"{stem}*.log") if not q.name.startswith("score_"))
    if not got:
        raise SystemExit(f"no {stem} log found in {F}; copy the worker logs in before extracting")
    return "\n".join(q.read_text(errors="ignore") for q in got)

txt = _logs("adcp")
events = sorted(
    [(m.start(), "hdr", m.group(1)) for m in re.finditer(r"=== ADCP (\S+) ", txt)] +
    [(m.start(), "dur", m.group(1)) for m in re.finditer(r"Docking performed in ([\d.]+) seconds", txt)]
)
cur = None
for _, kind, val in events:
    if kind == "hdr":
        cur = val
    elif cur is not None:
        rows.append(("ADCP", cur, float(val), "exact"))
        cur = None

uh = []
for _q in sorted(q for q in F.glob("unidock*.log") if not q.name.startswith("score_")):
    _seq = re.findall(r"=== Uni-Dock (\S+) .*?(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) ===",
                      _q.read_text(errors="ignore"))
    uh += list(zip(_seq, _seq[1:]))
for (cid, t0), (_, t1) in uh:
    a = dt.datetime.strptime(t0, "%Y-%m-%d %H:%M:%S")
    b = dt.datetime.strptime(t1, "%Y-%m-%d %H:%M:%S")
    sec = (b - a).total_seconds()
    if 0 < sec < 7200:
        rows.append(("Uni-Dock", cid, sec, "includes prep"))

_rdp = F / "rapidock.log"
rd = (_rdp if _rdp.exists() else F / "rapidock_driver.log").read_text(errors="ignore")
ts = re.findall(r"\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\] stage 4|\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\]\s+inference exit", rd)
flat = [x for pair in ts for x in pair if x]
if len(flat) >= 2:
    a = dt.datetime.strptime(flat[0], "%Y-%m-%d %H:%M:%S")
    b = dt.datetime.strptime(flat[1], "%Y-%m-%d %H:%M:%S")
    n = int(re.search(r"complexes with output: (\d+)", rd).group(1))
    rows.append(("RAPiDock", f"__batch_mean_of_{n}", (b - a).total_seconds() / n, "batch mean"))

out = F / "runtime_per_complex.csv"
with out.open("w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["program", "complex_id", "seconds", "quality"]); w.writerows(rows)

import collections
c = collections.Counter(r[0] for r in rows)
print(f"wrote {out}  ({len(rows)} rows)")
for k in ["HADDOCK3", "ADCP", "RAPiDock", "Uni-Dock"]:
    v = [r[2] for r in rows if r[0] == k]
    if v:
        v.sort()
        print(f"  {k:<9} n={c[k]:<4} median={v[len(v)//2]:8.1f}s  min={v[0]:.1f}  max={v[-1]:.1f}")
