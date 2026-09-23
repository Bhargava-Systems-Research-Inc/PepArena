import csv, json, os, shutil
from pathlib import Path

REPO = Path("/12TBDrive1/mega_pep_bench")
OUT = REPO / "runs/track_c/boltz2_affinity"
YML, SH = OUT / "yaml", OUT / "shards"

SLOTS = [
    ("yash_gpu0",   float(os.environ.get("R_YASH", 0.19)), 96, 10**9),
    ("a1_gpu0",     float(os.environ.get("R_A1", 0.05)),   16, 10**9),
    ("a1_gpu1",     float(os.environ.get("R_A1", 0.05)),   16, 10**9),
    ("spuser_gpu0", 0.13, 11, 350),
    ("spuser_gpu1", 0.15, 11, 350),
    ("spuser_gpu2", 0.11, 11, 350),
]

rows = {r["complex_id"]: r for r in csv.DictReader(open(REPO / "runs/track_c/affinity_targets.csv"))}
files = sorted(YML.glob("*.yaml"), key=lambda p: -len(rows.get(p.stem, {}).get("receptor_seq", "")))
print(f"{len(files)} yamls to place")

slots = [{"name": n, "thr": t, "vram": v, "cap": c, "files": []} for n, t, v, c in SLOTS]
total = sum(s["thr"] for s in slots)
T = len(files) / total
print(f"aggregate {total:.2f} it/s -> balanced structure pass ~{T/60:.0f} min\n")

for f in files:
    rl = len(rows.get(f.stem, {}).get("receptor_seq", ""))
    ok = [s for s in slots if rl <= s["cap"]]
    s = min(ok, key=lambda s: len(s["files"]) / s["thr"])
    s["files"].append(f)

if SH.exists():
    shutil.rmtree(SH)
manifest = {}
for s in slots:
    d = SH / s["name"]; d.mkdir(parents=True)
    for f in s["files"]:
        shutil.copy2(f, d / f.name)
    manifest[s["name"]] = len(s["files"])
    print(f"  {s['name']:12s} {len(s['files']):4d} yamls  thr={s['thr']:.2f}  "
          f"est_structure={len(s['files'])/s['thr']/60:.0f} min")
json.dump(manifest, open(SH / "manifest.json", "w"), indent=2)
print("\ntotal:", sum(manifest.values()))
