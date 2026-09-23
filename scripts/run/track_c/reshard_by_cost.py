import csv, json, os, shutil
from pathlib import Path

SCR = Path("/tmp")
REPO = Path("/12TBDrive1/mega_pep_bench")
OUT = REPO / "runs/track_c/boltz2_affinity"
YML, SH = OUT / "yaml", OUT / "shards2"

depth = {}
for l in open(SCR / "msa_depth.txt"):
    p = l.split()
    if len(p) == 2:
        depth[p[0]] = int(p[1])
rows = {r["complex_id"]: r for r in csv.DictReader(open(SCR / "targets.csv"))}
yash_has = {p.stem for p in (OUT / "shards/yash_gpu0").glob("*.yaml")} if (OUT / "shards/yash_gpu0").exists() else set()

def cost(cid):
    r = rows[cid]
    return len(r["receptor_seq"]) * depth.get(r["receptor_pdb_id"], 0)

remaining = [c for c in rows if c not in yash_has]
print(f"yash keeps {len(yash_has)} (running, 0 OOM); re-sharding {len(remaining)}")

CAP_SPUSER = 900_000
CAP_A1 = 1_400_000
slots = [
    {"name": "spuser_gpu0", "cap": CAP_SPUSER, "thr": 0.12, "files": []},
    {"name": "spuser_gpu1", "cap": CAP_SPUSER, "thr": 0.12, "files": []},
    {"name": "spuser_gpu2", "cap": CAP_SPUSER, "thr": 0.12, "files": []},
    {"name": "a1_gpu0",     "cap": CAP_A1,     "thr": 0.05, "files": []},
    {"name": "a1_gpu1",     "cap": CAP_A1,     "thr": 0.05, "files": []},
    {"name": "yash_gpu1",   "cap": 10**12,     "thr": 0.18, "files": []},
]
bands = {"spuser": 0, "a1": 0, "yash_tail": 0}
for cid in sorted(remaining, key=lambda c: -cost(c)):
    c = cost(cid)
    ok = [s for s in slots if c <= s["cap"]]
    s = min(ok, key=lambda s: len(s["files"]) / s["thr"])
    s["files"].append(cid)
    bands["spuser" if c <= CAP_SPUSER else ("a1" if c <= CAP_A1 else "yash_tail")] += 1

print(f"capacity bands: <=0.9M {bands['spuser']}   0.9-1.4M {bands['a1']}   >1.4M {bands['yash_tail']}")
if SH.exists():
    shutil.rmtree(SH)
manifest = {}
for s in slots:
    d = SH / s["name"]; d.mkdir(parents=True)
    for cid in s["files"]:
        shutil.copy2(YML / f"{cid}.yaml", d / f"{cid}.yaml")
    manifest[s["name"]] = len(s["files"])
    mx = max((cost(c) for c in s["files"]), default=0)
    print(f"  {s['name']:12s} {len(s['files']):4d}  max_cost={mx:,}")
json.dump(manifest, open(SH / "manifest.json", "w"), indent=2)
print("total:", sum(manifest.values()))
