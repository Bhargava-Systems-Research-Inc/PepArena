import json, pathlib, sys, time, urllib.request
import pandas as pd

R = pathlib.Path("/12TBDrive1/mega_pep_bench")
GQL = "https://data.rcsb.org/graphql"
Q = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
  rcsb_accession_info { deposit_date initial_release_date } } }"""

s = pd.read_parquet(R / "data/curated/structure_all.parquet")
ids = sorted({str(x).upper() for x in s.receptor_pdb_id.dropna() if len(str(x)) == 4})
print(f"{len(ids)} unique PDB ids")

rows, B = [], 300
for i in range(0, len(ids), B):
    chunk = ids[i:i + B]
    body = json.dumps({"query": Q, "variables": {"ids": chunk}}).encode()
    req = urllib.request.Request(GQL, data=body, headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            d = json.load(urllib.request.urlopen(req, timeout=90))
            break
        except Exception as e:
            if attempt == 3:
                print("FAILED chunk", i, e); d = {"data": {"entries": []}}
            time.sleep(2 * (attempt + 1))
    for e in (d.get("data", {}).get("entries") or []):
        if not e:
            continue
        a = e.get("rcsb_accession_info") or {}
        rows.append({"pdb_id": e["rcsb_id"].upper(),
                     "deposit_date": (a.get("deposit_date") or "")[:10] or None,
                     "release_date": (a.get("initial_release_date") or "")[:10] or None})
    print(f"  {min(i+B,len(ids))}/{len(ids)}", end="\r", flush=True)

df = pd.DataFrame(rows).drop_duplicates("pdb_id")
df.to_parquet(R / "data/curated/pdb_dates.parquet", index=False)
print(f"\nresolved {len(df)} / {len(ids)}")
miss = set(ids) - set(df.pdb_id)
if miss:
    print(f"unresolved ({len(miss)}): {sorted(miss)[:10]}")
d = df.dropna(subset=["deposit_date", "release_date"])
lag = (pd.to_datetime(d.release_date) - pd.to_datetime(d.deposit_date)).dt.days
print(f"release-minus-deposit lag: median {lag.median():.0f}d, 90th pct {lag.quantile(.9):.0f}d, max {lag.max():.0f}d")
