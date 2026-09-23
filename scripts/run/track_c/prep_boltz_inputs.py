import csv, os, sys, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from boltz.data.msa.mmseqs2 import run_mmseqs2

REPO = Path("/12TBDrive1/mega_pep_bench")
OUT = REPO / "runs/track_c/boltz2_affinity"
MSA, YML = OUT / "msa", OUT / "yaml"
for d in (MSA, YML):
    d.mkdir(parents=True, exist_ok=True)

CHUNK = int(os.environ.get("MSA_CHUNK", 8))
PAR = int(os.environ.get("MSA_PAR", 5))

rows = list(csv.DictReader(open(REPO / "runs/track_c/affinity_targets.csv")))
uniq = {}
for r in rows:
    uniq.setdefault(r["receptor_pdb_id"], r["receptor_seq"])
todo = [(p, s) for p, s in uniq.items() if not (MSA / f"{p}.a3m").exists()]
print(f"targets={len(rows)} unique_receptors={len(uniq)} have={len(uniq)-len(todo)} todo={len(todo)}", flush=True)

chunks = [todo[i:i + CHUNK] for i in range(0, len(todo), CHUNK)]
lock = threading.Lock()
done = 0

def fetch(idx_chunk):
    global done
    idx, chunk = idx_chunk
    t0 = time.time()
    try:
        a3ms = run_mmseqs2([s for _, s in chunk], str(MSA / f"_tmp_{idx}"),
                           use_env=True, use_pairing=False)
    except Exception as e:
        return f"chunk {idx}: FAILED {type(e).__name__}: {e}"
    for (pdb, _), a3m in zip(chunk, a3ms):
        (MSA / f"{pdb}.a3m").write_text(a3m)
    with lock:
        done += len(chunk)
        return f"chunk {idx}: +{len(chunk)} ({done}/{len(todo)}) in {time.time()-t0:.0f}s"

if chunks:
    with ThreadPoolExecutor(max_workers=PAR) as ex:
        futs = [ex.submit(fetch, (i, c)) for i, c in enumerate(chunks)]
        for f in as_completed(futs):
            print(" ", f.result(), flush=True)

have = {p.stem for p in MSA.glob("*.a3m") if p.stat().st_size > 0}
print(f"MSAs on disk: {len(have)}/{len(uniq)}", flush=True)

written = skipped = 0
for r in rows:
    if r["receptor_pdb_id"] not in have:
        skipped += 1; continue
    (YML / f"{r['complex_id'].replace('/', '_')}.yaml").write_text(
        "version: 1\nsequences:\n  - protein:\n      id: R\n"
        f"      sequence: {r['receptor_seq']}\n"
        f"      msa: MSA_DIR/{r['receptor_pdb_id']}.a3m\n"
        "  - ligand:\n      id: L\n"
        f"      smiles: '{r['smiles']}'\n"
        "properties:\n  - affinity:\n      binder: L\n")
    written += 1
print(f"YAMLs written: {written}  skipped(no MSA): {skipped} -> {YML}", flush=True)
