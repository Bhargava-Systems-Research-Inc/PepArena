from __future__ import annotations
import argparse, glob, hashlib, json, shutil, subprocess, sys
from pathlib import Path

ROOT = Path("/12TBDrive1/mega_pep_bench")
OUT = ROOT / "ncaa_ccd/_msa"
DBDIR = Path("/data/alphafast_db/mmseqs")
DBS = ["uniref90_padded", "small_bfd_padded"]
MM = str(Path.home() / "anaconda3/envs/pepgym-curate/bin/mmseqs")
TMP = Path("/12TBDrive1/tmp_mmseqs")

def sid(seq):
    return hashlib.sha1(seq.encode()).hexdigest()[:12]

def unique_receptors():
    seqmap, cmap = {}, {}
    for f in sorted(glob.glob(str(ROOT / "ncaa_ccd/protenix/*.json"))):
        cid = Path(f).stem
        d = json.load(open(f))
        recs = [s["proteinChain"]["sequence"] for s in d[0]["sequences"][:-1]]
        ids = []
        for s in recs:
            i = sid(s); seqmap[i] = s; ids.append(i)
        cmap[cid] = ids
    return seqmap, cmap

def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"FAIL: {' '.join(str(c) for c in cmd)}\n{r.stderr[-1500:]}\n")
        raise SystemExit(1)

def a3m_seqs(path):
    if not Path(path).exists():
        return []
    out, h, s = [], None, []
    for ln in open(path):
        if ln.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = ln.rstrip("\n"), []
        else:
            s.append(ln.rstrip("\n"))
    if h is not None:
        out.append((h, "".join(s)))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--search-only", action="store_true", help="GPU searches only (pause MD around this)")
    ap.add_argument("--post-only", action="store_true", help="CPU merge only (MD can run)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); TMP.mkdir(parents=True, exist_ok=True)
    seqmap, cmap = unique_receptors()
    ids = list(seqmap)[:args.limit] if args.limit else list(seqmap)
    json.dump(cmap, open(OUT / "receptor_map.json", "w"))
    work = TMP / "work"; work.mkdir(parents=True, exist_ok=True)

    todo = [i for i in ids if not (OUT / f"{i}.a3m").exists()]
    print(f"{len(seqmap)} unique receptors ({len(cmap)} complexes); {len(todo)} to MSA", flush=True)
    if not todo:
        return
    if not args.post_only:
        shutil.rmtree(work, ignore_errors=True); work.mkdir(parents=True)
        for d in TMP.glob("t_*"):
            shutil.rmtree(d, ignore_errors=True)
    qfa = work / "query.fasta"
    qfa.write_text("".join(f">{i}\n{seqmap[i]}\n" for i in todo))
    qdb = work / "qdb"
    if not args.post_only:
        run([MM, "createdb", str(qfa), str(qdb)])
        for db in DBS:
            print(f"  search {db} ...", flush=True)
            run([MM, "search", str(qdb), str(DBDIR / db), str(work / f"res_{db}"), str(TMP / f"t_{db}"),
                 "-a", "-s", "7.5", "-e", "1e-4", "--threads", "16", "--max-seqs", "5000", "--gpu", "1"])
            run([MM, "result2msa", str(qdb), str(DBDIR / db), str(work / f"res_{db}"),
                 str(work / f"msa_{db}"), "--msa-format-mode", "5"])
            od = work / f"unp_{db}"; od.mkdir(exist_ok=True)
            run([MM, "unpackdb", str(work / f"msa_{db}"), str(od), "--unpack-name-mode", "1", "--unpack-suffix", ".a3m"])
        if args.search_only:
            print("searches done (run --post-only next)", flush=True); return
    key2id = {ln.split("\t")[0]: ln.split("\t")[1] for ln in open(work / "qdb.lookup")}
    for key, i in key2id.items():
        merged, seen_q = [], False
        for db in DBS:
            rows = a3m_seqs(work / f"unp_{db}" / f"{key}.a3m")
            for j, (h, s) in enumerate(rows):
                if j == 0:
                    if not seen_q:
                        merged.append((f">{i}", seqmap[i])); seen_q = True
                    continue
                merged.append((h, s))
        if not seen_q:
            merged = [(f">{i}", seqmap[i])]
        (OUT / f"{i}.a3m").write_text("\n".join(f"{h}\n{s}" for h, s in merged) + "\n")
    todo = list(key2id.values())
    depths = [len(a3m_seqs(OUT / f"{i}.a3m")) for i in todo]
    print(f"wrote {len(todo)} a3m -> {OUT} | depth min/med/max "
          f"{min(depths)}/{sorted(depths)[len(depths)//2]}/{max(depths)}", flush=True)

if __name__ == "__main__":
    main()
