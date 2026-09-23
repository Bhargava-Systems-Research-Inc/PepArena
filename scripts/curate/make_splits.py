from __future__ import annotations
import hashlib, shutil, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd

MMSEQS = shutil.which("mmseqs") or str(Path(sys.executable).parent / "mmseqs")
if not Path(MMSEQS).exists():
    raise SystemExit(f"mmseqs not found (looked on PATH and at {MMSEQS}); "
                     f"activate pepgym-curate or set it on PATH")

def _bucket(key, mod=5):
    return int(hashlib.sha1(str(key).encode()).hexdigest(), 16) % mod

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "curate"))
import struct_utils

ALL = ROOT / "data/curated/structure_all.parquet"
CIF_DIR = ROOT / "data/raw/structures/cif"
OUT = ROOT / "data/splits/leakage.parquet"

CUTOFFS = {"after_af2_2018": "2018-04-30", "after_af3_2021": "2021-09-30",
           "after_2023": "2023-01-01", "after_2024": "2024-01-01"}

def receptor_seq(row):
    chains = str(row.get("receptor_chains") or "").replace(":", ",").split(",")
    chains = [c.strip() for c in chains if c.strip()]
    path = row.get("native_complex_path")
    has_path = isinstance(path, str) and path and Path(path).exists()
    try:
        if has_path:
            st = struct_utils.read_any(path)
        elif pd.notna(row.get("receptor_pdb_id")):
            cif = CIF_DIR / f"{str(row['receptor_pdb_id']).lower()}.cif.gz"
            if not cif.exists():
                return ""
            st = struct_utils.read_any(cif)
        else:
            return ""
    except Exception:
        return ""
    seqs = []
    for ch in chains:
        s, _, _ = struct_utils.chain_seq_and_ncaa(st, ch)
        if s:
            seqs.append(s)
    return "".join(seqs)

def mmseqs_clusters(id_to_seq, min_id=0.3, cov=0.8):
    items = {k: v for k, v in id_to_seq.items() if v}
    if not items:
        return {k: k for k in id_to_seq}
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        fa = td / "in.fasta"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in items.items()))
        subprocess.run([MMSEQS, "easy-cluster", str(fa), str(td / "res"), str(td / "tmp"),
                        "--min-seq-id", str(min_id), "-c", str(cov), "-v", "0"],
                       check=True, stdout=subprocess.DEVNULL)
        clu = {}
        for line in (td / "res_cluster.tsv").read_text().splitlines():
            rep, mem = line.split("\t")
            clu[mem] = rep
    for k in id_to_seq:
        clu.setdefault(k, k)
    return clu

def main():
    df = pd.read_parquet(ALL)
    print(f"extracting receptor sequences for {len(df)} complexes...")
    rec = {r.complex_id: receptor_seq(r) for _, r in df.iterrows()}
    pep = {r.complex_id: (r.peptide_seq or "") for _, r in df.iterrows()}
    n_rec = sum(1 for v in rec.values() if v)
    print(f"  got receptor seq for {n_rec}/{len(df)} (monomers/missing -> singletons)")

    rec_clu = mmseqs_clusters(rec)
    pep_clu = mmseqs_clusters(pep)

    out = pd.DataFrame({"complex_id": df.complex_id})
    out["receptor_cluster"] = out.complex_id.map(rec_clu)
    out["peptide_cluster"] = out.complex_id.map(pep_clu)
    dep = pd.to_datetime(df.set_index("complex_id").deposition_date, errors="coerce")
    for flag, cut in CUTOFFS.items():
        out[flag] = out.complex_id.map((dep > pd.Timestamp(cut)).to_dict()).fillna(False)

    out["split_random"] = out.complex_id.map(lambda c: "test" if _bucket(c) == 0 else "train")
    out["split_homology"] = out.receptor_cluster.map(
        lambda c: "test" if _bucket(c) == 0 else "train")
    pre2024 = set(out.loc[~out.after_2024, "receptor_cluster"])
    out["split_time"] = [
        "test" if (a24 and rc not in pre2024) else "train"
        for a24, rc in zip(out.after_2024, out.receptor_cluster)]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)

    print(f"\nreceptor clusters: {out.receptor_cluster.nunique()} (from {len(df)} complexes)")
    print(f"peptide clusters:  {out.peptide_cluster.nunique()}")
    print("time tiers (count after each cutoff):")
    for flag in CUTOFFS:
        print(f"  {flag:18s}: {int(out[flag].sum())}")
    print("split assignments (train/test):")
    for s in ("split_random", "split_homology", "split_time"):
        vc = out[s].value_counts()
        print(f"  {s:15s}: train {int(vc.get('train', 0))} / test {int(vc.get('test', 0))}")
    print(f"\nwrote {OUT}")

if __name__ == "__main__":
    main()
