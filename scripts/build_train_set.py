import os, subprocess, json
import pandas as pd, gemmi

ROOT = "/12TBDrive1/mega_pep_bench"
OUT = "/12TBDrive1/mega_pep_bench/runs/track_a/_protenix_ft"
os.makedirs(OUT, exist_ok=True)
MMSEQS = "/home/yash/anaconda3/envs/pepgym-curate/bin/mmseqs"

d = pd.read_parquet(f"{ROOT}/data/curated/structure_all.parquet")
d["dep"] = pd.to_datetime(d.deposition_date, errors="coerce")
train = d[(d.dep >= "2021-09-30") & (d.dep < "2023-06-01")].copy()
evalset = d[d.dep >= "2023-06-01"].copy()
print(f"train(pre-filter)={len(train)} eval={len(evalset)}")

def rec_seq(path, chains):
    if not isinstance(path, str) or not os.path.exists(path):
        return ""
    try:
        st = gemmi.read_structure(path)
    except Exception:
        return ""
    if not len(st):
        return ""
    want = set()
    for c in str(chains).replace(",", " ").split():
        want.add(c.strip())
    seqs = []
    for ch in st[0]:
        if ch.name in want:
            s = ch.get_polymer().make_one_letter_sequence()
            if s:
                seqs.append(s.replace("-", "").replace("(", "").replace(")", ""))
    return "".join(seqs)

def write_fasta(df, fp):
    n = 0
    with open(fp, "w") as f:
        for _, r in df.iterrows():
            s = rec_seq(r.native_complex_path, r.receptor_chains)
            if len(s) >= 10:
                f.write(f">{r.complex_id}\n{s}\n"); n += 1
    return n

tf = f"{OUT}/train_receptors.fasta"; ef = f"{OUT}/eval_receptors.fasta"
nt = write_fasta(train, tf); ne = write_fasta(evalset, ef)
print(f"train receptors fasta={nt}  eval receptors fasta={ne}")

m8 = f"{OUT}/homology.m8"; tmp = f"{OUT}/mmtmp"
subprocess.run([MMSEQS, "easy-search", tf, ef, m8, tmp,
                "--format-output", "query,target,fident", "-c", "0.5", "--cov-mode", "0",
                "-s", "6", "--max-seqs", "50", "--threads", "16"],
               check=True, capture_output=True)
maxid = {}
for line in open(m8):
    p = line.split("\t")
    if len(p) >= 3:
        q, fi = p[0], float(p[2])
        maxid[q] = max(maxid.get(q, 0.0), fi)

train["maxid_eval"] = train.complex_id.map(lambda c: maxid.get(c, 0.0))
kept = train[train.maxid_eval < 0.40].copy()
print(f"train after homology filter (<40% id to eval): {len(kept)} / {len(train)}")

kept[["complex_id", "receptor_pdb_id", "peptide_chain", "peptide_seq", "peptide_len",
      "is_cyclic", "has_ncaa", "deposition_date", "maxid_eval", "native_complex_path"]].to_csv(
    f"{OUT}/train_meta.csv", index=False)
pdbs = sorted(set(kept.receptor_pdb_id.str.lower()))
open(f"{OUT}/train_pdb_list.txt", "w").write("\n".join(pdbs) + "\n")
print(f"unique train PDBs: {len(pdbs)} -> {OUT}/train_pdb_list.txt")
print("cyclic in train:", int(kept.is_cyclic.sum()), " ncaa:", int(kept.has_ncaa.sum()))
