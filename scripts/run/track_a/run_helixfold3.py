from __future__ import annotations
import argparse, os, shutil, subprocess
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
HF3 = "/data/PaddleHelix/apps/protein_folding/helixfold3"
PY = str(Path.home() / "anaconda3/envs/helixfold3/bin/python")
B = str(Path.home() / "anaconda3/envs/helixfold3/bin")
INP = ROOT / "prediction_inputs/helixfold3_hf"
MSAROOT = ROOT / "runs/track_a/_hf3_msa"
PRED = ROOT / "runs/track_a/predictions/helixfold3"
DBS = ROOT / "runs/track_a/_hf3_dbs"
CKPT = f"{HF3}/init_models/HelixFold3-params-20250714/HelixFold3-20250714.pdparams"
CCD = f"{HF3}/ccd_preprocessed_etkdg.pkl.gz"
NATIVE = os.environ.get("HF3_NATIVE") == "1"
AF3DB = "/data/af3_db"

def outdir(cid):
    return PRED / cid

def predicted(cid):
    d = outdir(cid)
    return d.exists() and any(d.glob("**/*.cif"))

def run_one(cid):
    out = outdir(cid)
    out.mkdir(parents=True, exist_ok=True)
    dummy = str(DBS / "dummy.fasta")
    if not NATIVE:
        src = MSAROOT / cid / "msas"
        if src.exists():
            shutil.copytree(src, out / "msas", dirs_exist_ok=True)
    uniref90 = f"{AF3DB}/uniref90_2022_05.fa" if NATIVE else dummy
    mgnify = f"{AF3DB}/mgy_clusters_2022_05.fa" if NATIVE else dummy
    uniprot = f"{AF3DB}/uniprot_all_2021_04.fa" if NATIVE else dummy
    red_bfd = f"{AF3DB}/bfd-first_non_consensus_sequences.fasta" if NATIVE else dummy
    cmd = [PY, "inference.py",
           "--jackhmmer_binary_path", f"{B}/jackhmmer", "--hhblits_binary_path", f"{B}/hhblits",
           "--hhsearch_binary_path", f"{B}/hhsearch", "--kalign_binary_path", f"{B}/kalign",
           "--hmmsearch_binary_path", f"{B}/hmmsearch", "--hmmbuild_binary_path", f"{B}/hmmbuild",
           "--nhmmer_binary_path", f"{B}/nhmmer", "--preset", "reduced_dbs",
           "--reduced_bfd_database_path", red_bfd, "--uniprot_database_path", uniprot,
           "--pdb_seqres_database_path", "/data/af3_db/pdb_seqres_2022_09_28.fasta",
           "--uniref90_database_path", uniref90, "--mgnify_database_path", mgnify,
           "--template_mmcif_dir", "/data/af3_db/mmcif_files",
           "--obsolete_pdbs_path", str(DBS / "obsolete.dat"),
           "--ccd_preprocessed_path", CCD, "--rfam_database_path", dummy,
           "--max_template_date", "2021-09-30",
           "--input_json", str((INP / f"{cid}.json").resolve()),
           "--output_dir", str(PRED.resolve()),
           "--model_name", "allatom_demo", "--init_model", CKPT,
           "--infer_times", "1", "--precision", "fp32"]
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="0",
               FLAGS_fraction_of_gpu_memory_to_use="0.14",
               FLAGS_allocator_strategy="auto_growth")
    subprocess.run(cmd, cwd=HF3, env=env, check=False)
    return predicted(cid)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=["canonical", "all"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6, help="parallel workers (CPU template search overlaps)")
    args = ap.parse_args()
    man = pd.read_csv(ROOT / "prediction_inputs/inputs_manifest.csv")
    if args.tier == "canonical":
        man = man[~man.has_ncaa]
    todo = [c for c in man.complex_id if (INP / f"{c}.json").exists() and not predicted(c)]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(man)} complexes; {len(todo)} to fold. workers={args.workers}", flush=True)
    PRED.mkdir(parents=True, exist_ok=True)
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import time as _t
    done, t0 = 0, _t.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_one, cid): cid for cid in todo}
        for i, f in enumerate(as_completed(futs), 1):
            cid = futs[f]
            try:
                ok = f.result()
            except Exception:
                ok = False
            done += ok
            rate = i / (_t.time() - t0 + 1e-9)
            print(f"  [{i}/{len(todo)}] {cid}: {'ok' if ok else 'FAIL'} "
                  f"({rate*60:.1f}/min, ETA {(len(todo)-i)/rate/60:.0f} min)", flush=True)
    print(f"DONE: {sum(predicted(c) for c in man.complex_id)}/{len(man)} predicted", flush=True)

if __name__ == "__main__":
    main()
