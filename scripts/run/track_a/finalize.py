import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

R = Path("/12TBDrive1/mega_pep_bench")
SRC = {"canonical": R / "runs/track_a/scores_pep", "ncaa": R / "runs/track_a/scores_pep_ncaa"}
COPY = R / "runs/track_a/scores"
PROV = R / "runs/track_a/FINAL_PROVENANCE.json"
MODELS = ["protenix", "af3", "boltz2", "chai1", "highfold3", "helixfold3", "rfaa", "esmfold2"]

STEPS = [
    ("track_a/track_a_stats_v2.py", ["--tier", "canonical"]),
    ("track_a/track_a_stats_v2.py", ["--tier", "ncaa"]),
    ("track_a/strat_tables_v2.py", []),
    ("track_a/ncaa_chemistry.py", []),
    ("track_a/compare_models.py", []),
    ("track_a/restored54.py", []),
    ("track_a/tier_sensitivity.py", []),
    ("track_a/build_leaderboard.py", ["--tier", "canonical"]),
    ("track_a/build_leaderboard.py", ["--tier", "ncaa"]),
    ("track_a/paired_holm.py", []),
    ("track_a/length_model.py", []),
]
TABLES = ["leaderboard_v2_canonical.csv", "leaderboard_v2_ncaa.csv", "coverage_v2_canonical.csv",
          "coverage_v2_ncaa.csv", "paired_delta_v2_canonical.csv", "paired_delta_v2_ncaa.csv",
          "paired_holm_canonical.csv", "length_bins_v2_canonical.csv",
          "linear_cyclic_v2_canonical.csv", "length_multivariable.csv",
          "common_set_composition_v2_canonical.csv", "common_set_composition_v2_ncaa.csv",
          "leaderboard.csv", "leaderboard_ncaa.csv",
          "ncaa_chemistry_counts.csv", "ncaa_chemistry_by_model.csv",
          "MODEL_COMPARISON.md", "LEADERBOARD.md", "LEADERBOARD_ncaa.md",
          "restored_54_comparison.csv", "restored_54_comparison.json",
          "tier_sensitivity.csv", "tier_sensitivity.json"]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]

def hashes():
    src = {}
    for tier, d in SRC.items():
        for m in MODELS:
            if (d / f"{m}.parquet").exists():
                src[f"{tier}/{m}"] = sha(d / f"{m}.parquet")
    for m in MODELS:
        if (COPY / f"{m}.parquet").exists():
            src[f"promoted/{m}"] = sha(COPY / f"{m}.parquet")
    tab = {t: sha(R / "runs/track_a" / t) for t in TABLES if (R / "runs/track_a" / t).exists()}
    return src, tab

def agreement():
    bad = []
    pairs = [(SRC["canonical"], COPY, "canonical"),
             (SRC["ncaa"], R / "runs/track_a/scores_ncaa", "ncaa")]
    for src, dst, tier in pairs:
        for m in MODELS:
            bad += _agree(src, dst, m, tier)
    return bad

def _agree(src, dst, m, tier):
    import pandas as pd
    bad = []
    for _ in (0,):
        a, b = src / f"{m}.parquet", dst / f"{m}.parquet"
        if not (a.exists() and b.exists()):
            continue
        x = pd.read_parquet(a).set_index("complex_id").DockQ.sort_index()
        y = pd.read_parquet(b).set_index("complex_id").DockQ.sort_index()
        if len(x) != len(y) or not x.index.equals(y.index) or not x.round(6).equals(y.round(6)):
            bad.append(f"{tier}/{m}: the promoted copy differs from the peptide rescoring "
                       f"({len(x)} vs {len(y)} complexes)")
    return bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="verify the released tables are the ones this script produced")
    a = ap.parse_args()

    if a.check:
        if not PROV.exists():
            sys.exit("FINAL_PROVENANCE.json missing: run finalize.py before building")
        want = json.load(open(PROV))
        src, tab = hashes()
        bad = [f"score file {k} changed since the tables were built"
               for k in want["scores"] if src.get(k) != want["scores"][k]]
        bad += [f"table {k} is not the one finalize.py produced"
                for k in want["tables"] if tab.get(k) != want["tables"][k]]
        missing = [k for k in want["tables"] if k not in tab]
        bad += [f"table {k} is missing" for k in missing]
        bad += agreement()
        if bad:
            print("Track A provenance FAILED:")
            for b in bad:
                print("  " + b)
            sys.exit(1)
        print(f"Track A provenance OK: {len(want['tables'])} tables, all from the scoring run "
              f"finalized {want['finalized']}")
        return

    for script, args in STEPS:
        p = R / "scripts/run" / script
        print(f"=== {script} {' '.join(args)}", flush=True)
        r = subprocess.run([sys.executable, str(p), *args], cwd=R)
        if r.returncode != 0:
            sys.exit(f"{script} failed with rc={r.returncode}; tables not finalized")

    dis = agreement()
    if dis:
        print("the promoted score copy disagrees with the rescoring it came from:")
        for b in dis:
            print("  " + b)
        sys.exit(1)

    src, tab = hashes()
    import datetime
    json.dump({"finalized": datetime.datetime.now().isoformat(timespec="seconds"),
               "scores": src, "tables": tab,
               "note": "every table listed was regenerated from the score files listed, in one "
                       "pass, by scripts/run/track_a/finalize.py"},
              open(PROV, "w"), indent=1)
    print(f"\nfinalized {len(tab)} tables from {len(src)} score files -> {PROV}")

if __name__ == "__main__":
    main()
