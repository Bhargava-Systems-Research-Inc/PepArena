#!/usr/bin/env python3
import importlib.util
import shutil
import sys
from pathlib import Path

WITHDRAWN = "cycesmfold2"

def load_scorer(repo):
    spec = importlib.util.spec_from_file_location(
        "score_dockq", Path(repo) / "scripts/run/track_a/score_dockq.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def main(repo, dest, dry=False):
    repo, dest = Path(repo), Path(dest)
    sc = load_scorer(repo)
    import pandas as pd

    total = missing = 0
    for tier, scores_dir, out_dir in [("predictions", "scores", "predictions"),
                                      ("predictions_ncaa", "scores_ncaa", "predictions_ncaa")]:
        sc.PRED_BASE = repo / "runs/track_a" / tier
        for parquet in sorted((repo / "runs/track_a" / scores_dir).glob("*.parquet")):
            model = parquet.stem.lstrip("_")
            if model.startswith(WITHDRAWN):
                continue
            n = 0
            for cid in pd.read_parquet(parquet)["complex_id"]:
                src = sc.find_pred(model, cid)
                if not src or not Path(src).exists():
                    missing += 1
                    continue
                dst = dest / out_dir / model / f"{cid}{Path(src).suffix}"
                if not dry:
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, dst)
                n += 1
            total += n
            print(f"  {tier}/{model}: {n} best-ranked structures")
    print(f"  total {total} structures ({missing} scored complexes had no locatable prediction)")
    return 0

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sys.exit(main(*args, dry="--dry-run" in sys.argv))
