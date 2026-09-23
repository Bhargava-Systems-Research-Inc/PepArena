import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RUNS = ROOT / "runs/track_d"

def main():
    rows = []
    for mf in RUNS.glob("*/metrics.json"):
        m = json.loads(mf.read_text())
        t = m.get("test", {})
        rows.append({"method": m.get("model", "?"), "head": m.get("head", "?"),
                     "split": m.get("split", "?"), "n_train": m.get("n_train"),
                     **{k: t.get(k) for k in ("auroc", "aupr", "mcc", "f1", "accuracy", "n")}})
    lb = pd.DataFrame(rows).sort_values(["split", "auroc"], ascending=[True, False])
    lb.to_csv(RUNS / "leaderboard.csv", index=False)
    print(f"{len(lb)} runs across {lb['method'].nunique()} methods x {lb['head'].nunique()} heads "
          f"x {lb['split'].nunique()} splits -> {RUNS/'leaderboard.csv'}")
    for split in sorted(lb["split"].unique()):
        sub = lb[lb["split"] == split]
        best = sub.loc[sub.groupby("method")["auroc"].idxmax()].sort_values("auroc", ascending=False)
        print(f"\n=== best head per method — {split} ===")
        print(best[["method", "head", "auroc", "aupr", "mcc", "f1"]].to_string(index=False))

if __name__ == "__main__":
    main()
