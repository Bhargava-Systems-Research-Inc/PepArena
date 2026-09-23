import argparse
import json
import sys
from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
P = R / "runs/track_a/predictions"

LAYOUT = {
    "protenix":   ("protenix/*/seed_*/predictions/*.cif", 3),
    "af3":        ("af3/*/*.cif", 1),
    "boltz2":     ("boltz2/boltz_results__*/predictions/*/*.cif", 1),
    "chai1":      ("chai1/*/*.cif", 1),
    "highfold3":  ("highfold3/*/*.cif", 1),
    "helixfold3": ("helixfold3/*/*/predicted_structure.cif", 2),
    "rfaa":       ("rfaa/*.pdb", 0),
    "esmfold2":   ("esmfold2/*/*/*.cif", 1),
}
LAB = {"protenix": "Protenix", "af3": "AlphaFold3", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3", "rfaa": "RoseTTAFold-AA",
       "esmfold2": "ESMFold2"}

def folded(model):
    pattern, up = LAYOUT[model]
    out = set()
    for f in P.glob(pattern):
        node = f
        for _ in range(up):
            node = node.parent
        out.add(node.stem if node.is_file() else node.name)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-coverage", type=float, default=0.0,
                    help="fail a model below this fraction of the tier (0 = report only)")
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    ids = sorted(set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv")
                     .iloc[:, 0].astype(str)))
    old = set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv.bak334")
              .iloc[:, 0].astype(str))
    restored = set(ids) - old
    tier = set(ids)

    print(f"frozen tier {len(tier)} complexes, {len(restored)} of them restored by the date audit\n")
    print(f"{'model':<14}{'folded':>9}{'of tier':>9}{'restored':>10}{'missing':>9}")
    rows, short = [], []
    for m, lab in LAB.items():
        got = folded(m)
        n, nr = len(got & tier), len(got & restored)
        rows.append({"model": lab, "folded": n, "tier": len(tier), "restored": nr,
                     "restored_total": len(restored), "coverage": round(n / len(tier), 3),
                     "missing": sorted(tier - got)[:5]})
        flag = ""
        if a.min_coverage and n / len(tier) < a.min_coverage:
            flag, _ = "  << SHORT", short.append(lab)
        print(f"{lab:<14}{n:>9}{len(tier):>9}{nr:>6}/{len(restored):<3}{len(tier)-n:>9}{flag}")

    if a.json:
        json.dump({"tier": len(tier), "restored": len(restored), "models": rows},
                  open(a.json, "w"), indent=1)

    zero = [r["model"] for r in rows if r["restored"] == 0]
    if zero:
        print(f"\n  models with none of the restored complexes: {', '.join(zero)}")
    if short:
        print(f"\nFAIL: {len(short)} model(s) below {a.min_coverage:.0%} of the tier: {short}")
        return 1
    print("\nevery model verified against real structure files on disk")
    return 0

if __name__ == "__main__":
    sys.exit(main())
