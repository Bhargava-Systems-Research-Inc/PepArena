from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
D = R / "runs/track_a/pose_validity"
LAB = {"protenix": "Protenix", "af3": "AlphaFold3", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3", "rfaa": "RoseTTAFold-AA",
       "esmfold2": "ESMFold2"}
META = ("complex_id", "prepared", "error")
NATIVE_FLOOR = 0.95

def load(name):
    f = D / f"{name}.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d = d[d.prepared.fillna(False).astype(bool)]
    for c in d.columns:
        if c not in META:
            d[c] = d[c].map({True: True, False: False, "True": True, "False": False}).astype(bool)
    return d

def main():
    nat = load("native")
    if nat is None:
        raise SystemExit("native control missing - run pose_validity.py native first")
    checks = [c for c in nat.columns if c not in META]
    native_pass = nat[checks].mean()
    excluded = [c for c in checks if native_pass[c] < NATIVE_FLOOR]
    retained = [c for c in checks if c not in excluded]

    rows = []
    for key, lab in [("native", "native (control)")] + list(LAB.items()):
        d = nat if key == "native" else load(key)
        if d is None or d.empty:
            continue
        have = [c for c in retained if c in d.columns]
        r = {"model": lab, "n": len(d),
             "valid_pose_rate": round(float(d[have].all(axis=1).mean()), 3)}
        for c in checks:
            if c in d.columns:
                r[c] = round(float(d[c].mean()), 3)
        rows.append(r)
    S = pd.DataFrame(rows)
    S.to_csv(R / "runs/track_a/pose_validity_summary.csv", index=False)

    pd.DataFrame({"check": checks,
                  "native_pass": native_pass.round(3).values,
                  "retained": [c in retained for c in checks]}) \
      .sort_values("native_pass").to_csv(R / "runs/track_a/pose_validity_checks.csv", index=False)

    print(f"native control n={len(nat)}; {len(retained)} checks retained, "
          f"{len(excluded)} excluded (native pass < {NATIVE_FLOOR:.0%}): {excluded}\n")
    show = ["model", "n", "valid_pose_rate"] + [c for c in checks
                                                if S[c].min() < 1.0] if len(S) else []
    print(S[show].to_string(index=False))

if __name__ == "__main__":
    main()
