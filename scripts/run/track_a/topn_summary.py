from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
LAB = {"protenix": "Protenix", "boltz2": "Boltz-2", "chai1": "Chai-1",
       "highfold3": "HighFold3", "helixfold3": "HelixFold3"}

def main():
    rows = []
    for m, lab in LAB.items():
        f = R / f"runs/track_a/scores_topn/{m}.parquet"
        if not f.exists():
            continue
        d = pd.read_parquet(f)
        d = d[d.n_samples >= 2]
        if d.empty:
            continue
        rows.append({
            "model": lab, "n": len(d),
            "mean_samples": round(d.n_samples.mean(), 2),
            "top1_mean": round(d.top1.mean(), 3),
            "oracle_mean": round(d.oracle.mean(), 3),
            "gap_mean": round((d.oracle - d.top1).mean(), 3),
            "worst_mean": round(d.worst.mean(), 3),
            "spread_mean": round((d.oracle - d.worst).mean(), 3),
            "frac_any_near_zero": round(float((d.worst < 0.05).mean()), 3),
            "acc_top1": round((d.top1 >= 0.23).mean(), 3),
            "acc_oracle": round((d.oracle >= 0.23).mean(), 3),
            "med_top1": round((d.top1 >= 0.49).mean(), 3),
            "med_oracle": round((d.oracle >= 0.49).mean(), 3),
            "med_ranked_given_available": round(
                (d.top1 >= 0.49).sum() / max((d.oracle >= 0.49).sum(), 1), 3),
        })
    S = pd.DataFrame(rows).sort_values("oracle_mean", ascending=False)
    S.to_csv(R / "runs/track_a/topn_summary.csv", index=False)
    print(S.to_string(index=False))

if __name__ == "__main__":
    main()
