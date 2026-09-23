from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd, numpy as np

ROOT = Path(__file__).resolve().parents[3]
LABEL = {"protenix": "Protenix (v1.0.0)", "boltz2": "Boltz-2", "chai1": "Chai-1",
         "esmfold2": "ESMFold2 (ESMC-6B)",
         "af3": "AlphaFold3", "helixfold3": "HelixFold3",
         "highfold3": "HighFold3", "rfaa": "RoseTTAFold-AA"}
LBINS, LNAMES = [0, 5, 10, 15, 25, 60], ["<=5", "6-10", "11-15", "16-25", "26+"]

TIERS = {
    "canonical": dict(
        sc="runs/track_a/scores", out="runs/track_a/LEADERBOARD.md",
        csv="runs/track_a/leaderboard.csv",
        models=["protenix", "boltz2", "chai1", "esmfold2", "af3",
                "helixfold3", "highfold3", "rfaa"],
        title="canonical tier, post-2023-06-01, leakage-clean",
        subtitle="Canonical tier = 388 (375 linear + 13 cyclic).", chemistry=False),
    "ncaa": dict(
        sc="runs/track_a/scores_ncaa", out="runs/track_a/LEADERBOARD_ncaa.md",
        csv="runs/track_a/leaderboard_ncaa.csv",
        models=["protenix", "boltz2", "af3", "highfold3"],
        title="non-canonical (ncAA) tier, post-2023-06-01, leakage-clean",
        subtitle="ncAA tier = 225 ncAA peptide–protein complexes. Only the 4 CCD-capable "
                 "cofolders run here; ESMFold2/Chai-1/RoseTTAFold-AA/HelixFold3 cannot take "
                 "non-canonical residue input.", chemistry=True),
}

def stats(dq):
    dq = dq.dropna()
    return dict(n=len(dq), mean=round(dq.mean(), 3), median=round(dq.median(), 3),
                success_0p23=round((dq >= 0.23).mean(), 3),
                acceptable_0p49=round((dq >= 0.49).mean(), 3),
                high_0p80=round((dq >= 0.80).mean(), 3))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="canonical", choices=list(TIERS))
    args = ap.parse_args()
    cfg = TIERS[args.tier]
    SC = ROOT / cfg["sc"]

    loaded = {m: pd.read_parquet(SC / f"{m}.parquet").dropna(subset=["DockQ"])
              for m in cfg["models"] if (SC / f"{m}.parquet").exists()}
    if not loaded:
        print("no score parquets yet"); return
    common = set.intersection(*[set(d.complex_id) for d in loaded.values()])

    rows = []
    for m, d in loaded.items():
        base = {"model": LABEL.get(m, m)}
        s = stats(d.DockQ); base.update({f"all_{k}": v for k, v in s.items()})
        cyc = d.groupby("is_cyclic").DockQ.mean().round(3)
        base["all_linear_mean"] = cyc.get(False, np.nan)
        base["all_cyclic_mean"] = cyc.get(True, np.nan)
        dc = d[d.complex_id.isin(common)]
        base["common_n"] = len(dc)
        base["common_mean"] = round(dc.DockQ.mean(), 3)
        base["common_acceptable_0p49"] = round((dc.DockQ >= 0.49).mean(), 3)
        rows.append(base)
    lb = pd.DataFrame(rows).sort_values("common_mean", ascending=False)
    lb.to_csv(ROOT / cfg["csv"], index=False)

    L = [f"# Track A — cofolding leaderboard ({cfg['title']})", "",
         "Config: 1 seed x **5 diffusion samples**, **10 recycling / 200 sampling steps**. "
         f"Scored on the best-ranked sample, DockQ `--capri_peptide`. {cfg['subtitle']}", "",
         "## Overall (each model on all complexes it has predicted)", "",
         "| Model | n | mean DockQ | median | success@0.23 | acceptable@0.49 | high@0.80 | linear | cyclic |",
         "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for _, r in lb.iterrows():
        L.append(f"| {r.model} | {r.all_n} | {r.all_mean} | {r.all_median} | {r.all_success_0p23} "
                 f"| {r.all_acceptable_0p49} | {r.all_high_0p80} | {r.all_linear_mean} | {r.all_cyclic_mean} |")
    L += ["", f"## Common set — apples-to-apples on the {len(common)} complexes ALL models predicted", "",
          "| Model | mean DockQ | acceptable@0.49 |", "|---|--:|--:|"]
    for _, r in lb.iterrows():
        L.append(f"| {r.model} | {r.common_mean} | {r.common_acceptable_0p49} |")
    L += ["", "## Mean DockQ by peptide length", "",
          "| Model | " + " | ".join(LNAMES) + " |", "|---|" + "--:|" * len(LNAMES)]
    for m, d in loaded.items():
        d = d.assign(b=pd.cut(d.peptide_len, LBINS, labels=LNAMES))
        by = d.groupby("b", observed=True).DockQ.mean().round(2)
        L.append(f"| {LABEL.get(m, m)} | " + " | ".join(str(by.get(n, "-")) for n in LNAMES) + " |")
    if cfg["chemistry"]:
        cts = sorted({t for d in loaded.values() for t in d[d.is_cyclic].cyclization_type.dropna().unique()})
        L += ["", "## Mean DockQ by cyclization chemistry (cyclic subset)", "",
              "| Model | " + " | ".join(f"{c} (n)" for c in cts) + " |",
              "|---|" + "--:|" * len(cts)]
        for m, d in loaded.items():
            cells = []
            for c in cts:
                sub = d[(d.is_cyclic) & (d.cyclization_type == c)].DockQ
                cells.append(f"{sub.mean():.3f} ({len(sub)})" if len(sub) else "-")
            L.append(f"| {LABEL.get(m, m)} | " + " | ".join(cells) + " |")
    L += ["", f"_Models scored: {', '.join(LABEL.get(m, m) for m in loaded)}._"]
    L.append("")
    (ROOT / cfg["out"]).write_text("\n".join(L))
    print(lb.to_string(index=False))
    print(f"\n-> {cfg['csv']} + {cfg['out']} ({len(common)} common complexes)")

if __name__ == "__main__":
    main()
