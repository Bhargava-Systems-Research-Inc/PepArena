from pathlib import Path

import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
MODELS = {"protenix": "Protenix", "boltz2": "Boltz-2", "chai1": "Chai-1"}
ROWS = [
    ("pTM", "ptm", "all"),
    ("blended score (shipped)", "score", "all"),
    ("global ipTM", "iptm", "all"),
    ("peptide–receptor ipTM", "iptm_pep", "all"),
    ("interface PDE", "ipde", "boltz2"),
    ("ipSAE (d0res, default)", "ipSAE", "ipsae"),
    ("ipSAE (d0chn)", "ipSAE_d0chn", "ipsae"),
    ("ipTM (ipsae chain-pair)", "ipTM_af", "ipsae"),
    ("LIS", "LIS", "ipsae"),
    ("pDockQ", "pDockQ", "ipsae"),
    ("pDockQ2", "pDockQ2", "ipsae"),
]

def main():
    cc = pd.read_csv(R / "runs/track_a/confidence_components.csv")
    common = sorted(set.intersection(*[set(g.complex_id) for _, g in cc.groupby("model")]))
    cc = cc[cc.complex_id.isin(common)]
    ip = pd.read_csv(R / "runs/track_a/ipsae_boltz2.csv")
    ip = ip[ip.complex_id.isin(common)]

    def rho(d, c):
        return round(float(d[c].corr(d.DockQ, method="spearman")), 3)

    out = []
    for label, col, where in ROWS:
        r = {"metric": label, "n": len(common)}
        if where == "all":
            for k, lab in MODELS.items():
                r[lab] = rho(cc[cc.model == k], col)
            r["n_distinct_boltz2"] = int(cc[cc.model == "boltz2"][col].nunique())
        elif where == "boltz2":
            r["Boltz-2"] = rho(cc[cc.model == "boltz2"], col)
            r["n_distinct_boltz2"] = int(cc[cc.model == "boltz2"][col].nunique())
        else:
            r["Boltz-2"] = rho(ip, col)
            r["n_distinct_boltz2"] = int(ip[col].nunique())
        out.append(r)
    S = pd.DataFrame(out)[["metric", "n", "Protenix", "Boltz-2", "Chai-1", "n_distinct_boltz2"]]
    S.to_csv(R / "runs/track_a/confidence_table.csv", index=False)
    print(S.to_string(index=False, na_rep="—"))

if __name__ == "__main__":
    main()
