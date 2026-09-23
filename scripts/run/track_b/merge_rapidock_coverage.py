import pandas as pd
from pathlib import Path
F = Path("/12TBDrive1/mega_pep_bench/runs/track_b/fresh")
cov = pd.read_csv(F / "coverage.csv")
rd = pd.read_csv(F / "rapidock_coverage.csv")
print("rapidock_coverage columns:", list(rd.columns))
m = {str(r.complex_id): ("ok" if bool(r.ok) else
                         f"refused: {r.unsupported}" if str(r.unsupported) not in ("", "nan")
                         else "refused")
     for r in rd.itertuples()}
cov = cov.drop(columns=["rapidock"], errors="ignore")
cov["rapidock"] = cov.complex_id.astype(str).map(m).fillna("refused: not prepared")
cov.to_csv(F / "coverage.csv", index=False)
print({c: int((cov[c] == "ok").sum()) for c in ("haddock3", "adcp", "unidock", "rapidock")})
none = ~((cov.haddock3 == "ok") | (cov.adcp == "ok") | (cov.unidock == "ok") | (cov.rapidock == "ok"))
print("accepted by none of four:", int(none.sum()), round(float(none.mean()), 3))
