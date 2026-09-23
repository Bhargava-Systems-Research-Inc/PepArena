#!/usr/bin/env python3
import os
import json
from pathlib import Path

import gemmi
import pandas as pd

REPO = Path("/12TBDrive1/mega_pep_bench")
FRESH = Path(os.environ.get("FRESH", str(REPO / "runs/track_b/fresh")))

CANON = {"GLY": "G", "ALA": "A", "VAL": "V", "ILE": "I", "LEU": "L", "MET": "M", "PHE": "F",
         "TYR": "Y", "TRP": "W", "PRO": "P", "SER": "S", "THR": "T", "ASN": "N", "GLN": "Q",
         "ASP": "D", "GLU": "E", "CYS": "C", "ARG": "R", "HIS": "H", "LYS": "K"}
NCAA = """HYP SEP TYS ALY TPO PTR DAL MLE M3L DLE DLY AIB MSE DPR MVA NLE MLY SAR ABA FME
DAR ORN CGU DPN DTY DTR 4BF DGL DCY MK8 MP8 GHP ALC BMT MLZ DVA 3FG DAS 7ID DSN AR7 MEA PHI
MAA LPD KCR PCA DGN 2MR DHI ASA MLU YCP DSG DTH OMY FP9 DPP HCS SET DBB BTK DAM IIL 3MY SLL
PFF HRG DIL DNE MED D0C""".split()
SUPPORTED = set(CANON) | set(NCAA)

def peptide_residues(cif, chain):
    st = gemmi.read_structure(str(cif)); st.setup_entities(); st.remove_waters()
    for ch in st[0]:
        if ch.name == chain:
            return [r.name.strip().upper() for r in ch]
    best = min((c for c in st[0] if len(c)), key=len, default=None)
    return [r.name.strip().upper() for r in best] if best else []

def main():
    wl = pd.read_csv(FRESH / "worklist.csv")
    rows = []
    for _, w in wl.iterrows():
        res = peptide_residues(w.native_path, w.peptide_chain)
        if not res:
            rows.append({"complex_id": w.complex_id, "n_res": 0, "ok": False,
                         "unsupported": "could not read peptide chain", "seq": ""})
            continue
        CAPS = {"ACE", "NH2", "NME", "NHE"}
        bad = sorted({r for r in res if r not in SUPPORTED})
        bad_nocap = sorted({r for r in res if r not in SUPPORTED and r not in CAPS})
        seq = "".join(CANON[r] if r in CANON else f"[{r}]" for r in res)
        rows.append({"complex_id": w.complex_id, "n_res": len(res), "ok": not bad,
                     "ok_ignoring_caps": not bad_nocap,
                     "unsupported": ",".join(bad), "seq": seq,
                     "has_ncaa": bool(w.has_ncaa), "cyclization_type": w.cyclization_type})
    d = pd.DataFrame(rows)
    d.to_csv(FRESH / "rapidock_coverage.csv", index=False)
    n = len(d); ok = int(d.ok.sum()); okc = int(d.ok_ignoring_caps.sum())
    print(f"RAPiDock can express {ok}/{n} ({100*ok/n:.0f}%) strictly, "
          f"{okc}/{n} ({100*okc/n:.0f}%) if terminal caps (ACE/NH2/NME) are stripped")
    print(d.groupby("has_ncaa").ok.agg(["sum", "size"]).to_string())
    print("\nby cyclization:")
    print(d.groupby("cyclization_type").ok.agg(["sum", "size"]).to_string())
    unsup = d.loc[~d.ok, "unsupported"].str.split(",").explode().value_counts()
    print(f"\nmost common unsupported residues:\n{unsup.head(12).to_string()}")
    json.dump({"n": n, "expressible": ok, "expressible_ignoring_caps": okc},
              open(FRESH / "rapidock_coverage.json", "w"), indent=2)

if __name__ == "__main__":
    main()
