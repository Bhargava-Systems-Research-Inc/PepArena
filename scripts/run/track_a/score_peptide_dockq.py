import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import gemmi
import pandas as pd

R = Path("/12TBDrive1/mega_pep_bench")
DOCKQ = str(Path.home() / "anaconda3/envs/dockq/bin/DockQ")
MEM_CAP = "12G"

MAX_RESIDUES = 6000
DOCKQ_ATTEMPTS = 4

def guarded(cmd):
    return list(cmd)
NATIVE = R / "data/curated/structures/native"

def two_chain(path, pep_len, out_pdb, shift=None):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    model = st[0]
    water = {"HOH", "DOD", "WAT"}
    chains = []
    for ch in model:
        res = [r for r in ch if r.name not in water and len(r) > 1]
        if len(res) >= 2:
            chains.append((ch.name, res))
    if len(chains) < 2:
        return False
    pep_name, _ = min(chains, key=lambda x: (abs(len(x[1]) - pep_len), len(x[1])))

    new = gemmi.Structure()
    new.add_model(gemmi.Model("1"))
    pep = gemmi.Chain("P")
    rec = gemmi.Chain("R")
    n = 1
    chains = sorted(chains, key=lambda x: (x[0] != pep_name, x[0]))
    for name, aa in chains:
        for res in aa:
            r = gemmi.Residue()
            r.name = res.name
            r.seqid = gemmi.SeqId(n, " ")
            n += 1
            for at in res:
                if at.element == gemmi.Element("H"):
                    continue
                a = gemmi.Atom()
                a.name, a.element, a.pos = at.name, at.element, at.pos
                a.b_iso, a.occ = at.b_iso, 1.0
                if shift is not None and name == pep_name:
                    a.pos = gemmi.Position(at.pos.x + shift, at.pos.y, at.pos.z)
                r.add_atom(a)
            if len(r):
                (pep if name == pep_name else rec).add_residue(r)
        if name == pep_name:
            n = 1
    new[0].add_chain(pep)
    new[0].add_chain(rec)
    new.setup_entities()
    Path(out_pdb).write_text(new.make_pdb_string())
    return True

def too_large(path):
    try:
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        return sum(len(ch) for ch in st[0]) > MAX_RESIDUES
    except Exception:
        return False

def score(model_path, native_path, pep_len, shift=None):
    score.retries = getattr(score, 'retries', 0)
    if too_large(native_path) or too_large(model_path):
        return None
    with tempfile.TemporaryDirectory() as td:
        m, nat = f"{td}/m.pdb", f"{td}/n.pdb"
        if not two_chain(model_path, pep_len, m, shift=shift):
            return None
        if not two_chain(native_path, pep_len, nat):
            return None
        out = f"{td}/o.json"

        cmd = guarded([DOCKQ, m, nat, "--capri_peptide", "--allowed_mismatches", "10",
                       "--json", out, "--short"])
        d, last = None, None
        for attempt in range(DOCKQ_ATTEMPTS):
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
                d = json.loads(Path(out).read_text())
                if attempt:
                    score.retries += 1
                break
            except Exception as e:
                last = e
                err = (r.stderr or "") if "r" in dir() and r is not None else ""
                if "partially initialized module" not in err and "circular import" not in err:
                    break
                time.sleep(0.5 * (attempt + 1))
        if d is None:
            e = last or RuntimeError("DockQ produced no JSON")
            err = (r.stderr or "")[-400:] if "r" in dir() and r is not None else ""
            score.last_error = (f"{type(e).__name__}: {e}".strip()[:200]
                                + (f" | dockq stderr: {err.strip()[-200:]}" if err.strip() else "")
                                + f" | after {DOCKQ_ATTEMPTS} attempts")
            return None
    best = d.get("best_result") or {}
    if len(best) == 0:
        return {"DockQ": 0.0, "iRMSD": None, "LRMSD": None, "fnat": 0.0, "fnonnat": None,
                "F1": None, "clashes": None, "n_interfaces": 0}
    if len(best) != 1:
        return None
    iface = next(iter(best.values()))
    return {"DockQ": iface.get("DockQ"), "iRMSD": iface.get("iRMSD"), "LRMSD": iface.get("LRMSD"),
            "fnat": iface.get("fnat"), "fnonnat": iface.get("fnonnat"), "F1": iface.get("F1"),
            "clashes": iface.get("clashes"), "n_interfaces": len(best)}

def score_fallback(model_path, native_path, pep_len):
    with tempfile.TemporaryDirectory() as td:
        out = f"{td}/o.json"

        try:
            subprocess.run(guarded([DOCKQ, str(model_path), str(native_path), "--capri_peptide",
                                    "--allowed_mismatches", "10", "--json", out, "--short"]),
                           capture_output=True, text=True, timeout=900)
            d = json.loads(Path(out).read_text())
        except Exception:
            return None
    best = d.get("best_result") or {}
    if not best:
        return None
    cand = []
    for k, v in best.items():
        l1, l2 = v.get("len1"), v.get("len2")
        if l1 is None or l2 is None:
            continue
        cand.append((abs(min(l1, l2) - pep_len), -max(l1, l2), k, v))
    if not cand:
        return None
    cand.sort()
    gap, _, key, iface = cand[0]
    if gap > max(5, 0.5 * pep_len):
        return None
    return {"DockQ": iface.get("DockQ"), "iRMSD": iface.get("iRMSD"), "LRMSD": iface.get("LRMSD"),
            "fnat": iface.get("fnat"), "fnonnat": iface.get("fnonnat"), "F1": iface.get("F1"),
            "clashes": iface.get("clashes"), "n_interfaces": len(best), "route": "per_interface"}

def selftest():
    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    ok = True
    for cid in ["peppcbench__8tmz_A", "peppcbench__8jzw_E"]:
        nat = NATIVE / f"{cid}.cif"
        if not nat.exists():
            continue
        plen = int(sp.peptide_len.get(cid, 10))
        same = score(nat, nat, plen)
        moved = score(nat, nat, plen, shift=40.0)
        print(f"  {cid}: native-vs-native DockQ = {same['DockQ']:.3f} (expect 1.0); "
              f"peptide displaced 40 A = {moved['DockQ']:.3f} (expect ~0)")
        ok &= (same["DockQ"] > 0.99) and (moved["DockQ"] < 0.1)
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=None)
    ap.add_argument("--tier", default="canonical", choices=["canonical", "ncaa"])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    sys.path.insert(0, str(R / "scripts/run/track_a"))
    import importlib.util
    spec = importlib.util.spec_from_file_location("sd", R / "scripts/run/track_a/score_dockq.py")
    sd = importlib.util.module_from_spec(spec); spec.loader.exec_module(sd)
    sub = "predictions" if a.tier == "canonical" else "predictions_ncaa"
    sd.PRED_BASE = R / "runs/track_a" / sub
    sdir = R / "runs/track_a" / ("scores" if a.tier == "canonical" else "scores_ncaa")
    out_dir = R / "runs/track_a" / ("scores_pep" if a.tier == "canonical" else "scores_pep_ncaa")
    out_dir.mkdir(parents=True, exist_ok=True)

    sp = pd.read_parquet(R / "data/curated/structure_all.parquet").set_index("complex_id")
    if a.tier == "canonical":
        ids = sorted(set(pd.read_csv(R / "runs/track_a/esmfold2_inputs_canonical.csv")
                         .iloc[:, 0].astype(str)))
    else:
        ids = sorted(set(pd.read_csv(R / "ncaa_ccd/manifest.csv").complex_id.astype(str)))
    if a.limit:
        ids = ids[:a.limit]
    _oldf = sdir / f"{a.model}.parquet"
    old = (pd.read_parquet(_oldf) if _oldf.exists()
           else pd.DataFrame(columns=["complex_id", "DockQ"]))
    rows, unscorable = [], []
    for k, cid in enumerate(ids, 1):
        pred = sd.find_pred(a.model, cid)
        nat = NATIVE / f"{cid}.cif"
        if not pred or not nat.exists():
            unscorable.append({"complex_id": cid,
                               "reason": "no prediction file" if not pred else "no native"})
            continue
        plen = int(sp.peptide_len.get(cid, 10))
        s = score(pred, nat, plen)
        if s is None:
            s = score_fallback(pred, nat, plen)
        if s is None:
            why = getattr(score, "last_error", "") or "no interface after both routes"
            unscorable.append({"complex_id": cid, "reason": why})
            score.last_error = ""
            continue
        s.setdefault("route", "merged")
        rows.append({"complex_id": cid, **s})
        if k % 25 == 0:
            print(f"  {a.model}: {k}/{len(ids)}", flush=True)
    d = pd.DataFrame(rows)
    d.to_parquet(out_dir / f"{a.model}.parquet", index=False)
    unsc_path = out_dir / f"{a.model}_unscorable.csv"
    if unscorable:
        pd.DataFrame(unscorable).to_csv(unsc_path, index=False)
        print(f"  {len(unscorable)} complexes unscorable -> {unsc_path.name}")
    elif unsc_path.exists():
        unsc_path.unlink()
        print(f"  no unscorable complexes; removed stale {unsc_path.name}")
    if getattr(score, "retries", 0):
        print(f"  {score.retries} DockQ start(s) retried (transient multiprocessing import race)")
    j = d.merge(old[["complex_id", "DockQ"]], on="complex_id", suffixes=("", "_old"))
    n_new = len(d) - len(j)
    print(f"{a.model} [{a.tier}]: rescored {len(d)}; mean DockQ {d.DockQ.mean():.3f}")
    if len(j):
        print(f"  on the {len(j)} also scored before: {j.DockQ.mean():.3f} now "
              f"vs {j.DockQ_old.mean():.3f} before "
              f"(delta {j.DockQ.mean() - j.DockQ_old.mean():+.3f})")
        print(f"  rows changing by >0.05: "
              f"{(j.DockQ - j.DockQ_old).abs().gt(0.05).sum()}/{len(j)}")
    if n_new:
        print(f"  {n_new} complex(es) scored for the first time "
              f"(mean {d[~d.complex_id.isin(j.complex_id)].DockQ.mean():.3f})")

if __name__ == "__main__":
    sys.exit(main())
