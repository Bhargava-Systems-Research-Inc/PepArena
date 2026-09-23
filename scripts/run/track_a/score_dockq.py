from __future__ import annotations
import argparse, json, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
NATIVE = ROOT / "data/curated/structures/native"
DOCKQ = str(Path.home() / "anaconda3/envs/dockq/bin/DockQ")
MANIFEST = ROOT / "prediction_inputs/inputs_manifest.csv"

PRED_BASE = ROOT / "runs/track_a/predictions"

def find_pred(model, cid):
    base = PRED_BASE / model
    if model == "boltz2":
        hits = sorted(base.glob(f"boltz_results_*/predictions/{cid}/{cid}_model_0.cif"))
        return hits[0] if hits else None
    if model == "protenix":
        preds = sorted(base.glob(f"{cid}/seed_*/predictions/{cid}_sample_*.cif"))
        best, best_score = None, -1e18
        for p in preds:
            n = p.stem.split("_sample_")[-1]
            cj = p.parent / f"{cid}_summary_confidence_sample_{n}.json"
            score = -1e18
            if cj.exists():
                try:
                    score = float(json.loads(cj.read_text()).get("ranking_score", -1e18))
                except Exception:
                    pass
            if score > best_score:
                best, best_score = p, score
        return best
    if model == "chai1":
        import numpy as np
        best, best_score = None, -1e18
        for p in sorted((base / cid).glob("pred.model_idx_*.cif")):
            n = p.stem.split("model_idx_")[-1]
            npz = p.parent / f"scores.model_idx_{n}.npz"
            score = -1e18
            if npz.exists():
                try:
                    score = float(np.load(npz)["aggregate_score"].ravel()[0])
                except Exception:
                    pass
            if score > best_score:
                best, best_score = p, score
        return best
    if model == "rfaa":
        p = base / f"{cid}.pdb"
        return p if p.exists() else None
    if model == "helixfold3":
        p = base / cid / f"{cid}-rank1" / "predicted_structure.cif"
        return p if p.exists() else None
    if model in ("highfold3", "af3"):
        for name in (cid, cid.lower()):
            for fn in (f"{name}_model.cif", f"{name}_model_0.cif"):
                q = base / name / fn
                if q.exists():
                    return q
        return None
    for pat in (f"**/{cid}/{cid}_model_0.cif", f"**/{cid}*_model_0.cif",
                f"**/{cid}/*rank*1*.cif", f"**/{cid}*.pdb"):
        hits = sorted(base.glob(pat))
        if hits:
            return hits[0]
    return None

def cif_to_pdb(cif_path, pdb_path):
    cols, rows, in_loop = [], [], False
    for ln in open(cif_path):
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".")[1]); in_loop = True
        elif in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
            rows.append(s.split())
    idx = {c: i for i, c in enumerate(cols)}
    g = lambda r, k, d="": r[idx[k]] if k in idx and idx[k] < len(r) else d
    with open(pdb_path, "w") as f:
        for n, r in enumerate(rows, 1):
            atom = g(r, "auth_atom_id") or g(r, "label_atom_id")
            res = g(r, "auth_comp_id") or g(r, "label_comp_id")
            ch = (g(r, "auth_asym_id") or g(r, "label_asym_id") or "A")[0]
            seq = g(r, "auth_seq_id") or g(r, "label_seq_id") or "1"
            el = g(r, "type_symbol", atom[:1])
            an = atom if len(atom) >= 4 else f" {atom:<3s}"
            f.write(f"{g(r,'group_PDB','ATOM'):<6s}{n:5d} {an:<4s} {res:>3s} {ch}{int(seq):4d}    "
                    f"{float(g(r,'Cartn_x')):8.3f}{float(g(r,'Cartn_y')):8.3f}{float(g(r,'Cartn_z')):8.3f}"
                    f"{1.0:6.2f}{float(g(r,'B_iso_or_equiv','0')):6.2f}          {el:>2s}\n")
        f.write("END\n")

def _run_dockq(model_path, native_path):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
        out = tf.name
    try:
        subprocess.run([DOCKQ, str(model_path), str(native_path), "--capri_peptide",
                        "--allowed_mismatches", "10", "--json", out, "--short"],
                       capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return None
    try:
        return json.loads(Path(out).read_text())
    except Exception:
        return None

def dockq(model_path, native_path, convert_first=False):
    if convert_first and str(model_path).endswith(".cif"):
        try:
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tf:
                pdb = tf.name
            cif_to_pdb(model_path, pdb)
            model_path = pdb
        except Exception:
            pass
    d = _run_dockq(model_path, native_path)
    if d is None and str(model_path).endswith(".cif"):
        try:
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tf:
                pdb = tf.name
            cif_to_pdb(model_path, pdb)
            d = _run_dockq(pdb, native_path)
        except Exception:
            return None
    if d is None:
        return None
    best = d.get("best_result") or {}
    iface = next(iter(best.values()), {}) if isinstance(best, dict) else {}
    g = d.get("GlobalDockQ", iface.get("DockQ"))
    return {"DockQ": g, "iRMSD": iface.get("iRMSD"), "LRMSD": iface.get("LRMSD"),
            "fnat": iface.get("fnat"), "fnonnat": iface.get("fnonnat"),
            "F1": iface.get("F1"), "clashes": iface.get("clashes")}

def main():
    global PRED_BASE
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--tier", default="canonical", choices=["canonical", "ncaa"])
    args = ap.parse_args()
    ncaa = args.tier == "ncaa"
    if ncaa:
        PRED_BASE = ROOT / "runs/track_a/predictions_ncaa"
    out = ROOT / ("runs/track_a/scores_ncaa" if ncaa else "runs/track_a/scores") / f"{args.model}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if out.exists():
        done = {row["complex_id"]: row for row in pd.read_parquet(out).to_dict("records")}

    if ncaa:
        man = pd.read_csv(ROOT / "ncaa_ccd/manifest.csv")
        man["is_cyclic"] = man.cyclization_type != "linear"
    else:
        man = pd.read_csv(MANIFEST)
        man = man[~man.has_ncaa]
    rows, n_new, n_missing = list(done.values() if done else []), 0, 0
    scored = set(done)
    for _, m in man.iterrows():
        cid = m.complex_id
        if cid in scored:
            continue
        pred = find_pred(args.model, cid)
        nat = NATIVE / f"{cid}.cif"
        if not pred or not nat.exists():
            n_missing += 1; continue
        sc = dockq(pred, nat, convert_first=(args.model in ("highfold3", "af3", "helixfold3")))
        if not sc:
            n_missing += 1; continue
        row = {"complex_id": cid, "is_cyclic": bool(m.is_cyclic),
               "cyclization_type": m.cyclization_type, "peptide_len": int(m.peptide_len), **sc}
        if ncaa:
            row.update(n_ncaa=int(m.n_ncaa), n_bonds=int(m.n_bonds), tier=m.tier)
        rows.append(row)
        n_new += 1
        if n_new % 25 == 0:
            pd.DataFrame(rows).to_parquet(out, index=False)
            print(f"  scored {n_new} (+{len(scored)} cached)", flush=True)
    df = pd.DataFrame(rows)
    df.to_parquet(out, index=False)
    print(f"\n{args.model}: scored {len(df)} complexes ({n_new} new, {n_missing} missing pred) -> {out}")
    if len(df):
        d = df.dropna(subset=["DockQ"])
        print(f"  mean DockQ {d.DockQ.mean():.3f} | success@0.23 {(d.DockQ>=0.23).mean():.2f} | "
              f"acceptable@0.49 {(d.DockQ>=0.49).mean():.2f}")
        print("  by cyclic:", d.groupby("is_cyclic").DockQ.mean().round(3).to_dict())

if __name__ == "__main__":
    main()
