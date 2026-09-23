from __future__ import annotations
import json, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/12TBDrive1/mega_pep_bench")
SC = ROOT / "runs/track_a/scores"
P = ROOT / "runs/track_a/predictions"
PNCAA = ROOT / "runs/track_a/predictions_ncaa"
OUT = ROOT / "runs/track_a/figures"; OUT.mkdir(exist_ok=True)

LABEL = {"protenix": "Protenix (v1.0.0)", "af3": "AlphaFold3", "boltz2": "Boltz-2",
         "chai1": "Chai-1", "highfold3": "HighFold3", "helixfold3": "HelixFold3",
         "esmfold2": "ESMFold2 (ESMC-6B)",
         "rfaa": "RoseTTAFold-AA"}
DOCKQ_MODELS = ["protenix", "af3", "boltz2", "chai1", "highfold3", "helixfold3",
                "esmfold2", "rfaa"]

BLUE, INK, MUTED, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#e6e5e2"

CAPRI = [("Incorrect", 0.00, 0.23, "#f4a8a2"), ("Acceptable", 0.23, 0.49, "#f8c98a"),
         ("Medium", 0.49, 0.80, "#a9c9ec"), ("High", 0.80, 1.01, "#a6d8a8")]

def dockq(m):
    return pd.read_parquet(SC / f"{m}.parquet").dropna(subset=["DockQ"]).set_index("complex_id")

_esm_summary = None
def _esm():
    global _esm_summary
    if _esm_summary is None:
        _esm_summary = pd.read_csv(P / "esmfold2/summary.csv").set_index("pdb_id")["iptm"]
    return _esm_summary

def iptm(m, cid):
    try:
        if m == "protenix":
            best = None
            for f in glob.glob(f"{P}/protenix/{cid}/seed_*/predictions/{cid}_summary_confidence_sample_*.json"):
                d = json.load(open(f))
                if best is None or d["ranking_score"] > best[0]: best = (d["ranking_score"], d["iptm"])
            return best[1] if best else np.nan
        if m == "boltz2":
            f = f"{P}/boltz2/boltz_results__boltz5_batch/predictions/{cid}/confidence_{cid}_model_0.json"
            return json.load(open(f))["iptm"] if Path(f).exists() else np.nan
        if m == "chai1":
            best = None
            for f in glob.glob(f"{P}/chai1/{cid}/scores.model_idx_*.npz"):
                d = np.load(f); s = float(d["aggregate_score"].ravel()[0])
                if best is None or s > best[0]: best = (s, float(d["iptm"].ravel()[0]))
            return best[1] if best else np.nan
        if m == "highfold3":
            f = f"{P}/highfold3/{cid}/{cid}_summary_confidences.json"
            return json.load(open(f))["iptm"] if Path(f).exists() else np.nan
        if m == "helixfold3":
            f = f"{P}/helixfold3/{cid}/{cid}-rank1/all_results.json"
            return json.load(open(f))["iptm"] if Path(f).exists() else np.nan
        if m == "esmfold2":
            s = _esm(); return float(s.loc[cid]) if cid in s.index else np.nan
            return json.load(open(f))["iptm"] if Path(f).exists() else np.nan
    except Exception:
        return np.nan
    return np.nan

def barfig(labels, values, title, subtitle, fname, valfmt="{:.3f}", ylabel="mean DockQ"):
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(8, 0.95 * n + 2), 5.4))
    x = np.arange(n)
    bars = ax.bar(x, values, width=0.68, color=BLUE, zorder=3)
    for xi, v in zip(x, values):
        ax.text(xi, v + max(values) * 0.012, valfmt.format(v), ha="center", va="bottom",
                fontsize=10, color=INK, fontweight="medium")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=32, ha="right", fontsize=10, color=INK)
    ax.set_ylabel(ylabel, fontsize=11, color=MUTED)
    ax.set_ylim(0, max(values) * 1.16)
    ax.set_title(title, fontsize=14, fontweight="bold", color=INK, loc="left", pad=14)
    ax.text(0, 1.015, subtitle, transform=ax.transAxes, fontsize=9.5, color=MUTED, va="bottom")
    ax.yaxis.grid(True, color=GRID, lw=1, zorder=0); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED); ax.tick_params(length=0, colors=MUTED)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", OUT / fname)

dq = {m: dockq(m) for m in DOCKQ_MODELS if (SC / f"{m}.parquet").exists()}
rows = [(LABEL[m], d.DockQ.mean(), len(d)) for m, d in dq.items()]
rows.sort(key=lambda r: r[1], reverse=True)
labels = [f"{r[0]}\n(n={r[2]})" for r in rows]
barfig(labels, [r[1] for r in rows],
       "Mean DockQ — each model on its full post-cutoff set",
       "Every complex is deposited after 2023-06-01, i.e. past every model's training cutoff. "
       "n differs by how many complexes each model predicted.",
       "fig1_dockq_per_model_own_cutoff.png")

common = sorted(set.intersection(*[set(d.index) for d in dq.values()]))
rows2 = [(LABEL[m], d.loc[common].DockQ.mean()) for m, d in dq.items()]
rows2.sort(key=lambda r: r[1], reverse=True)
lbl2 = [r[0] for r in rows2]
barfig(lbl2, [r[1] for r in rows2],
       "Mean DockQ — shared common set (apples-to-apples)",
       f"n={len(common)} complexes predicted by ALL models, every one past every training cutoff.",
       "fig2_dockq_shared_clean_common.png")

def r2_of(dqframe, ip):
    j = dqframe.join(pd.Series(ip, name="iptm")).dropna(subset=["DockQ", "iptm"])
    if len(j) >= 30 and j.iptm.nunique() > 1:
        return pearsonr(j.iptm, j.DockQ)[0] ** 2, len(j)
    return None

r2rows = []
for m, d in dq.items():
    if m == "af3":
        continue
    res = r2_of(d, {cid: iptm(m, cid) for cid in d.index})
    if res: r2rows.append((LABEL[m], *res))

af3d = pd.read_parquet(ROOT / "runs/track_a/scores_ncaa/af3.parquet").dropna(subset=["DockQ"]).set_index("complex_id")
def af3_iptm(cid):
    f = PNCAA / "af3" / cid / f"{cid}_summary_confidences.json"
    return json.load(open(f))["iptm"] if f.exists() else np.nan
res = r2_of(af3d, {cid: af3_iptm(cid) for cid in af3d.index})
if res: r2rows.append((LABEL["af3"], *res))

r2rows.sort(key=lambda r: r[1], reverse=True)
lbl3 = [f"{r[0]}\n(n={r[2]})" for r in r2rows]
barfig(lbl3, [r[1] for r in r2rows],
       "Confidence calibration — R² of ipTM vs DockQ",
       "Pearson R²; higher = ipTM better predicts DockQ. AF3 from its ncAA-tier run "
       "(webserver run has no ipTM); RoseTTAFold-AA has none.",
       "fig3_r2_iptm_dockq.png", valfmt="{:.3f}",
       ylabel="R² (ipTM vs DockQ)")

def capri_counts(d):
    q = d.DockQ.values
    return [int(((q >= lo) & (q < hi)).sum()) for _, lo, hi, _ in CAPRI]

srows = [(LABEL[m], capri_counts(d), len(d)) for m, d in dq.items()]
srows.sort(key=lambda r: (r[1][3] + r[1][2]) / r[2])
fig, ax = plt.subplots(figsize=(10, 0.62 * len(srows) + 2.2))
y = np.arange(len(srows))
for i, (_, counts, tot) in enumerate(srows):
    left = 0.0
    for k, (name, lo, hi, col) in enumerate(CAPRI):
        pct = 100 * counts[k] / tot
        ax.barh(i, pct, left=left, color=col, zorder=3, height=0.66)
        if pct >= 5:
            ax.text(left + pct / 2, i, f"{counts[k]}\n{pct:.0f}%", ha="center", va="center",
                    fontsize=8.5, color=INK, fontweight="medium", zorder=4)
        left += pct
ax.set_yticks(y); ax.set_yticklabels([f"{r[0]}  (n={r[2]})" for r in srows], fontsize=10, color=INK)
ax.set_xlim(0, 100); ax.set_xlabel("% of predictions", fontsize=11, color=MUTED, labelpad=8)
ax.set_title("DockQ quality distribution per model (CAPRI bins)", fontsize=14,
             fontweight="bold", color=INK, loc="left", pad=30)
ax.text(0, 1.015, "Each model on its full post-cutoff set (every complex past every training cutoff). "
        "Bars: n and % per bin.", transform=ax.transAxes, fontsize=9.5, color=MUTED, va="bottom")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, _, c in CAPRI]
leg = [f"{n} (DockQ {lo:.2f}–{hi:.2f})" if hi < 1 else f"{n} (DockQ ≥ {lo:.2f})"
       for n, lo, hi, c in CAPRI]
ax.legend(handles, leg, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.13),
          frameon=False, fontsize=9, handlelength=1.1, columnspacing=1.6)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(MUTED); ax.tick_params(length=0, colors=MUTED)
fig.tight_layout()
fig.savefig(OUT / "fig4_dockq_capri_bins.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote", OUT / "fig4_dockq_capri_bins.png")

print("\ncommon set n =", len(common))
