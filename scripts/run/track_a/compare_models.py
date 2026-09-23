import pandas as pd, numpy as np, json, glob
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score
P=Path("runs/track_a/predictions"); SC=Path("runs/track_a/scores")
M=["protenix","boltz2","chai1"]; NAME={"protenix":"Protenix","boltz2":"Boltz-2","chai1":"Chai-1"}

def pconf(cid):
    b=None
    for f in glob.glob(f"{P}/protenix/{cid}/seed_*/predictions/{cid}_summary_confidence_sample_*.json"):
        d=json.load(open(f))
        if b is None or d["ranking_score"]>b["ranking_score"]: b=d
    return b and dict(conf=b["ranking_score"],iptm=b["iptm"],ptm=b["ptm"])
def cconf(cid):
    b=None
    for f in glob.glob(f"{P}/chai1/{cid}/scores.model_idx_*.npz"):
        d=np.load(f); s=float(d["aggregate_score"].ravel()[0])
        if b is None or s>b["conf"]: b=dict(conf=s,iptm=float(d["iptm"].ravel()[0]),ptm=float(d["ptm"].ravel()[0]))
    return b
def bconf(cid):
    f=f"{P}/boltz2/boltz_results__boltz5_batch/predictions/{cid}/confidence_{cid}_model_0.json"
    if not Path(f).exists(): return None
    d=json.load(open(f)); return dict(conf=d["confidence_score"],iptm=d["iptm"],ptm=d["ptm"])
EX={"protenix":pconf,"boltz2":bconf,"chai1":cconf}
dq={m:pd.read_parquet(SC/f"{m}.parquet").dropna(subset=["DockQ"]).set_index("complex_id") for m in M}
common=sorted(set.intersection(*[set(v.index) for v in dq.values()]))
D={}
for m in M:
    df=dq[m].loc[common].copy()
    cf=pd.DataFrame({cid:EX[m](cid) for cid in common if EX[m](cid)}).T
    df=df.join(cf); D[m]=df
mat=pd.DataFrame({m:D[m].DockQ for m in M})

def r(cells): return "| "+" | ".join(cells)+" |"
L=["# Track A — Protenix vs Boltz-2 vs Chai-1 (head-to-head)","",
   f"Leakage-clean canonical set, **common {len(common)} complexes** all three predicted. "
   "All WITH MSA, 1 seed x 5 diffusion samples / 10 recycle / 200 step, best-ranked sample scored "
   "(DockQ `--capri_peptide`). Confidence = each model's own ranking metric (Protenix ranking_score, "
   "Boltz-2 confidence_score, Chai-1 aggregate_score).","",
   "## Verdict","",
   "**Protenix wins on both accuracy and calibration.** Boltz-2 is the most *confident* but "
   "over-confident (highest self-score, middling accuracy, weakest calibration). Chai-1 is the most "
   "conservative and cleanest (fewest clashes, tight linear calibration) but lowest raw accuracy.",""]

L+=["## 1. Accuracy (DockQ)","", r(["Metric"]+[NAME[m] for m in M]), r(["---"]+["--:"]*3)]
def col(fn,f="{:.3f}"): return [f.format(fn(D[m])) for m in M]
for lab,fn in [("mean DockQ",lambda x:x.DockQ.mean()),("median DockQ",lambda x:x.DockQ.median()),
  ("std",lambda x:x.DockQ.std()),("IQR 25th",lambda x:x.DockQ.quantile(.25)),("IQR 75th",lambda x:x.DockQ.quantile(.75)),
  ("incorrect <0.23",lambda x:(x.DockQ<0.23).mean()),("acceptable >=0.23",lambda x:(x.DockQ>=0.23).mean()),
  ("medium >=0.49",lambda x:(x.DockQ>=0.49).mean()),("high >=0.80",lambda x:(x.DockQ>=0.80).mean())]:
    L.append(r([lab]+col(fn)))

L+=["","## 2. Interface quality (means)","", r(["Metric","Protenix","Boltz-2","Chai-1","better"]), r(["---","--:","--:","--:","---"])]
for lab,c,b in [("iRMSD (A)","iRMSD","lower"),("LRMSD (A)","LRMSD","lower"),("fnat","fnat","higher"),
  ("fnonnat","fnonnat","lower"),("F1","F1","higher"),("clashes","clashes","lower")]:
    L.append(r([lab]+["{:.3f}".format(D[m][c].mean()) for m in M]+[b]))

L+=["","## 3. Confidence calibration — does the model know when it's right?","",
    r(["Metric"]+[NAME[m] for m in M]+["better"]), r(["---","--:","--:","--:","---"])]
for lab,fn,b in [
  ("Spearman(conf, DockQ)",lambda x:spearmanr(x.conf,x.DockQ).statistic,"higher"),
  ("Pearson(conf, DockQ)",lambda x:pearsonr(x.conf,x.DockQ)[0],"higher"),
  ("Spearman(ipTM, DockQ)",lambda x:spearmanr(x.iptm,x.DockQ).statistic,"higher"),
  ("AUROC conf->success(>=0.23)",lambda x:roc_auc_score((x.DockQ>=0.23).astype(int),x.conf),"higher"),
  ("AUROC conf->high(>=0.80)",lambda x:roc_auc_score((x.DockQ>=0.80).astype(int),x.conf),"higher"),
  ("mean self-confidence",lambda x:x.conf.mean(),"—"),
  ("mean ipTM",lambda x:x.iptm.mean(),"—")]:
    L.append(r([lab]+["{:.3f}".format(fn(D[m])) for m in M]+[b]))

L+=["","## 4. Selective prediction — mean DockQ keeping only the most-confident","",
    r(["Keep"]+[NAME[m] for m in M]), r(["---","--:","--:","--:"])]
for fr in [1.0,0.75,0.5,0.25]:
    def sel(x,fr=fr): return x.sort_values("conf",ascending=False).head(int(len(x)*fr)).DockQ.mean()
    L.append(r([f"top {int(fr*100)}%"]+["{:.3f}".format(sel(D[m])) for m in M]))

L+=["","## 5. Mean DockQ by peptide length","", r(["Length"]+[NAME[m] for m in M]+["n"]), r(["---","--:","--:","--:","--:"])]
bins=[0,5,10,15,25,60]; names=["<=5","6-10","11-15","16-25","26+"]
for m in M: D[m]=D[m].assign(b=pd.cut(D[m].peptide_len,bins,labels=names))
for nm in names:
    n=int((D[M[0]].b==nm).sum())
    L.append(r([nm]+["{:.3f}".format(D[m][D[m].b==nm].DockQ.mean()) for m in M]+[str(n)]))

L+=["","## 6. Head-to-head (per-complex winner by DockQ)","", r(["","Protenix","Boltz-2","Chai-1"]), r(["---","--:","--:","--:"])]
L.append(r(["#1 (best) count"]+[str(int((mat.idxmax(axis=1)==m).sum())) for m in M]))
for a in M:
    L.append(r([f"{NAME[a]} beats ->"]+[("—" if a==bm else f"{int((mat[a]>mat[bm]).sum())} ({(mat[a]>mat[bm]).mean()*100:.0f}%)") for bm in M]))
L.append("")
Path("runs/track_a/MODEL_COMPARISON.md").write_text("\n".join(L))
print("\n".join(L))
