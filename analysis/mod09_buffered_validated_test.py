"""
R39 Module 9 -re-do the field-scale test with a BUFFERED manifestation join
(the 0 m point-in-polygon was the bug: polygons have gaps, caught only 30% of
liquefied sites). VALIDATE the susceptibility proxy before trusting the test.
"""
import numpy as np, pandas as pd, json
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score

P = Path(__file__).resolve().parents[1]
R39=P; OUT=R39/"outputs"
rep=[]; w=lambda s="":(rep.append(s),print(s))

# site table from mod07 (qc1Ncs) + mod08 (LSN)
d7=pd.read_csv(OUT/"mod07_site_table.csv")
d8=pd.read_csv(OUT/"mod08_site_table.csv")[["oid","LSN_S","LSN_F"]]
d=d7.merge(d8,on="oid",how="left")

def load(name):
    gj=json.load(open(R39/f"data/ecan_mapped_liquefaction_{name}.geojson"))
    pl=[shape(f["geometry"]) for f in gj["features"] if f.get("geometry")]
    return STRtree(pl),pl
tS,pS=load("sept2010"); tF,pF=load("feb2011")
def hit(tree,pl,E,N,buf):
    pt=Point(E,N); qg=pt.buffer(buf)
    for i in tree.query(qg):
        if pl[i].distance(pt)<=buf: return 1
    return 0
BUF=20
d["mS"]=[hit(tS,pS,E,N,BUF) for E,N in zip(d.E,d.N)]
d["mF"]=[hit(tF,pF,E,N,BUF) for E,N in zip(d.E,d.N)]
w(f"# R39 Module 9 -buffered ({BUF} m) manifestation join + validated test\n")
w(f"sites: {len(d)} | Sept manifest {int(d.mS.sum())} ({d.mS.mean():.2f}) | Feb manifest {int(d.mF.sum())} ({d.mF.mean():.2f})")

# -- VALIDATE susceptibility proxies now discriminate --
w("\n## VALIDATION (susceptibility must predict manifestation)")
for v,lab in [("qc1Ncs","qc1Ncs (lower=weaker)"),("LSN_F","LSN_Feb (higher=worse)")]:
    s=d.dropna(subset=[v])
    # qc1Ncs: lower -> more liq, so AUC on -qc1Ncs
    sign=-1 if v=="qc1Ncs" else 1
    try:
        auc=roc_auc_score(s["mF"], sign*s[v]); w(f"  {lab}: AUC={auc:.3f}")
    except Exception as e: w(f"  {lab}: {e}")
d["eband"]=pd.qcut(d["E"],4,labels=["W","CW","CE","E"])
w("  Feb manifest by Easting band (E floodplain should be highest):")
w(d.groupby("eband",observed=True)["mF"].agg(["mean","size"]).to_string())

valid = roc_auc_score(d["mF"], -d["qc1Ncs"])>0.58
w(f"\n  -> susceptibility {'VALIDATES (proceed)' if valid else 'STILL INVALID (do not trust test)'}")

# -- density-paradox test (controlling validated susceptibility) --
def run(sub,ctrl,label):
    sub=sub.dropna(subset=[ctrl]).copy()
    if sub["mF"].nunique()<2 or sub["mS"].sum()<5: w(f"\n[{label}] insufficient"); return
    sub["lnF"]=np.log(np.clip(sub["pgaF"],3,None))
    f="mF ~ "+ctrl+" + lnF + mS" if ctrl=="qc1Ncs" else "mF ~ "+ctrl+" + mS"
    m=smf.logit(f,data=sub).fit(disp=0)
    w(f"\n[{label}] n={len(sub)}  P(Feb|Sept)={sub[sub.mS==1]['mF'].mean():.2f} P(Feb|noSept)={sub[sub.mS==0]['mF'].mean():.2f}")
    for t in [x for x in [ctrl,"lnF","mS"] if x in m.params.index]:
        ci=m.conf_int().loc[t]
        w(f"    {t:7s} OR={np.exp(m.params[t]):.2f} (95%CI {np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f}) p={m.pvalues[t]:.3g}")
w("\n## DENSITY-PARADOX TEST (does prior liquefaction add risk beyond susceptibility-)")
for ctrl in ["qc1Ncs","LSN_F"]:
    run(d,ctrl,f"ALL, control={ctrl}")
    run(d[d.year<=2010],ctrl,f"pre-Feb(2010), control={ctrl}")
d.to_csv(OUT/"mod09_site_table.csv",index=False)
(OUT/"mod09_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod09_SUMMARY.md'}")
