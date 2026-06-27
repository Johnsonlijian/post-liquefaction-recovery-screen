"""
R39 Module 4 -DECISIVE identifiability test (responds to fatal review issue #1).

The within-coordinate slope vs calendar time is contaminated: post-soundings
straddle the Feb/Jun/Dec 2011 re-liquefaction events, so the slope is net
multi-event drift, not a clean single-process recovery rate. Here we ask the
honest question: are there EVENT-FREE segments (two soundings of the same
coordinate within ONE inter-event quiescent window) from which a clean recovery
rate is identifiable, and what is it-
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
BASE = Path(__file__).resolve().parents[1]
PROC=BASE/"data/processed/nzgd_profile_parse_2026-06-08"; OUT=BASE/"outputs"
rep=[]; w=lambda s="":(rep.append(s),print(s))

EVENTS=[pd.Timestamp("2010-09-04"),pd.Timestamp("2011-02-22"),
        pd.Timestamp("2011-06-13"),pd.Timestamp("2011-12-23")]
def window_id(d):  # which inter-event window a date falls in
    return sum(d>=e for e in EVENTS)

inv=pd.read_csv(PROC/"repeat_pair_qc1n_inventory.csv")
# reconstruct per-sounding qc1N (liquefiable critical-layer median) at each (coord, date)
grid=pd.read_csv(PROC/"repeat_pair_qc1n_depth_grid.csv.gz")
# build a per-sounding table from the pair endpoints (pre and post are soundings)
recs=[]
for _,r in inv.iterrows():
    recs.append((r["coord"],r["sounding_pre"],pd.Timestamp(r["date_pre"]),r["pre_NZGD_ID"],"pre",r["pair_id"]))
    recs.append((r["coord"],r["sounding_post"],pd.Timestamp(r["date_post"]),r["post_NZGD_ID"],"post",r["pair_id"]))
snd=pd.DataFrame(recs,columns=["coord","sounding","date","nzgd","role","pair_id"]).drop_duplicates(["coord","date","nzgd"])
# liquefiable critical-layer median qc1N per profile (nzgd id) below crust/WT
g=grid.copy()
glq=g[(g["pre_ic"]<=2.6)&(g["depth_m"]>=2.0)]
pre_q=glq.groupby("pre_NZGD_ID")["pre_qc1N"].median(); post_q=glq.groupby("post_NZGD_ID")["post_qc1N"].median()
def qof(nz):
    if nz in post_q.index: return post_q[nz]
    if nz in pre_q.index: return pre_q[nz]
    return np.nan
snd["qc1N"]=snd["nzgd"].map(qof)
snd["win"]=snd["date"].map(window_id)
snd=snd.dropna(subset=["qc1N"]).drop_duplicates(["coord","date"])

w("# R39 Module 4 -event-free identifiability test\n")
w("Inter-event windows: 0=pre-Darfield,1=Darfield-Feb,2=Feb-Jun,3=Jun-Dec,4=post-Dec")
# clean segments: same coord, same window, >=2 distinct dates
clean=[]
for c,s in snd.groupby("coord"):
    for win,sw in s.groupby("win"):
        sw=sw.sort_values("date")
        if sw["date"].nunique()>=2:
            for i in range(len(sw)):
                for j in range(i+1,len(sw)):
                    dtm=(sw["date"].iloc[j]-sw["date"].iloc[i]).days/30.44
                    if dtm>0:
                        clean.append(dict(coord=c,win=int(win),
                            d0=sw["date"].iloc[i].date(),d1=sw["date"].iloc[j].date(),
                            dt_months=dtm,dq=sw["qc1N"].iloc[j]-sw["qc1N"].iloc[i],
                            rate=(sw["qc1N"].iloc[j]-sw["qc1N"].iloc[i])/dtm))
cl=pd.DataFrame(clean)
w(f"\nclean event-free segments (same coord, same inter-event window, >=2 dates): {len(cl)}")
if len(cl):
    w(f"  distinct coordinates: {cl['coord'].nunique()}")
    w(f"  segment durations (months): min {cl['dt_months'].min():.2f}, median {cl['dt_months'].median():.2f}, max {cl['dt_months'].max():.2f}")
    w(f"  segments with duration >= 1 month: {(cl['dt_months']>=1).sum()}")
    w(f"  segments with duration >= 3 months: {(cl['dt_months']>=3).sum()}")
    w("\n  by window:")
    w(cl.groupby("win").agg(n=("rate","size"),coords=("coord","nunique"),
        med_dt=("dt_months","median"),med_rate=("rate","median")).to_string())
    w("\n  all clean segments:")
    w(cl.sort_values(["win","coord"])[["coord","win","d0","d1","dt_months","dq","rate"]].to_string(index=False,float_format=lambda x:f"{x:.2f}"))
    longish=cl[cl["dt_months"]>=1]
    if len(longish):
        w(f"\n  >=1-month clean segments: median rate {longish['rate'].median():.2f} qc1N/month "
          f"(IQR {longish['rate'].quantile(.25):.2f}..{longish['rate'].quantile(.75):.2f}); "
          f"{(longish['rate']>0).mean()*100:.0f}% positive")
    cl.to_csv(OUT/"mod04_clean_event_free_segments.csv",index=False)
w("\n## VERDICT")
n1=(cl['dt_months']>=1).sum() if len(cl) else 0
w(f"Clean event-free segments >=1 month: {n1}. If small, a clean single-process")
w("recovery RATE is NOT identifiable from these event-straddling repeats; the")
w("calendar-time slope must be reported as a NET multi-event sequence trajectory,")
w("not a recovery-kinetics rate. This is the honest, review-driven conclusion.")
(OUT/"mod04_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod04_SUMMARY.md'}")
