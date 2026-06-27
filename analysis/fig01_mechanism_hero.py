"""
R39 Figure 1 (hero/mechanism). Four panels:
 (a) two-process recovery model: reconsolidation step + log-time aging ramp,
     with the Canterbury inter-event window showing it samples only the slow early part.
 (b) within-coordinate trajectories (mean-reversion-immune evidence) at the 4
     rate-identifying coordinates.
 (c) depth profile of the recovery rate beta(z): concentration at 3-6 m.
 (d) demand-referenced recovery clock: rebuilt reserve vs required reserve for the
     next-event demand, vs the inter-event windows.
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy import stats

BASE = Path(__file__).resolve().parents[1]
PROC=BASE/"data/processed/nzgd_profile_parse_2026-06-08"; OUT=BASE
FIG=OUT/"figures"; FIG.mkdir(parents=True,exist_ok=True)

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.linewidth":0.8,
    "axes.spines.top":False,"axes.spines.right":False,"xtick.direction":"out","ytick.direction":"out",
    "axes.titlesize":9.5,"axes.titleweight":"bold","figure.dpi":150})
C={"sand":"#d9b35c","recon":"#2c7fb8","age":"#cc5500","clock":"#444444","req":"#b2182b"}
COORDC=["#1b9e77","#d95f02","#7570b3","#e7298a"]

# -- data --
inv=pd.read_csv(PROC/"repeat_pair_qc1n_inventory.csv")
inv["date_pre"]=pd.to_datetime(inv["date_pre"]); inv["date_post"]=pd.to_datetime(inv["date_post"])
inv["dt_months"]=(inv["date_post"]-inv["date_pre"]).dt.days/30.44
grid=pd.read_csv(PROC/"repeat_pair_qc1n_depth_grid.csv.gz")
g=grid.merge(inv[["pair_id","coord","dt_months","date_post","sounding_pre"]],on="pair_id",how="left").dropna(subset=["dqc1N","pre_ic"])
liq=g[g["pre_ic"]<=2.6]
po=(liq.groupby(["coord","pair_id","dt_months","sounding_pre"]).agg(dq=("dqc1N","median")).reset_index())
ndate=inv.groupby("coord")["date_post"].nunique(); c3=ndate[ndate>=3].index.tolist()
beta_z=pd.read_csv(OUT/"outputs/mod01_beta_depth_profile.csv")

fig=plt.figure(figsize=(7.2,6.4))
gs=fig.add_gridspec(2,2,hspace=0.42,wspace=0.32,left=0.09,right=0.97,top=0.93,bottom=0.09)

# (a) two-process mechanism
ax=fig.add_subplot(gs[0,0])
t=np.logspace(-3,2.3,400)  # months (1e-3 ~ 0.7 day)
step=20*(1-np.exp(-t/0.02))                 # fast reconsolidation step (days)
aging=8.5*np.log10(np.clip(t/0.1,1,None))   # log-time aging ramp
q=37.5+step+aging
ax.semilogx(t,q,color="k",lw=2,zorder=5)
ax.fill_between(t,37.5,q,where=(t<0.1),color=C["recon"],alpha=0.25)
ax.fill_between(t,37.5,q,where=(t>=0.1),color=C["age"],alpha=0.18)
ax.axvspan(3.6,6.3,color="0.5",alpha=0.18)
ax.text(0.01,62,"reconsolidation\n(drainage, t approx H^2/c_v)",color=C["recon"],fontsize=7,ha="left")
ax.text(8,52,"log-time\naging",color=C["age"],fontsize=7.5,ha="left")
ax.text(4.7,40,"Canterbury\ninter-event\nwindow",color="0.3",fontsize=6.6,ha="center")
ax.set_xlabel("time since liquefaction (months)"); ax.set_ylabel("critical-layer qc1Ncs")
ax.set_title("a  Two-process recovery model",loc="left")
ax.set_ylim(34,70)

# (b) within-coordinate trajectories
ax=fig.add_subplot(gs[0,1])
for i,c in enumerate(c3):
    s=po[po["coord"]==c].sort_values("dt_months")
    lab=s["sounding_pre"].iloc[0]
    sl,ic,*_=stats.linregress(s["dt_months"],s["dq"])
    recovers = sl>0 and (s["dq"].iloc[-1]>s["dq"].iloc[0])
    if recovers:
        ax.plot(s["dt_months"],s["dq"],"o-",color=COORDC[i],ms=4,lw=1.4,label=lab)
        xx=np.array([s["dt_months"].min(),s["dt_months"].max()])
        ax.plot(xx,ic+sl*xx,"-",color=COORDC[i],lw=0.8,alpha=0.5)
    else:  # non-recovering counterexample (heterogeneity)
        ax.plot(s["dt_months"],s["dq"],"o:",color="0.55",ms=4,mfc="white",lw=1.0,
                label=lab+" (no recovery)")
ax.axhline(0,color="0.7",lw=0.7)
ax.set_xlabel("months after baseline sounding"); ax.set_ylabel("within-coord Δqc1N (median)")
ax.set_title("b  Within-coord recovery (heterogeneous)",loc="left")
ax.legend(fontsize=6.2,frameon=False,loc="upper left",ncol=1)

# (c) depth profile beta(z)
ax=fig.add_subplot(gs[1,0])
bz=beta_z.dropna(subset=["z_mid"])
ax.plot(bz["beta_median"],bz["z_mid"],"o-",color=C["recon"],lw=1.5,ms=4)
ax.axvline(0,color="0.7",lw=0.7)
ax.fill_betweenx([2,6],-3,12,color=C["sand"],alpha=0.18)
ax.text(8.5,4,"liquefiable /\nreconsolidating\nzone",color="0.35",fontsize=6.6,va="center")
ax.set_ylim(10,1); ax.set_xlim(-3,11)
ax.set_xlabel("recovery rate β(z)  (qc1N/month)"); ax.set_ylabel("depth (m)")
ax.set_title("c  Recovery concentrates at mid-depth",loc="left")

# (d) demand-referenced clock
ax=fig.add_subplot(gs[1,1])
tt=np.linspace(0,12,100); q0=37.5; beta=1.95
ax.plot(tt,q0+beta*tt,color="k",lw=2,label="rebuilt reserve (β=1.95/mo)")
ax.fill_between(tt,q0+(beta-1.6)*tt,q0+(beta+1.6)*tt,color="0.6",alpha=0.2)
ax.axhline(153.8,color=C["req"],lw=1.5,ls="-")
ax.text(0.3,150,"required for FS=1 at Feb-2011 demand (~154)",color=C["req"],fontsize=6.6,va="top")
for x,lab in [(3.65,"Feb→Jun"),(6.34,"Jun→Dec")]:
    ax.axvline(x,color="0.5",lw=0.7,ls=":")
ax.annotate("",xy=(11,153.8),xytext=(11,q0+beta*11),
    arrowprops=dict(arrowstyle="<->",color=C["req"],lw=1))
ax.text(11.2,100,"deficit\n~118",color=C["req"],fontsize=7,va="center")
ax.set_xlim(0,13.5); ax.set_ylim(30,170)
ax.set_xlabel("months after event"); ax.set_ylabel("critical-layer qc1Ncs")
ax.set_title("d  Demand-referenced recovery clock",loc="left")
ax.legend(fontsize=6.4,frameon=False,loc="lower right")

fig.savefig(FIG/"Fig1_mechanism_hero.png",dpi=320,bbox_inches="tight")
fig.savefig(FIG/"Fig1_mechanism_hero.pdf",bbox_inches="tight")
print("[saved]",FIG/"Fig1_mechanism_hero.png")
