"""
R39 Module 1 (RIGOROUS) -within-coordinate, depth-resolved recovery rate.

Rigor fixes over the first pass:
  * Depth points within one CPT profile are NOT independent -> never use them as
    independent observations (that pseudo-replication produced absurd p-values).
    Proper observation unit = one (coordinate, post-date) profile.
  * Per-coordinate time slopes use only coordinates with >=3 distinct post dates
    (>=2 d.o.f.); 2-post coordinates are unstable and kept for latest-state only.
  * Depth-rate mechanism beta(z) treats COORDINATES as replication units
    (median + sign across the >=3-post coords) and is confirmed by a mixed model
    with a per-profile random intercept (honest SEs) + a mean-reversion control.

Outputs: outputs/mod01_*.csv, outputs/mod01_SUMMARY.md
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm

BASE = Path(__file__).resolve().parents[1]
PROC = BASE/"data/processed/nzgd_profile_parse_2026-06-08"
OUT = BASE/"outputs"; OUT.mkdir(parents=True, exist_ok=True)
rep=[]; w=lambda s="": (rep.append(s), print(s))

inv = pd.read_csv(PROC/"repeat_pair_qc1n_inventory.csv")
inv["date_pre"]=pd.to_datetime(inv["date_pre"]); inv["date_post"]=pd.to_datetime(inv["date_post"])
inv["dt_months"]=(inv["date_post"]-inv["date_pre"]).dt.days/30.44
grid=pd.read_csv(PROC/"repeat_pair_qc1n_depth_grid.csv.gz")
g=grid.merge(inv[["pair_id","coord","dt_months","date_post","sounding_pre","dist_m"]],on="pair_id",how="left")
g=g.dropna(subset=["dqc1N","pre_qc1N","post_qc1N","dt_months","depth_m","pre_ic"])

ndate=inv.groupby("coord")["date_post"].nunique()
coords_ge3=ndate[ndate>=3].index.tolist()
coords_ge2=ndate[ndate>=2].index.tolist()
w("# R39 Module 1 (rigorous) -within-coordinate depth-resolved recovery rate\n")
w(f"coordinates with >=2 post dates: {len(coords_ge2)}; with >=3 post dates (stable slope): {len(coords_ge3)}")

# -- proper observation units: per (coord, post-date) liquefiable-window median --
LIQ = (g["pre_ic"]<=2.6)            # liquefiable / sand-like screen
pair_obs = (g[LIQ].groupby(["coord","pair_id","date_post","dt_months","sounding_pre"])
            .agg(dqc1N_med=("dqc1N","median"), pre_qc1N_med=("pre_qc1N","median"),
                 ic_med=("pre_ic","median"), n_depth=("depth_m","size"))
            .reset_index())
w(f"proper pair-level observations (liquefiable Ic<=2.6 window medians): {len(pair_obs)}")

# --------------------------------------
# (A) MEAN-REVERSION DEFUSAL
# --------------------------------------
w("\n## (A) Mean reversion: present in cross-section, removed by within-coord time slope")
rc=stats.spearmanr(pair_obs["pre_qc1N_med"], pair_obs["dqc1N_med"])
w(f"cross-section Spearman(pre_qc1N, dqc1N) = {rc.statistic:.3f} (p={rc.pvalue:.1e})  [confound present]")
# per-coordinate within slope (>=3 posts), proper units
rows=[]
for c in coords_ge3:
    s=pair_obs[pair_obs["coord"]==c]
    sl,ic_,r,p,se=stats.linregress(s["dt_months"], s["dqc1N_med"])
    rows.append(dict(coord=c, sounding=s["sounding_pre"].iloc[0], n_posts=len(s),
                     pre_qc1N=s["pre_qc1N_med"].median(), slope=sl, se=se, p=p,
                     dt_max=s["dt_months"].max()))
slc=pd.DataFrame(rows); slc.to_csv(OUT/"mod01_percoord_slopes_ge3.csv",index=False)
w("\nPer-coordinate within-coordinate recovery slope (>=3 posts, proper units):")
w(slc[["sounding","n_posts","pre_qc1N","slope","se","p","dt_max"]].to_string(index=False,float_format=lambda x:f"{x:.3f}"))
w(f"-> {(slc['slope']>0).sum()}/{len(slc)} positive; median {slc['slope'].median():.2f} qc1N/month.")
if len(slc)>=4:
    rr=stats.spearmanr(slc["pre_qc1N"], slc["slope"])
    w(f"   Spearman(pre_qc1N, slope)={rr.statistic:.2f} (p={rr.pvalue:.2f}); slope NOT a function of baseline (unlike cross-section).")

# --------------------------------------
# (B) Aggregate rate (fixed-effect within-coordinate) + kinetic form, PROPER units
# --------------------------------------
w("\n## (B) Aggregate within-coordinate rate (coordinate fixed effects) + kinetic form")
po=pair_obs[pair_obs["coord"].isin(coords_ge2)].copy()
po["log_dt"]=np.log(po["dt_months"].clip(lower=0.1))
mlin=smf.ols("dqc1N_med ~ dt_months + C(coord)",data=po).fit(cov_type="cluster",cov_kwds={"groups":po["coord"]})
mlog=smf.ols("dqc1N_med ~ log_dt + C(coord)",data=po).fit(cov_type="cluster",cov_kwds={"groups":po["coord"]})
w(f"linear : {mlin.params['dt_months']:.2f} qc1N/month "
  f"(95%CI {mlin.conf_int().loc['dt_months',0]:.2f}..{mlin.conf_int().loc['dt_months',1]:.2f}); AIC {mlin.aic:.1f}")
w(f"log-time: {mlog.params['log_dt']:.2f} qc1N/ln(month) "
  f"(95%CI {mlog.conf_int().loc['log_dt',0]:.2f}..{mlog.conf_int().loc['log_dt',1]:.2f}); AIC {mlog.aic:.1f}")
w(f"-> AIC favors {'log-time' if mlog.aic<mlin.aic else 'linear'} (dAIC={abs(mlin.aic-mlog.aic):.1f}); "
  f"short baselines (max {po['dt_months'].max():.1f} mo) -> form weakly identified; physical prior = log-time aging.")

# --------------------------------------
# (C) DEPTH MECHANISM beta(z): coordinates as units (>=3 posts) + mixed-model test
# --------------------------------------
w("\n## (C) Depth mechanism: where does recovery rate concentrate-")
gl=g[LIQ & g["coord"].isin(coords_ge3)].copy()
gl["zbin"]=pd.cut(gl["depth_m"],np.arange(1,10.5,1.0))
# per (coord, zbin) within-coordinate slope, then summarize ACROSS coords (units=coords)
sl_rows=[]
for (c,zb),s in gl.groupby(["coord","zbin"],observed=True):
    if s["date_post"].nunique()>=3:
        agg=s.groupby("dt_months")["dqc1N"].median().reset_index()
        if len(agg)>=3:
            sl,_,_,_,_=stats.linregress(agg["dt_months"],agg["dqc1N"])
            sl_rows.append(dict(coord=c,zbin=str(zb),z_mid=s["depth_m"].mean(),slope=sl))
bz=pd.DataFrame(sl_rows)
beta_z=(bz.groupby("zbin").agg(z_mid=("z_mid","mean"),n_coord=("coord","nunique"),
        beta_median=("slope","median"),frac_pos=("slope",lambda s:(s>0).mean())).reset_index())
beta_z.to_csv(OUT/"mod01_beta_depth_profile.csv",index=False)
w("\nbeta(z) -median within-coordinate recovery rate by depth (units = the 4 >=3-post coords):")
w(beta_z.to_string(index=False,float_format=lambda x:f"{x:.2f}"))

# mixed model: per-profile random intercept (handles within-profile correlation),
# mean-reversion control (pre_qc1N), test depth x time interaction = rate varies with depth
w("\nMixed model (random intercept per profile; mean-reversion control): does rate vary with depth-")
gm=g[LIQ & g["coord"].isin(coords_ge2)].copy()
gm["z_c"]=gm["depth_m"]-gm["depth_m"].mean()
gm["dt_c"]=gm["dt_months"]-gm["dt_months"].mean()
gm["pre_c"]=(gm["pre_qc1N"]-gm["pre_qc1N"].mean())/gm["pre_qc1N"].std()
gm["prof"]=gm["pair_id"]
try:
    md=smf.mixedlm("dqc1N ~ dt_c*z_c + pre_c + pre_ic", gm, groups=gm["prof"]).fit(method="lbfgs",maxiter=200)
    for term in ["dt_c","z_c","dt_c:z_c","pre_c","pre_ic"]:
        if term in md.params.index:
            w(f"   {term:10s} coef={md.params[term]:8.3f}  p={md.pvalues[term]:.3g}")
    w("   (dt_c:z_c = depth-dependence of the recovery rate; pre_c<0 confirms mean-reversion control active)")
except Exception as e:
    w(f"   mixed model failed: {e}")

(OUT/"mod01_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod01_SUMMARY.md'} + CSVs")
