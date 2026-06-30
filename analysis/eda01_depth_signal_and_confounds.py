"""
R39 EDA 01 -Quantify the depth-resolved post-liquefaction recovery signal and
stress-test the main confounds BEFORE committing to a mechanistic model.

Inputs (read-only):
  data/processed/nzgd_profile_parse_2026-06-08/repeat_pair_qc1n_inventory.csv
  data/processed/nzgd_profile_parse_2026-06-08/repeat_pair_qc1n_depth_grid.csv.gz
  data/processed/nzgd_profile_parse_2026-06-08/pga_manifestation_join.csv

Questions:
  1. Sample structure: pairs, coordinates, multi-time coordinates, time gaps.
  2. Depth dependence of dqc1N (mechanism: drainage/reconsolidation -> depth).
  3. Soil-type (Ic) dependence of dqc1N (mechanism: fines control kinetics).
  4. Regression-to-the-mean confound: is dqc1N driven by low pre_qc1N-
     (Diagnostic concern: low initial state can create apparent "recovery".)
  5. Time dependence after controlling for depth/Ic/initial state.
  6. Within-coordinate (true repeat) time trajectories.
"""
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
PROC = BASE / "data/processed/nzgd_profile_parse_2026-06-08"
OUT = BASE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

inv = pd.read_csv(PROC / "repeat_pair_qc1n_inventory.csv")
grid = pd.read_csv(PROC / "repeat_pair_qc1n_depth_grid.csv.gz")
pga = pd.read_csv(PROC / "pga_manifestation_join.csv")

inv["date_pre"] = pd.to_datetime(inv["date_pre"])
inv["date_post"] = pd.to_datetime(inv["date_post"])
inv["dt_days"] = (inv["date_post"] - inv["date_pre"]).dt.days
inv["dt_months"] = inv["dt_days"] / 30.44

def line(s=""): print(s)

line("="*70)
line("1. SAMPLE STRUCTURE")
line("="*70)
line(f"pairs (inventory rows): {len(inv)}")
line(f"computable pairs: {(inv['profile_pair_status']=='computable').sum()}")
line(f"unique coordinates: {inv['coord'].nunique()}")
# multi-time coordinates = same coord, >=2 distinct post dates
mt = inv.groupby('coord')['date_post'].nunique()
line(f"coords with >=2 distinct post dates (rate-identifying): {(mt>=2).sum()}")
line(f"coords with >=3 distinct post dates: {(mt>=3).sum()}")
line(f"distance distribution (m): min={inv['dist_m'].min():.2f} "
     f"median={inv['dist_m'].median():.2f} max={inv['dist_m'].max():.2f}")
line(f"exact (0 m) pairs: {(inv['dist_m']==0).sum()}")
line("")
line("event_bracketed counts:")
line(inv['event_bracketed'].value_counts().to_string())
line("")
line("time gap (months) summary:")
line(inv['dt_months'].describe().to_string())

line("")
line("="*70)
line("2-4. DEPTH-RESOLVED SIGNAL AND CONFOUNDS (point-level)")
line("="*70)
g = grid.dropna(subset=["dqc1N","pre_qc1N","post_qc1N","depth_m"]).copy()
# attach pair-level time gap + event bracket
g = g.merge(inv[["pair_id","dt_months","event_bracketed","dist_m","coord"]],
            on="pair_id", how="left")
line(f"depth-point observations (clean): {len(g)}")
line(f"pairs represented: {g['pair_id'].nunique()}")
line(f"depth range: {g['depth_m'].min():.2f} - {g['depth_m'].max():.2f} m")

# 2. depth dependence
line("")
line("-- dqc1N by depth bin --")
g["zbin"] = pd.cut(g["depth_m"], [0,2,4,6,8,10,15,30])
line(g.groupby("zbin", observed=True)["dqc1N"].agg(["count","mean","median","std"]).to_string())

# 3. Ic dependence (use pre_ic as the soil type at baseline)
line("")
line("-- dqc1N by pre_ic (soil behaviour type) bin --")
g["icbin"] = pd.cut(g["pre_ic"], [0,1.31,2.05,2.6,2.95,4.0],
                    labels=["gravel/sand(<1.31)","clean-silty sand(1.31-2.05)",
                            "silt mix(2.05-2.6)","silty clay(2.6-2.95)","clay(>2.95)"])
line(g.groupby("icbin", observed=True)["dqc1N"].agg(["count","mean","median"]).to_string())

# 4. REGRESSION TO THE MEAN -the key confound
line("")
line("-- CONFOUND: dqc1N vs pre_qc1N (regression to mean test) --")
from scipy import stats
r_rtm = stats.spearmanr(g["pre_qc1N"], g["dqc1N"])
line(f"Spearman(pre_qc1N, dqc1N) = {r_rtm.statistic:.3f} (p={r_rtm.pvalue:.2e})")
line("  -> strong negative would mean low-initial points 'recover' spuriously.")
line("dqc1N by pre_qc1N quintile:")
g["q_pre"] = pd.qcut(g["pre_qc1N"], 5, duplicates="drop")
line(g.groupby("q_pre", observed=True)["dqc1N"].agg(["count","mean","median"]).to_string())
# also relative change to separate proportional growth from mean reversion
g["rel_change"] = g["dqc1N"] / g["pre_qc1N"].clip(lower=1)
line("")
line("relative change dqc1N/pre_qc1N by pre_qc1N quintile (median):")
line(g.groupby("q_pre", observed=True)["rel_change"].median().to_string())

line("")
line("="*70)
line("5. TIME DEPENDENCE controlling for depth/Ic/initial (OLS, point-level)")
line("="*70)
# Crude pooled OLS just to see partial signs; NOT the final model.
import statsmodels.formula.api as smf
gg = g.dropna(subset=["dt_months","pre_ic"]).copy()
gg["log_dt"] = np.log(gg["dt_months"].clip(lower=0.1))
m = smf.ols("dqc1N ~ dt_months + depth_m + pre_ic + pre_qc1N", data=gg).fit(
    cov_type="cluster", cov_kwds={"groups": gg["coord"]})
line("Pooled OLS dqc1N ~ dt_months + depth + pre_ic + pre_qc1N (cluster-robust by coord):")
line(m.params.to_string())
line("p-values:")
line(m.pvalues.to_string())
line("")
m2 = smf.ols("dqc1N ~ log_dt + depth_m + pre_ic + pre_qc1N", data=gg).fit(
    cov_type="cluster", cov_kwds={"groups": gg["coord"]})
line("With log(dt) instead (aging-style):")
line(m2.params.to_string())
line(f"R2 linear-dt={m.rsquared:.3f}  R2 log-dt={m2.rsquared:.3f}")

line("")
line("="*70)
line("6. WITHIN-COORDINATE TRUE-REPEAT TRAJECTORIES (rate-identifying)")
line("="*70)
mt_coords = mt[mt>=2].index.tolist()
line(f"rate-identifying coordinates: {len(mt_coords)}")
for c in mt_coords:
    sub = inv[inv["coord"]==c].sort_values("date_post")
    sp = sub["sounding_pre"].iloc[0]
    line(f"\ncoord {c} (pre={sp}), {len(sub)} pairs:")
    for _, r in sub.iterrows():
        line(f"   {r['date_pre'].date()} -> {r['date_post'].date()} "
             f"({r['dt_months']:.2f} mo, {r['dist_m']:.1f} m): "
             f"dqc1N(median)={r['dqc1N']:.1f}  pre_qc1N={r['qc1N_pre']:.0f}")

# Save merged depth-level dataset for modeling
g.to_csv(OUT / "depth_level_merged_R39.csv", index=False)
line("")
line(f"[saved] {OUT/'depth_level_merged_R39.csv'}  ({len(g)} rows)")
