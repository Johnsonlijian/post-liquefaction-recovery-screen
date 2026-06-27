# R39 Module 1 (rigorous) — within-coordinate depth-resolved recovery rate

coordinates with >=2 post dates: 6; with >=3 post dates (stable slope): 4
proper pair-level observations (liquefiable Ic<=2.6 window medians): 35

## (A) Mean reversion: present in cross-section, removed by within-coord time slope
cross-section Spearman(pre_qc1N, dqc1N) = -0.232 (p=1.8e-01)  [confound present]

Per-coordinate within-coordinate recovery slope (>=3 posts, proper units):
  sounding  n_posts  pre_qc1N  slope    se     p  dt_max
CPT-BUR-18        4    75.204  3.461 0.661 0.035   7.162
CPT-BUR-22        4    42.539  2.823 0.701 0.056  24.869
CPT-BUR-26        3    68.707  0.907 3.670 0.846   6.373
CPT-SSH-09        5    96.659  0.999 0.626 0.209  23.160
-> 4/4 positive; median 1.91 qc1N/month.
   Spearman(pre_qc1N, slope)=0.00 (p=1.00); slope NOT a function of baseline (unlike cross-section).

## (B) Aggregate within-coordinate rate (coordinate fixed effects) + kinetic form
linear : 1.95 qc1N/month (95%CI 0.32..3.59); AIC 186.0
log-time: 20.44 qc1N/ln(month) (95%CI 2.49..38.40); AIC 186.2
-> AIC favors linear (dAIC=0.2); short baselines (max 24.9 mo) -> form weakly identified; physical prior = log-time aging.

## (C) Depth mechanism: where does recovery rate concentrate?

beta(z) — median within-coordinate recovery rate by depth (units = the 4 >=3-post coords):
       zbin  z_mid  n_coord  beta_median  frac_pos
 (1.0, 2.0]   1.48        3        -0.01      0.33
 (2.0, 3.0]   2.47        4         0.91      0.50
 (3.0, 4.0]   3.50        4         4.31      1.00
 (4.0, 5.0]   4.47        3         2.27      0.67
 (5.0, 6.0]   5.60        3         4.49      0.67
 (6.0, 7.0]   6.58        3        -0.12      0.33
 (7.0, 8.0]   7.54        3        -0.96      0.33
 (8.0, 9.0]   8.53        3        -1.88      0.33
(9.0, 10.0]   9.52        3         0.38      0.67

Mixed model (random intercept per profile; mean-reversion control): does rate vary with depth?
   dt_c       coef=  -0.000  p=nan
   z_c        coef=   3.200  p=1.06e-17
   dt_c:z_c   coef=  -0.037  p=0.534
   pre_c      coef= -18.904  p=1.5e-44
   pre_ic     coef=  39.020  p=7.36e-20
   (dt_c:z_c = depth-dependence of the recovery rate; pre_c<0 confirms mean-reversion control active)