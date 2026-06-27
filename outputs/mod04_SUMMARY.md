# R39 Module 4 — event-free identifiability test

Inter-event windows: 0=pre-Darfield,1=Darfield-Feb,2=Feb-Jun,3=Jun-Dec,4=post-Dec

clean event-free segments (same coord, same inter-event window, >=2 dates): 6
  distinct coordinates: 6
  segment durations (months): min 0.30, median 0.43, max 0.46
  segments with duration >= 1 month: 0
  segments with duration >= 3 months: 0

  by window:
     n  coords   med_dt   med_rate
win                               
2    6       6  0.42707  25.619953

  all clean segments:
                coord  win         d0         d1  dt_months     dq   rate
1573171.14,5181557.65    2 2011-03-01 2011-03-14       0.43 100.19 234.60
1573678.78,5183942.60    2 2011-02-28 2011-03-09       0.30  -8.45 -28.58
1573681.08,5183528.51    2 2011-02-28 2011-03-14       0.46  17.40  37.84
1573725.52,5181312.16    2 2011-02-25 2011-03-10       0.43  18.14  42.47
1574310.74,5183287.96    2 2011-02-28 2011-03-14       0.46   6.16  13.40
1579212.13,5178865.59    2 2011-03-01 2011-03-10       0.30   0.95   3.21

## VERDICT
Clean event-free segments >=1 month: 0. If small, a clean single-process
recovery RATE is NOT identifiable from these event-straddling repeats; the
calendar-time slope must be reported as a NET multi-event sequence trajectory,
not a recovery-kinetics rate. This is the honest, review-driven conclusion.