# R39 Module 7 — field-scale density-paradox test (open data)

bridged CPTs with coordinates: 857
sites with computable critical-layer qc1Ncs + spatial join: 845
  qc1Ncs: median 91 (IQR 81-108)
  PGA(%g): Sept median 22, Feb median 69
  manifest: Sept 54, Feb 245

contingency (Sept x Feb):
mS  mF
0   0     571
    1     220
1   0      29
    1      25

[ALL bridged CPTs (mixed measurement year)] n=845  P(Feb|Sept)=0.46  P(Feb|noSept)=0.28
    q    coef=+0.005 OR=1.01 (95%CI 1.00-1.01) p=0.0798
    lnF  coef=-0.391 OR=0.68 (95%CI 0.47-0.98) p=0.0405
    mS   coef=+0.509 OR=1.66 (95%CI 0.86-3.20) p=0.128

[pre-Feb subset (2010 soundings)] n=300  P(Feb|Sept)=0.49  P(Feb|noSept)=0.16
    q    coef=+0.001 OR=1.00 (95%CI 0.99-1.01) p=0.931
    lnF  coef=-0.955 OR=0.38 (95%CI 0.20-0.73) p=0.00352
    mS   coef=+0.922 OR=2.51 (95%CI 1.15-5.52) p=0.0216

[<=2011 soundings] n=703  P(Feb|Sept)=0.48  P(Feb|noSept)=0.19
    q    coef=+0.000 OR=1.00 (95%CI 0.99-1.01) p=0.95
    lnF  coef=-0.494 OR=0.61 (95%CI 0.39-0.96) p=0.0338
    mS   coef=+0.905 OR=2.47 (95%CI 1.20-5.08) p=0.0139

-> mS (Sept_manifest) OR>1 controlling resistance+demand = field evidence for
   net weakening (prior liquefaction adds re-liquefaction risk); OR<1 = strengthening.