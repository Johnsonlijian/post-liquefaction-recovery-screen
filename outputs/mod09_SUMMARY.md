# R39 Module 9 — buffered (20 m) manifestation join + validated test

sites: 845 | Sept manifest 80 (0.09) | Feb manifest 675 (0.80)

## VALIDATION (susceptibility must predict manifestation)
  qc1Ncs (lower=weaker): AUC=0.378
  LSN_Feb (higher=worse): AUC=0.613
  Feb manifest by Easting band (E floodplain should be highest):
           mean  size
eband                
W      0.712264   212
CW     0.739336   211
CE     0.796209   211
E      0.947867   211

  -> susceptibility STILL INVALID (do not trust test)

## DENSITY-PARADOX TEST (does prior liquefaction add risk beyond susceptibility?)

[ALL, control=qc1Ncs] n=845  P(Feb|Sept)=0.64 P(Feb|noSept)=0.82
    qc1Ncs  OR=1.01 (95%CI 1.00-1.01) p=0.143
    lnF     OR=3.14 (95%CI 1.99-4.96) p=9.55e-07
    mS      OR=1.24 (95%CI 0.65-2.38) p=0.512

[pre-Feb(2010), control=qc1Ncs] n=300  P(Feb|Sept)=0.68 P(Feb|noSept)=0.71
    qc1Ncs  OR=1.01 (95%CI 1.00-1.02) p=0.12
    lnF     OR=6.77 (95%CI 3.39-13.54) p=6.28e-08
    mS      OR=4.57 (95%CI 1.97-10.62) p=0.000417

[ALL, control=LSN_F] n=845  P(Feb|Sept)=0.64 P(Feb|noSept)=0.82
    LSN_F   OR=1.06 (95%CI 1.03-1.10) p=7.29e-05
    mS      OR=0.50 (95%CI 0.30-0.82) p=0.00648

[pre-Feb(2010), control=LSN_F] n=300  P(Feb|Sept)=0.68 P(Feb|noSept)=0.71
    LSN_F   OR=1.07 (95%CI 1.02-1.12) p=0.00228
    mS      OR=1.18 (95%CI 0.64-2.18) p=0.597