# R39 Module 5 — field-scale re-liquefaction test

ShakeMaps loaded (Darfield 2010, Christchurch Feb 2011).
polygons: Sept2010 3000, Feb2011 6144

## (A) CPT-controlled sites in mapped domain: 74
manifest_S  manifest_F
0           0             52
            1             20
1           1              2
  P(Feb manifest | Sept manifest)   = 1.00
  P(Feb manifest | NO Sept manifest)= 0.28

  logit P(Feb manifest) ~ qc1Ncs + ln(Feb PGA) + Sept_manifest:
    q            coef=-0.013  OR=0.99  p=0.0536
    lnpgaF       coef=-3.244  OR=0.04  p=0.0111
    manifest_S   coef=+28.045  OR=1512365056185.69  p=1
    -> Sept_manifest coef>0 (controlling demand+resistance) = prior liquefaction
       ADDS re-liquefaction risk beyond static susceptibility (net weakening);
       <0 = net strengthening; ~0 = susceptibility persists, no extra event effect.