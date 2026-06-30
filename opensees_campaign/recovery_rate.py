# Convert simulated per-event reconsolidation densification -> resistance-gain (recovery) rate,
# using the calibrated CRR15(Dr) curve, and compare to the demand gap.  Honest about settlement source.
import json, os, numpy as np, pandas as pd
from pathlib import Path
SCR = Path(__file__).resolve().parent

# ---- CRR15(Dr) from the calibrated CSR-N sims (N=15 cycles), incl. directly-anchored low Dr ----
df = pd.concat([pd.read_csv(os.path.join(SCR, "csrn_low.csv")),
                pd.read_csv(os.path.join(SCR, "csrn_curves.csv"))], ignore_index=True)
def crr15(Dr):
    s = df[(df.Dr == Dr) & df.N_liq.notna() & (df.N_liq > 1)].sort_values('CSR')
    N = s.N_liq.values; C = s.CSR.values
    return float(np.interp(15.0, N[::-1], C[::-1]))
Drs = sorted(df.Dr.unique())
CRR = {d: crr15(d) for d in Drs}
print("CRR15 from calibrated PM4Sand CSR-N:", {round(d, 2): round(v, 3) for d, v in CRR.items()})
slope = (crr15(0.40) - crr15(0.30)) / 0.10          # LOCAL dCRR15/dDr at the in-situ (loose) density
print(f"local dCRR15/dDr at Dr~0.30-0.40 ~ {slope:.3f} per unit Dr  (flatter than high-Dr range)")

# ---- per-event reconsolidation densification from the independent runs ----
res = json.load(open(os.path.join(SCR, "events_independent.json"), encoding="utf-8"))
H_liq = 5.0          # liquefiable layer thickness (m)
e0 = 0.78            # in-situ void ratio of liquefiable sand (mat2)
emax, emin = 0.90, 0.50
Dr_insitu = (emax - e0) / (emax - emin)
print(f"\nin-situ Dr of liquefiable layer ~ {Dr_insitu:.2f} (e0={e0})")
print(f"\n{'event':20s} peak_ru  settle(mm)  eps_v(%)  dDr(%)  dCRR15")
dDr_list = []
for r in res:
    if r.get('settle_mm') is None: continue
    epsv = r['settle_mm'] / (H_liq * 1000)          # volumetric strain of liq layer
    dDr = epsv * (1 + e0) / (emax - emin)           # densification -> Dr increment
    dCRR = slope * dDr
    dDr_list.append(dDr)
    print(f"{r['name']:20s}  {r['peak_ru']:.2f}     {r['settle_mm']:.2f}      "
          f"{epsv*100:.3f}    {dDr*100:.2f}    {dCRR:.4f}")

dDr_ev = float(np.mean(dDr_list)); dCRR_ev = slope * dDr_ev
# demand gap: lift CRR15 from in-situ Dr up to the weakest liquefying demand (CSR ~ 0.18 here)
CRR_insitu = crr15(0.30)                               # DIRECTLY computed (not extrapolated)
CSR_demand = 0.18
gap = CSR_demand - CRR_insitu
print(f"\nrecovery per event: dDr ~ {dDr_ev*100:.2f}%, dCRR15 ~ {dCRR_ev:.4f}")
print(f"CRR15 at in-situ Dr ~ {CRR_insitu:.3f}; demand CSR ~ {CSR_demand:.2f}; gap ~ {gap:.3f}")
print(f"events needed to close the gap by densification alone: {gap/dCRR_ev:.0f}")
# sequence timescale: 4 major events in ~16 months (Sep2010 - Dec2011)
yrs = 16/12.0
rate_per_yr = dCRR_ev * (len([r for r in res if r.get('settle_mm')]) / yrs)
print(f"simulated recovery rate ~ {rate_per_yr:.4f} CRR15/yr; to close gap needs "
      f"~{gap/rate_per_yr:.0f} yr  vs sequence {yrs:.1f} yr  => recovery {gap/rate_per_yr/yrs:.0f}x too slow")
print("\nNOTE: PM4Sand/PDMY02 reconsolidation strain is modest; Ishihara-Yoshimine (1992) gives "
      "eps_v~2-3% for Dr~0.3 full liquefaction -> even the GENEROUS densification is far too slow.")
# Ishihara-Yoshimine generous bound
for epsv_gen in (0.02, 0.03):
    dDr_g = epsv_gen * (1 + e0) / (emax - emin); dCRR_g = slope * dDr_g
    print(f"  IY eps_v={epsv_gen*100:.0f}%: dDr~{dDr_g*100:.1f}%/event, dCRR~{dCRR_g:.3f}, "
          f"events-to-close~{gap/dCRR_g:.0f}")
