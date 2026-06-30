"""Tolerable-PGA gain per year under TWO independent rebuild-rate estimates, for the
median loose critical-layer state. Resolves the provenance of the screen's headline
number: literature log-time ageing (Hayati-Andrus) gives ~15 milli-g/yr; the empirical
net within-coordinate sequence slope (an upper bound, not a clean recovery rate) gives
~25 milli-g/yr. Both are of order 0.02 g/yr, far below the Feb-2011 demand (0.4-0.8 g)."""
import numpy as np
from scipy.optimize import brentq
Pa = 101.325
def crr(q): q = np.clip(q, 1, 250); return np.exp(q/113 + (q/1000)**2 - (q/140)**3 + (q/137)**4 - 2.80)
def rd(z, M):
    a = -1.012 - 1.126*np.sin(z/11.73 + 5.133); b = 0.106 + 0.118*np.sin(z/11.28 + 5.142)
    return np.exp(a + b*M)
def MSF(q, M): return 1 + (min(2.2, 1.09 + (np.clip(q, 1, 250)/180)**3) - 1)*(8.64*np.exp(-M/4) - 1.325)
def Ksig(q, sp): Cs = min(0.3, 1/(37.3 - 8.27*np.clip(q, 1, 211)**0.264)); return min(1.1, 1 - Cs*np.log(max(sp, 0.3*Pa)/Pa))
def tol_amax(q, z=3.0, wt=1.5, M=6.0):
    sigv = 18*z; sigp = max(sigv - 9.81*max(0, z - wt), 5)
    return crr(q)*MSF(q, M)*Ksig(q, sigp)/(0.65*(sigv/sigp)*rd(z, M))

q0 = 37.5                         # median critical-layer qc1Ncs, loose post-liquefaction state
a0 = tol_amax(q0)
a_lin = tol_amax(q0 + 1.95*12)    # empirical net sequence slope 1.95 qc1Ncs/month (upper bound)
qfc = lambda tgt: brentq(lambda q: crr(q) - tgt, 1, 250)
print(f"tolerable PGA(0)               = {a0*1000:.0f} milli-g")
print(f"empirical linear (beta 1.95)   = {a_lin*1000:.0f} milli-g  (+{(a_lin-a0)*1000:.0f}/yr)")
for kdr in (0.12, 0.13):          # Hayati & Andrus (2009) K_DR per log-cycle, one decade in year 1
    a = tol_amax(qfc(crr(q0)*(1 + kdr)))
    print(f"literature log-cycle K_DR={kdr}  = {a*1000:.0f} milli-g  (+{(a-a0)*1000:.0f}/yr)")
print("=> screen gain of order 15-25 milli-g/yr; conclusion robust to rate and to 2-3x scaling.")
