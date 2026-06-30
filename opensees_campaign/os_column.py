"""1D effective-stress soil column (Christchurch profile), PM4Sand SSPquadUP stack.
FIX vs v1: pin pore-pressure at ALL nodes during gravity (fully drained, exactly like the
validated single element), then free interior p for the dynamic stage. Validate sigma'_v(z)."""
import openseespy.opensees as op
import numpy as np
from pathlib import Path

g = 9.81; gw = 9.81
H = 12.0; dz = 0.5; nE = int(H/dz); nL = nE + 1
wt_depth = 1.5
layers = [(0.0, 2.0, 0.70), (2.0, 6.0, 0.45), (6.0, 12.0, 0.78)]
def Dr_at(depth):
    for z0, z1, Dr in layers:
        if z0 - 1e-6 <= depth < z1 + 1e-6: return Dr
    return layers[-1][2]

op.wipe(); op.model('basic', '-ndm', 2, '-ndf', 3)
for i in range(nL):
    op.node(2*i+1, 0.0, i*dz); op.node(2*i+2, 1.0, i*dz)
for i in range(nL):
    op.equalDOF(2*i+1, 2*i+2, 1, 2)               # periodic: 1D shear column
# base: fixed in u; pore pressure pinned everywhere during gravity (drained)
op.fix(1, 1, 1, 1); op.fix(2, 1, 1, 1)
for i in range(1, nL):
    op.fix(2*i+1, 0, 0, 1); op.fix(2*i+2, 0, 0, 1)
for e in range(1, nE+1):
    depth = H - (e-0.5)*dz; Dr = Dr_at(depth); rho = 1.42 + 0.3*Dr
    op.nDMaterial('PM4Sand', e, Dr, 476.0, 0.53, rho, 101.325, -1.0, 0.8, 0.5,
                  0.5, 0.1, -1.0, -1.0, 250.0, -1.0, 26.0, 0.3333, 2.0, -1.0, -1.0,
                  10.0, 1.5, 0.01, -1.0, -1.0)
    gamma_eff = rho*g - (gw if depth > wt_depth else 0.0)
    op.element('SSPquadUP', e, 2*(e-1)+1, 2*(e-1)+2, 2*e+2, 2*e+1, e,
               1.0, 2.2e6, 1.0, 1e-5, 1e-5, 0.6, 1.0e-9, 0.0, -gamma_eff)

damp, w1, w2 = 0.05, 0.2, 25.0; a1 = 2*damp/(w1+w2); a0 = a1*w1*w2
op.constraints('Transformation'); op.test('NormDispIncr', 1e-5, 50, 0)
op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('UmfPack')
op.integrator('Newmark', 0.6, 0.3025); op.rayleigh(a0, a1, 0.0, 0.0)
op.analysis('Transient')
for e in range(1, nE+1): op.updateMaterialStage('-material', e, '-stage', 0)
ok1 = op.analyze(30, 100.0)
for e in range(1, nE+1): op.updateMaterialStage('-material', e, '-stage', 1)
ok2 = op.analyze(50, 100.0)
print("gravity ok:", ok1 == 0 and ok2 == 0)

prof = sorted([(H-(e-0.5)*dz, -op.eleResponse(e, 'stress')[1]) for e in range(1, nE+1)])
prof = np.array(prof)
sv_exp = np.array([(g*0.45)*min(d, wt_depth) + (g*(1.7-1.0))*max(0, d-wt_depth) for d, _ in prof])
print("depth  sigma'v_sim  sigma'v_expected(gamma'~6.9)")
for (d, sv), se in list(zip(prof, sv_exp))[::3]:
    print(f"  {d:4.1f}   {sv:6.1f}      {se:6.1f}")
mono = bool(np.all(np.diff(prof[:, 1]) >= -2)) if not np.isnan(prof[:, 1]).any() else False
print(f"base sigma'_v = {prof[-1,1]:.0f} kPa | monotonic w/ depth: {mono} | NaN: {np.isnan(prof[:,1]).any()}")
np.savetxt(Path(__file__).resolve().parent / "col_profile.txt", prof)
