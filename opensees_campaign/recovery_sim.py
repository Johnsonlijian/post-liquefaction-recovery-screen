"""CORE physics: element-level liquefaction -> reconsolidation -> re-liquefaction.
Continuous PM4Sand state evolution. Measures N1 (first event) vs N2 (re-liquefaction after
one reconsolidation), and the reconsolidation volumetric strain (densification). Directly:
(a) simulated recovery from one event; (b) density-paradox test of susceptibility change."""
import openseespy.opensees as op
import numpy as np

def model(Dr=0.50, sig_v0_kpa=100.0):
    atm = -101.325; sig_v0 = sig_v0_kpa*atm/101.325
    op.wipe(); op.model('basic', '-ndm', 2, '-ndf', 3)
    op.node(1, 0., 0.); op.node(2, 1., 0.); op.node(3, 1., 1.); op.node(4, 0., 1.)
    op.fix(1, 1, 1, 1); op.fix(2, 1, 1, 1); op.fix(3, 0, 0, 1); op.fix(4, 0, 0, 1)
    op.equalDOF(3, 4, 1, 2)
    op.nDMaterial('PM4Sand', 1, Dr, 476.0, 0.53, 1.42, 101.325, -1.0, 0.8, 0.5,
                  0.5, 0.1, -1.0, -1.0, 250.0, -1.0, 26.0, 0.3333, 2.0, -1.0, -1.0,
                  10.0, 1.5, 0.01, -1.0, -1.0)
    op.element('SSPquadUP', 1, 1, 2, 3, 4, 1, 1.0, 2.2e6, 1.0, 1e-9, 1e-9, 0.6, 1.0e-9)
    damp, w1, w2 = 0.02, 0.2, 20.0; a1 = 2*damp/(w1+w2); a0 = a1*w1*w2
    op.constraints('Transformation'); op.test('NormDispIncr', 1e-5, 35, 0)
    op.algorithm('Newton'); op.numberer('RCM'); op.system('FullGeneral')
    op.integrator('Newmark', 5./6., 4./9.); op.rayleigh(a1, a0, 0., 0.); op.analysis('Transient')
    # consolidation, drained
    op.timeSeries('Path', 1, '-values', 0, 1, 1, '-time', 0., 100., 1e10)
    op.pattern('Plain', 1, 1, '-factor', 1.0); op.load(3, 0., sig_v0/2, 0.); op.load(4, 0., sig_v0/2, 0.)
    op.updateMaterialStage('-material', 1, '-stage', 0); op.analyze(100, 1.0)
    return sig_v0

def close_drain():
    vD = op.nodeDisp(3, 2)
    op.timeSeries('Path', 2, '-values', 1, 1, 1, '-time', 100., 8e4, 1e10)
    op.pattern('Plain', 2, 2, '-factor', 1.0); op.sp(3, 2, vD); op.sp(4, 2, vD)
    for i in range(4): op.remove('sp', i+1, 3)
    op.analyze(25, 1.0)
    op.updateMaterialStage('-material', 1, '-stage', 1)
    op.setParameter('-val', 0, '-ele', 1, 'FirstCall', '1'); op.analyze(25, 1.0)
    op.setParameter('-val', 0.3, '-ele', 1, 'poissonRatio', '1')

def cyclic(CSR, sig_v0, sigv0_eff, Cycle_max=40, strain_in=1.0e-5):
    devDisp = 0.03; controlDisp = 1.1*devDisp; tag = 9; numCycle = 0.0; Nliq = None
    b = op.eleResponse(1, 'stress')
    def seg(td, cond):
        nonlocal b
        h = op.nodeDisp(3, 1); ct = op.getTime(); tc = ct+(abs(td)+abs(h))/strain_in
        op.timeSeries('Path', tag, '-values', h, td, td, '-time', ct, tc, 1e10, '-factor', 1.0)
        op.pattern('Plain', tag, tag); op.sp(3, 1, 1.0); cnt = 0; liq = False
        while cond(b[2]) and cnt < 300000:
            if op.analyze(1, 1.0) != 0:
                op.algorithm('ModifiedNewton')
                if op.analyze(1, 1.0) != 0: op.algorithm('Newton'); break
                op.algorithm('Newton')
            b = op.eleResponse(1, 'stress'); cnt += 1
            if abs(op.nodeDisp(3, 1)) >= devDisp: liq = True; break
        op.remove('loadPattern', tag); op.remove('timeSeries', tag); op.remove('sp', 3, 1)
        return liq
    while numCycle <= Cycle_max:
        for td, cond, inc in ((controlDisp, lambda s: s <= CSR*sig_v0*-1, 0.25),
                              (-controlDisp, lambda s: s > CSR*sig_v0, 0.5),
                              (controlDisp, lambda s: s <= 0., 0.25)):
            seg(td, cond); numCycle += inc
            if Nliq is None and abs(op.nodeDisp(3, 1)) >= devDisp: Nliq = numCycle
        if Nliq is not None: break
    return Nliq

# ---------------- run: event1 -> reconsolidate -> event2 ----------------
Dr, CSR = 0.50, 0.16
sig_v0 = model(Dr); close_drain()
sigv0_eff = op.eleResponse(1, 'stress')[1]
uy0 = op.nodeDisp(3, 2)
N1 = cyclic(CSR, sig_v0, sigv0_eff)
print(f"EVENT 1: N_liq = {N1}  (Dr={Dr}, CSR={CSR})")

# reconsolidation: re-open drainage (pin pore pressure = 0), dissipate, densify
op.timeSeries('Constant', 20); op.pattern('Plain', 20, 20)
op.sp(3, 3, 0.0); op.sp(4, 3, 0.0)           # drain top (inside active pattern)
for k in range(120):
    if op.analyze(1, 50.0) != 0:
        op.algorithm('ModifiedNewton')
        if op.analyze(1, 50.0) != 0: op.algorithm('Newton'); break
        op.algorithm('Newton')
uy1 = op.nodeDisp(3, 2)
eps_v = (uy1 - uy0) / 1.0 * 100                # vertical (volumetric) strain, % (settlement)
sigv_recov = op.eleResponse(1, 'stress')[1]
print(f"RECONSOLIDATION: volumetric strain (densification) = {abs(eps_v):.2f}% | "
      f"sigma'_v recovered to {-sigv_recov:.1f} kPa")

# event 2: re-liquefy the reconsolidated (denser, evolved-fabric) element
op.remove('sp', 3, 3); op.remove('sp', 4, 3)
N2 = cyclic(CSR, sig_v0, sigv0_eff)
print(f"EVENT 2 (re-liquefaction): N_liq = {N2}")
if N1 and N2:
    print(f"\n=> re-liquefaction resistance change: N2/N1 = {N2/N1:.2f}  "
          f"({'MORE resistant (densification wins)' if N2 > N1*1.1 else 'LESS resistant (fabric)' if N2 < N1*0.9 else 'essentially UNCHANGED'})")
