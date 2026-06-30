"""Validated PM4Sand undrained cyclic DSS single element (faithful to OpenSeesPy docs
pm4sand_cyc_cal: 3-segment uniform-cycle loading). Returns N to liquefaction
(N at 3% single-amplitude shear strain, the standard criterion) and ru history.
Foundation for the Canterbury CSR-N + site-response campaign."""
import openseespy.opensees as op
import numpy as np

def cyclic_dss(Dr=0.55, CSR=0.16, sig_v0_kpa=100.0, Cycle_max=60, strain_in=1.0e-5, verbose=False):
    atm = -101.325
    sig_v0 = sig_v0_kpa * atm / 101.325                # compression negative
    devDisp = 0.03; perm = 1e-9                         # 3% strain = liquefaction criterion
    G0 = 476.0; hpo = 0.53; rho = 1.42; P_atm = 101.325
    e_max, e_min = 0.8, 0.5; e_ini = e_max - (e_max - e_min) * Dr
    op.wipe(); op.model('basic', '-ndm', 2, '-ndf', 3)
    op.node(1, 0., 0.); op.node(2, 1., 0.); op.node(3, 1., 1.); op.node(4, 0., 1.)
    op.fix(1, 1, 1, 1); op.fix(2, 1, 1, 1); op.fix(3, 0, 0, 1); op.fix(4, 0, 0, 1)
    op.equalDOF(3, 4, 1, 2)
    op.nDMaterial('PM4Sand', 1, Dr, G0, hpo, rho, P_atm, -1.0, e_max, e_min,
                  0.5, 0.1, -1.0, -1.0, 250.0, -1.0, 26.0, 0.3333, 2.0, -1.0, -1.0,
                  10.0, 1.5, 0.01, -1.0, -1.0)
    op.element('SSPquadUP', 1, 1, 2, 3, 4, 1, 1.0, 2.2e6, 1.0, perm, perm, e_ini, 1.0e-5)
    damp, w1, w2 = 0.02, 0.2, 20.0
    a1 = 2.0 * damp / (w1 + w2); a0 = a1 * w1 * w2
    op.constraints('Transformation'); op.test('NormDispIncr', 1e-5, 35, 0)
    op.algorithm('Newton'); op.numberer('RCM'); op.system('FullGeneral')
    op.integrator('Newmark', 5.0/6.0, 4.0/9.0); op.rayleigh(a1, a0, 0.0, 0.0)
    op.analysis('Transient')
    pN = sig_v0 / 2.0
    op.timeSeries('Path', 1, '-values', 0, 1, 1, '-time', 0.0, 100.0, 1e10)
    op.pattern('Plain', 1, 1, '-factor', 1.0); op.load(3, 0., pN, 0.); op.load(4, 0., pN, 0.)
    op.updateMaterialStage('-material', 1, '-stage', 0); op.analyze(100, 1.0)
    vD = op.nodeDisp(3, 2)
    op.timeSeries('Path', 2, '-values', 1, 1, 1, '-time', 100.0, 80000.0, 1e10)
    op.pattern('Plain', 2, 2, '-factor', 1.0); op.sp(3, 2, vD); op.sp(4, 2, vD)
    for i in range(4): op.remove('sp', i + 1, 3)
    op.analyze(25, 1.0)
    op.updateMaterialStage('-material', 1, '-stage', 1)
    op.setParameter('-val', 0, '-ele', 1, 'FirstCall', '1'); op.analyze(25, 1.0)
    op.setParameter('-val', 0.3, '-ele', 1, 'poissonRatio', '1')
    sigv0 = op.eleResponse(1, 'stress')[1]
    controlDisp = 1.1 * devDisp; numCycle = 0.0; tag = 3; CAP = 400000
    trace = []                                          # (cycle, ru, |gamma|)
    b = op.eleResponse(1, 'stress')

    def seg(target_disp, cond):
        nonlocal b
        h = op.nodeDisp(3, 1); ct = op.getTime()
        tc = ct + (abs(target_disp) + abs(h)) / strain_in
        op.timeSeries('Path', tag, '-values', h, target_disp, target_disp, '-time', ct, tc, 1e10, '-factor', 1.0)
        op.pattern('Plain', tag, tag); op.sp(3, 1, 1.0)
        cnt = 0; liq = False
        while cond(b[2]) and cnt < CAP:
            if op.analyze(1, 1.0) != 0:
                op.algorithm('ModifiedNewton')
                if op.analyze(1, 1.0) != 0: op.algorithm('Newton'); break
                op.algorithm('Newton')
            b = op.eleResponse(1, 'stress'); cnt += 1
            if abs(op.nodeDisp(3, 1)) >= devDisp: liq = True; break
        op.remove('loadPattern', tag); op.remove('timeSeries', tag); op.remove('sp', 3, 1)
        return liq

    Nliq = None
    while numCycle <= Cycle_max:
        for tgt, cond, inc in ((controlDisp, lambda s: s <= CSR * sig_v0 * -1.0, 0.25),
                               (-controlDisp, lambda s: s > CSR * sig_v0, 0.50),
                               (controlDisp, lambda s: s <= 0.0, 0.25)):
            liq = seg(tgt, cond)
            numCycle += inc
            ru = 1 - op.eleResponse(1, 'stress')[1] / sigv0
            g = abs(op.nodeDisp(3, 1))
            trace.append((numCycle, ru, g))
            if Nliq is None and g >= devDisp:           # 3% strain criterion
                Nliq = numCycle
        if Nliq is not None:
            break
    tr = np.array(trace) if trace else np.zeros((1, 3))
    return dict(sigv0=sigv0, ru_max=float(tr[:, 1].max()), ncyc=numCycle, Nliq=Nliq, trace=tr)


def crr_idriss_boulanger(N, qc1Ncs=None):
    """Idriss-Boulanger MSF-based CRR(N)/CRR(15) for reference (b=0.34 medium dense)."""
    b = 0.34
    return (15.0 / N) ** b


if __name__ == "__main__":
    import sys
    CSR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.16
    Dr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.55
    r = cyclic_dss(Dr=Dr, CSR=CSR)
    print(f"Dr={Dr} CSR={CSR}: sig'v0={r['sigv0']:.1f} | ru_max={r['ru_max']:.3f} | "
          f"N_liq(3% strain)={r['Nliq']} | cycles run={r['ncyc']:.2f}")
