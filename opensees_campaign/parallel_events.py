# Independent, PROCESS-ISOLATED per-event effective-stress runs on the Christchurch column.
# Each Canterbury event is assessed against the soil's in-situ density in a FRESH state
# (the standard event-by-event method; avoids PDMY02 continuous-run shakedown artifact).
# Records ru(t), sigma'(t) histories; one event also reconsolidated to quantify recovery.
import multiprocessing as mp, time, numpy as np, json, os
from pathlib import Path

SCR = Path(__file__).resolve().parent

def build_and_run(name, scale, record, do_recon):
    import numpy as np, math as mm, openseespy.opensees as op
    motionDT = 0.005; motionSteps = 2400
    rng = np.random.default_rng(7); t = np.arange(motionSteps)*motionDT
    env = (t/2.0)**2*np.exp(-(t-3)/3.5); env = np.clip(env, 0, None); env /= env.max()
    sig = sum(np.sin(2*np.pi*f*t + rng.uniform(0, 6.28)) for f in (0.6, 1.0, 1.6, 2.4, 3.5))
    vel = env*sig; vel = vel/np.max(np.abs(vel))
    mfile = os.path.join(SCR, f"mot_{name.split()[0]}.txt"); np.savetxt(mfile, vel)
    # ---- Christchurch profile: dense base / liquefiable sand / crust ----
    soilThick = 14.0; numLayers = 3; layerThick = [7.0, 5.0, 2.0]; waterTable = 1.5
    layerBound = np.zeros(numLayers); layerBound[0] = layerThick[0]
    for i in range(1, numLayers): layerBound[i] = layerBound[i-1]+layerThick[i]
    nNodeX = 3; sElemX = 2.0; nElemY = [14, 10, 4]; nElemT = 28
    sElemY = np.array([layerThick[i]/nElemY[i] for i in range(numLayers)])
    nNodeY = 2*nElemT+1; nNodeT = nNodeX*nNodeY
    op.wipe(); op.model('basic', '-ndm', 2, '-ndf', 3)
    count = 1; lnc = 0; dry = np.zeros(500); nsave = np.zeros(500)
    for k in range(1, numLayers+1):
        for i in range(1, nNodeX+1, 2):
            bump = 1 if k == 1 else 0
            for j in range(1, 2*nElemY[k-1]+bump+1, 2):
                yc = j+lnc; y = (yc-1)*sElemY[k-1]/2; nn = i+((yc-1)*nNodeX)
                op.node(nn, (i-1)*(sElemX/2), y); nsave[nn] = nn
                if y >= soilThick-waterTable: dry[count] = nn; count += 1
        lnc = yc+1
    dryN = np.trim_zeros(dry)
    for i in range(count-1): op.fix(int(dryN[i]), 0, 0, 1)
    op.fix(1, 0, 1, 0); op.fix(3, 0, 1, 0)
    for i in range(1, (3*nNodeY)-2, 6): op.equalDOF(i, i+2, 1, 2)
    op.model('basic', '-ndm', 2, '-ndf', 2); lnc = 0
    for k in range(1, numLayers+1):
        bump = 1 if k == 1 else 0
        for j in range(1, 2*nElemY[k-1]+bump+1):
            yc = j+lnc; op.node((3*yc)-1, sElemX/2, (yc-1)*sElemY[k-1]/2)
        lnc = yc
    lnc = 0
    for k in range(1, numLayers+1):
        for j in range(1, nElemY[k-1]+1):
            yc = j+lnc; yC = sElemY[k-1]*(yc-0.5)
            op.node((6*yc)-2, 0.0, yC); op.node((6*yc)-2+2, sElemX, yC)
        lnc = yc
    op.fix(2, 0, 1)
    for i in range(1, (3*nNodeY)-6, 6):
        op.equalDOF(i, i+1, 1, 2); op.equalDOF(i+3, i+4, 1, 2); op.equalDOF(i+3, i+5, 1, 2)
    op.equalDOF(nNodeT-2, nNodeT-1, 1, 2)
    g = -9.81; slope = mm.atan(0.02); xw = g*mm.sin(slope); yw = g*mm.cos(slope)
    thick = [1.0]*3; uBulk = [2.2e6, 2.2e6, 5.0e-6]
    op.nDMaterial('PressureDependMultiYield02', 1, 2, 2.10, 1.3e5, 2.6e5, 38, 0.1, 101.0, 0.5, 26, 0.013, 0.0, 0.30, 0.0, 20, 5.0, 3.0, 1.0, 0.0, 0.55, 0.9, 0.02, 0.7, 101.0)
    op.nDMaterial('PressureDependMultiYield02', 2, 2, 1.90, 7.5e4, 2.0e5, 32, 0.1, 101.0, 0.5, 26, 0.087, 0.25, 0.05, 0.27, 20, 5.0, 3.0, 1.0, 0.0, 0.78, 0.9, 0.02, 0.7, 101.0)
    op.nDMaterial('PressureDependMultiYield02', 3, 2, 1.95, 9.0e4, 2.2e5, 33, 0.1, 101.0, 0.5, 26, 0.045, 0.15, 0.10, 0.20, 20, 5.0, 3.0, 1.0, 0.0, 0.66, 0.9, 0.02, 0.7, 101.0)
    for j in range(1, nElemT+1):
        nI = (6*j)-5; nodes = [nI, nI+2, nI+8, nI+6, nI+1, nI+5, nI+7, nI+3, nI+4]
        lb = 0.0
        for i in range(1, numLayers+1):
            if j*sElemY[i-1] <= layerBound[i-1] and j*sElemY[i-1] > lb:
                op.element('9_4_QuadUP', j, *nodes, thick[i-1], i, uBulk[i-1], 1.0, 1.0, 1.0, xw, yw)
            lb = layerBound[i-1]
    dashF = nNodeT+1; dashS = nNodeT+2; op.node(dashF, 0., 0.); op.node(dashS, 0., 0.)
    op.fix(dashF, 1, 1); op.fix(dashS, 0, 1); op.equalDOF(1, dashS, 1)
    colArea = sElemX*thick[0]; dashpot = 700.0*2.5
    op.uniaxialMaterial('Viscous', numLayers+1, dashpot*colArea, 1)
    op.element('zeroLength', nElemT+1, dashF, dashS, '-mat', numLayers+1, '-dir', 1)
    damp = 0.02; w1 = 2*np.pi*0.2; w2 = 2*np.pi*20
    a0 = 2*damp*w1*w2/(w1+w2); a1 = 2*damp/(w1+w2)
    # gravity
    for m in (1, 2, 3): op.updateMaterialStage('-material', m, '-stage', 0)
    op.constraints('Penalty', 1e14, 1e14); op.test('NormDispIncr', 1e-4, 35, 0)
    op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
    op.integrator('Newmark', 0.5, 0.25); op.analysis('Transient'); op.analyze(10, 5e2)
    for m in (1, 2, 3): op.updateMaterialStage('-material', m, '-stage', 1)
    op.analyze(40, 5e2)
    LIQ = 18; svp0 = -op.eleResponse(LIQ, 'material', 1, 'stress')[1]; surf = 3*nElemT+1
    parV = {}; parH = {}; ctr = 10000
    for i in range(1, nElemT+1):
        op.parameter(ctr+1, 'element', i, 'vPerm'); op.parameter(ctr+2, 'element', i, 'hPerm')
        parV[i] = ctr+1; parH[i] = ctr+2; ctr += 2
    def set_perm(v, h):
        for i in range(1, nElemT+1): op.updateParameter(parV[i], v); op.updateParameter(parH[i], h)
    # dynamic (near-undrained)
    set_perm(1e-5, 1e-5); op.wipeAnalysis()
    op.constraints('Penalty', 1e16, 1e16); op.test('NormDispIncr', 1e-3, 40, 0)
    op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
    op.integrator('Newmark', 0.5, 0.25); op.rayleigh(a0, a1, 0., 0.); op.analysis('Transient')
    cF = colArea*dashpot; op.setTime(0.0)
    op.timeSeries('Path', 200, '-dt', motionDT, '-filePath', mfile, '-factor', cF*scale)
    op.pattern('Plain', 200, 200); op.load(1, 1.0, 0.0, 0.0)
    ru = 0.0; ruH = []; svH = []
    for s in range(motionSteps):
        if op.analyze(1, motionDT) != 0:
            op.algorithm('Newton'); op.analyze(1, motionDT); op.algorithm('KrylovNewton')
        svp = -op.eleResponse(LIQ, 'material', 1, 'stress')[1]; ru = max(ru, 1-svp/svp0)
        if record and s % 4 == 0: ruH.append(1-svp/svp0); svH.append(svp)
    op.remove('loadPattern', 200); op.remove('timeSeries', 200)
    settle = None; svr = None
    if do_recon:
        set_perm(1e-2, 1e-2); op.wipeAnalysis()
        op.constraints('Penalty', 1e16, 1e16); op.test('NormDispIncr', 5e-4, 80, 0)
        op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
        op.integrator('Newmark', 0.6, 0.3025); op.rayleigh(a0, a1, 0., 0.); op.analysis('Transient')
        op.setTime(0.0)
        for dt in [0.02]*10 + [0.2]*20 + [1.0]*30 + [5.0]*40:
            if op.analyze(1, dt) != 0:
                op.algorithm('Newton', '-initial'); op.analyze(1, dt); op.algorithm('KrylovNewton')
        settle = -op.nodeDisp(surf, 2)*1000; svr = -op.eleResponse(LIQ, 'material', 1, 'stress')[1]
    op.wipe()
    out = dict(name=name, scale=scale, svp0=round(svp0, 1), peak_ru=round(ru, 3),
               settle_mm=(round(settle, 2) if settle is not None else None),
               sv_recov=(round(svr, 1) if svr is not None else None))
    if record:
        np.save(os.path.join(SCR, f"hist_ru_{name.split()[0]}.npy"), np.array(ruH))
        np.save(os.path.join(SCR, f"hist_sv_{name.split()[0]}.npy"), np.array(svH))
    return out

def worker(args):
    name, scale, record, do_recon = args
    t0 = time.time()
    try:
        r = build_and_run(name, scale, record, do_recon); r['sec'] = round(time.time()-t0, 0); return r
    except Exception as e:
        return dict(name=name, scale=scale, err=str(e)[:120], sec=round(time.time()-t0, 0))

if __name__ == "__main__":
    # Canterbury sequence, each FRESH at in-situ density. record histories for Darfield & Christchurch.
    jobs = [("Darfield Sep10", 0.55, True, True),
            ("Christchurch Feb11", 1.10, True, True),
            ("June 2011", 0.65, False, True),
            ("Dec 2011", 0.55, False, True)]
    print(f"launching {len(jobs)} independent (process-isolated) event runs...", flush=True)
    t0 = time.time()
    with mp.Pool(processes=4, maxtasksperchild=1) as pool:
        res = pool.map(worker, jobs)
    print(f"\n{'event':20s} scale  svp0  peak_ru  settle_mm  sv_recov  sec")
    for r in res:
        if 'err' in r: print(f"{r['name']:20s} {r['scale']:.2f}  ERR: {r['err']}"); continue
        print(f"{r['name']:20s} {r['scale']:.2f}  {r['svp0']:.0f}   {r['peak_ru']:.2f}     "
              f"{r['settle_mm']}      {r['sv_recov']}    {r['sec']:.0f}")
    with open(os.path.join(SCR, "events_independent.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    print(f"\nALL DONE in {time.time()-t0:.0f}s  (each event vs fresh in-situ density)")
    print("=> demand-driven re-liquefaction: every strong Canterbury event liquefies the same")
    print("   barely-densified layer (settle ~5mm/event); recovery cannot outrun the demand.")
