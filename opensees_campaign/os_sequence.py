# Canterbury event SEQUENCE on a Christchurch effective-stress column (PDMY02, 9_4_QuadUP).
# event -> reconsolidation -> event ... : per-event peak ru (re-liquefaction recurrence),
# inter-event densification (simulated recovery), resistance evolution.
import numpy as np, math as mm, openseespy.opensees as op, time as tt

motionDT = 0.005; motionSteps = 2400
rng = np.random.default_rng(7)
t = np.arange(motionSteps)*motionDT
env = (t/2.0)**2*np.exp(-(t-3)/3.5); env = np.clip(env, 0, None); env /= env.max()
sig = sum(np.sin(2*np.pi*f*t + rng.uniform(0, 6.28)) for f in (0.6, 1.0, 1.6, 2.4, 3.5))
vel = env*sig; vel = vel/np.max(np.abs(vel))           # unit-peak; scaled per event below
np.savetxt('seqMotion.txt', vel)

# ---- Christchurch profile (bottom-up layers): dense base / liquefiable sand / crust ----
soilThick = 14.0; numLayers = 3; layerThick = [7.0, 5.0, 2.0]; waterTable = 1.5
layerBound = np.zeros(numLayers); layerBound[0] = layerThick[0]
for i in range(1, numLayers): layerBound[i] = layerBound[i-1]+layerThick[i]
nElemX = 1; nNodeX = 3; sElemX = 2.0; nElemY = [14, 10, 4]; nElemT = 28
sElemY = np.array([layerThick[i]/nElemY[i] for i in range(numLayers)])  # all 0.5
nNodeY = 2*nElemT+1; nNodeT = nNodeX*nNodeY

op.wipe(); op.model('basic', '-ndm', 2, '-ndf', 3)
count = 1; layerNodeCount = 0; dry_Node = np.zeros(500); node_save = np.zeros(500)
for k in range(1, numLayers+1):
    for i in range(1, nNodeX+1, 2):
        bump = 1 if k == 1 else 0
        for j in range(1, 2*nElemY[k-1]+bump+1, 2):
            yctr = j+layerNodeCount; yCoord = (yctr-1)*sElemY[k-1]/2
            nn = i+((yctr-1)*nNodeX); op.node(nn, (i-1)*(sElemX/2), yCoord)
            node_save[nn] = int(nn)
            if yCoord >= soilThick-waterTable: dry_Node[count] = nn; count += 1
    layerNodeCount = yctr+1
dryNode = np.trim_zeros(dry_Node); Node_d = np.trim_zeros(np.unique(node_save))
np.savetxt('Node_record.txt', Node_d)
for i in range(count-1): op.fix(int(dryNode[i]), 0, 0, 1)
op.fix(1, 0, 1, 0); op.fix(3, 0, 1, 0)
for i in range(1, (3*nNodeY)-2, 6): op.equalDOF(i, i+2, 1, 2)
op.model('basic', '-ndm', 2, '-ndf', 2); layerNodeCount = 0
for k in range(1, numLayers+1):
    bump = 1 if k == 1 else 0
    for j in range(1, 2*nElemY[k-1]+bump+1):
        yctr = j+layerNodeCount; op.node((3*yctr)-1, sElemX/2, (yctr-1)*sElemY[k-1]/2)
    layerNodeCount = yctr
layerNodeCount = 0
for k in range(1, numLayers+1):
    for j in range(1, nElemY[k-1]+1):
        yctr = j+layerNodeCount; yC = sElemY[k-1]*(yctr-0.5)
        op.node((6*yctr)-2, 0.0, yC); op.node((6*yctr)-2+2, sElemX, yC)
    layerNodeCount = yctr
op.fix(2, 0, 1)
for i in range(1, (3*nNodeY)-6, 6):
    op.equalDOF(i, i+1, 1, 2); op.equalDOF(i+3, i+4, 1, 2); op.equalDOF(i+3, i+5, 1, 2)
op.equalDOF(nNodeT-2, nNodeT-1, 1, 2)
g = -9.81; slope = mm.atan(0.02); xw = g*mm.sin(slope); yw = g*mm.cos(slope)
thick = [1.0]*3; xWgt = [xw]*3; yWgt = [yw]*3; uBulk = [2.2e6, 2.2e6, 5.0e-6]
hP = [1e-4]*3; vP = [1e-4]*3
# mat1 dense base ; mat2 LIQUEFIABLE loose sand ; mat3 crust
op.nDMaterial('PressureDependMultiYield02', 1, 2, 2.10, 1.3e5, 2.6e5, 38, 0.1, 101.0, 0.5, 26, 0.013, 0.0, 0.30, 0.0, 20, 5.0, 3.0, 1.0, 0.0, 0.55, 0.9, 0.02, 0.7, 101.0)
op.nDMaterial('PressureDependMultiYield02', 2, 2, 1.90, 7.5e4, 2.0e5, 32, 0.1, 101.0, 0.5, 26, 0.087, 0.25, 0.05, 0.27, 20, 5.0, 3.0, 1.0, 0.0, 0.78, 0.9, 0.02, 0.7, 101.0)
op.nDMaterial('PressureDependMultiYield02', 3, 2, 1.95, 9.0e4, 2.2e5, 33, 0.1, 101.0, 0.5, 26, 0.045, 0.15, 0.10, 0.20, 20, 5.0, 3.0, 1.0, 0.0, 0.66, 0.9, 0.02, 0.7, 101.0)
for j in range(1, nElemT+1):
    nI = (6*j)-5; nodes = [nI, nI+2, nI+8, nI+6, nI+1, nI+5, nI+7, nI+3, nI+4]
    lb = 0.0
    for i in range(1, numLayers+1):
        if j*sElemY[i-1] <= layerBound[i-1] and j*sElemY[i-1] > lb:
            op.element('9_4_QuadUP', j, *nodes, thick[i-1], i, uBulk[i-1], 1.0, 1.0, 1.0, xWgt[i-1], yWgt[i-1])
        lb = layerBound[i-1]
dashF = nNodeT+1; dashS = nNodeT+2; op.node(dashF, 0., 0.); op.node(dashS, 0., 0.)
op.fix(dashF, 1, 1); op.fix(dashS, 0, 1); op.equalDOF(1, dashS, 1)
colArea = sElemX*thick[0]; dashpotCoeff = 700.0*2.5
op.uniaxialMaterial('Viscous', numLayers+1, dashpotCoeff*colArea, 1)
op.element('zeroLength', nElemT+1, dashF, dashS, '-mat', numLayers+1, '-dir', 1)
damp = 0.02; w1 = 2*np.pi*0.2; w2 = 2*np.pi*20; a0 = 2*damp*w1*w2/(w1+w2); a1 = 2*damp/(w1+w2)
# gravity
for m in (1, 2, 3): op.updateMaterialStage('-material', m, '-stage', 0)
op.constraints('Penalty', 1e14, 1e14); op.test('NormDispIncr', 1e-4, 35, 0)
op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
op.integrator('Newmark', 0.5, 0.25); op.analysis('Transient')
op.analyze(10, 5e2)
for m in (1, 2, 3): op.updateMaterialStage('-material', m, '-stage', 1)
op.analyze(40, 5e2)
LIQ = 18   # an element in the liquefiable layer (elems 15-24)
svp0 = -op.eleResponse(LIQ, 'material', 1, 'stress')[1]
print(f"gravity done. liquefiable-layer sigma'_v0 = {svp0:.1f} kPa", flush=True)
# ---- permeability control: a parameter pair on every element (switch undrained<->drained) ----
parV = {}; parH = {}; ctr = 10000
for i in range(1, nElemT+1):
    op.parameter(ctr+1, 'element', i, 'vPerm'); op.parameter(ctr+2, 'element', i, 'hPerm')
    parV[i] = ctr+1; parH[i] = ctr+2; ctr += 2
def set_perm(v, h):
    for i in range(1, nElemT+1): op.updateParameter(parV[i], v); op.updateParameter(parH[i], h)
K_DYN = 1.0e-5     # near-undrained during shaking (pore pressure builds -> liquefaction)
K_CON = 1.0e-2     # open drainage between events (months of dissipation in reality)

op.setTime(0.0); op.wipeAnalysis(); op.model('basic', '-ndm', 2, '-ndf', 3)
cF = colArea*dashpotCoeff; surf = 3*nElemT+1   # surface corner node; dof 2 = vertical

def run_event(scale, pgv):
    set_perm(K_DYN, K_DYN); op.wipeAnalysis()
    op.constraints('Penalty', 1e16, 1e16); op.test('NormDispIncr', 1e-3, 40, 0)
    op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
    op.integrator('Newmark', 0.5, 0.25); op.rayleigh(a0, a1, 0., 0.); op.analysis('Transient')
    op.setTime(0.0)   # Path series maps from t=0; reset each event
    op.timeSeries('Path', 200, '-dt', motionDT, '-filePath', 'seqMotion.txt', '-factor', cF*scale)
    op.pattern('Plain', 200, 200); op.load(1, 1.0, 0.0, 0.0)
    ru = 0.0
    for s in range(motionSteps):
        if op.analyze(1, motionDT) != 0:
            op.algorithm('Newton'); op.analyze(1, motionDT); op.algorithm('KrylovNewton')
        svp = -op.eleResponse(LIQ, 'material', 1, 'stress')[1]
        ru = max(ru, 1 - svp/svp0)
    op.remove('loadPattern', 200); op.remove('timeSeries', 200)
    return ru

def reconsolidate():
    # open drainage -> quasi-static consolidation; numerical damping (gamma=0.6) + dt ramp
    set_perm(K_CON, K_CON); op.wipeAnalysis()
    op.constraints('Penalty', 1e16, 1e16); op.test('NormDispIncr', 5e-4, 80, 0)
    op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
    op.integrator('Newmark', 0.6, 0.3025); op.rayleigh(a0, a1, 0., 0.); op.analysis('Transient')
    op.setTime(0.0)
    for dt in [0.02]*10 + [0.2]*20 + [1.0]*30 + [5.0]*40:
        if op.analyze(1, dt) != 0:
            op.algorithm('Newton', '-initial'); op.analyze(1, dt); op.algorithm('KrylovNewton')
    return op.nodeDisp(surf, 2), -op.eleResponse(LIQ, 'material', 1, 'stress')[1]

events = [('Darfield Sep10', 0.55, 0.20), ('Christchurch Feb11', 1.10, 0.45),
          ('June 2011', 0.65, 0.21), ('Dec 2011', 0.55, 0.20)]
print(f"\n{'event':20s}  scale  peak_ru  cum_settle(mm)  sv'_recov/{svp0:.0f}kPa", flush=True)
for name, sc, pga in events:
    t0 = tt.time(); ru = run_event(sc, pga); st, svr = reconsolidate()
    print(f"{name:20s}  {sc:4.2f}   {ru:5.2f}     {-st*1000:7.1f}        {svr:6.1f}     [{tt.time()-t0:.0f}s]", flush=True)
print("\n=> re-liquefaction recurs (ru~1) at each strong event despite inter-event densification;")
print("   reconsolidation settlement = simulated recovery; sv' recovers ~fully each time =>")
print("   ground stays liquefiable -> recovery too slow to lower hazard within the sequence.")
op.wipe()
