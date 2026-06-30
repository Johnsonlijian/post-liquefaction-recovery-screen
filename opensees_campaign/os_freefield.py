# Faithful reproduction of OpenSeesPy freeFieldEffective (McGann/Arduino/UW),
# fixes: openseespy import, np.float->float, np.int->int. Synthetic motion to validate.
import numpy as np
import math as mm
import openseespy.opensees as op
import time as tt

# ---- generate a synthetic rock-outcrop VELOCITY motion (validation) ----
motionDT = 0.005
motionSteps = 3000
rng = np.random.default_rng(11)
t = np.arange(motionSteps)*motionDT
env = (t/2.0)**2*np.exp(-(t-3)/4.0); env = np.clip(env, 0, None); env /= env.max()
sig = np.zeros(motionSteps)
for f in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0):
    sig += np.sin(2*np.pi*f*t + rng.uniform(0, 6.28))
vel = env*sig; vel = 0.35*vel/np.max(np.abs(vel))         # ~0.35 m/s peak (strong)
np.savetxt('velocityHistory.txt', vel)

op.wipe()
nodes_dict = dict()
soilThick = 30.0; numLayers = 3; layerThick = [20.0, 8.0, 2.0]; waterTable = 2.0
layerBound = np.zeros(numLayers); layerBound[0] = layerThick[0]
for i in range(1, numLayers): layerBound[i] = layerBound[i-1]+layerThick[i]
nElemX = 1; nNodeX = 2*nElemX+1; sElemX = 2.0
nElemY = [40, 16, 4]; nElemT = 60
sElemY = np.zeros(numLayers)
for i in range(numLayers): sElemY[i] = layerThick[i-1]/nElemY[i-1]
nNodeY = 2*nElemT+1; nNodeT = nNodeX*nNodeY

op.model('basic', '-ndm', 2, '-ndf', 3)
count = 1; layerNodeCount = 0
dry_Node = np.zeros(500); node_save = np.zeros(500)
for k in range(1, numLayers+1):
    for i in range(1, nNodeX+1, 2):
        bump = 1 if k == 1 else 0
        j_end = 2*nElemY[k-1]+bump
        for j in range(1, j_end+1, 2):
            xCoord = (i-1)*(sElemX/2); yctr = j+layerNodeCount
            yCoord = (yctr-1)*(float(sElemY[k-1]))/2
            nodeNum = i+((yctr-1)*nNodeX)
            op.node(nodeNum, xCoord, yCoord)
            nodes_dict[nodeNum] = (nodeNum, xCoord, yCoord); node_save[nodeNum] = int(nodeNum)
            waterHeight = soilThick-waterTable
            if yCoord >= waterHeight: dry_Node[count] = nodeNum; count += 1
    layerNodeCount = yctr+1
dryNode = np.trim_zeros(dry_Node)
Node_d = np.trim_zeros(np.unique(node_save)); np.savetxt('Node_record.txt', Node_d)
for i in range(count-1):
    op.fix(int(dryNode[i]), 0, 0, 1)
op.fix(1, 0, 1, 0); op.fix(3, 0, 1, 0)
for i in range(1, ((3*nNodeY)-2), 6): op.equalDOF(i, i+2, 1, 2)

op.model('basic', '-ndm', 2, '-ndf', 2)
xCoord = float(sElemX/2); layerNodeCount = 0
for k in range(1, numLayers+1):
    bump = 1 if k == 1 else 0
    j_end = 2*nElemY[k-1]+bump
    for j in range(1, j_end+1, 1):
        yctr = j+layerNodeCount; yCoord = (yctr-1)*(float(sElemY[k-1]))/2
        op.node((3*yctr)-1, xCoord, yCoord)
    layerNodeCount = yctr
layerNodeCount = 0
for k in range(1, numLayers+1):
    for j in range(1, nElemY[k-1]+1):
        yctr = j+layerNodeCount; yCoord = float(sElemY[k-1])*(yctr-0.5)
        nL = (6*yctr)-2; nR = nL+2
        op.node(nL, 0.0, yCoord); op.node(nR, sElemX, yCoord)
    layerNodeCount = yctr
op.fix(2, 0, 1)
for i in range(1, ((3*nNodeY)-6), 6):
    op.equalDOF(i, i+1, 1, 2); op.equalDOF(i+3, i+4, 1, 2); op.equalDOF(i+3, i+5, 1, 2)
op.equalDOF(nNodeT-2, nNodeT-1, 1, 2)

grade = 2.0; slope = mm.atan(grade/100.0); g = -9.81
xw = g*mm.sin(slope); yw = g*mm.cos(slope)
thick = [1.0, 1.0, 1.0]; xWgt = [xw]*3; yWgt = [yw]*3
uBulk = [6.88E6, 5.06E6, 5.0E-6]; hPerm = [1.0E-4]*3; vPerm = [1.0E-4]*3
op.nDMaterial('PressureDependMultiYield02', 3, 2, 1.8, 9.0e4, 2.2e5, 32, 0.1, 101.0, 0.5, 26, 0.067, 0.23, 0.06, 0.27, 20, 5.0, 3.0, 1.0, 0.0, 0.77, 0.9, 0.02, 0.7, 101.0)
op.nDMaterial('PressureDependMultiYield02', 2, 2, 2.24, 9.0e4, 2.2e5, 32, 0.1, 101.0, 0.5, 26, 0.067, 0.23, 0.06, 0.27, 20, 5.0, 3.0, 1.0, 0.0, 0.77, 0.9, 0.02, 0.7, 101.0)
op.nDMaterial('PressureDependMultiYield02', 1, 2, 2.45, 1.3e5, 2.6e5, 39, 0.1, 101.0, 0.5, 26, 0.010, 0.0, 0.35, 0.0, 20, 5.0, 3.0, 1.0, 0.0, 0.47, 0.9, 0.02, 0.7, 101.0)
for j in range(1, nElemT+1):
    nI = (6*j)-5; nJ = nI+2; nK = nI+8; nL = nI+6; nM = nI+1; nN = nI+5; nP = nI+7; nQ = nI+3; nR = nI+4
    lowerBound = 0.0
    for i in range(1, numLayers+1):
        if j*sElemY[i-1] <= layerBound[i-1] and j*sElemY[i-1] > lowerBound:
            op.element('9_4_QuadUP', j, nI, nJ, nK, nL, nM, nN, nP, nQ, nR, thick[i-1], i, uBulk[i-1], 1.0, 1.0, 1.0, xWgt[i-1], yWgt[i-1])
        lowerBound = layerBound[i-1]

dashF = nNodeT+1; dashS = nNodeT+2
op.node(dashF, 0.0, 0.0); op.node(dashS, 0.0, 0.0)
op.fix(dashF, 1, 1); op.fix(dashS, 0, 1); op.equalDOF(1, dashS, 1)
colArea = sElemX*thick[0]; rockVS = 700.0; rockDen = 2.5; dashpotCoeff = rockVS*rockDen
op.uniaxialMaterial('Viscous', numLayers+1, dashpotCoeff*colArea, 1)
op.element('zeroLength', nElemT+1, dashF, dashS, '-mat', numLayers+1, '-dir', 1)

damp = 0.02; omega1 = 2*np.pi*0.2; omega2 = 2*np.pi*20
a0 = 2*damp*omega1*omega2/(omega1+omega2); a1 = 2*damp/(omega1+omega2)
gamma, beta = 0.5, 0.25
# ---- gravity ----
for m in (1, 2, 3): op.updateMaterialStage('-material', m, '-stage', 0)
op.constraints('Penalty', 1.0E14, 1.0E14); op.test('NormDispIncr', 1e-4, 35, 0)
op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
op.integrator('Newmark', gamma, beta); op.analysis('Transient')
ok = op.analyze(10, 5.0E2); print('elastic gravity ok:', ok, flush=True)
for m in (1, 2, 3): op.updateMaterialStage('-material', m, '-stage', 1)
ok = op.analyze(40, 5.0e2); print('plastic gravity ok:', ok, flush=True)
# effective stress profile check (corner element bottom layer vs surface)
sv_base = op.eleResponse(1, 'material', 1, 'stress')[1]; sv_top = op.eleResponse(nElemT, 'material', 1, 'stress')[1]
print(f"sigma'_yy: base elem={-sv_base:.1f} kPa, top elem={-sv_top:.1f} kPa (base should be >> top)", flush=True)
# ---- switch permeability low -> dynamic ----
ctr = 10000.0
for i in range(1, nElemT+1):
    op.parameter(int(ctr+1.0), 'element', i, 'vPerm'); op.parameter(int(ctr+2.0), 'element', i, 'hPerm'); ctr += 2.0
ctr = 10000.0
for j in range(1, nElemT+1):
    lowerBound = 0.0
    for i in range(1, numLayers+1):
        if j*sElemY[i-1] <= layerBound[i-1] and j*sElemY[i-1] > lowerBound:
            op.updateParameter(int(ctr+1.0), vPerm[i-1]); op.updateParameter(int(ctr+2.0), hPerm[i-1])
        lowerBound = layerBound[i-1]
    ctr += 2.0
op.setTime(0.0); op.wipeAnalysis()
# ---- dynamic ----
op.model('basic', '-ndm', 2, '-ndf', 3)
cFactor = colArea*dashpotCoeff
op.timeSeries('Path', 2, '-dt', motionDT, '-filePath', 'velocityHistory.txt', '-factor', cFactor)
op.pattern('Plain', 10, 2); op.load(1, 1.0, 0.0, 0.0)
op.constraints('Penalty', 1.0E16, 1.0E16); op.test('NormDispIncr', 1e-3, 35, 0)
op.algorithm('KrylovNewton'); op.numberer('RCM'); op.system('ProfileSPD')
op.integrator('Newmark', gamma, beta); op.rayleigh(a0, a1, 0.0, 0.0); op.analysis('Transient')
# track pore pressure in the liquefiable layer (a mid node) during shaking
liq_node = int(6*30)   # a corner pore-pressure node ~ mid-depth
ru_peak = 0.0; svp0 = -op.eleResponse(30, 'material', 1, 'stress')[1]
startT = tt.time()
nS = motionSteps; dT = motionDT
for s in range(nS):
    if op.analyze(1, dT) != 0:
        op.algorithm('Newton')
        if op.analyze(1, dT) != 0: print('dyn stop @ step', s, flush=True); break
        op.algorithm('KrylovNewton')
    svp = -op.eleResponse(30, 'material', 1, 'stress')[1]
    ru_peak = max(ru_peak, 1 - svp/svp0)
print(f"DYNAMIC done: peak ru in liquefiable layer (elem30) = {ru_peak:.2f}  [{tt.time()-startT:.0f}s]", flush=True)
op.wipe()
