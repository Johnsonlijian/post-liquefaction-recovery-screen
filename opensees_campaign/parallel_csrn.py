"""Parallel CSR-N liquefaction-resistance sweep (PM4Sand undrained cyclic DSS).
Fans out one OpenSees process per (Dr, CSR) case across the 96-core box.
Validates the engine (N must rise as CSR falls) and gives the CRR-N input to the screen."""
import multiprocessing as mp, itertools, time, numpy as np
from pathlib import Path
import pandas as pd

def worker(args):
    Dr, CSR = args
    import importlib, pm4_dss; importlib.reload(pm4_dss)
    t0 = time.time()
    try:
        r = pm4_dss.cyclic_dss(Dr=Dr, CSR=CSR, sig_v0_kpa=100.0, Cycle_max=30, strain_in=1.5e-5)
        # N to ru=0.95 from the per-step history
        ru = r['ru']; N = None
        if ru.size and ru.max() > 0.95:
            N = r['ncyc'] * (int(np.argmax(ru > 0.95)) / max(ru.size, 1))
        return dict(Dr=Dr, CSR=CSR, ru_max=round(float(ru.max()) if ru.size else 0, 3),
                    N_liq=round(N, 2) if N else None, half_cycles=r['ncyc'],
                    sec=round(time.time()-t0, 1))
    except Exception as e:
        return dict(Dr=Dr, CSR=CSR, ru_max=np.nan, N_liq=None, err=str(e)[:60],
                    sec=round(time.time()-t0, 1))

if __name__ == "__main__":
    Drs = [0.45, 0.55]
    CSRs = [0.10, 0.13, 0.16, 0.20, 0.25, 0.30]
    grid = list(itertools.product(Drs, CSRs))
    print(f"launching {len(grid)} OpenSees cases on {min(len(grid), 12)} parallel processes...")
    t0 = time.time()
    with mp.Pool(processes=min(len(grid), 12), maxtasksperchild=1) as pool:
        res = pool.map(worker, grid)
    df = pd.DataFrame(res).sort_values(['Dr', 'CSR'])
    df.to_csv(Path(__file__).resolve().parent / "csrn_results.csv", index=False)
    print(df.to_string(index=False))
    print(f"\ntotal wall time {time.time()-t0:.0f}s for {len(grid)} nonlinear effective-stress simulations")
    # quick validation: N should increase as CSR decreases
    for Dr in Drs:
        sub = df[(df.Dr == Dr) & df.N_liq.notna()].sort_values('CSR')
        if len(sub) >= 2:
            mono = all(np.diff(sub['N_liq'].values) <= 0)
            print(f"  Dr={Dr}: N(CSR) monotonic-decreasing = {mono}  | {list(zip(sub.CSR, sub.N_liq))}")
