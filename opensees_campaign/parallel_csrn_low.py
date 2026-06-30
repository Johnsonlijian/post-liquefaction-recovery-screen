"""Anchor the LOW-density resistance directly (no extrapolation): CSR-N at Dr=0.30/0.35/0.40."""
import multiprocessing as mp, itertools, time, numpy as np, pandas as pd
from pathlib import Path
def worker(args):
    Dr, CSR = args
    import pm4_dss, time
    t0 = time.time()
    try:
        r = pm4_dss.cyclic_dss(Dr=Dr, CSR=CSR, sig_v0_kpa=65.0, Cycle_max=60, strain_in=1.5e-5)
        return dict(Dr=Dr, CSR=CSR, N_liq=r['Nliq'], ru_max=round(r['ru_max'], 3), sec=round(time.time()-t0, 1))
    except Exception as e:
        return dict(Dr=Dr, CSR=CSR, N_liq=None, ru_max=np.nan, err=str(e)[:50], sec=round(time.time()-t0, 1))
if __name__ == "__main__":
    Drs = [0.30, 0.35, 0.40]; CSRs = [0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
    grid = list(itertools.product(Drs, CSRs))
    print(f"launching {len(grid)} low-Dr CSR-N cases (sig_v0=65 kPa)...", flush=True)
    t0 = time.time()
    with mp.Pool(processes=min(len(grid), 24), maxtasksperchild=1) as pool:
        res = pool.map(worker, grid)
    df = pd.DataFrame(res).sort_values(['Dr', 'CSR'])
    out = Path(__file__).resolve().parent / "csrn_low.csv"
    df.to_csv(out, index=False); print(df.to_string(index=False))
    for Dr in Drs:
        s = df[(df.Dr == Dr) & df.N_liq.notna() & (df.N_liq > 1)]
        if len(s) >= 2:
            crr15 = float(np.interp(15.0, s.N_liq.values[::-1], s.CSR.values[::-1]))
            print(f"  Dr={Dr}: CRR15 = {crr15:.3f}")
    print(f"{len(grid)} sims in {time.time()-t0:.0f}s")
