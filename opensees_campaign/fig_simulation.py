"""Effective-stress simulation figure for the EG manuscript.

Panels:
  A. Calibrated CSR-N engine versus Idriss-Boulanger style scaling.
  B. Effective-stress liquefaction cycle for one Canterbury event.
  C. Demand-recovery race showing why densification is too slow.
"""
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCR = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.0,
    "axes.linewidth": 0.8,
    "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C_LIQ = "#c0392b"
C_SIM = "#2c3e50"
C_DEM = "#e67e22"
C_REC = "#27ae60"
C_BASE = "#2980b9"

# Load data.
hi = pd.read_csv(SCR / "csrn_curves.csv")
lo_path = SCR / "csrn_low.csv"
lo = pd.read_csv(lo_path) if lo_path.exists() else pd.DataFrame(columns=hi.columns)
allc = pd.concat([lo, hi], ignore_index=True)


def crr15_of(df, dr):
    s = df[(df.Dr == dr) & df.N_liq.notna() & (df.N_liq > 1)].sort_values("CSR")
    if len(s) < 2:
        return None
    return float(np.interp(15.0, s.N_liq.values[::-1], s.CSR.values[::-1]))


Drs = sorted(allc.Dr.unique())
CRR = {dr: crr15_of(allc, dr) for dr in Drs}
CRR = {dr: value for dr, value in CRR.items() if value}
dr_arr = np.array(sorted(CRR))
crr_arr = np.array([CRR[dr] for dr in dr_arr])
b, a = np.polyfit(dr_arr, crr_arr, 1)  # CRR15 ~ a + b*Dr

with open(SCR / "events_independent.json", encoding="utf-8") as f:
    res = json.load(f)
ru_d = np.load(SCR / "hist_ru_Darfield.npy")
sv_d = np.load(SCR / "hist_sv_Darfield.npy")
tt = np.arange(len(ru_d)) * 0.02

fig = plt.figure(figsize=(7.3, 4.55))
gs = fig.add_gridspec(
    2,
    2,
    left=0.080,
    right=0.970,
    bottom=0.090,
    top=0.885,
    hspace=0.560,
    wspace=0.420,
    height_ratios=[1.0, 1.05],
)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, :])

# Panel A: CSR-N engine versus Idriss-Boulanger.
cmap = {0.30: "#9b59b6", 0.35: "#8e44ad", 0.40: "#5dade2", 0.45: "#2980b9", 0.55: "#16a085", 0.65: "#27ae60"}
Nline = np.linspace(1.5, 60, 100)
for dr in Drs:
    if dr not in CRR:
        continue
    s = allc[(allc.Dr == dr) & allc.N_liq.notna() & (allc.N_liq > 1)]
    col = cmap.get(dr, "#555555")
    axA.scatter(s.N_liq, s.CSR, s=24, color=col, zorder=3, edgecolor="w", linewidth=0.4)
    axA.plot(Nline, CRR[dr] * (15.0 / Nline) ** 0.34, color=col, lw=1.3, label=f"$D_r$={dr:.2f}")
axA.axvline(15, color="grey", ls=":", lw=0.9)
axA.text(15, 0.052, "$N$=15", rotation=90, va="bottom", ha="right", fontsize=7.5, color="grey")
axA.set_xscale("log")
axA.set_xlim(1.5, 62)
axA.set_ylim(0.05, 0.30)
axA.set_xlabel("Cycles to liquefaction, $N$")
axA.set_ylabel("Cyclic stress ratio, CSR")
axA.set_title("(A)  Calibrated effective-stress engine", fontsize=7.0, loc="left", fontweight="bold")
axA.legend(fontsize=4.8, frameon=False, ncol=2, loc="upper right", handlelength=1.3, columnspacing=0.7)
axA.text(
    0.03,
    0.04,
    "points: PM4Sand sims\nlines: $(15/N)^{0.34}$ (Idriss-Boulanger)",
    transform=axA.transAxes,
    fontsize=7,
    va="bottom",
    color="#444444",
)

# Panel B: liquefaction cycle.
axB.plot(tt, ru_d, color=C_LIQ, lw=1.4, label="$r_u$ (excess PWP ratio)")
axB.axhline(1.0, color=C_LIQ, ls=":", lw=0.8)
axB.set_ylim(0, 1.18)
axB.set_xlim(0, tt[-1])
axB.set_xlabel("Time (s)")
axB.set_ylabel(r"$r_u = 1-\sigma'/\sigma'_0$", color=C_LIQ)
axB.tick_params(axis="y", colors=C_LIQ)
axB.text(tt[-1] * 0.97, 1.03, r"liquefaction ($r_u\!\to\!1$)", ha="right", va="bottom", color=C_LIQ, fontsize=8)
axB2 = axB.twinx()
axB2.plot(tt, sv_d, color=C_BASE, lw=1.2, label=r"$\sigma'_v$")
axB2.set_ylabel("")
axB2.tick_params(axis="y", colors=C_BASE, right=False, labelright=False)
axB2.set_ylim(-3, 70)
axB2.spines["right"].set_visible(False)
axB.text(0.92, 0.08, r"$\sigma'_v$ (kPa)", transform=axB.transAxes, ha="right", va="bottom", color=C_BASE, fontsize=6.2)
axB.set_title("(B)  Effective-stress collapse under one event", fontsize=7.0, loc="left", fontweight="bold")
axB.annotate(
    "post-event reconsolidation:\n$\\sigma'_v$ recovers to ~94%, but\ndensification only ~0.1% / event",
    xy=(tt[-1] * 0.62, 0.06),
    xytext=(tt[-1] * 0.30, 0.30),
    fontsize=5.6,
    color="#333333",
    bbox=dict(boxstyle="round,pad=0.3", fc="#fdf6ec", ec=C_DEM, lw=0.8),
)

# Panel C: demand-recovery race.
drx = np.linspace(0.26, 0.72, 100)
axC.plot(drx, a + b * drx, color=C_SIM, lw=1.6, zorder=2, label="CRR$_{15}$($D_r$) (calibrated)")
axC.scatter(dr_arr, crr_arr, s=26, color=C_SIM, zorder=3, edgecolor="w", linewidth=0.4)
Dr0 = 0.30
CRR0 = a + b * Dr0
CSR_dem = 0.18
Dr_need = (CSR_dem - a) / b
axC.axhline(CSR_dem, color=C_DEM, ls="--", lw=1.3)
axC.text(0.555, CSR_dem + 0.004, "demand CSR ~= 0.18", ha="right", color=C_DEM, fontsize=7.8)
axC.scatter([Dr0], [CRR0], s=60, color=C_LIQ, zorder=5, edgecolor="w", linewidth=0.6)
axC.annotate(
    "in-situ $D_r$ ~= .30\nCRR$_{15}$ ~= %.2f" % CRR0,
    xy=(Dr0, CRR0),
    xytext=(0.305, 0.115),
    fontsize=7.5,
    color=C_LIQ,
    arrowprops=dict(arrowstyle="-", color=C_LIQ, lw=0.7),
)
dDr = 0.0051
xc = Dr0
for _ in range(7):
    axC.annotate("", xy=(xc + dDr, a + b * (xc + dDr)), xytext=(xc, a + b * xc), arrowprops=dict(arrowstyle="->", color=C_REC, lw=1.2))
    xc += dDr
axC.annotate(
    "densification recovery\n~0.5% $D_r$/event -> tens of events\n(~decades) to reach the demand",
    xy=(xc, a + b * xc),
    xytext=(0.405, 0.072),
    fontsize=7.4,
    color=C_REC,
    arrowprops=dict(arrowstyle="->", color=C_REC, lw=0.8),
)
axC.axvline(Dr_need, color=C_DEM, ls=":", lw=0.8)
axC.text(
    0.283,
    0.233,
    "recovery ~20x too slow\n(~27 yr vs 1.3-yr seq.)",
    ha="left",
    va="top",
    fontsize=7.6,
    color=C_LIQ,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", fc="#fdecea", ec=C_LIQ, lw=0.8),
)
axC.set_xlim(0.26, 0.72)
axC.set_ylim(0.04, 0.24)
axC.set_xlabel("Relative density, $D_r$")
axC.set_ylabel("Liquefaction resistance / demand, CRR$_{15}$")
axC.set_title("(C)  The demand-recovery race", fontsize=7.0, loc="left", fontweight="bold")

# Inset: every event liquefies.
axins = axC.inset_axes([0.64, 0.64, 0.32, 0.30])
names = [r["name"].split()[0] for r in res]
rus = [r["peak_ru"] for r in res]
axins.bar(range(len(rus)), rus, color=C_LIQ, width=0.62)
axins.axhline(1.0, color="grey", ls=":", lw=0.7)
axins.set_xticks(range(len(names)))
axins.set_xticklabels(["Dar", "Chc", "Jun", "Dec"], fontsize=6.2)
axins.set_ylim(0, 1.15)
axins.set_yticks([0, 1])
axins.tick_params(labelsize=6.2)
axins.set_title("peak $r_u$, every event", fontsize=6.6)

fig.suptitle(
    "Effective-stress simulation: the Canterbury layer re-liquefies under each event's demand;\n"
    "recovery (densification) is far too slow to intervene",
    fontsize=7.3,
    y=0.970,
    fontweight="bold",
)
for ext in ("png", "pdf"):
    fig.savefig(SCR / f"Fig_simulation.{ext}", dpi=600 if ext == "png" else 300)
print(
    "saved Fig_simulation.png/pdf | CRR15:",
    {round(k, 2): round(v, 3) for k, v in CRR.items()},
    "| fit a=%.3f b=%.3f | Dr_need=%.2f" % (a, b, Dr_need),
)
