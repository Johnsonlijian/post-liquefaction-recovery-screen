"""Converged Fig 1: mechanism schematic + recovery screen + driver/null forest + identifiability."""
import numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R39 = Path(__file__).resolve().parents[1]
FIG = R39 / "figures"; FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 9.5, "axes.titleweight": "bold"})
C = {"dem": "#b2182b", "sus": "#2166ac", "ok": "#1a9850", "gray": "#555555"}
fig = plt.figure(figsize=(7.4, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.34, left=0.10, right=0.97, top=0.93, bottom=0.09)

# -------------------------------- (a) mechanism schematic
ax = fig.add_subplot(gs[0, 0])
evx = [0, 5.6, 9.3, 15.6]
evnm = ["Darfield\nSep 2010", "Christchurch\nFeb 2011", "Jun 2011", "Dec 2011"]
for x in evx:
    ax.axvline(x, color=C["gray"], lw=0.9, alpha=0.5, zorder=1)
for x, nm in zip(evx, evnm):
    ax.text(x, 1.045, nm, fontsize=5.7, ha="center", va="bottom", color=C["gray"])
demand = 0.92
ax.axhline(demand, color=C["dem"], ls="-", lw=1.3, zorder=3)
ax.text(16.9, demand + 0.015, "seismic demand", color=C["dem"], fontsize=6.0, ha="right", va="bottom")
t = np.linspace(0, 16.4, 600)
floor, rate = 0.42, 0.030
res = np.array([floor + rate * (x - max([e for e in evx if e <= x] + [0])) for x in t])
ax.plot(t, res, color=C["sus"], lw=1.8, zorder=4)
ax.fill_between(t, res, demand, color=C["dem"], alpha=0.06, zorder=2)
ax.text(11.7, 0.745, "gap never closes", color=C["dem"], fontsize=6.0, ha="center")
ax.annotate("H1: slow recovery\n(rebuild too slow)", xy=(8.0, floor + rate * (8.0 - 5.6)),
            xytext=(9.3, 0.255), fontsize=6.0, color=C["sus"], ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=C["sus"], lw=0.8))
ax.annotate("H2: prior-manifestation\nsignal? (paradox)", xy=(5.6, floor),
            xytext=(2.5, 0.66), fontsize=6.0, color=C["dem"], ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=C["dem"], lw=0.8))
ax.text(8.2, 0.10, "finding: recovery too slow + no stable\nprior-manifestation signal -> demand sets re-liquefaction",
        fontsize=6.0, ha="center", va="center", color=C["ok"], weight="bold",
        bbox=dict(boxstyle="round,pad=0.3", fc="#eaf6ea", ec=C["ok"], lw=0.8))
ax.set_xlim(-1, 17.4); ax.set_ylim(0.0, 1.16); ax.set_yticks([])
ax.set_xlabel("months into sequence"); ax.set_ylabel("resistance reserve (schematic)")
ax.set_title("a  Two hypotheses, one sequence", loc="left")

# -------------------------------- (b) recovery screen
ax = fig.add_subplot(gs[0, 1])
t = np.linspace(0, 24, 100)
tol = 0.110 + (0.135 - 0.110) * (t / 12)
ax.plot(t, tol, color=C["sus"], lw=2, label="tolerable next-event PGA")
ax.axhline(0.6, color=C["dem"], ls="-", lw=1.4)
ax.text(0.5, 0.62, "Feb-2011 demand (0.4-0.8 g)", color=C["dem"], fontsize=6.6, va="bottom")
ax.axhline(0.057, color=C["gray"], ls=":", lw=1)
ax.text(12, 0.062, "re-trigger threshold (Quigley)", color=C["gray"], fontsize=6.2, va="bottom")
ax.fill_between(t, tol, 0.6, color=C["dem"], alpha=0.07)
ax.set_xlim(0, 24); ax.set_ylim(0, 0.7)
ax.set_xlabel("months of recovery"); ax.set_ylabel("PGA (g)")
ax.set_title("b  Recovery screen: ~25 milli-g/yr", loc="left")
ax.text(13, 0.30, "deficit never\nclosed within\na sequence", color=C["dem"], fontsize=6.8, ha="center")

# -------------------------------- (c) forest plot
ax = fig.add_subplot(gs[1, 0])
rows = [("PGA  (robust driver)", 4.6, 3.0, 7.0, C["dem"]),
        ("LSN severity (weak, AUC 0.61)", 1.06, 1.03, 1.10, C["sus"]),
        ("prior-liq: pre-event, LSN", 1.18, 0.64, 2.18, C["gray"]),
        ("prior-liq: all, LSN", 0.50, 0.30, 0.82, C["gray"]),
        ("prior-liq: all, qc1Ncs", 1.24, 0.65, 2.38, C["gray"]),
        ("prior-liq: pre-event, qc1Ncs", 4.57, 1.97, 10.62, C["gray"])]
labs = []
for i, (lab, o, lo, hi, c) in enumerate(rows):
    y = len(rows) - i; labs.append((y, lab))
    ax.plot([lo, hi], [y, y], color=c, lw=1.4); ax.plot(o, y, "o", color=c, ms=5)
ax.axvline(1, color="0.6", lw=0.8, ls="-")
ax.set_xscale("log"); ax.set_xlim(0.25, 16); ax.set_ylim(0.4, len(rows) + 0.6)
ax.set_yticks([y for y, _ in labs]); ax.set_yticklabels([l for _, l in labs], fontsize=6.0)
ax.set_xlabel("odds ratio for Feb-2011 re-liquefaction")
ax.set_title("c  Demand drives it; prior-liq not resolvable", loc="left")

# -------------------------------- (d) identifiability
ax = fig.add_subplot(gs[1, 1])
ax.scatter([0.3, 0.33, 0.43, 0.43, 0.46, 0.46], [1] * 6, s=30, color=C["ok"], zorder=5,
           label="event-free (n=6, 9-14 d)")
strad_dt = [3.4, 3.6, 3.9, 4.8, 5.2, 6.4, 6.6, 7.2, 8.3, 23, 24.9]
ax.scatter(strad_dt, [0] * len(strad_dt), s=24, color=C["dem"], marker="x",
           label="straddle re-liq event(s)")
ax.set_yticks([0, 1]); ax.set_yticklabels(["straddle\nre-liq events", "event-\nfree"], fontsize=7)
ax.set_xlabel("repeat-segment duration (months)"); ax.set_xlim(-1, 27); ax.set_ylim(-0.7, 2.05)
ax.set_title("d  No clean recovery rate is identifiable", loc="left")
ax.text(13, 1.62, "a clean rate needs event-free, multi-month gaps:\nnone exist (the 6 event-free spans are all 9-14 d)",
        fontsize=5.7, ha="center", va="center", color=C["gray"])

fig.savefig(FIG / "Fig1_converged.png", dpi=320, bbox_inches="tight")
fig.savefig(FIG / "Fig1_converged.pdf", bbox_inches="tight")
print("[saved]", FIG / "Fig1_converged.png", "and .pdf")
