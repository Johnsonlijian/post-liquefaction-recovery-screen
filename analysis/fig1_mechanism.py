"""Figure 1 - physically concrete mechanism figure for the EG paper.

Narrative left-to-right: liquefiable ground + CPT -> earthquake -> two
processes (H1 recovery / H2 susceptibility) -> demand-referenced recovery
screen (rebuildable vs required reserve; the gap never closes).
"""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image

R39 = Path(__file__).resolve().parents[1]
OUT = R39 / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C = dict(
    crust="#d8c7a2",
    sand="#eab94e",
    base="#a9a195",
    water="#2f7fb8",
    qc="#2b2b2b",
    dem="#b2182b",
    rec="#1f5fa6",
    ok="#1a8a43",
    ink="#1d1d1d",
    gray="#6b6b6b",
)

fig = plt.figure(figsize=(7.3, 4.4), dpi=300)
fig.patch.set_facecolor("white")
row_y, row_h = 0.360, 0.440

# Left: soil stratigraphy.
axs = fig.add_axes([0.060, row_y, 0.080, row_h])
for z0, z1, col, hatch in [
    (0.0, 1.2, C["crust"], ""),
    (1.2, 3.0, "#dcc98f", ""),
    (3.0, 6.0, C["sand"], "////"),
    (6.0, 8.0, C["base"], ""),
]:
    axs.add_patch(
        Rectangle((0, z0), 1, z1 - z0, facecolor=col, edgecolor="#7a7a7a", lw=0.6, hatch=hatch)
    )
axs.add_patch(Rectangle((0, 3.0), 1, 3.0, fill=False, edgecolor=C["dem"], lw=1.6))
axs.plot([0, 1], [1.2, 1.2], color=C["water"], lw=1.3)
axs.plot(0.5, 1.2, marker="v", color=C["water"], ms=7, mec="white", mew=0.5)
axs.set_xlim(0, 1)
axs.set_ylim(8, 0)
axs.set_xticks([])
axs.set_ylabel("depth (m)", fontsize=7.5)
axs.tick_params(labelsize=7)
axs.set_title("ground", fontsize=8.2, weight="bold", color=C["ink"], pad=3)

# CPT trace.
axq = fig.add_axes([0.160, row_y, 0.110, row_h])
z = np.linspace(0, 8, 400)
qc = (
    4.2 * np.exp(-((z - 0.6) / 0.8) ** 2)
    + 6.0 * np.exp(-((z - 2.1) / 0.9) ** 2)
    + 1.8 * (np.abs(z - 4.5) < 1.5)
    + 0.9
    + 12.0 / (1 + np.exp(-(z - 6.6) * 2.2))
)
qc = np.clip(qc + np.sin(z * 7.5) * 0.22, 0.4, None)
axq.axhspan(3.0, 6.0, color=C["sand"], alpha=0.25)
axq.plot(qc, z, color=C["qc"], lw=1.3)
axq.set_xlim(0, 16)
axq.set_ylim(8, 0)
axq.set_yticks([])
axq.set_xlabel(r"$q_c$ (MPa)", fontsize=7.5)
axq.tick_params(labelsize=7)
axq.set_title("CPT", fontsize=8.2, weight="bold", color=C["ink"], pad=3)
axq.annotate(
    "liquefiable\ncritical layer",
    xy=(2.2, 4.5),
    xytext=(7.0, 4.6),
    fontsize=6.5,
    color=C["dem"],
    va="center",
    weight="bold",
    arrowprops=dict(arrowstyle="->", color=C["dem"], lw=0.9),
)

# Earthquake load.
axE = fig.add_axes([0.055, 0.825, 0.230, 0.135])
axE.axis("off")
axE.set_xlim(0, 1)
axE.set_ylim(0, 1)
axE.plot([0.30, 0.40, 0.34, 0.46], [0.05, 0.55, 0.55, 1.0], color=C["dem"], lw=2.0, solid_capstyle="round")
for x0 in (0.16, 0.30, 0.44, 0.58):
    axE.annotate("", xy=(x0, 0.02), xytext=(x0, 0.42), arrowprops=dict(arrowstyle="->", color=C["dem"], lw=1.0))
axE.text(0.66, 0.45, "seismic demand\n(PGA)", fontsize=7.0, color=C["dem"], weight="bold", va="center")

# Middle: the two processes.
axM = fig.add_axes([0.300, row_y, 0.165, row_h])
axM.axis("off")
axM.set_xlim(0, 1)
axM.set_ylim(0, 1)
axM.annotate("", xy=(0.97, 0.60), xytext=(0.02, 0.60), arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=2.0))
axM.text(0.50, 0.68, "after liquefaction", fontsize=7.0, color=C["ink"], ha="center", style="italic")
axM.annotate("", xy=(0.64, 0.92), xytext=(0.44, 0.64), arrowprops=dict(arrowstyle="->", color=C["rec"], lw=1.6))
axM.text(0.50, 0.99, "H1  slow recovery\n(ageing, reconsolidation)", fontsize=6.4, color=C["rec"], ha="center", va="top", weight="bold")
axM.annotate("", xy=(0.64, 0.28), xytext=(0.44, 0.56), arrowprops=dict(arrowstyle="->", color=C["gray"], lw=1.6))
axM.text(0.50, 0.20, "H2  susceptibility\nchange?  (paradox,\nnot resolvable)", fontsize=6.4, color=C["gray"], ha="center", va="top")

# Right: demand-referenced recovery screen.
axr = fig.add_axes([0.520, row_y, 0.430, row_h])
months = np.linspace(0, 24, 200)
req, valid = 185.0, 160.0
rebuild = 85 + (95 - 85) * (months / 24)
axr.axhspan(valid, 205, color="#bdbdbd", alpha=0.32, hatch="xx", ec="#9a9a9a", lw=0)
axr.axhline(req, color=C["dem"], ls="--", lw=1.7)
axr.text(
    0.5,
    187,
    "required for FS = 1 at Feb-2011 demand",
    fontsize=6.7,
    color=C["dem"],
    va="bottom",
    weight="bold",
    bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.72),
)
axr.fill_between(months, rebuild, req, color=C["dem"], alpha=0.08)
axr.plot(months, rebuild, color=C["rec"], lw=2.6)
axr.annotate(
    "rebuildable reserve\n$\\approx$ 15-25 milli-$g$/yr",
    xy=(16.5, rebuild[137]),
    xytext=(9.6, 120),
    fontsize=6.6,
    color=C["rec"],
    weight="bold",
    ha="center",
    arrowprops=dict(arrowstyle="->", color=C["rec"], lw=0.9),
)
axr.text(18.5, 150, "deficit never\ncloses within\na sequence", fontsize=7.6, color=C["dem"], ha="center", va="center", weight="bold")
axr.axvline(4, color=C["ink"], ls=":", lw=1.1)
axr.text(4.4, 69, "next strong event", fontsize=6.4, color=C["ink"], va="center")
axr.set_xlim(0, 24)
axr.set_ylim(60, 205)
axr.set_xlabel("time since event (months)", fontsize=7.6)
axr.set_ylabel(r"clean-sand resistance  $q_{c1Ncs}$", fontsize=7.6)
axr.tick_params(labelsize=7)
for side in ["top", "right"]:
    axr.spines[side].set_visible(False)
axr.set_title("demand-referenced recovery screen", fontsize=8.6, weight="bold", color=C["ink"], pad=5)

# Bottom: Canterbury timeline.
axt = fig.add_axes([0.060, 0.165, 0.890, 0.085])
axt.axis("off")
axt.set_xlim(0, 1)
axt.set_ylim(0, 1)
axt.text(0.50, 0.97, "inter-event windows: weeks to months  (too short for recovery)", fontsize=6.5, ha="center", va="top", color=C["gray"], style="italic")
axt.annotate("", xy=(0.99, 0.60), xytext=(0.0, 0.60), arrowprops=dict(arrowstyle="-|>", color=C["gray"], lw=1.3))
for x, nm in [(0.04, "Darfield\nSep 2010"), (0.36, "Christchurch\nFeb 2011"), (0.55, "Jun 2011"), (0.92, "Dec 2011")]:
    axt.plot(x, 0.60, marker="o", ms=6.5, color=C["dem"], mec="white", mew=0.6, zorder=3)
    axt.text(x, 0.40, nm, fontsize=6.2, ha="center", va="top", color=C["gray"])

fig.add_artist(
    FancyBboxPatch(
        (0.060, 0.028),
        0.890,
        0.088,
        transform=fig.transFigure,
        boxstyle="round,pad=0.003",
        facecolor="#eaf6ea",
        edgecolor=C["ok"],
        lw=1.0,
    )
)
fig.text(0.505, 0.073, "The deficit never closes within a sequence, and no susceptibility change is resolvable", fontsize=7.6, ha="center", va="center", color=C["ok"], weight="bold")
fig.text(0.505, 0.045, "$\\rightarrow$  re-liquefaction is governed by seismic demand", fontsize=7.6, ha="center", va="center", color=C["ok"], weight="bold")

fig.savefig(OUT / "Fig1_mechanism.png", dpi=600, facecolor="white")
fig.savefig(OUT / "Fig1_mechanism.pdf", facecolor="white")
print("size px:", Image.open(OUT / "Fig1_mechanism.png").size)
