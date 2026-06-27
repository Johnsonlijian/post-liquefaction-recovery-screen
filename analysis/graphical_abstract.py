"""Graphical abstract for the Engineering Geology submission (lead image).
Elsevier spec: >=1328x531 px, ratio 2.5:1, >=300 dpi, legible at 200 px high."""
import numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R39 = Path(__file__).resolve().parents[1]
OUT = R39 / "figures"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
C = {"dem": "#b2182b", "sus": "#1f5fa6", "ok": "#1a8a43", "ink": "#1d1d1d", "gray": "#6a6a6a"}

fig = plt.figure(figsize=(7.5, 3.0), dpi=300)   # 2250 x 900 px = 2.5:1
fig.patch.set_facecolor("white")

# -- left: mechanism chart ------------------------
ax = fig.add_axes([0.075, 0.21, 0.40, 0.55])
evx = [0, 5.6, 9.3, 15.6]; evnm = ["Sep'10", "Feb'11", "Jun'11", "Dec'11"]
for x in evx:
    ax.axvline(x, color="#c2c2c2", lw=0.8, zorder=1)
for x, nm in zip(evx, evnm):
    ax.text(x, 1.03, nm, fontsize=7.0, ha="center", va="bottom", color=C["gray"])
demand = 0.92
ax.axhline(demand, color=C["dem"], ls="-", lw=1.8, zorder=3)
ax.text(16.5, demand + 0.03, "seismic demand", color=C["dem"], fontsize=7.6, ha="right",
        va="bottom", weight="bold")
t = np.linspace(0, 16.3, 500); floor, rate = 0.40, 0.030
res = np.array([floor + rate * (x - max([e for e in evx if e <= x] + [0])) for x in t])
ax.plot(t, res, color=C["sus"], lw=2.8, zorder=4)
ax.fill_between(t, res, demand, color=C["dem"], alpha=0.10, zorder=2)
ax.text(8.0, 0.72, "gap never closes", color=C["dem"], fontsize=8.2, ha="center", style="italic")
ax.annotate("recovery\n+25 milli-g/yr", xy=(8.0, floor + rate * (8.0 - 5.6)), xytext=(12.3, 0.31),
            fontsize=7.0, color=C["sus"], ha="center", va="center", weight="bold",
            arrowprops=dict(arrowstyle="->", color=C["sus"], lw=1.1))
ax.set_xlim(-1, 17.4); ax.set_ylim(0.18, 1.15); ax.set_yticks([]); ax.set_xticks([])
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
ax.set_xlabel("time into earthquake sequence", fontsize=7.6, color=C["ink"])
ax.set_ylabel("re-liquefaction\nresistance reserve", fontsize=7.6, color=C["ink"])
ax.set_title("Reserve rebuilds far too slowly",
             fontsize=8.8, weight="bold", color=C["ink"], loc="left", pad=12)

# -- divider -------------------------------
fig.add_artist(plt.Line2D([0.505, 0.505], [0.13, 0.87], color="#dcdcdc", lw=1.1,
                          transform=fig.transFigure))

# -- right: takeaways ---------------------------
ax2 = fig.add_axes([0.53, 0.0, 0.45, 1.0]); ax2.axis("off")
ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)
ax2.text(0.0, 0.93, "Does liquefied ground recover\nbefore the next earthquake-",
         fontsize=11.0, weight="bold", color=C["ink"], va="top", linespacing=1.25)
ax2.text(0.0, 0.635, "No: far too slow to matter", fontsize=13.5, weight="bold",
         color=C["dem"], va="top")
bullets = [("Recovery adds only ~25 milli-g per year", C["sus"]),
           ("No resolvable susceptibility change", C["ink"]),
           ("Re-liquefaction driven by seismic demand", C["ok"])]
y = 0.45
for txt, col in bullets:
    ax2.plot(0.022, y - 0.01, "s", color=col, ms=6)
    ax2.text(0.075, y, txt, fontsize=9.0, color=C["ink"], va="top")
    y -= 0.125
ax2.text(0.0, 0.05, "Open CPT pipeline  |  845 sites  |  2010-2011 Canterbury, NZ",
         fontsize=7.2, color=C["gray"], va="bottom", style="italic")

fig.savefig(OUT / "GraphicalAbstract.png", dpi=300, facecolor="white")
fig.savefig(OUT / "GraphicalAbstract.pdf", facecolor="white")
from PIL import Image
im = Image.open(OUT / "GraphicalAbstract.png")
print("size px:", im.size, "ratio %.3f" % (im.size[0] / im.size[1]), "(need >=1328x531, ratio 2.5)")
