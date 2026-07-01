"""Figure 2 - the empirical evidence (3 panels): recovery-screen quantification,
field-scale odds-ratio test, and the identifiability limit."""
import numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
R39 = Path(__file__).resolve().parents[1]
OUT = R39 / "figures"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 8.5,
    "axes.titleweight": "bold",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})
C = dict(dem="#b2182b", rec="#1f5fa6", ok="#1a8a43", gray="#555555", ink="#1d1d1d")

fig = plt.figure(figsize=(7.3, 2.7), dpi=300)

# ----------------------------------------------- (a) recovery screen (PGA space)
ax = fig.add_axes([0.060, 0.21, 0.205, 0.63])
t = np.linspace(0, 24, 100)
tol = 0.110 + (0.135 - 0.110) * (t / 12)
ax.axhspan(0.40, 0.80, color=C["dem"], alpha=0.10)
ax.axhline(0.6, color=C["dem"], ls="--", lw=1.4)
ax.text(12, 0.66, "Feb-2011 demand", color=C["dem"], fontsize=6.4, ha="center", weight="bold")
ax.plot(t, tol, color=C["rec"], lw=2.4)
ax.annotate("tolerable PGA\n" r"$\approx$ 15-25 milli-$g$/yr", xy=(17, tol[70]), xytext=(12, 0.30),
            fontsize=6.4, color=C["rec"], ha="center", weight="bold",
            arrowprops=dict(arrowstyle="->", color=C["rec"], lw=0.9))
ax.axhline(0.057, color=C["gray"], ls=":", lw=1)
ax.text(12, 0.082, "re-trigger threshold", color=C["gray"], fontsize=5.8, ha="center")
ax.set_xlim(0, 24); ax.set_ylim(0, 0.82)
ax.set_xlabel("months since event", fontsize=7.2)
ax.set_ylabel("PGA (g)", fontsize=7.2); ax.tick_params(labelsize=6.6)
ax.set_title("a  Recovery screen", loc="left")

# ----------------------------------------------- (b) field-scale odds-ratio forest plot
ax = fig.add_axes([0.420, 0.21, 0.300, 0.63])
rows = [("Demand  (all eras, $q_c$ ctrl)", 3.14, 1.99, 4.96, C["dem"], "845"),
        ("LSN (AUC 0.61)", 1.06, 1.03, 1.10, C["rec"], "845"),
        ("prior-manif.: pre, LSN", 1.18, 0.64, 2.18, C["gray"], "300"),
        ("prior-manif.: all, LSN", 0.50, 0.30, 0.82, C["gray"], "845"),
        ("prior-manif.: all, $q_c$", 1.24, 0.65, 2.38, C["gray"], "845"),
        ("prior-manif.: pre, $q_c$", 4.57, 1.97, 10.62, C["gray"], "300")]
ax.axvspan(0.7, 1.45, color="0.87", alpha=0.7, zorder=0)
labs = []
for i, (lab, o, lo, hi, c, n) in enumerate(rows):
    y = len(rows) - i; labs.append((y, lab))
    ax.plot([lo, hi], [y, y], color=c, lw=1.5, zorder=3); ax.plot(o, y, "o", color=c, ms=4.6, zorder=4)
    ax.text(20, y, f"n={n}", fontsize=5.4, color=C["gray"], va="center", ha="right")
ax.axvline(1, color="0.4", lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.25, 21); ax.set_ylim(0.4, len(rows) + 0.6)
ax.set_yticks([y for y, _ in labs]); ax.set_yticklabels([l for _, l in labs], fontsize=5.9)
ax.set_xlabel("odds ratio for Feb-2011 re-liquefaction", fontsize=7.2); ax.tick_params(axis="x", labelsize=6.6)
ax.set_title("b  Demand controls", loc="left")

# ----------------------------------------------- (c) identifiability
ax = fig.add_axes([0.790, 0.21, 0.198, 0.63])
ax.scatter([0.3, 0.33, 0.43, 0.43, 0.46, 0.46], [1] * 6, s=34, color=C["ok"], zorder=5,
           edgecolor="white", linewidth=0.5)
strad = [3.4, 3.6, 3.9, 4.8, 5.2, 6.4, 6.6, 7.2, 8.3, 23, 24.9]
ax.scatter(strad, [0] * len(strad), s=28, color=C["dem"], marker="x")
ax.set_yticks([0, 1]); ax.set_yticklabels(["straddle\nre-liq", "event-\nfree (n=6)"], fontsize=6.4)
ax.set_xlabel("repeat-segment months", fontsize=7.2); ax.set_xlim(-1, 27); ax.set_ylim(-0.7, 1.85)
ax.tick_params(axis="x", labelsize=6.6)
ax.set_title("c  No clean rate", loc="left")
ax.text(13, 1.52, "6 event-free spans\nare all 9$-$14 days", fontsize=6.0,
        ha="center", va="center", color=C["gray"], style="italic")

fig.savefig(OUT / "Fig2_evidence.png", dpi=600, facecolor="white")
fig.savefig(OUT / "Fig2_evidence.pdf", facecolor="white")
from PIL import Image
print("size px:", Image.open(OUT / "Fig2_evidence.png").size)

