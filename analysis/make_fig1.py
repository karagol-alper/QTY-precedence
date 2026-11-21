"""Figure 1 for the manuscript, drawn from the run's own output tables.
Usage: python3 make_fig1.py <run_dir> <out.png>"""
import sys, pandas as pd, numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator
from matplotlib.lines import Line2D

RUN, OUT = sys.argv[1], sys.argv[2]
C_ENR, C_DEP, C_NEU = "#2a78d6", "#e34948", "#c6c5be"
mpl.rcParams.update({"figure.dpi": 130, "savefig.dpi": 400, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False, "savefig.bbox": "tight",
                     "font.family": "DejaVu Sans"})

pool = pd.read_csv(f"{RUN}/tables/pooled_by_swap.csv")
per  = pd.read_csv(f"{RUN}/tables/by_protein.csv")
pool["swap"] = pool.native + "→" + pool.target
per["swap"]  = per.native + "→" + per.target
POOLED = ["EAA1 / SLC1A3", "VPAT / SLC18B1"]
order = ["L→Q", "I→T", "V→T", "F→Y"]
pool = pool.set_index("swap").loc[order].reset_index()

n = len(pool)
fig, ax = plt.subplots(figsize=(7.0, 0.95 * n + 1.5))
ax.axvline(1, color="#8a8984", ls="--", lw=1.1, zorder=1)

for i, r in pool.iterrows():
    y = n - 1 - i
    sig = r.p < 0.05
    col = (C_DEP if r.ratio < 1 else C_ENR) if sig else C_NEU
    # 95% CI whisker
    ax.plot([r.ci_lo, r.ci_hi], [y, y], "-", color=col, lw=2.2, alpha=.55,
            solid_capstyle="butt", zorder=2)
    for x in (r.ci_lo, r.ci_hi):
        ax.plot([x, x], [y - .09, y + .09], "-", color=col, lw=2.2, alpha=.55, zorder=2)
    ax.plot(r.ratio, y, "D", ms=10, color=col, zorder=3)
    # per-protein estimates below
    vs = [float(per[(per.protein == p) & (per.swap == r.swap)].ratio.iloc[0])
          for p in POOLED if len(per[(per.protein == p) & (per.swap == r.swap)])]
    for v in vs:
        ax.plot(v, y - 0.27, "o", ms=4.5, color=C_NEU, mec="white", mew=.9, zorder=3)
    lab = f"{r.ratio:.2f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}]"
    ax.annotate(lab, (r.ci_hi, y), xytext=(7, 0), textcoords="offset points",
                va="center", ha="left", fontsize=9, fontweight="bold")

ax.set_yticks(range(n))
ax.set_yticklabels([f"{r.swap}\n$n$={int(r.n_sites)}" for _, r in pool[::-1].iterrows()],
                   fontsize=10)
ax.tick_params(axis="y", length=0)
lo, hi = 0.08, 7.2
ax.set_xscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(-0.75, n - 0.45)
ticks = [0.125, 0.25, 0.5, 1, 2, 4]
ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_minor_locator(NullLocator())
ax.set_xticklabels([f"{t:g}×" for t in ticks])
ax.set_xlabel("observed rate ÷ within-protein background rate   (log scale)", fontsize=9.5)
ax.text(1 / 1.09, n - 0.52, "depleted", ha="right", va="bottom", fontsize=8.5,
        style="italic", color="#8a8984")
ax.text(1 * 1.09, n - 0.52, "already favoured", ha="left", va="bottom", fontsize=8.5,
        style="italic", color="#8a8984")

ax.legend(handles=[
    Line2D([], [], marker="D", ls="-", ms=8, color=C_ENR, lw=2.2, alpha=.8,
           label="pooled estimate, 95% CI"),
    Line2D([], [], marker="D", ls="", ms=8, color=C_NEU, label="not significant"),
    Line2D([], [], marker="o", ls="", ms=4.5, color=C_NEU, label="each transporter")],
    loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, fontsize=8.5)
fig.savefig(OUT)
print("wrote", OUT)
for _, r in pool.iterrows():
    print(f"  {r.swap}  n={int(r.n_sites):>3}  {r.ratio:.2f} [{r.ci_lo:.2f}, {r.ci_hi:.2f}]  p={r.p:.1e}")
