"""Figure 2: what the alignment columns at EAA1 leucine positions contain.

Ranks all twenty residues by the share of EAA1/SLC1A3 leucine positions whose
ConSurf alignment column already contains that residue, and marks where
glutamine — the residue the QTY code installs at every one of them — falls.

Usage: python3 make_fig2.py <out.png|out.eps>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import core

OUT = sys.argv[1] if len(sys.argv) > 1 else "figure2.png"
PROTEIN = "EAA1 / SLC1A3"
NATIVE, TARGET = "L", "Q"

C_HI, C_NEU, C_INK, C_GREY = "#2a78d6", "#c6c5be", "#2b2a26", "#8a8984"
mpl.rcParams.update({"figure.dpi": 130, "savefig.dpi": 400, "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "font.family": "DejaVu Sans"})

per = pd.read_csv(f"{core.RUN}/tables/per_site.csv")
sites = per[(per["protein"] == PROTEIN) & (per["aa"] == NATIVE)]
alphabets = sites["variety"].fillna("").apply(
    lambda s: {x for x in str(s).split("|") if x and x != "-"})
n = len(alphabets)
if n == 0:
    sys.exit(f"no {NATIVE} positions found for {PROTEIN}")

pct = {a: 100 * sum(a in v for v in alphabets) / n for a in core.AA20}
order = sorted(pct, key=lambda a: -pct[a])
rank = order.index(TARGET) + 1

fig, ax = plt.subplots(figsize=(8.2, 6.0))
ys = range(len(order))
for y, a in zip(ys, order):
    hi = a == TARGET
    ax.barh(y, pct[a], color=C_HI if hi else C_NEU, height=0.74, zorder=2)
    ax.annotate(f"{pct[a]:.0f}%", (pct[a], y), xytext=(6, 0),
                textcoords="offset points", va="center",
                fontsize=11.5 if hi else 11,
                color=C_INK if hi else C_GREY,
                fontweight="bold" if hi else "normal")

ax.set_yticks(list(ys))
ax.set_yticklabels(order, fontsize=13)
ax.invert_yaxis()
ax.set_xlim(0, 108)
ax.set_xlabel(f"% of the {n} leucine positions whose alignment column "
              "contains the residue", labelpad=8)
ax.tick_params(axis="y", length=0)

y_q = order.index(TARGET)
ax.annotate(f"glutamine: rank {rank} of {len(order)}",
            xy=(pct[TARGET] + 1.0, y_q + 0.42), xytext=(44, y_q + 2.9),
            fontsize=12, color=C_HI, va="center",
            arrowprops=dict(arrowstyle="->", color=C_HI, lw=1.6,
                            connectionstyle="arc3,rad=-0.35"))

fig.savefig(OUT)
print(f"wrote {OUT}   ({NATIVE}→{TARGET} rank {rank} of {len(order)}, "
      f"{pct[TARGET]:.0f}% of {n} positions)")
