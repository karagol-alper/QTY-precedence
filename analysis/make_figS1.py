"""Supplementary Figure S1, drawn from the extended-analysis tables.

Usage: python3 make_figS1.py <out.png|out.eps>

Colour carries significance and marker shape carries the series, so neither
has to be inferred from position. A ratio of exactly zero cannot be drawn on
a log axis, so it is clamped to the floor of the scale and marked with a
downward triangle rather than a round point.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import core
from matplotlib.ticker import FixedLocator, NullFormatter

OUT = sys.argv[1]
TAB = core.EXTENDED

C_ENR, C_DEP, C_NEU = "#2a78d6", "#e34948", "#c6c5be"
C_INK, C_RULE, C_EDGE = "#2b2a26", "#8a8984", "#9a9992"
mpl.rcParams.update({"figure.dpi": 130, "savefig.dpi": 400, "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "font.family": "DejaVu Sans"})

ORDER = ["L→Q", "I→T", "V→T", "F→Y"]
SHAPE = dict(zip(ORDER, ["o", "s", "^", "D"]))
FLOOR = 0.085


def colour(ratio, p):
    if not np.isfinite(p) or p >= 0.05:
        return C_NEU
    return C_DEP if ratio < 1 else C_ENR


def logaxis(ax, ticks, lo, hi):
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticklabels(["0" if t == FLOOR else
                        (f"{t:g}" if t < 1 else f"{t:g}") for t in ticks])


def point(ax, x, r, shape):
    """One estimate with its interval; returns nothing."""
    c = colour(r["ratio"], r["p"])
    lo = max(r["ci_lo"], FLOOR) if np.isfinite(r["ci_lo"]) else FLOOR
    hi = r["ci_hi"] if np.isfinite(r["ci_hi"]) else FLOOR
    if hi > lo:
        ax.plot([x, x], [lo, hi], "-", color=c, lw=1.8, alpha=.55,
                solid_capstyle="butt", zorder=2)
    zero = (not np.isfinite(r["ratio"])) or r["ratio"] <= 0
    ax.plot(x, FLOOR if zero else r["ratio"], "v" if zero else shape,
            ms=6.2 if zero else 5.6, color=c, zorder=3, mec="white", mew=.7)


fig = plt.figure(figsize=(7.2, 8.0))
gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.05, 1.2], hspace=0.55)

# ---- (a) the four background definitions ---------------------------------
ax = fig.add_subplot(gs[0])
s7 = pd.read_csv(f"{TAB}/S7_background_definitions.csv")
s7["p"] = s7["p"].astype(float)
defs = list(dict.fromkeys(s7["definition"]))
short = ["as published", "self-exclusive", "alphabet-size matched",
         "alphabet and conservation matched"]
dshape = ["o", "s", "^", "D"]
ax.axhline(1, color=C_RULE, ls="--", lw=1.0, zorder=1)
for j, d in enumerate(defs):
    sub = s7[s7["definition"] == d].set_index("swap").loc[ORDER]
    for i, (_, r) in enumerate(sub.iterrows()):
        point(ax, i + (j - 1.5) * 0.19, r, dshape[j])
logaxis(ax, [0.1, 0.25, 0.5, 1, 2, 4], 0.085, 6.0)
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(ORDER)
ax.set_xlim(-0.55, len(ORDER) - 0.45)
ax.set_ylabel("precedence ratio")
ax.set_title("a   Robustness to the background definition", loc="left",
             fontsize=9.5, fontweight="bold", pad=36)
ax.legend(handles=[Line2D([], [], marker=m, ls="", color=C_RULE, ms=5.4,
                          label=lab) for m, lab in zip(dshape, short)],
          loc="lower left", bbox_to_anchor=(0, 1.005), frameon=False,
          fontsize=7, handletextpad=.3, ncol=2, columnspacing=1.4)

# ---- (b) conservation grade strata ---------------------------------------
ax = fig.add_subplot(gs[1])
s8 = pd.read_csv(f"{TAB}/S8b_grade_strata.csv")
s8["p"] = pd.to_numeric(s8["p"], errors="coerce")
bands = ["1-3 variable", "4-6 intermediate", "7-8 conserved",
         "9 most conserved"]
ax.axhline(1, color=C_RULE, ls="--", lw=1.0, zorder=1)
for j, sw in enumerate(ORDER):
    sub = s8[s8["swap"] == sw].set_index("stratum")
    for i, b in enumerate(bands):
        if b not in sub.index:
            continue
        r = sub.loc[b]
        if not np.isfinite(r["ratio"]) and not np.isfinite(r["ci_hi"]):
            continue          # stratum with no site at all: nothing to show
        point(ax, i + (j - 1.5) * 0.19, r, SHAPE[sw])
logaxis(ax, [0.085, 0.5, 1, 2, 5, 10], 0.07, 14)
ax.set_xticks(range(len(bands)))
ax.set_xticklabels([b.replace(" ", "\n", 1) for b in bands], fontsize=7.4)
ax.set_xlim(-0.55, len(bands) - 0.45)
ax.set_ylabel("precedence ratio")
ax.set_xlabel("ConSurf conservation grade", labelpad=2)
ax.set_title("b   Precedence within conservation strata", loc="left",
             fontsize=9.5, fontweight="bold", pad=22)
ax.legend(handles=[Line2D([], [], marker=SHAPE[s], ls="", color=C_RULE,
                          ms=5.4, label=s) for s in ORDER]
          + [Line2D([], [], marker="o", ls="", color=C_ENR, ms=5.4,
                    label="p < 0.05"),
             Line2D([], [], marker="v", ls="", color=C_NEU, ms=5.8,
                    label="ratio of zero")],
          loc="lower left", bbox_to_anchor=(0, 1.005), frameon=False,
          fontsize=7, handletextpad=.3, ncol=6, columnspacing=.9)

# ---- (c) concordance with general exchangeability ------------------------
ax = fig.add_subplot(gs[2])
d = pd.read_csv(f"{TAB}/S10c_all_pairs_with_scores.csv")
corr = pd.read_csv(f"{TAB}/S10a_matrix_concordance.csv").set_index("matrix")
d["y"] = np.maximum(d["ratio"].astype(float), FLOOR)
d["zero"] = d["ratio"].astype(float) <= 0
rng = np.random.default_rng(0)
d["xj"] = d["BLOSUM62"] + (rng.random(len(d)) - 0.5) * 0.26
ax.axhline(1, color=C_RULE, ls="--", lw=1.0, zorder=1)
oth = d[~d.is_qty]
for zero, mark, size in [(False, "o", 20), (True, "v", 30)]:
    s = oth[oth.zero == zero]
    ax.scatter(s["xj"], s["y"], s=size, marker=mark, c=C_NEU, lw=.5,
               edgecolors=C_EDGE, zorder=2)
q = d[d.is_qty]
ax.scatter(q["BLOSUM62"], q["y"], s=70, marker="D",
           c=[colour(r, p) for r, p in zip(q["ratio"], q["p"])],
           edgecolors="white", lw=.9, zorder=4)
for _, r in q.iterrows():
    ax.annotate(r["swap"], (r["BLOSUM62"], r["y"]), xytext=(9, 2),
                textcoords="offset points", fontsize=7.6, color=C_INK,
                fontweight="bold")
logaxis(ax, [FLOOR, 0.25, 0.5, 1, 2, 4], 0.07, 5.0)
ax.set_xlabel("BLOSUM62 exchangeability score", labelpad=2)
ax.set_ylabel("precedence ratio")
ax.set_title("c   All 76 possible substitutions at L, I, V and F positions",
             loc="left", fontsize=9.5, fontweight="bold", pad=6)
ax.annotate(f"Spearman $\\rho$ = {corr.loc['BLOSUM62', 'spearman_rho']:.2f}"
            f"   (GONNET1992 {corr.loc['GONNET1992', 'spearman_rho']:.2f},"
            f" JTT {corr.loc['JONES', 'spearman_rho']:.2f})",
            (0.985, 0.05), xycoords="axes fraction", ha="right",
            fontsize=7.2, color=C_INK)
ax.legend(handles=[
    Line2D([], [], marker="o", ls="", color=C_NEU, mec=C_EDGE, ms=4.8,
           label="other targets (n = 72)"),
    Line2D([], [], marker="D", ls="", color=C_RULE, ms=5.6,
           label="QTY substitutions")],
    loc="upper left", frameon=False, fontsize=7, handletextpad=.3)

fig.savefig(OUT)
print("wrote", OUT)
