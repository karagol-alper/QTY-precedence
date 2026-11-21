"""Check every numeric claim added to the manuscript against the analysis CSVs.

Each check names the sentence it guards; a FAIL means the text and the tables
disagree and the text is wrong.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pathlib import Path

import numpy as np
import pandas as pd

import core

TAB = core.EXTENDED
# The manuscript sources are optional: present in a working copy, absent in a
# data-and-code release. The numeric checks always run; the text checks are
# skipped with a notice when the manuscript is not alongside.
MS = Path(os.environ.get("QTY_MANUSCRIPT", core.ROOT / "manuscript"))
TEXT = ((MS / "proteins_submission.md").read_text(encoding="utf-8")
        if (MS / "proteins_submission.md").exists() else None)

fails = []


def check(name, cond, detail=""):
    print(f"{'ok  ' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        fails.append(name)


def approx(a, b, tol=0.005):
    return abs(float(a) - float(b)) <= tol


R = {n: pd.read_csv(f"{TAB}/{n}") for n in os.listdir(TAB) if n.endswith(".csv")}

# --- S7: the matched-background claims ------------------------------------
s7 = R["S7_background_definitions.csv"]
full = s7[s7.definition.str.startswith("alphabet-size matched, conservation")]
pub = s7[s7.definition.str.startswith("as published")]
fy = float(full[full.swap.str.startswith("F")].ratio.iloc[0])
vt = float(full[full.swap.str.startswith("V")].ratio.iloc[0])
fy0 = float(pub[pub.swap.str.startswith("F")].ratio.iloc[0])
vt0 = float(pub[pub.swap.str.startswith("V")].ratio.iloc[0])
check("F->Y 2.97 -> 3.31 under fullest matching",
      approx(fy0, 2.97) and approx(fy, 3.31), f"{fy0:.2f} -> {fy:.2f}")
check("V->T 1.48 -> 1.34 under fullest matching",
      approx(vt0, 1.48) and approx(vt, 1.34), f"{vt0:.2f} -> {vt:.2f}")
check("no null result significant under any definition",
      bool((s7[s7.swap.str.startswith(("L", "I"))].p_adj.astype(float)
            >= 0.05).all()))
check("V->T significant by three orders of magnitude under full matching",
      float(full[full.swap.str.startswith("V")].p_adj.iloc[0]) < 1e-3)

# --- S8b: conservation strata ---------------------------------------------
g = R["S8b_grade_strata.csv"]
g["p_adj"] = pd.to_numeric(g.p_adj, errors="coerce")
g9 = g[(g.swap.str.startswith("F")) & (g.stratum.str.startswith("9"))].iloc[0]
check("F->Y grade 9 ratio 6.43", approx(g9.ratio, 6.43), f"{g9.ratio}")
check("F->Y grade 9 observed 5 of 7, expected 0.78",
      int(g9.observed) == 5 and int(g9.n_sites) == 7 and approx(g9.expected, 0.78, 0.01))
fsig = (g[g.swap.str.startswith("F")].p_adj < 0.05).sum()
check("F->Y significant in three of four grade strata", fsig == 3, f"{fsig}")
vsig = g[(g.swap.str.startswith("V")) & (g.p_adj < 0.05)]
check("V->T significant in the conserved stratum alone",
      len(vsig) == 1 and vsig.iloc[0].stratum.startswith("7-8")
      and approx(vsig.iloc[0].ratio, 1.66, 0.01), f"{list(vsig.stratum)}")

# --- S12/S13: leave-one-out and depth -------------------------------------
l = R["S11_leave_one_segment_out.csv"]
check("22 segment drops; V->T and F->Y significant in all, nulls in none",
      bool((l.n_drops == 22).all())
      and set(l[l.swap.str.startswith(("V", "F"))].n_p_below_05) == {22}
      and set(l[l.swap.str.startswith(("L", "I"))].n_p_below_05) == {0})
d = R["S12_bilayer_depth.csv"].set_index(["swap", "zone"])
core, intf = "core (inner half)", "interfacial (outer half)"
vals = {(s[0], z): float(d.loc[(sw, z), "ratio"])
        for sw in d.index.get_level_values(0).unique() for s in [sw]
        for z in [core, intf]}
check("depth F->Y 3.75 core / 2.64 interface",
      approx(vals[("F", core)], 3.75) and approx(vals[("F", intf)], 2.64))
check("depth V->T 1.73 core / 1.46 interface",
      approx(vals[("V", core)], 1.73) and approx(vals[("V", intf)], 1.46))

# --- S10: ranks among all nineteen alternatives ---------------------------
a = R["S9_all_targets.csv"]
q = a[a.is_qty].set_index("native")
check("Y ranks 1 of 19 at F", int(q.loc["F", "rank_in_native"]) == 1)
check("T ranks 2 of 19 at V", int(q.loc["V", "rank_in_native"]) == 2)
check("T ranks 5 of 19 at I", int(q.loc["I", "rank_in_native"]) == 5)
check("Q ranks 12 of 19 at L", int(q.loc["L", "rank_in_native"]) == 12)
fr = a[a.native == "F"].sort_values("ratio", ascending=False)
enr = fr[(fr.p_adj < 0.05) & (fr.ratio > 1)]
check("W is second at F, and Y and W are the only enriched targets there",
      fr.iloc[1].target == "W" and set(enr.target) == {"Y", "W"},
      f"enriched: {sorted(enr.target)}")
dep = fr[(fr.p_adj < 0.05) & (fr.ratio < 1)]
check("G, A and S are significantly depleted at F",
      set(dep.target) == {"G", "A", "S"}, f"depleted: {sorted(dep.target)}")
lr = a[(a.native == "L")].sort_values("ratio", ascending=False)
check("four best targets at L are hydrophobic",
      set(lr.head(4).target) <= set("MIFVLWAC"), f"{list(lr.head(4).target)}")

# --- S11: matrix concordance ----------------------------------------------
c = R["S10a_matrix_concordance.csv"].set_index("matrix")
check("Spearman 0.83 BLOSUM62 / 0.88 GONNET1992 / 0.81 JTT",
      approx(c.loc["BLOSUM62", "spearman_rho"], 0.83)
      and approx(c.loc["GONNET1992", "spearman_rho"], 0.88)
      and approx(c.loc["JONES", "spearman_rho"], 0.81))
pc = R["S10b_qty_vs_matrices.csv"].set_index("swap")
gaps = {k[0]: int(v) for k, v in pc["gap"].items()}
check("L->Q 11 points below, V->T 17 above",
      gaps["L"] == -11 and gaps["V"] == 17, f"{gaps}")

# --- S15: reciprocal -------------------------------------------------------
rc = R["S14_reciprocal.csv"]
yf = rc[rc.forward.str.startswith("F")].iloc[0]
ql = rc[rc.forward.str.startswith("L")].iloc[0]
check("Y->F 1.86, adjusted p 0.034",
      approx(yf.reciprocal_ratio, 1.86) and approx(float(yf.p_adj), 0.034, 0.001))
check("Q->L rests on 4 sites, minimum detectable 1.70",
      int(ql.n_sites) == 4 and approx(ql.min_detectable, 1.70))

# --- the manuscript text itself -------------------------------------------
if TEXT is None:
    print("\nmanuscript not found at", MS, "- skipping the text checks")
    print("\n" + ("ALL NUMERIC CHECKS PASSED" if not fails
                   else f"{len(fails)} FAILED: {fails}"))
    sys.exit(1 if fails else 0)

for frag in ["2.97 to 3.31", "1.48 to 1.34", "6.43", "five of seven",
             "−0.78 and −0.80", "3.75 against 2.64",
             "1.73 against 1.46", "twenty-two", "ranks twelfth",
             "0.83 against BLOSUM62", "0.88 against GONNET1992",
             "0.81 against JTT", "eleven percentile points",
             "seventeen points above", "1.86, adjusted *p* = 0.034",
             "below 1.70", "Table S7", "Table S9", "Table S10", "Table S11",
             "Table S12", "Table S13", "Table S14", "Table S15",
             "Figure S1a", "Figure S1b", "Figure S1c"]:
    check(f"text contains {frag!r}", frag in TEXT)

# every supplementary table referenced must exist in the supplementary file
sup = (MS / "supplementary.md").read_text(encoding="utf-8")
refs = sorted({int(m) for m in re.findall(r"Table S(\d+)", TEXT)})
have = sorted({int(m) for m in re.findall(r"^## Table S(\d+)\.", sup, re.M)})
al = R["S8a_alphabet_strata.csv"]
al["p_adj"] = pd.to_numeric(al.p_adj, errors="coerce")
nfy = int((al[al.swap.str.startswith("F")].p_adj < 0.05).sum())
check("F->Y enriched in three of five alphabet-size bands", nfy == 3, f"{nfy}")
check("text contains 'three of five bands'", "three of five bands" in TEXT)

check("every supplementary table is referenced at least once",
      set(have) == set(refs), f"referenced {refs}, present {have}")
check("supplementary tables numbered without gaps",
      have == list(range(1, len(have) + 1)), f"{have}")

print("\n" + ("ALL CHECKS PASSED" if not fails
               else f"{len(fails)} FAILED: {fails}"))
sys.exit(1 if fails else 0)
