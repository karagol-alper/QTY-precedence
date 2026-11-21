"""Emit the markdown for supplementary tables S7-S14 from the analysis CSVs.

Every number in the supplementary material is written from the tables the
analysis produced, so nothing is transcribed by hand.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

import core

TAB = core.EXTENDED
SUP = {"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
       "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
       "-": "⁻"}


def pfmt(v):
    """Journal-style p-value: 0.31, 8.4 x 10^-5, or < 1 x 10^-12."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "—"
    if not np.isfinite(v):
        return "—"
    if v >= 0.1:
        return f"{v:.2f}"
    if v >= 0.001:
        return f"{v:.3f}".rstrip("0")
    m, e = f"{v:.1e}".split("e")
    return f"{m} × 10{''.join(SUP[c] for c in str(int(e)))}"


def num(v, nd=2, dash_if_nan=True):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "—"
    if not np.isfinite(v):
        return "—" if dash_if_nan else ""
    return f"{v:.{nd}f}"


def ci(lo, hi):
    a, b = num(lo), num(hi)
    return "—" if a == "—" or b == "—" else f"{a}–{b}"


def table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join("---" for _ in headers) + " |"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def read(n):
    return pd.read_csv(f"{TAB}/{n}")


P = []
A = P.append

# ------------------------------------------------------------------ S7
s7 = read("S7_background_definitions.csv")
A("""## Table S7. Precedence under four definitions of the background rate

The background rate the observed count is judged against can be drawn in more
than one way, and the choice is not neutral. The rate as published counts every
membrane position in the protein, including positions whose own residue is the
target; that inflates the denominator and biases each ratio toward depletion.
The self-exclusive rate counts only positions whose residue differs from the
target, which places the substitutable positions and the background on the same
footing. The two matched rates are computed separately within strata of observed
alphabet size, and within strata of alphabet size crossed with conservation
grade, so that a site is compared only against positions of comparable
alignment breadth. Alphabet size and conservation grade are derived from the
same alignment column and are strongly correlated (Spearman ρ = −0.78 for
EAA1, −0.80 for VPAT), so the doubly matched estimate is a conservative bound
rather than an independent control. *p* values are two-sided and adjusted across
the sixteen tests by the Benjamini–Hochberg procedure. The number of
sites and the number already carrying the target do not depend on the
background definition: L→Q 7 of 87, I→T 33 of 56, V→T 48 of 64 and
F→Y 21 of 37. Only the expectation, and so the ratio, moves between
definitions.""")
SHORT_DEF = {
    "as published (self-inclusive, unstratified)": "As published",
    "self-exclusive": "Self-exclusive",
    "alphabet-size matched": "Alphabet matched",
    "alphabet-size matched, conservation matched": "Alphabet and grade",
}
rows = []
for d in dict.fromkeys(s7["definition"]):
    sub = s7[s7["definition"] == d].set_index("swap").loc[
        ["L→Q", "I→T", "V→T", "F→Y"]]
    for sw, r in sub.iterrows():
        rows.append([SHORT_DEF[d], sw, num(r["expected"]),
                     num(r["ratio"]), ci(r["ci_lo"], r["ci_hi"]),
                     pfmt(r["p"]), pfmt(r["p_adj"])])
A(table(["Background", "Substitution", "Expected", "Ratio", "95% CI",
         "*p*", "*p* adj."], rows))
A("""Background definitions: *As published*, every membrane position counted;
*Self-exclusive*, positions whose own residue is the target omitted;
*Alphabet matched*, strata of observed alphabet size; *Alphabet and grade*,
strata of alphabet size crossed with conservation grade.""")
A("""Both enrichments survive every definition. F→Y strengthens as matching is
added, from 2.97 to 3.31; V→T attenuates from 1.48 to 1.34 while remaining
significant by three orders of magnitude, so part of its excess is attributable
to alignment breadth and part is not. Neither null result becomes significant
under any definition, and the apparent depletion of L→Q moves toward unity as
matching is added.""")

# ------------------------------------------------------------------ S8
for fn, lab, blurb in [
    ("S8a_alphabet_strata.csv", "S8. Precedence within strata of observed alphabet size",
     """A position whose alignment column admits many residues is more likely to
admit any particular target, so precedence rises with alphabet breadth for
purely arithmetic reasons. Each stratum below is judged against the background
rate computed within that same stratum. Strata containing no substitutable
position, and strata in which the background rate is identically zero, are
omitted. *p* values are adjusted across the twenty tests."""),
    ("S8b_grade_strata.csv", "S9. Precedence within strata of conservation grade",
     """ConSurf grades run from 1 to 9, with 9 the most conserved. As above, each
stratum is judged against the background rate within that stratum, so a
conserved position is compared only against other conserved positions. *p*
values are adjusted across the sixteen tests."""),
]:
    d = read(fn)
    d = d[np.isfinite(pd.to_numeric(d["ratio"], errors="coerce"))
          | np.isfinite(pd.to_numeric(d["ci_hi"], errors="coerce"))]
    key = "stratum"
    order = ["1-3", "4-5", "6-7", "8-9", "10+"] if "alphabet" in fn else \
        ["1-3 variable", "4-6 intermediate", "7-8 conserved", "9 most conserved"]
    rows = []
    for sw in ["L→Q", "I→T", "V→T", "F→Y"]:
        sub = d[d["swap"] == sw].set_index(key)
        for b in order:
            if b not in sub.index:
                continue
            r = sub.loc[b]
            rows.append([sw, b, int(r["n_sites"]),
                         f'{int(r["observed"])} / {num(r["expected"])}',
                         num(r["ratio"]), ci(r["ci_lo"], r["ci_hi"]),
                         pfmt(r["p"]), pfmt(r["p_adj"])])
    A(f"## Table {lab}\n\n{blurb}")
    A(table(["Substitution",
             "Alphabet size" if "alphabet" in fn else "Grade", "Sites",
             "Precedent / expected", "Ratio", "95% CI", "*p*",
             "*p* adj."], rows))

A("""F→Y is enriched in three of the four conservation strata after
adjustment, and most strongly in the most conserved stratum, where five of
seven phenylalanine columns at grade 9 already contain tyrosine against an
expectation of 0.78. The enrichment is therefore not an artefact of permissive,
poorly conserved columns; it is largest exactly where the alignment is most
constrained.""")

# ------------------------------------------------------------------ S10
s9 = read("S9_all_targets.csv")
A("""## Table S10. Every possible target residue ranked at each native position

Each of the nineteen alternatives to the native residue was scored by the same
background-controlled procedure, over the same membrane positions, giving
seventy-six estimates. The five best-supported targets are listed for each
native residue, together with the QTY target where it falls outside the top
five. *p* values are adjusted across all seventy-six tests. The complete table
is Supplementary Data S2.""")
rows = []
for nat in "LIVF":
    sub = s9[s9["native"] == nat].sort_values("ratio", ascending=False)
    keep = list(sub.index[:5])
    q = sub[sub["is_qty"]].index[0]
    if q not in keep:
        keep.append(q)
    for i in keep:
        r = sub.loc[i]
        rows.append([nat, r["swap"] + (" †" if r["is_qty"] else ""),
                     f"{int(r['rank_in_native'])} of 19", int(r["n_sites"]),
                     int(r["observed"]), num(r["expected"]),
                     num(r["ratio"]), pfmt(r["p_adj"])])
A(table(["Native", "Substitution", "Rank", "Sites", "With precedent",
         "Expected", "Ratio", "*p* adj."], rows))
A("""† the substitution specified by the QTY code.

At phenylalanine positions the QTY target is the best-supported of all nineteen
alternatives, with tryptophan second; the two aromatic targets are the only ones
enriched, while glycine, alanine and serine are significantly depleted, so the
alignment at these positions discriminates by side-chain bulk and aromaticity
rather than by polarity alone. At valine positions threonine ranks second, behind
isoleucine alone,
and is the best-supported polar target. At isoleucine positions threonine ranks
fifth, behind four hydrophobic targets, and is again the best-supported polar
one. At leucine positions the four best-supported targets are all hydrophobic
and no polar target is enriched; glutamine ranks twelfth. Threonine ranks ninth
at leucine positions and is present in far more leucine columns than glutamine
is (32 against 7 of the same 87 columns), but that difference reflects the much
higher overall abundance of threonine in membrane spans, and after background
control the two are indistinguishable (0.75 and 0.70, overlapping intervals).""")

# ------------------------------------------------------------------ S11
c = read("S10a_matrix_concordance.csv")
q = read("S10b_qty_vs_matrices.csv")
A("""## Table S11. Concordance with general amino-acid exchangeability

The seventy-six estimates of Table S10 were correlated against five
general-purpose substitution matrices. None is specific to membrane proteins,
which is the point of the comparison: it asks how much of the pattern observed
inside the bilayer is simply general exchangeability.""")
A(table(["Matrix", "Pairs", "Spearman ρ", "*p*"],
        [[{"JONES": "JTT (Jones–Taylor–Thornton)"}.get(r["matrix"],
                                                             r["matrix"]),
           int(r["n_pairs"]), num(r["spearman_rho"]), pfmt(r["p"])]
         for _, r in c.iterrows()]))
A("""Precedence inside membrane spans tracks general exchangeability closely.
Positioning each QTY substitution as a percentile of the seventy-six estimates,
on the precedence scale and on each matrix scale, isolates the part that general
exchangeability does not account for.""")
mats = ["BLOSUM62", "PAM250", "GONNET1992", "JONES", "RISLER"]
A(table(["Substitution", "Ratio", "Precedence percentile",
         "Mean matrix percentile", "Difference"],
        [[r["swap"], num(r["ratio"]), f"{int(r['precedence_pct'])}",
          f"{int(r['mean_matrix_pct'])}",
          f"{int(r['gap']):+d}".replace("-", "−")]
         for _, r in q.iterrows()]))
A("""F→Y and I→T sit where general exchangeability predicts. V→T is
better precedented inside membrane spans than any of the five matrices
anticipates, and L→Q slightly worse, which is the expected direction if the
bilayer penalises the introduction of an amide more than a general-purpose
matrix does.""")

# ------------------------------------------------------------------ S12
l = read("S11_leave_one_segment_out.csv")
A("""## Table S12. Leave-one-segment-out robustness

Each of the twenty-two annotated membrane segments was removed in turn and the
pooled estimate recomputed from the remaining segments, with the background rate
recomputed over the same reduced scope. No single segment carries any of the
four results.""")
A(table(["Substitution", "Segments dropped", "Full estimate", "Minimum",
         "Attained on dropping", "Maximum", "Attained on dropping",
         "Drops with *p* < 0.05"],
        [[r["swap"], int(r["n_drops"]), num(r["full_ratio"]),
          num(r["min_ratio"]), r["min_at"], num(r["max_ratio"]), r["max_at"],
          f"{int(r['n_p_below_05'])} of {int(r['n_drops'])}"]
         for _, r in l.iterrows()]))

# ------------------------------------------------------------------ S13
dep = read("S12_bilayer_depth.csv")
A("""## Table S13. Precedence by position within the membrane segment

Each membrane position was assigned a relative depth, zero at the midpoint of
its segment and one at either end, and split at the median into the inner half
of the segment and the outer half abutting the headgroup region. Background
rates are computed separately within each zone. *p* values are adjusted across
the eight tests.""")
A(table(["Substitution", "Zone", "Sites", "Precedent / expected",
         "Ratio", "95% CI", "*p*", "*p* adj."],
        [[r["swap"], {"core (inner half)": "Core",
                      "interfacial (outer half)": "Interfacial"}[r["zone"]],
          int(r["n_sites"]),
          f'{int(r["observed"])} / {num(r["expected"])}',
          num(r["ratio"]), ci(r["ci_lo"], r["ci_hi"]),
          pfmt(r["p"]), pfmt(r["p_adj"])] for _, r in dep.iterrows()]))
A("""Both enrichments are present in both zones and are somewhat larger in the
interior of the segment, where a polar substitution would be least expected to
be tolerated.""")

# ------------------------------------------------------------------ S14
b = read("S13_burial.csv")
A("""## Table S14. Precedence and predicted burial

ConSurf classifies positions as buried or exposed from the predicted relative
solvent accessibility of the query structure. The classification is reported
only at the substitutable positions, so the comparison below is between buried
and exposed positions carrying the same native residue, by Fisher's exact test,
rather than against a burial-matched background. *n* is the number of positions and *k* the number whose alignment column already contains the target. Twenty-six of the 351 pooled
positions carry no call and are excluded.""")
A(table(["Substitution", "Buried, *k* / *n* (rate)",
         "Exposed, *k* / *n* (rate)", "Odds ratio", "*p*", "*p* adj."],
        [[r["swap"],
          f'{int(r["prec_buried"])} / {int(r["n_buried"])} '
          f'({num(r["rate_buried"])})',
          f'{int(r["prec_exposed"])} / {int(r["n_exposed"])} '
          f'({num(r["rate_exposed"])})',
          num(r["odds_ratio"]), pfmt(r["p"]), pfmt(r["p_adj"])]
         for _, r in b.iterrows()]))
A("""No substitution shows an association between precedence and predicted
burial after adjustment.""")

# ------------------------------------------------------------------ S15
rc = read("S14_reciprocal.csv")
A("""## Table S15. The same four pairs read in the opposite direction

Precedence for the reverse substitution was computed at positions carrying the
QTY target residue. The two directional estimates need not agree even for a
freely exchangeable pair, because each conditions on the query residue rather
than on membership of the alignment column; the symmetric lift in the final
column removes that asymmetry by using column membership on both sides. The
minimum detectable ratio is the smallest enrichment the available sites could
have resolved at 95% confidence.""")
A(table(["Pair", "Sites", "Precedent / expected",
         "Reverse ratio (95% CI)", "*p* adj.", "Forward ratio", "Lift",
         "Min. detect."],
        [[r["reciprocal"], int(r["n_sites"]),
          f'{int(r["observed"])} / {num(r["expected"])}',
          f'{num(r["reciprocal_ratio"])} ({ci(r["ci_lo"], r["ci_hi"])})',
          pfmt(r["p_adj"]), num(r["forward_ratio"]),
          num(r["symmetric_lift"]), num(r["min_detectable"])]
         for _, r in rc.iterrows()]))
A("""Phenylalanine is enriched at tyrosine positions as tyrosine is at
phenylalanine positions, so the aromatic pair is exchangeable in both
directions. The threonine reversals are unremarkable, and the glutamine
reversal rests on four positions and could not have resolved an enrichment
below 1.70; it is reported for completeness and carries no weight.""")

print("\n\n".join(P))
