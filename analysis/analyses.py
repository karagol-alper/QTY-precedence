"""The extended analyses: supplementary tables S7-S14.

Each writes a CSV to out/. Nothing here changes the primary estimates; these
are controls, robustness checks and an external concordance test. The primary
definition (background inclusive of positions whose residue is the target) is
retained wherever a number has to line up with the main text; the
self-excluding definition is carried alongside as a sensitivity.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy import stats

import core

OUT = core.EXTENDED
os.makedirs(OUT, exist_ok=True)
T = core.load_tables()
LAB = {("L", "Q"): "L→Q", ("I", "T"): "I→T",
       ("V", "T"): "V→T", ("F", "Y"): "F→Y"}


def _w(df, name):
    path = os.path.join(OUT, name)
    df.to_csv(path, index=False)
    print(f"\n=== {name} ===")
    print(df.to_string(index=False))
    return df


# ---------------------------------------------------------------- S7
def s7_confounds():
    """Four background definitions, crossed: self-inclusion x alphabet matching."""
    rows = []
    defs = [
        ("as published (self-inclusive, unstratified)", False, None),
        ("self-exclusive", True, None),
        ("alphabet-size matched", True, core.alphabet_band),
        ("alphabet-size matched, conservation matched", True,
         lambda f: core.alphabet_band(f).astype(str) + "|"
         + core.grade_band(f).astype(str)),
    ]
    for label, excl, strat in defs:
        for native, target in core.SWAPS:
            r = core.ratio_block(T, native, target, exclude_self=excl,
                                 strata=strat)
            rows.append(dict(definition=label, swap=LAB[(native, target)],
                             n_sites=r["n"], observed=r["k"],
                             expected=round(r["exp"], 2),
                             ratio=round(r["ratio"], 3),
                             ci_lo=round(r["ci_lo"], 3),
                             ci_hi=round(r["ci_hi"], 3), p=r["p"]))
    df = pd.DataFrame(rows)
    df["p_adj"] = core.fdr(df["p"])
    df["p"] = df["p"].map(lambda v: f"{v:.3g}")
    df["p_adj"] = df["p_adj"].map(lambda v: f"{v:.3g}")
    return _w(df, "S7_background_definitions.csv")


# ---------------------------------------------------------------- S8
def s8_strata(band_fn, name, fname):
    """Ratios within strata of a covariate, each against its own background."""
    rows = []
    for native, target in core.SWAPS:
        bands = sorted({str(b) for d in T.values()
                        for b in band_fn(d[d.in_tm & (d.aa == native)]).dropna()})
        for b in bands:
            r = core.ratio_block(
                T, native, target, exclude_self=True,
                strata=band_fn,
                site_filter=lambda f, b=b: band_fn(f).astype(str) == b)
            if r["n"] == 0:
                continue
            rows.append(dict(swap=LAB[(native, target)], stratum=b,
                             n_sites=r["n"], observed=r["k"],
                             expected=round(r["exp"], 2),
                             ratio=round(r["ratio"], 3),
                             ci_lo=round(r["ci_lo"], 3),
                             ci_hi=round(r["ci_hi"], 3), p=r["p"]))
    df = pd.DataFrame(rows)
    df["p_adj"] = core.fdr(df["p"])
    for c in ("p", "p_adj"):
        df[c] = df[c].map(lambda v: f"{v:.3g}")
    df.attrs["name"] = name
    return _w(df, fname)


# ---------------------------------------------------------------- S9
def s9_all_targets():
    """Every possible target at each native residue, ranked."""
    rows = []
    for native in "LIVF":
        for target in core.AA20:
            if target == native:
                continue
            r = core.ratio_block(T, native, target, exclude_self=False)
            rows.append(dict(native=native, target=target,
                             swap=f"{native}→{target}",
                             n_sites=r["n"], observed=r["k"],
                             expected=round(r["exp"], 3),
                             ratio=r["ratio"], p=r["p"],
                             is_qty=(core.QTY[native] == target)))
    df = pd.DataFrame(rows)
    df["p_adj"] = core.fdr(df["p"])
    df["rank_in_native"] = df.groupby("native")["ratio"].rank(
        ascending=False, method="min")
    df["n_targets_ranked"] = df.groupby("native")["ratio"].transform(
        lambda s: s.notna().sum())
    df = df.sort_values(["native", "ratio"], ascending=[True, False])
    _w(df.assign(ratio=df["ratio"].round(3),
                 p=df["p"].map(lambda v: f"{v:.3g}"),
                 p_adj=df["p_adj"].map(lambda v: f"{v:.3g}")),
       "S9_all_targets.csv")
    print("\n--- where the QTY target ranks among all targets at its residue ---")
    print(df[df.is_qty][["swap", "ratio", "rank_in_native",
                         "n_targets_ranked", "p_adj"]].to_string(index=False))
    return df


# ---------------------------------------------------------------- S10
def s10_matrices(s9):
    """Does observed precedence track general amino-acid exchangeability?"""
    from Bio.Align import substitution_matrices as sm
    names = ["BLOSUM62", "PAM250", "GONNET1992", "JONES", "RISLER"]
    mats = {n: sm.load(n) for n in names}
    d = s9[np.isfinite(s9["ratio"])].copy()
    rows = []
    for n, m in mats.items():
        d[n] = [float(m[a, b]) for a, b in zip(d["native"], d["target"])]
        rho, p = stats.spearmanr(d["ratio"], d[n])
        rows.append(dict(matrix=n, n_pairs=len(d), spearman_rho=round(rho, 3),
                         p=f"{p:.3g}"))
    corr = pd.DataFrame(rows)
    _w(corr, "S10a_matrix_concordance.csv")

    # Where each QTY swap sits on each scale, as a percentile among all pairs.
    pct = []
    d["ratio_pct"] = 100 * d["ratio"].rank(pct=True)
    for n in names:
        d[n + "_pct"] = 100 * d[n].rank(pct=True)
    for _, r in d[d.is_qty].iterrows():
        row = dict(swap=r["swap"], ratio=round(r["ratio"], 3),
                   precedence_pct=round(r["ratio_pct"]))
        for n in names:
            row[n] = r[n]
            row[n + "_pct"] = round(r[n + "_pct"])
        pct.append(row)
    pctd = pd.DataFrame(pct)
    pctd["mean_matrix_pct"] = pctd[[n + "_pct" for n in names]].mean(axis=1).round()
    pctd["gap"] = (pctd["precedence_pct"] - pctd["mean_matrix_pct"]).round()
    _w(pctd, "S10b_qty_vs_matrices.csv")
    d.to_csv(os.path.join(OUT, "S10c_all_pairs_with_scores.csv"), index=False)
    return corr, pctd


# ---------------------------------------------------------------- S11
def s11_loso():
    """Leave one membrane segment out, one protein at a time."""
    m = core.manifest()
    rows = []
    for native, target in core.SWAPS:
        full = core.ratio_block(T, native, target, exclude_self=False)["ratio"]
        vals = []
        for pname in T:
            for j in range(len(m["tm_regions"][pname])):
                r = core.ratio_block(
                    T, native, target, exclude_self=False,
                    mask_by=lambda d, pn=pname, jj=j: ~(
                        (d["protein"] == pn) & (d["seg"] == jj)))
                if np.isfinite(r["ratio"]):
                    vals.append((core.SHORT[pname], j + 1, r["ratio"], r["p"]))
        rat = [v[2] for v in vals]
        lo, hi = min(vals, key=lambda v: v[2]), max(vals, key=lambda v: v[2])
        sig = sum(v[3] < 0.05 for v in vals)
        rows.append(dict(swap=LAB[(native, target)], n_drops=len(vals),
                         full_ratio=round(full, 3),
                         min_ratio=round(min(rat), 3),
                         min_at=f"{lo[0]} TM{lo[1]}",
                         max_ratio=round(max(rat), 3),
                         max_at=f"{hi[0]} TM{hi[1]}",
                         n_p_below_05=sig,
                         conclusion_stable=("yes" if sig in (0, len(vals))
                                            else "no")))
    return _w(pd.DataFrame(rows), "S11_leave_one_segment_out.csv")


# ---------------------------------------------------------------- S12
def s12_depth():
    """Helix core against the bilayer interface, each with its own background."""
    def zone(f):
        return pd.cut(f["depth"], bins=[-0.01, 0.5, 1.01],
                      labels=["core (inner half)", "interfacial (outer half)"])
    rows = []
    for native, target in core.SWAPS:
        for z in ["core (inner half)", "interfacial (outer half)"]:
            r = core.ratio_block(
                T, native, target, exclude_self=True, strata=zone,
                site_filter=lambda f, z=z: zone(f).astype(str) == z)
            rows.append(dict(swap=LAB[(native, target)], zone=z,
                             n_sites=r["n"], observed=r["k"],
                             expected=round(r["exp"], 2),
                             ratio=round(r["ratio"], 3),
                             ci_lo=round(r["ci_lo"], 3),
                             ci_hi=round(r["ci_hi"], 3), p=r["p"]))
    df = pd.DataFrame(rows)
    df["p_adj"] = core.fdr(df["p"])
    for c in ("p", "p_adj"):
        df[c] = df[c].map(lambda v: f"{v:.3g}")
    return _w(df, "S12_bilayer_depth.csv")


# ---------------------------------------------------------------- S13
def s13_burial():
    """Is precedence associated with ConSurf's buried/exposed call?"""
    rows = []
    for native, target in core.SWAPS:
        tab = np.zeros((2, 2), dtype=int)   # [buried/exposed] x [prec/not]
        for d in T.values():
            s = core.sites(d, native)
            s = s[s["buried"].isin(["b", "e"])]
            prec = s["vset"].apply(lambda v: target in v)
            for i, b in enumerate(["b", "e"]):
                sel = s["buried"] == b
                tab[i, 0] += int((sel & prec).sum())
                tab[i, 1] += int((sel & ~prec).sum())
        odds, p = stats.fisher_exact(tab)
        n_b, n_e = tab[0].sum(), tab[1].sum()
        rows.append(dict(
            swap=LAB[(native, target)],
            n_buried=int(n_b), prec_buried=int(tab[0, 0]),
            rate_buried=round(tab[0, 0] / n_b, 3) if n_b else float("nan"),
            n_exposed=int(n_e), prec_exposed=int(tab[1, 0]),
            rate_exposed=round(tab[1, 0] / n_e, 3) if n_e else float("nan"),
            odds_ratio=round(odds, 3) if np.isfinite(odds) else float("nan"),
            p=p))
    df = pd.DataFrame(rows)
    df["p_adj"] = core.fdr(df["p"])
    for c in ("p", "p_adj"):
        df[c] = df[c].map(lambda v: f"{v:.3g}")
    return _w(df, "S13_burial.csv")


# ---------------------------------------------------------------- S14
def _symmetric_lift(native, target):
    """Co-occurrence lift with set membership on both sides, so symmetric.

    The two directional ratios need not agree even when a pair is freely
    exchangeable, because each conditions on the query residue rather than on
    membership of the column's alphabet. This statistic removes that asymmetry
    and so says how much of any forward/reverse gap is structural.
    """
    num = da = db = n = 0
    for d in T.values():
        tm = d[d.in_tm]
        a = tm["vset"].apply(lambda v: native in v)
        b = tm["vset"].apply(lambda v: target in v)
        num += int((a & b).sum())
        da += int(a.sum())
        db += int(b.sum())
        n += len(tm)
    return (num / n) / ((da / n) * (db / n))


def s14_reciprocal():
    """The same four pairs read in the opposite direction."""
    rows = []
    for native, target in core.SWAPS:
        r = core.ratio_block(T, target, native, exclude_self=False)
        fwd = core.ratio_block(T, native, target, exclude_self=False)
        rows.append(dict(
            forward=LAB[(native, target)],
            reciprocal=f"{target}→{native}",
            n_sites=r["n"], observed=r["k"], expected=round(r["exp"], 2),
            reciprocal_ratio=round(r["ratio"], 3) if np.isfinite(r["ratio"]) else "",
            ci_lo=round(r["ci_lo"], 3) if np.isfinite(r["ci_lo"]) else "",
            ci_hi=round(r["ci_hi"], 3) if np.isfinite(r["ci_hi"]) else "",
            p=r["p"], forward_ratio=round(fwd["ratio"], 3),
            symmetric_lift=round(_symmetric_lift(native, target), 3),
            min_detectable=(round(1 + 1.96 * r["sd"] / r["exp"], 2)
                            if r["exp"] > 0 else "")))
    df = pd.DataFrame(rows)
    df["p_adj"] = core.fdr(df["p"])
    for c in ("p", "p_adj"):
        df[c] = df[c].map(lambda v: f"{v:.3g}" if np.isfinite(v) else "")
    return _w(df, "S14_reciprocal.csv")


if __name__ == "__main__":
    s7_confounds()
    s8_strata(core.alphabet_band, "alphabet size", "S8a_alphabet_strata.csv")
    s8_strata(core.grade_band, "conservation grade", "S8b_grade_strata.csv")
    s9 = s9_all_targets()
    s10_matrices(s9)
    s11_loso()
    s12_depth()
    s13_burial()
    s14_reciprocal()
    print("\nwrote to", OUT)
