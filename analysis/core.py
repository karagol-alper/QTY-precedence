"""Shared machinery for the extended analyses.

Everything is rebuilt from the run's own ConSurf tables and its recorded TM
boundaries, then checked against the run's published pooled ratios
(check_repro.py). If that check fails, nothing downstream is trustworthy.

Two background definitions are supported. The published one counts every
in-scope position, including positions whose own residue is the target; that
inflates the denominator a little and biases every ratio toward depletion.
The self-excluding one counts only positions whose residue differs from the
target, putting numerator sites (whose residue is by construction not the
target) and background positions on the same footing. The self-excluding
definition is also what makes alphabet-size stratification coherent: a
position with a single observed residue can never carry a different target,
so its stratum contributes zero expected and zero observed rather than a
spurious expectation.
"""
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
# Overridable so a different run or a different ConSurf set can be analysed
# without editing the source.
RUN = os.environ.get("QTY_RUN", str(ROOT / "results" / "run_2026-11-04_161147"))
DATA = os.environ.get("QTY_DATA", str(ROOT / "data"))
EXTENDED = os.environ.get("QTY_OUT", str(ROOT / "results" / "extended"))

QTY = {"L": "Q", "I": "T", "V": "T", "F": "Y"}
SWAPS = [("L", "Q"), ("I", "T"), ("V", "T"), ("F", "Y")]
AA20 = list("ACDEFGHIKLMNPQRSTVWY")

FILES = {
    "EAA1 / SLC1A3": "eaa1_residue_variety.csv",
    "VPAT / SLC18B1": "vpat_residue_variety.csv",
}

# The run's EAA1 input carries no conservation grade. Taking grades from the
# per-site table instead would restrict them to L/I/V/F positions, which would
# make every grade-matched background hydrophobic-only and the matching
# worthless. This file covers SEQRES 43-499 at the same MSA depth (150), is
# identical to the run input over that span residue for residue, and contains
# the whole EAA1 membrane span (48-494); its grades agree exactly with an
# independent transcription of the same ConSurf job.
GRADE_FILES = {"EAA1 / SLC1A3": "eaa1_residue_variety_std.csv"}
SHORT = {"EAA1 / SLC1A3": "EAA1", "VPAT / SLC18B1": "VPAT"}


def manifest():
    with open(f"{RUN}/run_manifest.json") as fh:
        return json.load(fh)


def load_tables():
    """{protein: DataFrame} with vset, in_tm, grade, seg, depth, buried."""
    m = manifest()
    per_site = pd.read_csv(f"{RUN}/tables/per_site.csv")
    out = {}
    for name, fn in FILES.items():
        d = pd.read_csv(f"{DATA}/{fn}")
        d["vset"] = d["variety"].fillna("").apply(
            lambda s: frozenset(x for x in str(s).split("|") if x and x != "-"))
        regions = m["tm_regions"][name]
        pos = d["seqres_pos"].to_numpy()
        in_tm = np.zeros(len(d), dtype=bool)
        seg = np.full(len(d), -1, dtype=int)
        depth = np.full(len(d), np.nan)
        for j, (a, b) in enumerate(sorted(regions)):
            sel = (pos >= a) & (pos <= b)
            in_tm |= sel
            seg[sel] = j
            mid = (a + b) / 2.0
            half = max((b - a) / 2.0, 1e-9)
            # 0 at the helix midpoint, 1 at either end of the segment
            depth[sel] = np.abs(pos[sel] - mid) / half
        d["in_tm"], d["seg"], d["depth"] = in_tm, seg, depth

        if "grade" not in d.columns:
            g = pd.read_csv(f"{DATA}/{GRADE_FILES[name]}")
            chk = d.merge(g, on="seqres_pos", suffixes=("", "_g"))
            assert (chk["aa"] == chk["aa_g"]).all(), f"{name}: grade file disagrees"
            assert (chk["n_variety"] == chk["n_variety_g"]).all(), \
                f"{name}: grade file disagrees on alphabet size"
            d = d.merge(g[["seqres_pos", "grade", "score"]], on="seqres_pos",
                        how="left")
        # buried/exposed is only published at the substitutable positions
        ann = (per_site[per_site["protein"] == name][["seqres_pos", "buried"]]
               .drop_duplicates("seqres_pos"))
        d = d.merge(ann, on="seqres_pos", how="left")
        miss = int(d["in_tm"].sum() - d.loc[d["in_tm"], "grade"].notna().sum())
        assert miss == 0, f"{name}: {miss} membrane positions lack a grade"
        d["protein"] = name
        out[name] = d
    return out


def _scope(d, mask):
    sub = d if mask is None else d[mask]
    return sub[sub["in_tm"]]


def background(d, target, mask=None, exclude_self=True):
    """Share of in-scope positions whose observed alphabet contains `target`."""
    sub = _scope(d, mask)
    if exclude_self:
        sub = sub[sub["aa"] != target]
    if len(sub) == 0:
        return float("nan"), 0
    return float(sub["vset"].apply(lambda s: target in s).mean()), len(sub)


def sites(d, native, mask=None):
    sub = _scope(d, mask)
    return sub[sub["aa"] == native]


def _summarise(k, ps):
    """Poisson-binomial summary of k successes against per-site rates ps."""
    if not len(ps):
        nan = float("nan")
        return dict(n=0, k=k, exp=0.0, sd=0.0, ratio=nan, z=nan, p=nan,
                    ci_lo=nan, ci_hi=nan)
    p = np.asarray(ps, dtype=float)
    exp = float(p.sum())
    sd = math.sqrt(float((p * (1 - p)).sum()))
    ratio = k / exp if exp > 0 else float("nan")
    z = (k - exp) / sd if sd > 0 else float("nan")
    pv = math.erfc(abs(z) / math.sqrt(2)) if np.isfinite(z) else float("nan")
    return dict(n=len(p), k=int(k), exp=exp, sd=sd, ratio=ratio, z=z, p=pv,
                ci_lo=max(0.0, (k - 1.96 * sd) / exp) if exp > 0 else float("nan"),
                ci_hi=(k + 1.96 * sd) / exp if exp > 0 else float("nan"))


def ratio_block(tables, native, target, mask_by=None, exclude_self=True,
                strata=None, site_filter=None):
    """Pooled k / expected for native->target.

    Each protein contributes its own within-scope background, so the pooling is
    over sites rather than over ratios. `strata` is a callable returning a
    stratum label per row; when given, the background is computed separately
    within each stratum of each protein, which is what controls for a
    covariate (observed alphabet size, conservation grade) rather than merely
    describing it.
    """
    k_tot, ps, per = 0, [], {}
    for name, d in tables.items():
        mask = None if mask_by is None else mask_by(d)
        s = sites(d, native, mask)
        if site_filter is not None:
            s = s[site_filter(s)]
        if len(s) == 0:
            per[name] = dict(n=0, k=0, ratio=float("nan"))
            continue
        if strata is None:
            bg, n_scope = background(d, target, mask, exclude_self)
            rates = [bg] * len(s)
            kk = int(s["vset"].apply(lambda v: target in v).sum())
            per[name] = dict(n=len(s), k=kk, bg=bg, n_scope=n_scope,
                             ratio=(kk / (len(s) * bg)) if bg else float("nan"))
        else:
            scope = _scope(d, mask)
            lab_scope = strata(scope)
            lab_site = strata(s)
            rates, kk, exp_p = [], 0, 0.0
            for lab in pd.unique(lab_site):
                in_lab = scope[lab_scope == lab]
                if exclude_self:
                    in_lab = in_lab[in_lab["aa"] != target]
                bg = (float(in_lab["vset"].apply(lambda v: target in v).mean())
                      if len(in_lab) else float("nan"))
                grp = s[lab_site == lab]
                if not np.isfinite(bg):
                    continue
                rates.extend([bg] * len(grp))
                kk += int(grp["vset"].apply(lambda v: target in v).sum())
                exp_p += bg * len(grp)
            per[name] = dict(n=len(rates), k=kk,
                             ratio=(kk / exp_p) if exp_p > 0 else float("nan"))
        k_tot += per[name]["k"]
        ps.extend(rates)
    out = _summarise(k_tot, ps)
    out["per"] = per
    return out


def fdr(pvals):
    """Benjamini-Hochberg adjusted p-values, input order preserved."""
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    v = p[ok]
    if v.size == 0:
        return out
    order = np.argsort(v)
    n = v.size
    adj = v[order] * n / np.arange(1, n + 1)
    adj = np.minimum(np.minimum.accumulate(adj[::-1])[::-1], 1.0)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def grade_band(frame):
    """ConSurf grade collapsed into four bands; 9 is the most conserved."""
    g = frame["grade"]
    return pd.cut(g, bins=[0, 3, 6, 8, 9],
                  labels=["1-3 variable", "4-6 intermediate",
                          "7-8 conserved", "9 most conserved"])


def alphabet_band(frame):
    """Observed alphabet size collapsed into five bands."""
    return pd.cut(frame["n_variety"], bins=[0, 3, 5, 7, 9, 100],
                  labels=["1-3", "4-5", "6-7", "8-9", "10+"])
