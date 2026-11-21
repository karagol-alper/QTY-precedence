import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import core

T = core.load_tables()

print("=== collinearity of the two matching covariates (TM positions) ===")
for name, d in T.items():
    tm = d[d.in_tm].dropna(subset=["grade"])
    rho, p = stats.spearmanr(tm["n_variety"], tm["grade"])
    print(f"{core.SHORT[name]}: spearman(n_variety, grade) = {rho:+.3f}  p={p:.2g}  n={len(tm)}")

print("\n=== stratum sizes in the doubly matched analysis ===")
def joint(f):
    return core.alphabet_band(f).astype(str) + "|" + core.grade_band(f).astype(str)
for native, target in core.SWAPS:
    small = []
    for name, d in T.items():
        scope = d[d.in_tm]
        ls, lt = joint(scope), joint(d[d.in_tm & (d.aa == native)])
        for lab in pd.unique(lt):
            bgset = scope[(ls == lab) & (scope.aa != target)]
            nsite = int((lt == lab).sum())
            if len(bgset) < 10:
                small.append((core.SHORT[name], lab, len(bgset), nsite))
    tot = sum(s[3] for s in small)
    print(f"{native}->{target}: {len(small)} strata with <10 background positions, "
          f"covering {tot} sites")
    for s in small[:4]:
        print("      ", s)

print("\n=== symmetric co-occurrence lift (set membership both sides) ===")
rows = []
for native, target in core.SWAPS:
    num = den_a = den_b = n = 0
    for d in T.values():
        tm = d[d.in_tm]
        a = tm["vset"].apply(lambda v: native in v)
        b = tm["vset"].apply(lambda v: target in v)
        num += int((a & b).sum()); den_a += int(a.sum()); den_b += int(b.sum()); n += len(tm)
    lift = (num / n) / ((den_a / n) * (den_b / n))
    rows.append(dict(pair=f"{native}/{target}", n_tm=n, both=num,
                     p_native=round(den_a / n, 3), p_target=round(den_b / n, 3),
                     symmetric_lift=round(lift, 3)))
print(pd.DataFrame(rows).to_string(index=False))

print("\n=== is L->T better precedented than L->Q? ===")
rq = core.ratio_block(T, "L", "Q", exclude_self=False)
rt = core.ratio_block(T, "L", "T", exclude_self=False)
kq, kt = rq["k"], rt["k"]
print(f"L->Q {kq}/{rq['exp']:.2f} = {rq['ratio']:.3f}   L->T {kt}/{rt['exp']:.2f} = {rt['ratio']:.3f}")
# paired across the same 87 sites: McNemar on discordant pairs
disc_q = disc_t = 0
for d in T.values():
    s = core.sites(d, "L")
    q = s["vset"].apply(lambda v: "Q" in v); t = s["vset"].apply(lambda v: "T" in v)
    disc_q += int((q & ~t).sum()); disc_t += int((t & ~q).sum())
res = stats.binomtest(disc_q, disc_q + disc_t, 0.5) if (disc_q + disc_t) else None
print(f"same 87 leucine sites: Q-only={disc_q}, T-only={disc_t}, "
      f"exact p={res.pvalue:.3g}" if res else "no discordant pairs")
