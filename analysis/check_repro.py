"""Does rebuilding from the raw tables reproduce the run's published ratios?"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import core

t = core.load_tables()
pub = pd.read_csv(core.RUN + "/tables/pooled_by_swap.csv")
print("published:")
print(pub.to_string(index=False))
print("\nrebuilt:")
rows = []
for native, target in core.QTY.items():
    r = core.ratio_block(t, native, target, exclude_self=False)
    rows.append(dict(swap=f"{native}->{target}", n=r["n"], k=r["k"],
                     exp=round(r["exp"], 2), ratio=round(r["ratio"], 3),
                     p=f"{r['p']:.3g}",
                     ci=f"{r['ci_lo']:.2f}-{r['ci_hi']:.2f}",
                     **{k2.split()[0]: round(v["ratio"], 3)
                        for k2, v in r["per"].items()}))
print(pd.DataFrame(rows).to_string(index=False))
