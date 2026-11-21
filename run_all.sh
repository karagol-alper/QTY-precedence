#!/usr/bin/env bash
# Reproduce every number and figure in the paper from the ConSurf tables.
set -euo pipefail
cd "$(dirname "$0")"

RUN="${QTY_RUN:-results/run_2026-11-04_161147}"

echo "== 1. rebuilding the primary estimates from the raw tables =="
# If these do not match the published run, nothing downstream is trustworthy.
python3 analysis/check_repro.py

echo
echo "== 2. extended analyses -> results/extended =="
python3 analysis/analyses.py > results/extended/analyses_log.txt
echo "   wrote results/extended/ and analyses_log.txt"

echo
echo "== 3. figures =="
python3 analysis/make_fig1.py "$RUN" supplementary/figure1.png
python3 analysis/make_fig1.py "$RUN" supplementary/figure1.eps
python3 analysis/make_fig2.py supplementary/figure2.png
python3 analysis/make_fig2.py supplementary/figure2.eps
python3 analysis/make_figS1.py supplementary/figureS1.png
python3 analysis/make_figS1.py supplementary/figureS1.eps

echo
echo "== 4. supplementary tables S7-S15 (markdown) =="
python3 analysis/gen_supp.py > results/extended/supplementary_tables.md
echo "   wrote results/extended/supplementary_tables.md"

echo
echo "== 5. verification =="
python3 analysis/verify.py
