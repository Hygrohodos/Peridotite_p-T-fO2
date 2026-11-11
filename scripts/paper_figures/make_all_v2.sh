#!/usr/bin/env bash
# Revision v2 figure set in the manuscript's own style (Figs 1-11).
# Output: $PAPER_FIG_OUT (default scripts/paper_figures/out). Run from anywhere.
set -euo pipefail
cd "$(dirname "$0")"
export PAPER_FIG_OUT="${PAPER_FIG_OUT:-$(pwd)/out}"; mkdir -p "$PAPER_FIG_OUT"
python3 ../anchor_level.py >/dev/null
python3 fig1.py; python3 fig2.py; python3 fig3.py           # unchanged figures 1-3
python3 v2_fig04_all_xenoliths.py                           # new Fig. 4 (Paolo, line 255)
python3 ../make_fig04_weightcircles.py --outdir "$PAPER_FIG_OUT/_fig5tmp" >/dev/null   # Fig. 5 = experiments (was Fig. 4)
cp "$PAPER_FIG_OUT/_fig5tmp/Fig04_TA98_WZ07_measured_Fe3.png" "$PAPER_FIG_OUT/Fig5_final.png"
cp "$PAPER_FIG_OUT/_fig5tmp/Fig04_weights.csv" "$PAPER_FIG_OUT/Fig5_weights.csv"
python3 v2_fig06_fe3_three_levels.py; python3 v2_fig07_level_sensitivity.py; python3 v2_fig08_fe3_adopted_level.py
python3 v2_fig09_validation.py; python3 v2_fig10_fo2_measured.py; python3 v2_fig11_fo2_estimated.py
