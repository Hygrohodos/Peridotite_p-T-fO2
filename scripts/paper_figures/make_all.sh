#!/usr/bin/env bash
# Figures 1-10, Table 1 and the quoted numbers (see make_paper_figures.py; the match level is set there).
# Output: $PAPER_FIG_OUT (default scripts/paper_figures/out, git-ignored); nothing is written into data/.
set -euo pipefail
cd "$(dirname "$0")"
python3 make_paper_figures.py
