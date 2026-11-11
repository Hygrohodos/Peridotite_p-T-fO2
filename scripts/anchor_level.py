#!/usr/bin/env python3
"""Compute the anchored flat WZ07 match level (revision v3, option C).

Writes data/derived/level_anchor.json by default, or --out PATH (used by
scripts/paper_figures/make_all.sh so figure builds do not write into the repo).
--check exits non-zero if params.dt_match_c differs from the anchor by > 0.05 C.
"""
import argparse, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
import pandas as pd
import make_figures as mf
from calculations.pt_fo2.level_anchor import anchor_level
from calculations.pt_fo2.params import PARAMS

ap = argparse.ArgumentParser()
ap.add_argument("--out", default=str(ROOT / "data" / "derived" / "level_anchor.json"))
ap.add_argument("--check", action="store_true")
a_ = ap.parse_args()
df = pd.read_csv(mf.FULL, encoding="utf-8-sig", low_memory=False)
df = df.rename(columns={df.columns[0]: "SampleID"})
a = anchor_level(df)
a["params_dt_match_c"] = PARAMS.fe3.dt_match_c
Path(a_.out).write_text(json.dumps(a, indent=1))
print(json.dumps({k: v for k, v in a.items() if k != "sample_ids"}, indent=1))
if a_.check and abs(a["dt_level_C"] - PARAMS.fe3.dt_match_c) > 0.05:
    sys.exit(f"params.dt_match_c {PARAMS.fe3.dt_match_c} != anchor {a['dt_level_C']:.2f}")
