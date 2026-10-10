#!/usr/bin/env python3
"""
Derive the locked Eq. 5 coefficients and the ΔT match threshold.

Usage (from the project root, after ``python run_pipeline.py``):

    python scripts/fit_calibration.py                  # report + check params.py
    python scripts/fit_calibration.py --update-params  # write values to params.py

Steps (see calculations/pt_fo2/calibration.py):
  1. Eq. 5 = least-squares line of (measured − estimated) Fe3+/ΣFe on
     T_TA98 at ΔT = 25 °C, for the Slave samples inside the dashed
     envelope of Fig. 5b (error groups removed, T ≥ 1130 °C).
  2. Fig. 6 sweep of ΔT with Eq. 5 held fixed, scored on those samples
     only; the locked ΔT has the highest ±0.02 rate (ties: smallest |bias|).

Outputs:
  data/derived/emp_correction_fit.json
  data/derived/dt_match_sweep_corr.csv
If the coefficients are changed, re-run run_pipeline.py and make_figures.py.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calculations.pt_fo2.params import PARAMS  # noqa: E402
from calculations.pt_fo2 import calibration as cal  # noqa: E402

DERIVED = ROOT / "data" / "derived"
PARAMS_PY = ROOT / "calculations" / "pt_fo2" / "params.py"


def load_screened():
    df = pd.read_csv(DERIVED / "screened_PT.csv", encoding="utf-8-sig",
                     low_memory=False)
    return df.rename(columns={df.columns[0]: "SampleID"})


def run(df):
    fit = cal.fit_empirical_correction(df)
    sweep = cal.dt_sweep(df, fit["slope"], fit["intercept"],
                         mask=fit["masks"]["reg"])
    dt = cal.select_dt(sweep)
    return fit, sweep, dt


def _update_params(slope, intercept, dt, nreg):
    s = PARAMS_PY.read_text(encoding="utf-8")
    s = re.sub(r"(emp_slope: float = )[-0-9.eE+]+", rf"\g<1>{slope:.4e}", s)
    s = re.sub(r"(emp_intercept: float = )[-0-9.eE+]+", rf"\g<1>{intercept:.5f}", s)
    s = re.sub(r"(dt_match_c: float = )[-0-9.eE+]+", rf"\g<1>{float(dt):.1f}", s)
    s = re.sub(r"(envelope, n = )\d+", rf"\g<1>{nreg}", s)
    PARAMS_PY.write_text(s, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--update-params", action="store_true")
    args = ap.parse_args()

    df = load_screened()
    fit, sweep, dt = run(df)
    sweep.to_csv(DERIVED / "dt_match_sweep_corr.csv", index=False)
    rec = {k: v for k, v in fit.items() if k not in ("est_fit", "masks")}
    rec["selected_dt_match_c"] = dt
    rec["selected_row"] = sweep.loc[sweep["dt"] == dt].iloc[0].to_dict()
    (DERIVED / "emp_correction_fit.json").write_text(json.dumps(rec, indent=2))

    print(f"Eq. 5 (ΔT = {fit['fit_dt_c']:.0f} °C, T ≥ {fit['t_min_c']:.0f} °C, "
          f"n = {fit['n_regression']} of {fit['n_calibration']}; "
          f"blue {fit['n_blue']}, red {fit['n_red']}):")
    print(f"  XFe3_corr = XFe3_iter + {fit['slope']:.4e}·T + ({fit['intercept']:.5f})")
    r = rec["selected_row"]
    print(f"Fig. 6 selected ΔT = {dt} °C  (n = {r['n']}, ±0.02 = {r['p02']:.1%}, "
          f"±0.05 = {r['p05']:.1%}, mean = {r['mean']:+.4f})")

    fp = PARAMS.fe3
    same = (abs(fp.emp_slope - float(f"{fit['slope']:.4e}")) < 1e-12
            and abs(fp.emp_intercept - round(fit["intercept"], 5)) < 1e-9
            and abs(fp.dt_match_c - dt) < 1e-9)
    if args.update_params and not same:
        _update_params(fit["slope"], fit["intercept"], dt, fit["n_regression"])
        print("params.py updated; re-run run_pipeline.py and make_figures.py.")
    elif same:
        print("params.py matches the derived values.")
    else:
        print(f"WARNING: params.py has emp_slope={fp.emp_slope}, "
              f"emp_intercept={fp.emp_intercept}, dt_match_c={fp.dt_match_c}; "
              "run with --update-params.")


if __name__ == "__main__":
    main()
