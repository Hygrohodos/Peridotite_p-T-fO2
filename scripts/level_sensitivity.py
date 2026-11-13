#!/usr/bin/env python3
"""Sensitivity of the results to the match level (not used to choose it).

Levels, all expressed at the pivot: the adopted level (Slave median of the slope-removed values r) and
the Slave mean, Slave mean with the error groups excluded and validation median of r, and the Slave
median of TA98 - WZ07 before slope removal. Each is applied with the WZ07 correction in PARAMS.fe3;
make_paper_figures.py calls level_sensitivity() after setting it."""
import sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
import make_figures as mf
from calculations.pt_fo2 import calibration as calib
from calculations.pt_fo2.level_anchor import anchor_level


def level_sensitivity(df, level, to_pipeline):
    """One row per level: +-0.02 / +-0.05 success and 1 sigma on the Slave calibration set (cal), the Slave set
    without the error groups (rem) and the validation set (val), and the share of the screened samples left at
    the 0.01 starting value. to_pipeline converts a level at the pivot to the pipeline stopping value."""
    a = anchor_level(df); s = a["sensitivity"]; shift = -to_pipeline(0.0)
    levels = [("adopted level (Slave median of r)", level),
              ("Slave mean of r, error groups excluded", s["slave_mean_excl_error_groups_C"] + shift),
              ("Slave mean of r", s["slave_mean_C"] + shift), ("validation median of r", s["validation_median_C"] + shift),
              ("Slave median before slope removal", a["median_D_C"])]
    base, meas, blue, red, _, _ = calib.error_groups(df)
    base = base & meas.notna(); val = calib.validation_mask(df) & meas.notna()
    rem = base & ~(blue | red)
    screened = df["screen_pass"].astype(str).str.lower().isin(["true", "1", "1.0"])
    rows = []
    for label, lv in levels:
        est = mf._fe3_at_dt(df, to_pipeline(lv))
        r = dict(level=label, dt_C=round(float(lv), 2))
        for key, m in (("cal", base), ("rem", rem), ("val", val)):
            d = (est - meas)[m & est.notna()]
            r[f"{key}_n"] = int(len(d))
            r[f"{key}_p02"] = float((d.abs() <= 0.02 + 1e-9).mean())
            r[f"{key}_p05"] = float((d.abs() <= 0.05 + 1e-9).mean())
            r[f"{key}_sd"] = float(d.std(ddof=1))
        r["floor_share_screened"] = float((est[screened] <= 0.0105).mean()); r["n_screened"] = int(screened.sum())
        rows.append(r)
    return pd.DataFrame(rows)
