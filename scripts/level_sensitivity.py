#!/usr/bin/env python3
"""Sensitivity of the results to the match level (not used to choose it).

Levels, all expressed at the pivot: the adopted level (Slave-30 median of the slope-removed values r) and
the Slave-30 mean, Slave-25 mean (error groups excluded) and validation median of r, and the Slave-30
median of TA98 - WZ07 before slope removal. Each is applied with the WZ07 correction in PARAMS.fe3; make_paper_figures.py calls
level_sensitivity() after setting it."""
import sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
import make_figures as mf
from calculations.pt_fo2 import calibration as calib
from calculations.pt_fo2.level_anchor import anchor_level, BLUE_IDS, RED_IDS


def level_sensitivity(df, level, to_pipeline):
    """One row per level: +-0.02 / +-0.05 success and 1 sigma on the Slave-30, Slave-25 and validation sets,
    and the share of the screened samples left at the 0.01 starting value. to_pipeline converts a level at the
    pivot to the pipeline stopping value."""
    a = anchor_level(df); s = a["sensitivity"]; shift = -to_pipeline(0.0)
    levels = [("adopted level (Slave-30 median of r)", level),
              ("Slave-25 mean of r (error groups excluded)", s["slave25_mean_excl_error_groups_C"] + shift),
              ("Slave-30 mean of r", s["slave30_mean_C"] + shift), ("validation-27 median of r", s["validation27_median_C"] + shift),
              ("Slave-30 median before slope removal", a["median_D_C"])]
    meas = calib.measured_fe3(df)
    base = calib.calibration_mask(df) & meas.notna(); val = calib.validation_mask(df) & meas.notna()
    rem = base & ~df.SampleID.astype(str).str.strip().isin(BLUE_IDS + RED_IDS)
    screened = df["screen_pass"].astype(str).str.lower().isin(["true", "1", "1.0"])
    rows = []
    for label, lv in levels:
        est = mf._fe3_at_dt(df, to_pipeline(lv))
        r = dict(level=label, dt_C=round(float(lv), 2))
        for key, m in (("cal30", base), ("rem25", rem), ("val27", val)):
            d = (est - meas)[m & est.notna()]
            r[f"{key}_p02"] = float((d.abs() <= 0.02 + 1e-9).mean())
            r[f"{key}_p05"] = float((d.abs() <= 0.05 + 1e-9).mean())
            r[f"{key}_sd"] = float(d.std(ddof=1))
        r["floor_share_screened"] = float((est[screened] <= 0.0105).mean()); r["n_screened"] = int(screened.sum())
        rows.append(r)
    return pd.DataFrame(rows)
