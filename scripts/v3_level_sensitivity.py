#!/usr/bin/env python3
"""Sensitivity of the rev. v3 results to the choice of flat anchor (not used to choose it).
Levels: adopted Slave-30 median, Slave-25 mean (error groups excluded), Slave-30 mean, validation median.
Writes CSV to --out (default stdout)."""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
import make_figures as mf
from calculations.pt_fo2 import calibration as calib
from calculations.pt_fo2.level_anchor import anchor_level
ap = argparse.ArgumentParser(); ap.add_argument("--out"); a_ = ap.parse_args()
df = pd.read_csv(mf.FULL, encoding="utf-8-sig", low_memory=False); df = df.rename(columns={df.columns[0]: "SampleID"})
a = anchor_level(df); s = a["sensitivity"]
levels = [("adopted: Slave-30 median", a["dt_level_C"]), ("Slave-25 mean (error groups excluded)", s["slave25_mean_excl_error_groups_C"]),
          ("Slave-30 mean", s["slave30_mean_C"]), ("validation-27 median (for reference only)", s["validation27_median_C"])]
meas = calib.measured_fe3(df); base = calib.calibration_mask(df) & meas.notna(); val = calib.validation_mask(df) & meas.notna()
sid = df.SampleID.astype(str).str.strip(); rem = base & ~sid.isin(("14-107", "MX044", "MX088", "MX104", "MX5000"))
scr = df["screen_pass"].astype(str).str.lower().isin(["true", "1", "1.0"])
rows = []
for lab, lv in levels:
    est = mf._fe3_at_dt(df, lv)
    r = dict(level=lab, dt_C=round(lv, 2))
    for k, m in [("cal30", base), ("rem25", rem), ("val27", val)]:
        d = (est - meas)[m & est.notna()]
        r[f"{k}_p02"] = float((d.abs() <= 0.02 + 1e-9).mean()); r[f"{k}_p05"] = float((d.abs() <= 0.05 + 1e-9).mean()); r[f"{k}_sd"] = float(d.std(ddof=1))
    e = est[scr]; r["floor_share_450"] = float((e <= 0.0105).mean()); r["n_screened"] = int(scr.sum())
    rows.append(r)
t = pd.DataFrame(rows)
(t.to_csv(a_.out, index=False) if a_.out else print(t.round(3).to_string(index=False)))
