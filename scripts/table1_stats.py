#!/usr/bin/env python3
"""Write Table 1 (data/derived/Table1_calibration_stats.csv) from screened_PT.csv.
Run after run_pipeline.py; then scripts/fo2_uncertainty.py adds the log-unit column."""
import sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from calculations.pt_fo2 import calibration as calib
FULL = ROOT / "data" / "derived" / "screened_PT.csv"
df = pd.read_csv(FULL, encoding="utf-8-sig", low_memory=False); df = df.rename(columns={df.columns[0]: "SampleID"})
tab = calib.calibration_stats(df); out = ROOT / "data" / "derived" / "Table1_calibration_stats.csv"
tab.to_csv(out, index=False); print(tab.drop(columns="notes").round(4).to_string(index=False)); print("Wrote", out)
