"""Revision v2 helpers shared by fig4-7 / figX v2 scripts.

Calibration (Paolo Nimis, lines 265-280 and 320): WZ07 is corrected with the weighted
13-experiment regression from the stand-alone module scripts/weighted_thermometer_regression.py,
TA98 - WZ07 = b0 + b1 (T - 1350), held flat outside the experiments' TA98 range; Eq. 5 is dropped;
the match level dt_match_c is chosen on the 25 Slave samples outside the error groups.
This module checks that the pipeline parameters in params.py are the module's fit (fails loudly if not).
"""
import os, sys, json
from unified import *
import weighted_thermometer_regression as wtr

FP = PARAMS.fe3
TLO, THI = FP.wz07_offset_tmin_c, FP.wz07_offset_tmax_c
assert not FP.apply_empirical_correction, 'Eq. 5 must be off for revision v2'

RUNS = wtr.load_runs()
WT, FIT = wtr.fit(RUNS)
# pipeline: T_WZ07corr = T_WZ07raw - (c + s*T)  <=>  TA98 - WZ07 curve = b0 + b1 (T-1350) = -(c + s*T)
_c = -(FIT.b0 - FIT.b1 * FIT.t_ref); _s = -FIT.b1
assert abs(_c - FP.wz07_offset_c) < 0.01 and abs(_s - FP.wz07_offset_slope) < 1e-5, (_c, _s, FP.wz07_offset_c, FP.wz07_offset_slope)
assert abs(TLO - RUNS.T_TA98.min()) < 0.01 and abs(THI - RUNS.T_TA98.max()) < 0.01, (TLO, THI)

def curve(T):
    """Calibration curve TA98 - WZ07 (deg C), held flat outside [TLO, THI]."""
    return FIT.predict(np.clip(np.asarray(T, float), TLO, THI))

STATS = {}
def log(key, val):
    STATS[key] = val
    out = os.path.join(OUTD, 'figure_stats_v2.json')
    old = json.load(open(out)) if os.path.exists(out) else {}
    old[key] = val; json.dump(old, open(out, 'w'), indent=1, default=float)

def fe3_at(df, dt):
    return mf._fe3_at_dt(df, dt)          # Eq. 5 off, offset slope and T limits kept

def pct(d, tol):
    return float((d.abs() <= tol + 1e-9).mean())

# Display convention (revision v2): the match level is reported as the Slave calibration mean of the stopping
# target, i.e. params.dt_match_c + mean curve over the 25 anchor samples (-157.6 + 65.7 = -91.9 C).
import json as _json
MEANC=_json.load(open(os.path.join(REPO,'data','derived','level_anchor.json')))['mean_curve_C']
