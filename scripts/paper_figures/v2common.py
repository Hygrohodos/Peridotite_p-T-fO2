"""Revision v3 helpers shared by the v2_fig04-11 scripts.

Calibration (revision v3, option C): no experimental WZ07 correction is applied. The stopping rule is
TA98 - WZ07(Fe2+) <= dt_match_c, a single flat level equal to the median TA98 - WZ07 of the 30 Slave
calibration samples at their measured garnet Fe3+ (calculations/pt_fo2/level_anchor.py; no validation
samples, no outcome-based exclusions). The weighted 13-experiment regression (scripts/weighted_thermometer_regression.py)
is still loaded, but only to draw the experiments in Figs 4b and 5 for contrast.
"""
import os, sys, json
from unified import *
import weighted_thermometer_regression as wtr
from calculations.pt_fo2.level_anchor import anchor_level

FP = PARAMS.fe3
assert not FP.apply_empirical_correction, 'Eq. 5 must be off (revision v2/v3)'
assert FP.wz07_offset_c == 0.0 and FP.wz07_offset_slope == 0.0, 'revision v3: no WZ07 offset or slope'

RUNS = wtr.load_runs()
WT, FIT = wtr.fit(RUNS)
TLO, THI = float(RUNS.T_TA98.min()), float(RUNS.T_TA98.max())

def exp_curve(T):
    """Experimental fit TA98 - WZ07 (deg C), for display only (not applied)."""
    return FIT.predict(np.asarray(T, float))

_ANCHOR = None
def anchor():
    """Anchor computed from the pipeline output (also checks params.py)."""
    global _ANCHOR
    if _ANCHOR is None:
        _ANCHOR = anchor_level(load())
        assert abs(_ANCHOR['dt_level_C'] - FP.dt_match_c) < 0.05, (_ANCHOR['dt_level_C'], FP.dt_match_c)
        json.dump(_ANCHOR, open(os.path.join(OUTD, 'level_anchor.json'), 'w'), indent=1, default=float)
    return _ANCHOR

STATS = {}
def log(key, val):
    STATS[key] = val
    out = os.path.join(OUTD, 'figure_stats_v3.json')
    old = json.load(open(out)) if os.path.exists(out) else {}
    old[key] = val; json.dump(old, open(out, 'w'), indent=1, default=float)

def fe3_at(df, dt):
    return mf._fe3_at_dt(df, dt)          # Eq. 5 off, offset slope and T limits kept

def pct(d, tol):
    return float((d.abs() <= tol + 1e-9).mean())

