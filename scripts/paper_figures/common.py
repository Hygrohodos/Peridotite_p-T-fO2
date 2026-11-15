"""Helpers shared by the scripts of Figs 4 and 6-10.

WZ07 is corrected by the slope of the weighted experimental fit c(T) of Fig. 5, applied about a pivot
temperature, and garnet Fe3+/SumFe is raised until (Eq. 4)
    T_TA98 - T_WZ07(Fe2+) <= dT_match + c(T_TA98) - c(T_pivot).
The level dT_match and the pivot are set once, in make_paper_figures.py, through set_match_level().
"""
import os, json
from unified import *
from calculations.pt_fo2.level_anchor import curve

FP = PARAMS.fe3
_LEVEL = {}


def set_match_level(level_c, pivot_c):
    """Set dT_match (deg C, at the pivot) and the pivot temperature, and write them into PARAMS.fe3.

    The pipeline stops at T_TA98 - T_WZ07 - c(T) <= FP.dt_match_c, so FP.dt_match_c = level - c(pivot).
    The error groups (Fig. 7) are identified at the same level.
    """
    _LEVEL.update(level=float(level_c), pivot=float(pivot_c))
    FP.dt_match_c = FP.emp_fit_dt_c = to_pipeline(level_c)


def match_level():
    """dT_match of Eq. 4 at the pivot (deg C)."""
    return _LEVEL['level']


def pivot_t():
    """Pivot temperature of the experimental slope (deg C)."""
    return _LEVEL['pivot']


def to_pipeline(level_c):
    """Pipeline stopping value FP.dt_match_c for a level expressed at the pivot."""
    return float(level_c) - float(curve(_LEVEL['pivot']))


def slope_term(T, cr2o3=None):
    """c(T) - c(pivot): the slope of the correction about the pivot (deg C). Revision v7_cr: the TA98 and garnet
    Cr2O3 terms of calculations/pt_fo2/cr_correction.py (cr2o3 per sample, wt%); v7: the experimental slope."""
    return curve(T, cr2o3=cr2o3) - curve(_LEVEL['pivot'])


def log(key, val):
    """Record one figure's statistics in figure_stats.json in the output folder."""
    out = os.path.join(OUTD, 'figure_stats.json')
    stats = json.load(open(out)) if os.path.exists(out) else {}
    stats[key] = val
    json.dump(stats, open(out, 'w'), indent=1, default=float)


def fe3_at(df, level_c):
    """Iterated garnet Fe3+/SumFe at match level level_c (deg C, at the pivot)."""
    return mf._fe3_at_dt(df, to_pipeline(level_c))


def pct(d, tol):
    """Share of the differences d within +-tol."""
    return float((d.abs() <= tol + 1e-9).mean())
