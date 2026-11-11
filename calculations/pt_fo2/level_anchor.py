"""
level_anchor.py
---------------
Revision v2 (option B, Stephanie Greene 11 Nov 2026): the level of the WZ07
calibration correction is ANCHORED to the natural Slave calibration samples,
not tuned on the Fe3+ success score.

Two parts make up the correction applied to raw WZ07 (Eq. 4 of the revision):

1. Temperature dependence (slope) from the experiments: the weighted
   13-experiment regression of scripts/weighted_thermometer_regression.py,
       c(T) = b0 + b1 (T - 1350),  held flat outside 1162-1620 deg C
   (pipeline: params.Fe3Params.wz07_offset_c / wz07_offset_slope / T limits).

2. Level from the natural samples: the 25 Slave calibration samples outside
   the blue and red error groups, with their MEASURED garnet Fe3+ removed
   from the olivine-garnet calculation, give
       D_i = T_TA98,i - T_WZ07,i(measured Fe3+)    (mean -91.9 deg C)
   The match level is the mean residual of D about the experimental curve,
       dT_level = mean_i [ D_i - c(T_i) ]
   so that, on average, the iteration stops at the measured Fe3+ of the
   anchor set. Nothing is optimised against the +-0.02 / +-0.05 scores, and
   the Siberia + Kaapvaal validation samples are not used.

Equivalently: subtract the natural mean discrepancy (-91.9 deg C at the
anchor set's mean TA98) and let the experimental slope carry the
temperature dependence.

    python scripts/anchor_level.py      # prints and writes data/derived/level_anchor.json
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .params import PARAMS
from . import calibration as calib
from .fe3_correction import _wz07corr_fe2only

# Error-group membership is fixed by sample ID (the five samples of the
# published Fig. 5b / Section 3.7), so the anchor set does not depend on the
# level being computed.
BLUE_IDS = ("14-107", "MX044", "MX088")
RED_IDS = ("MX104", "MX5000")


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def curve(T, p=None):
    """Experimental TA98 - WZ07 curve c(T) (deg C) as applied by the pipeline."""
    fp = (p or PARAMS).fe3
    T = np.asarray(T, float)
    lo = -np.inf if fp.wz07_offset_tmin_c is None else fp.wz07_offset_tmin_c
    hi = np.inf if fp.wz07_offset_tmax_c is None else fp.wz07_offset_tmax_c
    Tc = np.clip(T, lo, hi)
    return -(fp.wz07_offset_c + fp.wz07_offset_slope * Tc)


def anchor_set(df):
    """Boolean mask of the 25 Slave calibration samples outside the error groups."""
    base = calib.calibration_mask(df)
    sid = df["SampleID"].astype(str).str.strip()
    return base & ~sid.isin(BLUE_IDS + RED_IDS)


def discrepancy_measured(df):
    """D = TA98 - WZ07 with the measured garnet Fe3+/SumFe removed (deg C)."""
    meas = calib.measured_fe3(df)
    gFe = _num(df["garnet_Fe"]).to_numpy(float)
    T_wz = _wz07corr_fe2only(
        _num(df["olivine_Mg"]).to_numpy(float), _num(df["olivine_Fe"]).to_numpy(float),
        _num(df["garnet_Mg"]).to_numpy(float), gFe * (1.0 - meas.fillna(0).to_numpy(float)),
        _num(df["garnet_Ca"]).to_numpy(float), _num(df["garnet_Mn"]).to_numpy(float),
        _num(df["P_CA91"]).to_numpy(float), offset_c=0.0,
    )
    D = _num(df["T_TA98"]).to_numpy(float) - T_wz
    return pd.Series(np.where(meas.notna(), D, np.nan), index=df.index)


def anchor_level(df, p=None):
    """Return a dict with the anchored match level and its ingredients."""
    m = anchor_set(df)
    D = discrepancy_measured(df)[m]
    T = _num(df["T_TA98"])[m]
    c = pd.Series(curve(T, p), index=T.index)
    r = D - c
    return dict(
        n=int(m.sum()),
        mean_D_C=float(D.mean()), sd_D_C=float(D.std(ddof=1)),
        mean_T_TA98_C=float(T.mean()),
        mean_curve_C=float(c.mean()),
        dt_level_C=float(r.mean()), se_dt_level_C=float(r.std(ddof=1) / np.sqrt(len(r))),
        median_residual_C=float(r.median()),
        sample_ids=df.loc[m, "SampleID"].astype(str).tolist(),
    )
