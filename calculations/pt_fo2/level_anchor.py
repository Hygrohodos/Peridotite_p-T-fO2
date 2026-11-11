"""
level_anchor.py
---------------
Revision v3 (option C): the WZ07 calibration correction is a FLAT level with
no temperature slope, ANCHORED to natural xenoliths (not tuned on the Fe3+
success score).

For each sample with measured garnet Fe3+/SumFe, the measured Fe3+ is removed
from the olivine-garnet calculation and the remaining two-pyroxene minus
olivine-garnet difference is

    D_i = T_TA98,i - T_WZ07,i(measured Fe3+)

The adopted match level is the MEDIAN of D over all 30 screened Slave
calibration samples:

* calibration samples only, so the Siberia + Kaapvaal validation samples play
  no part in setting it;
* no exclusions: the five "error-group" samples were identified from Fe3+
  residuals, so excluding them would make the anchor depend on the outcome;
* the median rather than the mean, because three of those five samples are
  thermometric outliers (D = -188 to -304 deg C) that would pull a mean.

Sensitivity values (not used) are reported alongside: the Slave-25 mean
(error groups excluded, the rev. v2 anchor), the Slave-30 mean, and the
validation median.

The iteration then stops when TA98 - WZ07(Fe2+) <= level (params
dt_match_c), with no offset and no slope applied to WZ07.

    python scripts/anchor_level.py [--out PATH]
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .params import PARAMS
from . import calibration as calib
from .fe3_correction import _wz07corr_fe2only

# The five error-group samples of Section 3.7 (fixed by ID; reported only for
# the sensitivity statistics, NOT excluded from the anchor).
BLUE_IDS = ("14-107", "MX044", "MX088")
RED_IDS = ("MX104", "MX5000")


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def curve(T, p=None):
    """TA98 - WZ07 correction c(T) applied by the pipeline (zero in rev. v3)."""
    fp = (p or PARAMS).fe3
    T = np.asarray(T, float)
    lo = -np.inf if fp.wz07_offset_tmin_c is None else fp.wz07_offset_tmin_c
    hi = np.inf if fp.wz07_offset_tmax_c is None else fp.wz07_offset_tmax_c
    Tc = np.clip(T, lo, hi)
    return -(fp.wz07_offset_c + fp.wz07_offset_slope * Tc)


def anchor_set(df):
    """Boolean mask of the anchor samples: all 30 screened Slave calibration samples."""
    return calib.calibration_mask(df) & calib.measured_fe3(df).notna()


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


def _boot_median_ci(x, n=4000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    b = np.median(rng.choice(x, size=(n, len(x)), replace=True), axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def anchor_level(df, p=None):
    """Return a dict with the anchored flat match level and sensitivity values."""
    m = anchor_set(df)
    D = discrepancy_measured(df)
    T = _num(df["T_TA98"])
    sid = df["SampleID"].astype(str).str.strip()
    eg = sid.isin(BLUE_IDS + RED_IDS)
    val = calib.validation_mask(df) & D.notna()
    d = D[m]
    c = pd.Series(curve(T[m], p), index=d.index)
    r = d - c
    lo, hi = _boot_median_ci(r)
    return dict(
        rule="median of TA98 - WZ07(measured Fe3+) over all screened Slave calibration samples",
        n=int(m.sum()),
        dt_level_C=float(r.median()), ci95_bootstrap_C=[lo, hi],
        mean_curve_C=float(c.mean()),
        mean_D_C=float(d.mean()), sd_D_C=float(d.std(ddof=1)), median_D_C=float(d.median()),
        mean_T_TA98_C=float(T[m].mean()),
        sensitivity=dict(
            slave25_mean_excl_error_groups_C=float(D[m & ~eg].mean()),
            slave30_mean_C=float(d.mean()),
            validation27_median_C=float(D[val].median()),
            all52_mean_excl_error_groups_C=float(D[(m | val) & ~eg].mean()),
        ),
        sample_ids=df.loc[m, "SampleID"].astype(str).tolist(),
    )
