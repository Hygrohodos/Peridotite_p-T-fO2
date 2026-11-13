"""
level_anchor.py
---------------
Match level of Eq. 4, anchored to natural xenoliths (not tuned on the Fe3+ success score).

For each sample with measured garnet Fe3+/SumFe, the measured Fe3+ is removed from the olivine-garnet
calculation and D_i = T_TA98,i - T_WZ07,i(measured Fe3+). The weighted experimental slope of c(T)
(Fig. 4; params wz07_offset_*) is removed about a pivot, r_i = D_i - [c(T_TA98,i) - c(pivot)]. Over all
30 screened Slave calibration samples

* the pivot is the MEDIAN of T_TA98 (1208.5 deg C), and
* the match level is the MEDIAN of r (-99.1 deg C),

so the iteration stops at T_TA98 - T_WZ07 - [c(T_TA98) - c(pivot)] <= level. Calibration samples only, no
exclusions (the five error-group samples were identified from Fe3+ residuals), and medians because three
of those five are thermometric outliers. median_r_C and the sensitivity values are on D - c(T); add
c(pivot) to express them at the pivot.
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
    """TA98 - WZ07 correction c(T) applied by the pipeline (the weighted experimental fit)."""
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
    """Return the medians of D, of D - c(T) and of T_TA98 over the 30 Slave samples, with sensitivity values
    on D - c(T) (add c(pivot) to express them at the pivot)."""
    m = anchor_set(df)
    D = discrepancy_measured(df)
    T = _num(df["T_TA98"])
    sid = df["SampleID"].astype(str).str.strip()
    eg = sid.isin(BLUE_IDS + RED_IDS)
    val = calib.validation_mask(df) & D.notna()
    R = D - pd.Series(curve(T, p), index=D.index)
    d = D[m]
    c = pd.Series(curve(T[m], p), index=d.index)
    r = R[m]
    lo, hi = _boot_median_ci(r)
    return dict(
        rule="pivot: median TA98; level: median of TA98 - WZ07(measured Fe3+) with the slope removed; all screened Slave calibration samples",
        n=int(m.sum()),
        median_r_C=float(r.median()), ci95_bootstrap_r_C=[lo, hi], ci95_bootstrap_D_C=list(_boot_median_ci(d)),
        mean_curve_C=float(c.mean()),
        mean_D_C=float(d.mean()), sd_D_C=float(d.std(ddof=1)), median_D_C=float(d.median()),
        mean_T_TA98_C=float(T[m].mean()), median_T_TA98_C=float(T[m].median()),
        sensitivity=dict(
            slave25_mean_excl_error_groups_C=float(R[m & ~eg].mean()),
            slave30_mean_C=float(r.mean()),
            validation27_median_C=float(R[val].median()),
            all52_mean_excl_error_groups_C=float(R[(m | val) & ~eg].mean()),
        ),
        sample_ids=df.loc[m, "SampleID"].astype(str).tolist(),
    )
