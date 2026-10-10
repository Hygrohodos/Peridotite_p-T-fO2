"""
calibration.py
--------------
Reproducible derivation of the locked comparative-geothermometry
parameters (manuscript Eq. 5, Fig. 5, Fig. 6 and Table 1).

Workflow (``python scripts/fit_calibration.py`` runs all of it):

1. Calibration set: Slave craton samples that pass every screen and have
   measured garnet Fe3+/ΣFe with 0 < Fe3+/ΣFe < 0.5 (n = 31).
2. Uncorrected estimates at ΔT = ``emp_fit_dt_c`` (25 °C; Fig. 5b).
3. Error groups at that ΔT (Fig. 5b):
     blue  measured − estimated > ``emp_blue_threshold`` (0.05)
     red   sample IDs in ``emp_red_ids``: the samples at T > 1250 °C with
           estimated > measured by > 0.05 (21-6, +0.054 at 1136 °C, meets
           the 0.05 criterion but lies on the high-T trend and is kept in
           the regression set)
4. Regression samples = calibration set minus blue and red, with
   T_TA98 ≥ ``t_emp_min_c`` (1130 °C).  These are the samples inside the
   dashed envelope of Fig. 5b.
5. Eq. 5 = ordinary least-squares line through
   (measured − estimated) versus T_TA98 for the regression samples.  This
   is exactly the solid line drawn in Fig. 5b (no level shift):
       XFe3_corr = XFe3_iter + emp_slope · T + emp_intercept   (T ≥ 1130 °C)
6. Fig. 6: with Eq. 5 held fixed, sweep the ΔT match threshold
   (0–75 °C) and score only the regression samples.  The locked ΔT is
   the value with the highest ±0.02 success rate; ties are broken by the
   smallest |mean(estimated − measured)|.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .params import PARAMS, CalculationParams, Fe3Params


VALIDATION_CRATONS = ("Siberian", "Kaapvaal")


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def measured_fe3(df):
    return _num(df["garnet_Fe3+"])


def _has_measured(df, screened=True):
    meas = measured_fe3(df)
    m = meas.notna() & (meas > 0) & (meas < 0.5)
    if screened:
        m = m & (df["screen_pass"] == True)
    return m


def calibration_mask(df, screened=True):
    """Slave samples with measured garnet Fe3+/ΣFe (screened by default)."""
    slave = df["Craton"].fillna("").astype(str).str.contains("Slave")
    return _has_measured(df, screened) & slave


def validation_mask(df, screened=True):
    """Siberian + Kaapvaal samples with non-zero measured garnet Fe3+/ΣFe."""
    cr = df["Craton"].fillna("").astype(str)
    val = cr.str.contains("|".join(VALIDATION_CRATONS))
    return _has_measured(df, screened) & val


def fe3_uncorrected(df, dt, params=None):
    """Iterated garnet Fe3+/ΣFe at match threshold ``dt`` (no Eq. 5)."""
    from dataclasses import replace
    from .fe3_correction import iterate_garnet_fe3

    p = params or PARAMS
    fp = replace(p.fe3, dt_match_c=float(dt), apply_empirical_correction=False)
    out = iterate_garnet_fe3(
        df, verbose=False,
        params=CalculationParams(screening=p.screening, fe3=fp, export=p.export),
    )
    return _num(out["XFe3_gt_iter"])


def apply_eq5(est, T, slope, intercept, t_min=None):
    """Eq. 5 applied to a Series of uncorrected estimates (clipped 0–1)."""
    t_min = PARAMS.fe3.t_emp_min_c if t_min is None else t_min
    est = _num(est)
    T = _num(T)
    out = np.where(T >= t_min, est + slope * T + intercept, est)
    return pd.Series(np.clip(out, 0.0, 1.0), index=est.index)


def error_groups(df, params=None, est_fit=None):
    """Return (base, meas, blue, red, remainder, regression) boolean masks."""
    p = params or PARAMS
    base = calibration_mask(df)
    meas = measured_fe3(df)
    if est_fit is None:
        est_fit = fe3_uncorrected(df, p.fe3.emp_fit_dt_c, p)
    T = _num(df["T_TA98"])
    blue = base & ((meas - est_fit) > p.fe3.emp_blue_threshold)
    sid = df["SampleID"].astype(str).str.strip()
    red = base & sid.isin(set(p.fe3.emp_red_ids))
    rem = base & ~blue & ~red
    reg = rem & T.notna() & (T >= p.fe3.t_emp_min_c)
    return base, meas, blue, red, rem, reg


def fit_empirical_correction(df, params=None):
    """Fit Eq. 5 exactly as drawn in Fig. 5b.  Returns a dict."""
    p = params or PARAMS
    est = fe3_uncorrected(df, p.fe3.emp_fit_dt_c, p)
    base, meas, blue, red, rem, reg = error_groups(df, p, est_fit=est)
    T = _num(df["T_TA98"])
    x = T[reg].to_numpy(dtype=float)
    y = (meas - est)[reg].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    return dict(
        fit_dt_c=float(p.fe3.emp_fit_dt_c),
        t_min_c=float(p.fe3.t_emp_min_c),
        slope=float(slope),
        intercept=float(intercept),
        n_calibration=int(base.sum()),
        n_blue=int(blue.sum()),
        n_red=int(red.sum()),
        n_regression=int(reg.sum()),
        regression_ids=sorted(df.loc[reg, "SampleID"].astype(str).str.strip()),
        est_fit=est,
        masks=dict(base=base, blue=blue, red=red, rem=rem, reg=reg),
    )


def dt_sweep(df, slope, intercept, dts=range(0, 76), params=None, mask=None):
    """Fig. 6 sweep: Eq. 5 held fixed, scored on ``mask`` (regression set)."""
    p = params or PARAMS
    meas = measured_fe3(df)
    T = _num(df["T_TA98"])
    if mask is None:
        mask = error_groups(df, p)[5]
    rows = []
    for dt in dts:
        est = apply_eq5(fe3_uncorrected(df, dt, p), T, slope, intercept,
                        p.fe3.t_emp_min_c)
        m = mask & est.notna()
        d = est[m] - meas[m]
        rows.append(dict(
            dt=int(dt), n=int(m.sum()), slope=slope, intercept=intercept,
            mean=float(d.mean()), absmean=float(d.abs().mean()),
            p02=float((d.abs() <= 0.02 + 1e-12).mean()),
            p05=float((d.abs() <= 0.05 + 1e-12).mean()),
        ))
    return pd.DataFrame(rows)


def select_dt(sweep):
    """Highest ±0.02 rate; ties broken by smallest |mean|, then ±0.05."""
    s = sweep.copy()
    s["_absbias"] = s["mean"].abs().round(6)
    s = s.sort_values(["p02", "_absbias", "p05"], ascending=[False, True, False])
    return int(s.iloc[0]["dt"])


def _row(name, d, notes=""):
    d = d.dropna()
    return dict(
        dataset=name, n=int(len(d)),
        p02=float((d.abs() <= 0.02 + 1e-12).mean()),
        p05=float((d.abs() <= 0.05 + 1e-12).mean()),
        sigma=float(d.std(ddof=1)),
        mean=float(d.mean()),
        notes=notes,
    )


def calibration_stats(df, est=None, params=None, fo2_per_002=None):
    """Table 1 rows.  ``est`` defaults to the pipeline XFe3_gt_corr column.

    Residual 1σ is the sample standard deviation (ddof = 1) of
    estimated − measured.  ``sigma_log`` converts 1σ (rounded to 3 dp)
    to log units with ``fo2_per_002`` log units per 0.02 Fe3+/ΣFe
    (0.44 in the manuscript).  Rates count |est − meas| ≤ 0.02 / 0.05.
    """
    p = params or PARAMS
    meas = measured_fe3(df)
    if est is None:
        est = _num(df["XFe3_gt_corr"])
    base, _, blue, red, rem, reg = error_groups(df, p)
    val = validation_mask(df)
    groups = blue | red
    rows = [
        _row("Slave calibration", (est - meas)[base]),
        _row("Slave, error groups excluded", (est - meas)[base & ~groups]),
        _row("Siberia + Kaapvaal validation", (est - meas)[val]),
        _row("Calibration + validation", (est - meas)[base | val]),
        _row("Combined, error groups excluded", (est - meas)[(base | val) & ~groups]),
        _row("Eq. 5 regression samples", (est - meas)[reg]),
    ]
    all_slave = calibration_mask(df, screened=False)
    for col, nm in (("XFe3_quinn", "Quinn et al. (2016) stoichiometry"),
                    ("XFe3_droop", "Droop (1987) stoichiometry")):
        if col in df.columns:
            s = _num(df[col])
            rows.append(_row(nm, (s - meas)[all_slave]))
            rows.append(_row(nm + ", screened", (s - meas)[base]))
    out = pd.DataFrame(rows)
    if fo2_per_002 is not None:
        # Converted from 1σ rounded to 3 dp (as reported), so quoted
        # pairs such as "0.028 (0.62 log units)" are reproducible.
        out["sigma_log"] = (out["sigma"].round(3) * fo2_per_002 / 0.02).round(2)
    return out
