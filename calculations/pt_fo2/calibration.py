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
7. fO2 uncertainty (Table 1 log column, Fig. 10 bars): for each sample
   with measured Fe3+/ΣFe, ΔFMQ is calculated (Stagno et al. 2013, same
   P, T and chemistry as the pipeline) from the estimated and from the
   measured garnet Fe3+/ΣFe; 1σ is the standard deviation of
   ΔFMQ(estimated) − ΔFMQ(measured).  For Fig. 10 the 58 calibration +
   validation samples are split at their median estimated ΔFMQ into a
   reduced and an oxidised half (``fo2_uncertainty``).
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


def fo2_from_xfe3(df, xfe3):
    """Stagno et al. (2013) ΔFMQ for ``df`` with garnet Fe3+/ΣFe = ``xfe3``.

    Total garnet Fe, P, T and all other inputs are those of the pipeline;
    only the Fe3+/Fe2+ split is replaced.
    """
    from .fo2 import calculate_fo2

    tmp = df.copy()
    fe_tot = _num(tmp.get("garnet_Fe"))
    if "Fe3_gt_corr" in tmp.columns and "Fe2_gt_corr" in tmp.columns:
        fe_tot = fe_tot.fillna(_num(tmp["Fe3_gt_corr"]) + _num(tmp["Fe2_gt_corr"]))
    x = _num(pd.Series(np.asarray(xfe3, dtype=float), index=tmp.index))
    tmp["Fe3_gt_corr"] = x * fe_tot
    tmp["Fe2_gt_corr"] = (1.0 - x) * fe_tot
    return _num(calculate_fo2(tmp, verbose=False)["delta_logfO2_FMQ"])


def fo2_residuals(df, est=None, params=None):
    """Per-sample fO2 residuals for the calibration + validation samples.

    Returns a DataFrame (index = df index) with measured and estimated
    Fe3+/ΣFe, ΔFMQ from each, ``res_log`` = ΔFMQ(est) − ΔFMQ(meas), the
    local sensitivity (log units per 0.02 Fe3+/ΣFe at the measured value),
    and flags for the dataset and the Fig. 5 error groups.
    """
    p = params or PARAMS
    base, meas, blue, red, _, _ = error_groups(df, p)
    val = validation_mask(df)
    pool = base | val
    if est is None:
        est = _num(df["XFe3_gt_corr"])
    u = df[pool]
    m, e = meas[pool], _num(est)[pool]
    f_m, f_e = fo2_from_xfe3(u, m), fo2_from_xfe3(u, e)
    h = 0.005
    sens = (fo2_from_xfe3(u, m + h) - fo2_from_xfe3(u, (m - h).clip(lower=1e-4))) / (2 * h) * 0.02
    return pd.DataFrame({
        "SampleID": u["SampleID"], "Craton": u["Craton"],
        "dataset": np.where(base[pool], "Slave calibration", "Validation"),
        "error_group": (blue | red)[pool],
        "XFe3_meas": m, "XFe3_est": e,
        "dFMQ_meas": f_m, "dFMQ_est": f_e, "res_log": f_e - f_m,
        "sens_per_002": sens,
    })


def fo2_uncertainty(res, split=None):
    """fO2 1σ (log units) overall and for reduced / oxidised halves.

    ``res`` is the output of ``fo2_residuals``.  All samples are used
    (error groups retained, the conservative case).  The split defaults
    to the median estimated ΔFMQ; reduced = ΔFMQ(est) < split.
    """
    r = res.dropna(subset=["res_log", "dFMQ_est"])
    cut = float(r["dFMQ_est"].median()) if split is None else float(split)
    lo = r["dFMQ_est"] < cut
    def _s(d):
        return dict(n=int(len(d)), sigma_log=float(d.std(ddof=1)),
                    rms_log=float(np.sqrt((d ** 2).mean())), bias_log=float(d.mean()))
    return dict(split=cut, all=_s(r["res_log"]),
                reduced=_s(r.loc[lo, "res_log"]), oxidised=_s(r.loc[~lo, "res_log"]))


def calibration_stats(df, est=None, params=None, fo2=True):
    """Table 1 rows.  ``est`` defaults to the pipeline XFe3_gt_corr column.

    Residual 1σ is the sample standard deviation (ddof = 1) of
    estimated − measured Fe3+/ΣFe.  With ``fo2=True``, ``sigma_log`` is
    the 1σ of ΔFMQ(estimated) − ΔFMQ(measured) for the same samples
    (``fo2_residuals``); it is left blank for the charge-balance rows.
    Rates count |est − meas| ≤ 0.02 / 0.05.
    """
    p = params or PARAMS
    meas = measured_fe3(df)
    if est is None:
        est = _num(df["XFe3_gt_corr"])
    base, _, blue, red, rem, reg = error_groups(df, p)
    val = validation_mask(df)
    groups = blue | red
    sets = [
        ("Slave calibration", base),
        ("Slave, error groups excluded", base & ~groups),
        ("Siberia + Kaapvaal validation", val),
        ("Calibration + validation", base | val),
        ("Combined, error groups excluded", (base | val) & ~groups),
        ("Eq. 5 regression samples", reg),
    ]
    rows = [_row(nm, (est - meas)[m]) for nm, m in sets]
    if fo2:
        res = fo2_residuals(df, est, p)["res_log"]
        for row, (_, m) in zip(rows, sets):
            d = res.reindex(df.index[m]).dropna()
            row["sigma_log"] = float(d.std(ddof=1)) if len(d) > 1 else np.nan
    # Charge-balance stoichiometry (Droop 1987 = Quinn et al. 2016; see
    # stoichiometry.py).  The two columns are identical to rounding, so one
    # row per sample set is reported: screened 31 first, then all 57 Slave
    # samples with measured Fe3+/ΣFe before screening.
    if "XFe3_quinn" in df.columns:
        cb = _num(df["XFe3_quinn"])
        if "XFe3_droop" in df.columns:
            dr = _num(df["XFe3_droop"])
            both = cb.notna() & dr.notna()
            if not np.allclose(cb[both], dr[both], rtol=0, atol=1e-9):
                raise ValueError("XFe3_droop and XFe3_quinn differ; "
                                 "check stoichiometry.py")
        all_slave = calibration_mask(df, screened=False)
        nm = "Charge balance (Droop 1987; Quinn et al. 2016)"
        rows.append(_row(nm + ", screened Slave", (cb - meas)[base],
                         f"Same {int(base.sum())} screened Slave samples as the calibration row"))
        rows.append(_row(nm + ", all Slave", (cb - meas)[all_slave],
                         f"All {int(all_slave.sum())} Slave samples with measured Fe3+/ΣFe, before screening"))
    return pd.DataFrame(rows)
