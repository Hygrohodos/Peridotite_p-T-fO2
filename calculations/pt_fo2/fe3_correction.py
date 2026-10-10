"""
fe3_correction.py
-----------------
Garnet Fe³⁺ estimation by comparative geothermometry, translated from
the "Comparative Geothermometry" sheet of SD1 - P-T-fO2.xlsm.

Procedure
---------
1. Assume olivine is entirely Fe²⁺.
2. Start with garnet XFe³⁺ = params.fe3.xfe3_start (0.01; XFe²⁺ = 0.99).
3. Recalculate the Wu & Zhao (2007) olivine–garnet thermometer using
   only Fe²⁺ in both phases and subtract the flat calibration offset
   (WZ07corr = WZ07raw − wz07_offset_c, 118 °C from experiment F1300o).
4. While (T_TA98 − T_WZ07corr) > dt_match_c (29 °C), increase XFe³⁺ by
   dxfe3 (0.001) and recompute.
5. Once converged, apply the empirical high-temperature correction
   (manuscript Eq. 5):
       if T < t_emp_min_c (1130 °C):  XFe³⁺_corr = XFe³⁺_iterated
       else:  XFe³⁺_corr = XFe³⁺_iterated + emp_slope·T + emp_intercept
   with emp_slope / emp_intercept from params.py (3.8319e-4, −0.46673),
   derived by calibration.fit_empirical_correction (Fig. 5b regression at
   ΔT = 25 °C); clipped to 0–1.
6. Convert the corrected fraction back to cation units.

All numbers above are defaults; params.py is authoritative.

Only samples that possess the necessary cation columns and a valid
TA98/CA91 solution are iterated.  Final corrected values are reported
for all such samples; downstream fO₂ code can further restrict to
screen_pass == True if desired.
"""

import numpy as np
import pandas as pd

from .params import PARAMS


# ---------------------------------------------------------------------------
# Internal helper – WZ07corr using Fe²⁺-only site fractions
# ---------------------------------------------------------------------------
def _wz07corr_fe2only(Mg_ol, Fe2_ol, Mg_gt, Fe2_gt, Ca_gt, Mn_gt, P_kbar, offset_c=None):
    """
    Wu & Zhao (2007) thermometer evaluated with Fe²⁺ only, minus the
    flat calibration offset (params.fe3.wz07_offset_c, 118 °C).

    All arguments are 1-D numpy arrays of the same length.
    Returns T in °C.
    """
    eps = 1e-16
    XOlMg = Mg_ol / (Mg_ol + Fe2_ol + eps)
    XOlFe = Fe2_ol / (Mg_ol + Fe2_ol + eps)

    sum_gt = Mg_gt + Fe2_gt + Ca_gt + Mn_gt + eps
    XGtMg = Mg_gt / sum_gt
    XGtCa = Ca_gt / sum_gt
    XGtFe = Fe2_gt / sum_gt
    XGtMn = Mn_gt / sum_gt

    R = 8.3143

    # Fea (BD)
    AG = (
        5.993 * XGtMg**2 - 9.67 * XGtCa**2
        - 11.006 * XGtFe * XGtMg + 0.674 * XGtFe * XGtCa
        - 9.5725 * XGtMg * XGtCa - 22.765 * XGtMg * XGtMn
        - 4.6665 * XGtCa * XGtMn
        + 11.006 * XGtFe**2 * XGtMg - 0.674 * XGtFe**2 * XGtCa
        - 11.986 * XGtMg**2 * XGtFe + 37.966 * XGtMg**2 * XGtCa
        + 46.02 * XGtMg**2 * XGtMn
        + 19.34 * XGtCa**2 * XGtFe + 25.746 * XGtCa**2 * XGtMg
        + 46.02 * XGtMn**2 * XGtMg
        + 19.145 * XGtFe * XGtMg * XGtCa + 45.53 * XGtFe * XGtMg * XGtMn
        + 9.333 * XGtFe * XGtCa * XGtMn + 77.876 * XGtMg * XGtCa * XGtMn
    )
    # Feb (BE)
    AH = (
        -0.034 * XGtMg**2 + 0.135 * XGtCa**2 + 0.024 * XGtMn**2
        + 0.1 * XGtFe * XGtMg + 0.08 * XGtFe * XGtCa + 0.048 * XGtFe * XGtMn
        - 0.0515 * XGtMg * XGtCa + 0.07 * XGtMg * XGtMn + 0.1765 * XGtCa * XGtMn
        - 0.1 * XGtFe**2 * XGtMg - 0.08 * XGtFe**2 * XGtCa - 0.048 * XGtFe**2 * XGtMn
        + 0.068 * XGtMg**2 * XGtFe - 0.136 * XGtMg**2 * XGtCa - 0.076 * XGtMg**2 * XGtMn
        - 0.27 * XGtCa**2 * XGtFe - 0.28 * XGtCa**2 * XGtMg - 0.13 * XGtCa**2 * XGtMn
        - 0.048 * XGtMn**2 * XGtFe - 0.076 * XGtMn**2 * XGtMg - 0.13 * XGtMn**2 * XGtCa
        + 0.103 * XGtFe * XGtMg * XGtCa - 0.14 * XGtFe * XGtMg * XGtMn
        - 0.353 * XGtFe * XGtCa * XGtMn - 0.414 * XGtMg * XGtCa * XGtMn
    )
    # Fec (BF)
    AI = (
        -5672 * XGtMg**2 + 19932 * XGtCa**2 + 1617 * XGtMn**2
        + 23244 * XGtFe * XGtMg - 2608 * XGtFe * XGtCa + 3234 * XGtFe * XGtMn
        + 31326.5 * XGtMg * XGtCa + 45841 * XGtMg * XGtMn + 12356 * XGtCa * XGtMn
        - 23244 * XGtFe**2 * XGtMg + 2608 * XGtFe**2 * XGtCa - 3234 * XGtFe**2 * XGtMn
        + 11344 * XGtMg**2 * XGtFe - 132228 * XGtMg**2 * XGtCa - 82498 * XGtMg**2 * XGtMn
        - 39864 * XGtCa**2 * XGtFe - 51518 * XGtCa**2 * XGtMg - 2850 * XGtCa**2 * XGtMn
        - 3234 * XGtMn**2 * XGtFe - 82498 * XGtMn**2 * XGtMg - 2850 * XGtMn**2 * XGtCa
        - 62653 * XGtFe * XGtMg * XGtCa - 91682 * XGtFe * XGtMg * XGtMn
        - 24712 * XGtFe * XGtCa * XGtMn - 177221 * XGtMg * XGtCa * XGtMn
    )
    # Mga (BG)
    AJ = (
        -5.503 * XGtFe**2 - 12.873 * XGtCa**2 - 23.01 * XGtMn**2
        + 11.986 * XGtFe * XGtMg - 9.5725 * XGtFe * XGtCa - 22.765 * XGtFe * XGtMn
        - 37.966 * XGtMg * XGtCa - 46.02 * XGtMg * XGtMn - 38.938 * XGtCa * XGtMn
        + 11.006 * XGtFe**2 * XGtMg - 0.674 * XGtFe**2 * XGtCa
        - 11.986 * XGtMg**2 * XGtFe + 37.966 * XGtMg**2 * XGtCa + 46.02 * XGtMg**2 * XGtMn
        + 19.34 * XGtCa**2 * XGtFe + 25.746 * XGtCa**2 * XGtMg
        + 46.02 * XGtMn**2 * XGtMg
        + 19.145 * XGtFe * XGtMg * XGtCa + 45.53 * XGtFe * XGtMg * XGtMn
        + 9.333 * XGtFe * XGtCa * XGtMn + 77.876 * XGtMg * XGtCa * XGtMn
    )
    # Mgb (BH)
    AK = (
        0.05 * XGtFe**2 + 0.14 * XGtCa**2 + 0.038 * XGtMn**2
        - 0.068 * XGtFe * XGtMg - 0.0515 * XGtFe * XGtCa + 0.07 * XGtFe * XGtMn
        + 0.136 * XGtMg * XGtCa + 0.076 * XGtMg * XGtMn + 0.207 * XGtCa * XGtMn
        - 0.1 * XGtFe**2 * XGtMg - 0.08 * XGtFe**2 * XGtCa - 0.048 * XGtFe**2 * XGtMn
        + 0.068 * XGtMg**2 * XGtFe - 0.136 * XGtMg**2 * XGtCa - 0.076 * XGtMg**2 * XGtMn
        - 0.27 * XGtCa**2 * XGtFe - 0.28 * XGtCa**2 * XGtMg - 0.13 * XGtCa**2 * XGtMn
        - 0.048 * XGtMn**2 * XGtFe - 0.076 * XGtMn**2 * XGtMg - 0.13 * XGtMn**2 * XGtCa
        + 0.103 * XGtFe * XGtMg * XGtCa - 0.14 * XGtFe * XGtMg * XGtMn
        - 0.353 * XGtFe * XGtCa * XGtMn - 0.414 * XGtMg * XGtCa * XGtMn
    )
    # Mgc (BI)
    AL = (
        11622 * XGtFe**2 + 25759 * XGtCa**2 + 41249 * XGtMn**2
        - 11344 * XGtFe * XGtMg + 31326.5 * XGtFe * XGtCa + 45841 * XGtFe * XGtMn
        + 132228 * XGtMg * XGtCa + 82498 * XGtMg * XGtMn + 88610.5 * XGtCa * XGtMn
        - 23244 * XGtFe**2 * XGtMg + 2608 * XGtFe**2 * XGtCa - 3234 * XGtFe**2 * XGtMn
        + 11344 * XGtMg**2 * XGtFe - 132228 * XGtMg**2 * XGtCa - 82498 * XGtMg**2 * XGtMn
        - 39864 * XGtCa**2 * XGtFe - 51518 * XGtCa**2 * XGtMg - 2850 * XGtCa**2 * XGtMn
        - 3234 * XGtMn**2 * XGtFe - 82498 * XGtMn**2 * XGtMg - 2850 * XGtMn**2 * XGtCa
        - 62653 * XGtFe * XGtMg * XGtCa - 91682 * XGtFe * XGtMg * XGtMn
        - 24712 * XGtFe * XGtCa * XGtMn - 177221 * XGtMg * XGtCa * XGtMn
    )

    P = np.asarray(P_kbar, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_term = np.log((XGtFe / (XGtMg + eps)) / (XOlFe / (XOlMg + eps)))
    ln_term = np.nan_to_num(ln_term, nan=0.0, posinf=0.0, neginf=0.0)

    num = (
        58136.3
        + (P * 1000.0) * (0.406 - 2 * AH + 2 * AK)
        - 2 * AI
        + 2 * AL
        - 23653.8 * (XOlFe - XOlMg)
    )
    den = 38.743 + 2 * AG - 2 * AJ + 6 * R * ln_term
    # Extra offset that distinguishes WZ07corr from the ordinary WZ07.
    # offset_c may be a scalar or an array (T-sensitive calibration).
    if offset_c is None:
        offset_c = PARAMS.fe3.wz07_offset_c
    T_C = num / np.clip(den, 1e-8, None) - 273.0 - np.asarray(offset_c, dtype=float)
    return T_C


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def iterate_garnet_fe3(
    df,
    max_steps=None,
    dXFe2=None,
    tol_C=None,
    XFe2_start=None,
    verbose=True,
    params=None,
):
    """
    Iterate garnet XFe²⁺ downward by `dXFe2` until
    (T_TA98 − T_WZ07corr) ≤ `tol_C`.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain olivine_Mg, olivine_Fe, garnet_Mg, garnet_Fe,
        garnet_Ca, garnet_Mn, T_TA98, P_CA91.
    max_steps : int
        Safety limit on the number of 0.001 decrements.
    dXFe2 : float
        Step size subtracted from XFe²⁺ each iteration (default 0.001).
    tol_C : float
        Temperature agreement tolerance in °C (default 25).
    XFe2_start : float
        Initial XFe²⁺ (default 0.99 → XFe³⁺ = 0.01).
    verbose : bool
        Print progress summary.

    Returns
    -------
    pandas.DataFrame
        Copy of df with added columns:
          XFe2_gt_iter, XFe3_gt_iter, T_WZ07corr,
          n_fe3_steps, Fe2_gt_iter, Fe3_gt_iter
    """
    fp = params.fe3 if params is not None else PARAMS.fe3
    if max_steps is None:
        max_steps = fp.max_steps
    if dXFe2 is None:
        dXFe2 = fp.dxfe3
    if tol_C is None:
        tol_C = fp.dt_match_c
    if XFe2_start is None:
        XFe2_start = 1.0 - fp.xfe3_start

    result = df.copy()
    n = len(result)

    # Output columns
    result["XFe2_gt_iter"] = np.nan
    result["XFe3_gt_iter"] = np.nan
    result["T_WZ07corr"] = np.nan
    result["n_fe3_steps"] = 0
    result["Fe2_gt_iter"] = np.nan
    result["Fe3_gt_iter"] = np.nan

    required = [
        "olivine_Mg", "olivine_Fe",
        "garnet_Mg", "garnet_Fe", "garnet_Ca", "garnet_Mn",
        "T_TA98", "P_CA91",
    ]
    missing = [c for c in required if c not in result.columns]
    if missing:
        if verbose:
            print(f"\n[fe3_correction] Missing columns {missing} – skipping iteration.")
        return result

    # Candidates: samples with positive Fe in both phases and finite T/P
    Mg_ol = result["olivine_Mg"].to_numpy(dtype=float)
    Fe_ol = result["olivine_Fe"].to_numpy(dtype=float)
    Mg_gt = result["garnet_Mg"].to_numpy(dtype=float)
    Fe_gt = result["garnet_Fe"].to_numpy(dtype=float)
    Ca_gt = result["garnet_Ca"].to_numpy(dtype=float)
    Mn_gt = result["garnet_Mn"].to_numpy(dtype=float)
    T_act = result["T_TA98"].to_numpy(dtype=float)
    P_kbar = result["P_CA91"].to_numpy(dtype=float)

    valid = (
        np.isfinite(Mg_ol) & np.isfinite(Fe_ol) & (Fe_ol > 0)
        & np.isfinite(Mg_gt) & np.isfinite(Fe_gt) & (Fe_gt > 0)
        & np.isfinite(Ca_gt) & np.isfinite(Mn_gt)
        & np.isfinite(T_act) & np.isfinite(P_kbar)
    )
    idx = np.where(valid)[0]
    n_valid = len(idx)

    if n_valid == 0:
        if verbose:
            print("\n[fe3_correction] No samples with complete mineral data for Fe³⁺ iteration.")
        return result

    if verbose:
        print(f"\n=== Garnet Fe³⁺ iteration (Comparative Geothermometry) ===")
        print(f"  Candidates with complete data: {n_valid}")
        print(f"  Start XFe²⁺ = {XFe2_start}, step = −{dXFe2}, tolerance = {tol_C} °C")

    # Work arrays for the valid subset
    XFe2 = np.full(n_valid, XFe2_start, dtype=float)
    steps = np.zeros(n_valid, dtype=int)
    T_wz = np.full(n_valid, np.nan)
    still_active = np.ones(n_valid, dtype=bool)

    Mg_ol_v = Mg_ol[idx]
    Fe_ol_v = Fe_ol[idx]          # assumed all Fe²⁺
    Mg_gt_v = Mg_gt[idx]
    Fe_gt_v = Fe_gt[idx]
    Ca_gt_v = Ca_gt[idx]
    Mn_gt_v = Mn_gt[idx]
    T_act_v = T_act[idx]
    P_v = P_kbar[idx]

    for it in range(max_steps):
        if not still_active.any():
            break

        Fe2_gt = Fe_gt_v * XFe2
        fp_off = params.fe3 if params is not None else PARAMS.fe3
        offset_v = fp_off.wz07_offset_c + fp_off.wz07_offset_slope * T_act_v
        T_wz = _wz07corr_fe2only(
            Mg_ol_v, Fe_ol_v, Mg_gt_v, Fe2_gt, Ca_gt_v, Mn_gt_v, P_v,
            offset_c=offset_v,
        )
        dT = T_act_v - T_wz
        needs_more = still_active & (dT > tol_C) & (XFe2 > dXFe2)
        if not needs_more.any():
            still_active[:] = False
            break

        XFe2[needs_more] -= dXFe2
        steps[needs_more] += 1
        # Samples that no longer need adjustment become inactive
        still_active = needs_more

    # Final evaluation at the converged XFe2
    Fe2_gt = Fe_gt_v * XFe2
    fp_off = params.fe3 if params is not None else PARAMS.fe3
    offset_v = fp_off.wz07_offset_c + fp_off.wz07_offset_slope * T_act_v
    T_wz = _wz07corr_fe2only(
        Mg_ol_v, Fe_ol_v, Mg_gt_v, Fe2_gt, Ca_gt_v, Mn_gt_v, P_v,
        offset_c=offset_v,
    )
    XFe3 = 1.0 - XFe2

    # Write results back
    result.loc[result.index[idx], "XFe2_gt_iter"] = XFe2
    result.loc[result.index[idx], "XFe3_gt_iter"] = XFe3
    result.loc[result.index[idx], "T_WZ07corr"] = T_wz
    result.loc[result.index[idx], "n_fe3_steps"] = steps
    result.loc[result.index[idx], "Fe2_gt_iter"] = Fe2_gt
    result.loc[result.index[idx], "Fe3_gt_iter"] = Fe_gt_v * XFe3

    n_converged = int(np.sum(np.isfinite(XFe3) & ((T_act_v - T_wz) <= tol_C)))
    if verbose:
        print(f"  Converged within {tol_C} °C: {n_converged} of {n_valid}")
        print(f"  Mean steps taken: {steps.mean():.1f}  (max {steps.max()})")

    return result


def apply_empirical_fe3_correction(df, verbose=True, params=None):
    """
    Apply the temperature-dependent empirical correction (manuscript Eq. 5).

        if T_TA98 < t_emp_min_c:
            XFe3_corr = XFe3_iterated
        else:
            XFe3_corr = XFe3_iterated + emp_slope * T + emp_intercept

    emp_slope and emp_intercept are the least-squares line of
    (measured − estimated) on T for the Fig. 5b regression samples
    (Slave, error groups removed, T ≥ t_emp_min_c) at ΔT = emp_fit_dt_c;
    see calibration.fit_empirical_correction.  No level shift is applied.

    Then convert to cation units using total garnet Fe.

    Adds columns:
      XFe3_gt_corr, Fe3_gt_corr, Fe2_gt_corr
    """
    fp = params.fe3 if params is not None else PARAMS.fe3
    if not fp.apply_empirical_correction:
        result = df.copy()
        if "XFe3_gt_corr" not in result.columns and "XFe3_gt_iter" in result.columns:
            result["XFe3_gt_corr"] = result["XFe3_gt_iter"]
        if verbose:
            print("[fe3_correction] Empirical correction disabled.")
        return result
    result = df.copy()

    if "XFe3_gt_iter" not in result.columns:
        if verbose:
            print("\n[fe3_correction] XFe3_gt_iter missing – empirical correction skipped.")
        return result

    T = result["T_TA98"].to_numpy(dtype=float)
    XFe3_iter = result["XFe3_gt_iter"].to_numpy(dtype=float)
    Fe_gt = result["garnet_Fe"].to_numpy(dtype=float)

    XFe3_corr = np.where(
        T < fp.t_emp_min_c,
        XFe3_iter,
        XFe3_iter + (fp.emp_slope * T + fp.emp_intercept),
    )
    # Guard against unphysical values
    XFe3_corr = np.clip(XFe3_corr, 0.0, 1.0)

    result["XFe3_gt_corr"] = XFe3_corr
    result["Fe3_gt_corr"] = XFe3_corr * Fe_gt
    result["Fe2_gt_corr"] = Fe_gt - result["Fe3_gt_corr"]

    n_corr = int(np.isfinite(XFe3_corr).sum())
    if verbose:
        print(f"\n[fe3_correction] Empirical correction applied to {n_corr} samples.")
        n_above = int(np.sum(np.isfinite(T) & (T >= fp.t_emp_min_c) & np.isfinite(XFe3_iter)))
        print(f"  of which {n_above} received the high-T adjustment (T ≥ {fp.t_emp_min_c:.0f} °C).")

    return result
