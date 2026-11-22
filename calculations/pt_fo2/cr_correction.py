"""
cr_correction.py
----------------
Calibration correction of TA98 - WZ07 that depends on TA98 temperature and garnet Cr2O3 (Eq. 4; v2 method,
manuscript s35_lowT, adopted 20 Nov 2026). This replaces the v7 Eq. 4 (weighted experimental slope 0.555 about
the Slave pivot, level -103.8 C) and the v7_cr level of -102 C.

The correction
==============
For a sample with two-pyroxene temperature T_TA98 (deg C) and garnet Cr2O3 (wt%), the calibration part of
the thermometer discrepancy is

    T_TA98 >= 1124 C:  C(T, Cr) = a + b*t + c*(Cr2O3 - Cr_ref)*t,      t = min(T_TA98, 1394) - T_pivot
    T_TA98 <  1124 C:  C(T, Cr) = D_low                                (flat, no slope or Cr2O3 term)

    a      = -92.77 +- 2.92 deg C  (CR_A_FIT, the joint WLS intercept at the pivot; the level is folded into Eq. 4)
    b      = -0.2850 +- 0.0472 deg C/deg C  (slope at Cr_ref)
    c      = +0.0673 +- 0.0199 deg C/deg C per wt% Cr2O3   (p = 0.0007)
    Cr_ref = 5.3643 wt%   (mean garnet Cr2O3 of the 37 natural fit samples)
    T_pivot = 1207.02 deg C (median TA98 of the 23 screened Slave calibration samples, as in v7)
    D_low  = -113.24 deg C (CR_D_LOW: mean TA98 - T_WZ07(measured Fe3+) of the 8 screened Slave, Siberia and
             Kaapvaal samples with TA98 < 1124 C, outside the fit exclusions; scripts/fit_cr_correction.py)

i.e.  T_WZ07corr = T_WZ07 + C(T, Cr) = T_WZ07 - 92.77 - 0.285*t + 0.067*(Cr2O3 - 5.36)*t   (TA98 >= 1124 C)
      T_WZ07corr = T_WZ07 - 113.24                                                       (TA98 <  1124 C)

Garnet Fe3+/SumFe is raised from the 0.01 starting value in 0.001 steps until

    T_TA98 - T_WZ07corr <= dT_match = +10 deg C

The pipeline subtracts -C(T, Cr) from WZ07 (fe3_correction.wz07_offset) and stops at dt_match_c = +10
(params.Fe3Params). Low-Cr garnets get a steeper (more negative) temperature slope than high-Cr garnets.
Before 20 Nov 2026 (fe606c1) the level was carried by dt_match_c = -102 C (the Fig. 6 +-0.02 optimum) and
the slope was held at 1124 C below 1124 C instead of the flat D_low.

Where it is applied
===================
Downstream of the screens: pipeline.run_full_pipeline computes cations and P-T, then applies the analytical-quality
screen (Nimis and Grutter 2010), the four-phase requirement, T_TA98 > 900 deg C and the three equilibrium screens
(screening.py; all on uncorrected temperatures, so they do not depend on this correction), and only then iterates
garnet Fe3+/SumFe with this correction (fe3_correction.iterate_garnet_fe3). Every screened sample uses its own
garnet Cr2O3. Pressures and TA98 temperatures are unchanged; only the iterated garnet Fe3+/SumFe and fO2 change.

Validity range and clamping
===========================
The fit is constrained over garnet Cr2O3 1.3-9.8 wt% (8 Matjuschkin runs 1.28-1.44 wt%; natural samples
1.77-9.80 wt%) and T_TA98 1124-1394 deg C (natural samples 1124-1366 deg C; the two 1400 deg C Matjuschkin runs
1394.1-1394.2 deg C). Outside that box the inputs are HELD AT THE LIMIT (never extrapolated):

    t  = min(T_TA98, 1394) - T_pivot   (below 1124 C the flat D_low applies instead)
    Cr = clip(Cr2O3, 1.3, 9.8);  samples without garnet Cr2O3 use Cr_ref (dCr = 0)

and the sample is flagged (flags(); cr_correction_flags.csv written by scripts/paper_figures/make_paper_figures.py,
one row per flagged iterated sample with SampleID, Craton, screen_pass, T_TA98, garnet_Cr2O3 and the boolean
columns no_Cr, Cr_low, Cr_high, T_low, T_high). Of the 424 fully screened xenoliths, 231 have T_TA98 < 1124 deg C and 6 have T_TA98 > 1394 deg C;
2 have Cr2O3 < 1.3 wt% and 2 have Cr2O3 > 9.8 wt%; none lacks garnet Cr2O3 (239 flagged in all).

Derivation (scripts/fit_cr_correction.py reproduces every coefficient from the pipeline table)
==========================================================================================
1. Natural samples: the screened Slave calibration (23) and Siberia + Kaapvaal validation (26) xenoliths
   (after the analytical-quality, four-phase, T > 900 deg C and equilibrium screens) with
   D = T_TA98 - T_WZ07(measured garnet Fe3+), level_anchor.discrepancy_measured.
   Kept: T_TA98 >= 1120 deg C, excluding the four error-group samples (14-107, MX044, MX088, MX104) and the
   Udachnaya outlier Uv-50-04  ->  n = 37 (15 Slave, 13 Kaapvaal, 9 Siberia).
   Below ~1120 deg C the natural discrepancies bend below the trend and are not used.
2. Experiments: the 8 Matjuschkin et al. (2014) flank-method runs with garnet Cr2O3 (Rzehak et al. 2020
   report no garnet Cr2O3), each with its paper sigma_i (scripts/weighted_thermometer_regression.py).
3. Joint weighted least squares (statsmodels WLS, weights 1/sigma^2, fixed scale):
       y = a*nat + e*exp + b*t + c*dCr*t
   natural sigma = 17.43 deg C (residual SD of the natural-only fit); the experiments get their own level e
   (+62 +- 31 deg C), so they constrain only the slopes. The WLS intercept a is -92.77 +- 2.92 deg C (CR_A_FIT);
   b and c are the coefficients above (chi2 = 40.9 for 45 points; experiments chi2 = 5.9 for 8 runs).
   The WLS intercept is the adopted level (folded into Eq. 4); dT_match = +10 C is the Fig. 5 threshold
   (79 % of the 19 screened Slave samples outside the error groups within +-0.02, 100 % within +-0.05).
4. Below 1124 C (the lower end of the fit range) the natural discrepancies bend below the trend, so a flat
   level D_low = mean D of the 8 screened natural samples with TA98 < 1124 C (outside FIT_EXCLUDED) is used.
   The natural-only fit gives the same slopes within error (a = -92.82, b = -0.283, c = 0.0652 +- 0.0205).
   Robustness of c (natural-only): Huber 0.066; leave-one-out 0.057-0.078 (all p <= 0.011); leave-one-craton-out
   0.052 (no Siberia) to 0.101 (no Kaapvaal); bootstrap 95 % interval 0.030-0.106, c > 0 in 99.8 % of 5000 resamples.

Slope check against the experiments
===================================
At the mean garnet Cr2O3 of the 8 Matjuschkin runs (1.331 wt%) the correction predicts a TA98 - WZ07 slope of
b + c*(1.331 - 5.364) = -0.5566 deg C/deg C, against the measured weighted slope of the 13 runs, -0.5547 +- 0.216
(Fig. 5; scripts/weighted_thermometer_regression.py). The experiments enter the fit only through b and c (they
have their own level), so this is a consistency check of the slope, not of the level. The natural high-Cr garnets alone have a much flatter slope
(-0.14 +- 0.06 for Cr2O3 >= 4.9 wt%, n = 23) than the low-Cr garnets (-0.44 +- 0.07, n = 15; difference p = 0.003).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Coefficients. b, c, CR_REF, the WLS intercept CR_A_FIT and the low-T level CR_D_LOW are from
# scripts/fit_cr_correction.py (joint natural n = 37 + 8 Matjuschkin runs; 8 low-T natural samples).
CR_A_FIT = -92.77342661716227       # deg C, joint WLS intercept at the pivot
CR_A_FIT_SE = 2.9164312998263755    # 1 sigma of the WLS intercept
CR_A = CR_A_FIT                     # deg C, adopted level at the pivot (folded into Eq. 4)
CR_A_SE = CR_A_FIT_SE
CR_B = -0.285024042750883           # deg C per deg C at Cr_ref
CR_B_SE = 0.04722822442158119
CR_C = 0.06733635988444961          # deg C per deg C per wt% Cr2O3
CR_C_SE = 0.019923184067913657
CR_REF = 5.364313828837838          # wt% garnet Cr2O3
CR_PIVOT = 1207.0242233332644       # deg C
CR_MIN, CR_MAX = 1.3, 9.8           # wt%, validity range (held at the limit outside)
T_MIN, T_MAX = 1124.0, 1394.0       # deg C, validity range (held at the limit above; flat D_low below)
CR_D_LOW = -113.24224691422943      # deg C, flat C(T, Cr) below T_LOW (mean of the 8 low-T natural samples)
T_LOW = 1124.0                      # deg C, flat level below this TA98
DT_MATCH = 10.0                     # deg C, adopted dT_match on the Eq. 4 scale (level folded in)
NATURAL_T_MIN_FIT = 1120.0          # deg C, natural samples used in the fit
FIT_EXCLUDED = ("14-107", "MX044", "MX088", "MX104", "Uv-50-04")


def _coef(fp=None):
    if fp is None:
        return CR_A, CR_B, CR_C, CR_REF, CR_PIVOT, CR_MIN, CR_MAX, T_MIN, T_MAX
    return (fp.cr_a, fp.cr_b, fp.cr_c, fp.cr_ref, fp.cr_pivot_c,
            fp.cr_min_wt, fp.cr_max_wt, fp.cr_tmin_c, fp.cr_tmax_c)


def _low(fp=None):
    """(D_low, T_low): the flat low-T level and the TA98 below which it applies (D_low None = no flat branch)."""
    if fp is None:
        return CR_D_LOW, T_LOW
    return getattr(fp, "cr_d_low_c", None), getattr(fp, "cr_low_t_c", T_LOW)


def cr_term(T_ta98, cr2o3=None, fp=None):
    """C(T, Cr) - a (deg C): b*t + c*dCr*t with T and Cr held within the validity range for TA98 >= T_low, and
    D_low - a (flat) below T_low.

    ``cr2o3`` may be None or NaN (no garnet Cr2O3): Cr_ref is used (dCr = 0)."""
    a, b, c, ref, piv, clo, chi, tlo, thi = _coef(fp)
    d_low, t_low = _low(fp)
    T = np.asarray(T_ta98, dtype=float)
    cr = np.full(T.shape, np.nan) if cr2o3 is None else np.broadcast_to(np.asarray(cr2o3, dtype=float), T.shape)
    t = np.clip(T, tlo, thi) - piv
    dcr = np.where(np.isfinite(cr), np.clip(cr, clo, chi), ref) - ref
    out = b * t + c * dcr * t
    if d_low is not None:
        out = np.where(T < t_low, d_low - a, out)
    return out if out.ndim else float(out)


def predicted_discrepancy(T_ta98, cr2o3=None, fp=None):
    """C(T, Cr): TA98 - WZ07 expected from calibration alone (deg C); a + b*t + c*dCr*t, or D_low below T_low."""
    return _coef(fp)[0] + cr_term(T_ta98, cr2o3, fp)


def corrected_wz07(T_wz07, T_ta98, cr2o3=None, fp=None):
    """T_WZ07corr = T_WZ07 + C(T_TA98, Cr2O3) (deg C); Fe3+ is raised until T_TA98 - T_WZ07corr <= dT_match."""
    return np.asarray(T_wz07, dtype=float) + predicted_discrepancy(T_ta98, cr2o3, fp)


def flags(df, fp=None):
    """Per-sample range flags for the correction (True = input held at the limit / Cr_ref used)."""
    _, _, _, _, _, clo, chi, tlo, thi = _coef(fp)
    T = pd.to_numeric(df["T_TA98"], errors="coerce")
    cr = pd.to_numeric(df["garnet_Cr2O3"], errors="coerce") if "garnet_Cr2O3" in df else pd.Series(np.nan, index=df.index)
    f = pd.DataFrame({"SampleID": df["SampleID"] if "SampleID" in df else df.index,
                      "Craton": df["Craton"] if "Craton" in df else None,
                      "screen_pass": df["screen_pass"] if "screen_pass" in df else None,
                      "T_TA98": T, "garnet_Cr2O3": cr})
    f["no_Cr"] = cr.isna()
    f["Cr_low"] = cr < clo
    f["Cr_high"] = cr > chi
    f["T_low"] = T < tlo
    f["T_high"] = T > thi
    f["flag"] = f[["no_Cr", "Cr_low", "Cr_high", "T_low", "T_high"]].any(axis=1)
    return f


def write_flags(df, path, fp=None):
    """Write the flagged rows (iterated samples only) to ``path``; return the counts per flag."""
    it = df["XFe3_gt_iter"].notna() if "XFe3_gt_iter" in df else pd.to_numeric(df["T_TA98"], errors="coerce").notna()
    f = flags(df, fp)[it.values]
    f[f.flag].to_csv(path, index=False)
    return {k: int(f[k].sum()) for k in ["no_Cr", "Cr_low", "Cr_high", "T_low", "T_high", "flag"]} | {"n_iterated": int(len(f))}
