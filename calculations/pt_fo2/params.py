"""
Editable calculation parameters
===============================

Change the numbers in this file, then re-run:

    python run_pipeline.py
    python scripts/fit_calibration.py
    python scripts/make_figures.py

You do not need to edit the calculation modules themselves.
Values below match the published chapter defaults (Greene, in prep.).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass
class ScreeningParams:
    # Minimum TA98 temperature (°C). Samples below this are excluded
    # because TA98 and WZ07 are not calibrated below ~900 °C.
    t_min_c: float = 900.0

    # |TA98 − NG10| threshold (°C) for opx–cpx equilibrium
    dt_opx_cpx_c: float = 80.0

    # |HA84 − quadratic in TA98| threshold (°C) for pyroxene–garnet equilibrium.
    # Replaces the KR00−HA84 test (2026-10-09). Fit is the median HA84 in
    # 50 °C TA98 bins for samples that passed the opx–cpx screen.
    dt_pyx_gt_c: float = 60.0
    # HA84_pred = ha84_poly_a + ha84_poly_b * TA98 + ha84_poly_c * TA98²
    ha84_poly_a: float = 1203.7
    ha84_poly_b: float = -0.9252
    ha84_poly_c: float = 7.195802e-4

    # |WZ07-predicted from HA84 − WZ07| threshold (°C) for garnet–olivine
    dt_gt_ol_c: float = 50.0

    # Oxide wt% SiO2 above which a phase is considered present
    sio2_present_wt: float = 1.0

    # Analytical quality (Nimis and Grütter 2010, p. 413): opx, cpx and garnet
    # oxide totals within 98.5–101.5 wt%; opx and cpx cation sums ≥ 3.990 apfu
    # on 6 O; garnet Cr2O3 > 1 wt% or Mg/(Mg + Fe) > 0.6.
    oxide_total_min_wt: float = 98.5
    oxide_total_max_wt: float = 101.5
    pyx_cation_sum_min: float = 3.990
    garnet_cr2o3_min_wt: float = 1.0
    garnet_mg_number_min: float = 0.6


@dataclass
class Fe3Params:
    # Initial assumed garnet Fe3+/ΣFe.
    # The iterator decreases XFe2 from (1 − xfe3_start), so 0.01 → XFe2 = 0.99.
    # 0.01 is below the typical mantle range (0.05–0.15) so that highly
    # depleted samples can be captured (a 0.05 start pinned shallow
    # estimated ΔFMQ to the BSE isopleth).  The iteration only increases
    # Fe3+/ΣFe, so estimates of 0.01 below t_emp_min_c are floor-limited.
    xfe3_start: float = 0.01

    # Increment in Fe3+/ΣFe (equivalently decrement in XFe2) per step
    dxfe3: float = 0.001

    # Temperature agreement level (°C): stop when
    #   TA98 − WZ07corr  ≤  dt_match_c
    # v2 method (manuscript s35_lowT, adopted 20 Nov 2026): +10 °C on the Eq. 4
    # scale, i.e. with the level folded into WZ07corr (cr_a below, and the flat
    # cr_d_low_c below cr_low_t_c).  Fig. 5 threshold: 79 % of the 19 screened
    # Slave samples outside the error groups within ±0.02, 100 % within ±0.05.
    # History: 29.0 with the −118 °C offset and Eq. 5 (original); −164.0 with
    # the experimental curve, tuned (rev. v1); −157.6 with the experimental
    # curve, anchored to the Slave-25 mean (rev. v2); −96.0, Slave-30 median
    # (rev. v3; calculations/pt_fo2/level_anchor.py); −102.0 with the level NOT
    # folded into WZ07corr (rev. v7_cr, fe606c1).
    dt_match_c: float = 10.0

    # Eq. 4 (calculations/pt_fo2/cr_correction.py has the full method):
    #   T_WZ07corr = T_WZ07 + C(TA98, Cr2O3),
    #   C = cr_a + cr_b*t + cr_c*(Cr2O3 - cr_ref)*t   for TA98 >= cr_low_t_c,
    #       t = clip(TA98, cr_tmin_c, cr_tmax_c) - cr_pivot_c,  Cr2O3 clipped to
    #       [cr_min_wt, cr_max_wt] (missing Cr2O3 -> cr_ref);
    #   C = cr_d_low_c (flat)                          for TA98 <  cr_low_t_c.
    # cr_d_low_c = None drops the flat branch (the slope is then held at cr_tmin_c).
    # "cr2o3" (default) or "experimental" (v7 Eq. 4: wz07_offset_* only, level in dt_match_c).
    wz07_correction: str = "cr2o3"
    cr_a: float = -92.77342661716227   # WLS intercept CR_A_FIT, folded into Eq. 4 (fe606c1: -102 via dt_match_c)
    cr_b: float = -0.285024042750883
    cr_c: float = 0.06733635988444961
    cr_ref: float = 5.364313828837838
    cr_pivot_c: float = 1207.0242233332644
    cr_min_wt: float = 1.3
    cr_max_wt: float = 9.8
    cr_tmin_c: float = 1124.0
    cr_tmax_c: float = 1394.0
    cr_d_low_c: float | None = -113.24224691422943   # mean TA98 - WZ07(measured) of the 8 low-T natural samples
    cr_low_t_c: float = 1124.0

    # Calibration offset applied to raw WZ07 (°C):
    #   T_WZ07corr = T_WZ07raw − (wz07_offset_c + wz07_offset_slope * T_TA98)
    # Slope = 0 is a flat offset. Non-zero slope is a T-sensitive
    # experimental calibration (TA98−WZ07 = a·T + b  ⇒  offset = −a·T − b).
    # Revision v3: no offset and no slope (the level is carried by dt_match_c).
    # Since v7_cr the Eq. 4 correction (cr_*) is added on top of this offset.
    # The weighted 13-experiment fit (scripts/weighted_thermometer_regression.py;
    # rev. v1/v2: wz07_offset_c −734.5397, slope 0.554746, held flat outside
    # 1162.03–1619.80 °C) is shown in the paper (Fig. 5) but is not applied,
    # because the natural samples do not reproduce its level or slope (Fig. 4).
    wz07_offset_c: float = 0.0   # was 118.0 (original), −734.5397 (rev. v1/v2)
    wz07_offset_slope: float = 0.0   # was 0.554746 (rev. v1/v2)
    # Optional TA98 range (°C) over which the T-sensitive offset is
    # evaluated; outside it the offset is held at its end value instead of
    # being extrapolated linearly (None = no limit, the original behaviour).
    # Set to the TA98 range of the calibrating experiments when the offset
    # comes from calculations/pt_fo2/exp_calibration.py.
    wz07_offset_tmin_c: float | None = None   # rev. v1/v2: 1162.03 (irrelevant with slope 0)
    wz07_offset_tmax_c: float | None = None   # rev. v1/v2: 1619.80

    # Safety cap on iteration steps
    max_steps: int = 200

    # Empirical high-T correction (manuscript Eq. 5), T ≥ t_emp_min_c only:
    #   XFe3_corr = XFe3_iter + emp_slope·T + emp_intercept
    # Fitted exactly as drawn in original-submission Fig. 5b (calibration.fit_empirical_correction):
    #   uncorrected estimates at ΔT = emp_fit_dt_c (25 °C); Slave calibration
    #   samples minus the blue (meas − est > emp_blue_threshold) and red
    #   (emp_red_ids) error groups; T_TA98 ≥ t_emp_min_c (the dashed
    #   envelope, n = 21); ordinary least squares of (meas − est) on T,
    #   no level shift.  Re-derive with:  python scripts/fit_calibration.py
    apply_empirical_correction: bool = False   # revision v1: Eq. 5 dropped
    t_emp_min_c: float = 1130.0
    emp_fit_dt_c: float = 10.0   # error groups evaluated at the adopted level (v7_cr: -102.0; v3: -96.0; was 25.0)
    emp_blue_threshold: float = 0.05
    emp_red_ids: tuple = ("14-78", "MX5011", "8-1", "MX104", "MX5000", "40-38")
    emp_offset: float = 0.0   # unused; kept for backward compatibility
    emp_slope: float = 3.8345e-04
    emp_intercept: float = -0.46732


@dataclass
class ExportParams:
    p_decimals: int = 1
    t_decimals: int = 0
    xfe3_decimals: int = 3
    fo2_decimals: int = 2


@dataclass
class CalculationParams:
    screening: ScreeningParams
    fe3: Fe3Params
    export: ExportParams

    def to_dict(self):
        return asdict(self)


# The single object imported by the rest of the code.
# Edit the nested fields above (or replace this instance) to change behaviour.
PARAMS = CalculationParams(
    screening=ScreeningParams(),
    fe3=Fe3Params(),
    export=ExportParams(),
)
