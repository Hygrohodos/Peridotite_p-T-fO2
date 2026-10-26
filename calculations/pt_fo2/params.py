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

    # Temperature agreement tolerance (°C): stop when
    #   TA98 − WZ07corr  ≤  dt_match_c
    # Chosen in Fig. 6 (calibration.dt_sweep / select_dt): Eq. 5 held
    # fixed, scored on the Eq. 5 regression samples only; highest ±0.02
    # rate, ties broken by smallest |mean bias|.
    # Re-derive with:  python scripts/fit_calibration.py
    dt_match_c: float = 29.0

    # Calibration offset applied to raw WZ07 (°C):
    #   T_WZ07corr = T_WZ07raw − (wz07_offset_c + wz07_offset_slope * T_TA98)
    # Slope = 0 is a flat offset. Non-zero slope is a T-sensitive
    # experimental calibration (TA98−WZ07 = a·T + b  ⇒  offset = −a·T − b).
    wz07_offset_c: float = 118.0
    wz07_offset_slope: float = 0.0

    # Safety cap on iteration steps
    max_steps: int = 200

    # Empirical high-T correction (manuscript Eq. 5), T ≥ t_emp_min_c only:
    #   XFe3_corr = XFe3_iter + emp_slope·T + emp_intercept
    # Fitted exactly as drawn in Fig. 5b (calibration.fit_empirical_correction):
    #   uncorrected estimates at ΔT = emp_fit_dt_c (25 °C); Slave calibration
    #   samples minus the blue (meas − est > emp_blue_threshold) and red
    #   (emp_red_ids) error groups; T_TA98 ≥ t_emp_min_c (the dashed
    #   envelope, n = 21); ordinary least squares of (meas − est) on T,
    #   no level shift.  Re-derive with:  python scripts/fit_calibration.py
    apply_empirical_correction: bool = True
    t_emp_min_c: float = 1130.0
    emp_fit_dt_c: float = 25.0
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
