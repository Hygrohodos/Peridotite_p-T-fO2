"""
screening.py
------------
Sample screening: analytical quality, then the equilibrium filters translated
from the Screening sheet of SD1.

Samples must pass successive tests:
  1. Analytical quality (Nimis and Grütter 2010, p. 413; see analytical_quality)
  2. Four required phases present (cpx, opx, olivine, garnet)
  3. TA98 temperature > 900 °C
  4. Opx–Cpx equilibrium (|TA98 – NG10| < 80 °C)
  5. Pyroxene–Garnet equilibrium (|HA84 − quadratic in TA98| < 60 °C)
  6. Garnet–Olivine equilibrium (ΔT < 50 °C after linear correction)

  screen_quality and screen_4phases are independent flags. Later screens are
  cumulative and start from the intersection of those two.

Note: measured garnet Fe³⁺ is NOT required for the calculation path.

The final column `screen_pass` is True only for samples that survive all filters.
"""

import pandas as pd
import numpy as np

from .params import PARAMS
from .cations import MINERAL_SETTINGS, oxides_to_apfu


def oxide_total(df, phase):
    """Oxide total (wt%) of one phase: the reported total, or the sum of the oxides (all Fe as FeO) where
    no total is reported. Reported values outside 50–150 wt% are not oxide totals (the 0 and ~500 wt%
    entries in the compilation) and are replaced by the oxide sum."""
    cols = [c for c in MINERAL_SETTINGS[phase]["oxide_cols"] if c in df.columns]
    oxide_sum = df[cols].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=1)
    if f"{phase}_Total" not in df.columns:
        return oxide_sum
    reported = pd.to_numeric(df[f"{phase}_Total"], errors="coerce")
    return reported.where(reported.between(50.0, 150.0), oxide_sum)


def analytical_quality(df, params=None):
    """True where the opx, cpx and garnet analyses meet the quality criteria of Nimis and Grütter (2010, p. 413):
    oxide totals within 98.5–101.5 wt% (oxide_total); opx and cpx cation sums ≥ 3.990 apfu on 6 O (unrounded,
    all Fe as Fe2+); garnet Cr2O3 > 1 wt% or Mg/(Mg + Fe) > 0.6. Thresholds are in params.screening."""
    sp = params.screening if params is not None else PARAMS.screening
    ok = pd.Series(True, index=df.index)
    for phase in ("opx", "cpx", "garnet"):
        ok &= oxide_total(df, phase).between(sp.oxide_total_min_wt, sp.oxide_total_max_wt)
    for phase in ("opx", "cpx"):
        cfg = MINERAL_SETTINGS[phase]
        cations = oxides_to_apfu(df, cfg["oxide_cols"], ideal_oxygens=cfg["ideal_oxygens"], round_digits=None)
        ok &= cations["Sum_Cations"] >= sp.pyx_cation_sum_min
    cr2o3 = pd.to_numeric(df["garnet_Cr2O3"], errors="coerce")
    mg, fe = (pd.to_numeric(df[f"garnet_{c}"], errors="coerce") for c in ("Mg", "Fe"))
    ok &= (cr2o3 > sp.garnet_cr2o3_min_wt) | (mg / (mg + fe) > sp.garnet_mg_number_min)
    return ok


def apply_equilibrium_screening(df, verbose=True, params=None):
    """
    Apply the full sequence of screens (analytical quality first, then four-phase, then equilibrium).

    Parameters
    ----------
    df : pandas.DataFrame
        Must already contain cation columns and the thermometer results
        (T_TA98, T_NG10, T_HA84, T_KR00, T_WZ07, …).
    verbose : bool
        Print how many samples pass each successive filter.

    Returns
    -------
    pandas.DataFrame
        Input DataFrame with added screen_* columns and diagnostic ΔT columns.
    dict
        Summary statistics useful for the website upload page.
    """
    result = df.copy()
    n_total = len(result)
    sp = params.screening if params is not None else PARAMS.screening

    if verbose:
        print("\n=== Equilibrium screening ===")

    # ------------------------------------------------------------------
    # Screen 1: analytical quality of the opx, cpx and garnet analyses
    # Applied first, ahead of the four-phase screen. The flag is not
    # cumulative: a sample can pass quality and still lack a phase.
    # ------------------------------------------------------------------
    result["screen_quality"] = analytical_quality(result, params)
    n_quality = int(result["screen_quality"].sum())
    if verbose:
        print(f"Analyses meet the quality criteria in {n_quality} of {n_total} samples.")

    # ------------------------------------------------------------------
    # Screen 2: four required phases present
    # Excel checks oxide wt% SiO2 > 1 (CU-*!B), not cation Si.
    # Using cation Si > 1 fails for olivine (4-O basis, Si ≈ 0.99–1.0).
    # The flag is the phase-presence test only (not cumulative with quality).
    # ------------------------------------------------------------------
    def _has_phase(prefix):
        # Prefer oxide column if present; else fall back to cation Si > 0.5
        ox = f"{prefix}_SiO2"
        cat = f"{prefix}_Si"
        if ox in result.columns:
            return result[ox].fillna(0) > sp.sio2_present_wt
        return result[cat].fillna(0) > 0.5

    has_opx = _has_phase("opx")
    has_cpx = _has_phase("cpx")
    has_ol  = _has_phase("olivine")
    has_gt  = _has_phase("garnet")

    result["screen_4phases"] = (has_opx & has_cpx & has_ol & has_gt)
    n_4ph = int(result["screen_4phases"].sum())
    n_quality_4ph = int((result["screen_quality"] & result["screen_4phases"]).sum())
    if verbose:
        print(f"Four required phases are present in {n_4ph} of {n_total} samples "
              f"({n_quality_4ph} after the quality screen).")

    # ------------------------------------------------------------------
    # Screen 3: TA98 temperature > 900 °C
    # ------------------------------------------------------------------
    T_threshold = sp.t_min_c
    result["screen_T900"] = (
        result["screen_quality"] & result["screen_4phases"] & (result["T_TA98"] > T_threshold)
    )
    n_T900 = int(result["screen_T900"].sum())
    if verbose:
        print(f"Temperature > {T_threshold:.0f} °C in {n_T900} of {n_quality_4ph} samples.")

    # ------------------------------------------------------------------
    # Screen 4: Opx–Cpx equilibrium (|TA98 – NG10| < 80 °C)
    # ------------------------------------------------------------------
    dT_opx_cpx = np.abs(result["T_TA98"] - result["T_NG10"])
    result["dT_opx_cpx"] = dT_opx_cpx
    result["screen_opx_cpx"] = result["screen_T900"] & (dT_opx_cpx < sp.dt_opx_cpx_c)
    n_opx_cpx = int(result["screen_opx_cpx"].sum())
    if verbose:
        print(f"Cpx and opx in equilibrium for {n_opx_cpx} of {n_T900} samples.")

    # Screen 5: Pyroxene–garnet equilibrium
    # |HA84 − quadratic in TA98| < 60 °C. The curve is the median HA84
    # in 50 °C TA98 bins for samples that passed the opx–cpx screen.
    TA = result["T_TA98"]
    HA_pred = sp.ha84_poly_a + sp.ha84_poly_b * TA + sp.ha84_poly_c * TA ** 2
    dT_pyx_gt = np.abs(result["T_HA84"] - HA_pred)
    result["dT_pyx_gt"] = dT_pyx_gt
    result["screen_pyx_gt"] = result["screen_opx_cpx"] & (dT_pyx_gt < sp.dt_pyx_gt_c)
    n_pyx_gt = int(result["screen_pyx_gt"].sum())
    if verbose:
        print(f"Pyroxene and garnet in equilibrium for {n_pyx_gt} of {n_opx_cpx} samples.")

    # ------------------------------------------------------------------
    # Screen 6: Garnet–Olivine equilibrium
    # ------------------------------------------------------------------
    WZ_pred = 88.99741613 + 0.9864930064 * result["T_HA84"]
    dT_gt_ol = np.abs(WZ_pred - result["T_WZ07"])
    result["dT_gt_ol"] = dT_gt_ol
    result["screen_gt_ol"] = result["screen_pyx_gt"] & (dT_gt_ol < sp.dt_gt_ol_c)
    n_gt_ol = int(result["screen_gt_ol"].sum())
    if verbose:
        print(f"Garnet and olivine in equilibrium for {n_gt_ol} of {n_pyx_gt} samples.")

    # Final flag
    result["screen_pass"] = result["screen_gt_ol"]
    n_pass = int(result["screen_pass"].sum())
    if verbose:
        print(f"\n>>> Fully screened (all filters passed): {n_pass} of {n_total} samples.")

    stats = {
        "n_total": n_total,
        "n_quality": n_quality,
        "n_4phases": n_4ph,
        "n_quality_4phases": n_quality_4ph,
        "n_T900": n_T900,
        "n_opx_cpx": n_opx_cpx,
        "n_pyx_gt": n_pyx_gt,
        "n_gt_ol": n_gt_ol,
        "n_pass": n_pass,
    }
    return result, stats
