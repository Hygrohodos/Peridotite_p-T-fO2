"""
screening.py
------------
Equilibrium screening filters translated from the Screening sheet of SD1.

Samples must pass successive tests:
  1. Four required phases present (cpx, opx, olivine, garnet)
  2. TA98 temperature > 900 °C
  3. Opx–Cpx equilibrium (|TA98 – NG10| < 80 °C)
  4. Pyroxene–Garnet equilibrium (|HA84 − quadratic in TA98| < 60 °C)
  5. Garnet–Olivine equilibrium (ΔT < 50 °C after linear correction)

Note: measured garnet Fe³⁺ is NOT required for the calculation path.

The final column `screen_pass` is True only for samples that survive all filters.
"""

import pandas as pd
import numpy as np

from .params import PARAMS


def apply_equilibrium_screening(df, verbose=True, params=None):
    """
    Apply the full sequence of equilibrium screens.

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
    # Screen 1: four required phases present
    # Excel checks oxide wt% SiO2 > 1 (CU-*!B), not cation Si.
    # Using cation Si > 1 fails for olivine (4-O basis, Si ≈ 0.99–1.0).
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
    if verbose:
        print(f"Four required phases are present in {n_4ph} of {n_total} samples.")

    # ------------------------------------------------------------------
    # Screen 2: TA98 temperature > 900 °C
    # ------------------------------------------------------------------
    T_threshold = sp.t_min_c
    result["screen_T900"] = result["screen_4phases"] & (result["T_TA98"] > T_threshold)
    n_T900 = int(result["screen_T900"].sum())
    if verbose:
        print(f"Temperature > {T_threshold:.0f} °C in {n_T900} of {n_4ph} samples.")

    # ------------------------------------------------------------------
    # Screen 3: Opx–Cpx equilibrium (|TA98 – NG10| < 80 °C)
    # ------------------------------------------------------------------
    dT_opx_cpx = np.abs(result["T_TA98"] - result["T_NG10"])
    result["dT_opx_cpx"] = dT_opx_cpx
    result["screen_opx_cpx"] = result["screen_T900"] & (dT_opx_cpx < sp.dt_opx_cpx_c)
    n_opx_cpx = int(result["screen_opx_cpx"].sum())
    if verbose:
        print(f"Cpx and opx in equilibrium for {n_opx_cpx} of {n_T900} samples.")

    # Screen 4: Pyroxene–garnet equilibrium
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
    # Screen 5: Garnet–Olivine equilibrium
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
        print(f"\n>>> Fully screened (all equilibrium filters passed): {n_pass} of {n_total} samples.")

    stats = {
        "n_total": n_total,
        "n_4phases": n_4ph,
        "n_T900": n_T900,
        "n_opx_cpx": n_opx_cpx,
        "n_pyx_gt": n_pyx_gt,
        "n_gt_ol": n_gt_ol,
        "n_pass": n_pass,
    }
    return result, stats
