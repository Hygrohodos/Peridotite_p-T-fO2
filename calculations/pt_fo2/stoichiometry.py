"""
stoichiometry.py
----------------
Garnet Fe3+/ΣFe from EMP oxides, translated from the Stoichiometry
sheet of SD1 - P-T-fO2.xlsm.

Two estimators are implemented:

1. Droop (1987)
   Recast the analysis on a 12-oxygen basis.  If the cation sum S > 8,
   Fe3+ cations = 2 × 12 × (1 − 8/S).

2. Quinn et al. (2016)
   Normalise cations to 8, compute the cation valence sum, and assign
   any deficit relative to 24 (the ideal anion charge of 12 oxygens)
   to Fe3+.

Oxide order and molecular weights follow the workbook
(CU-Garnet / Stoichiometry sheet), including the slightly non-standard
Na2O oxygen mass used there (15.99 rather than 15.999).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# Molecular weights exactly as Q3:Z3 on the Stoichiometry sheet
MW = {
    "SiO2": 28.085 + 15.999 * 2,          # 60.083
    "TiO2": 47.867 + 15.999 * 2,          # 79.865
    "Al2O3": 26.982 * 2 + 15.999 * 3,     # 101.943
    "Cr2O3": 51.996 * 2 + 15.999 * 3,     # 151.983
    "FeO": 55.845 + 15.999,               # 71.844
    "MnO": 54.938 + 15.999,               # 70.937
    "MgO": 24.305 + 15.999,               # 40.304
    "CaO": 40.078 + 15.999,               # 56.077
    "Na2O": 22.99 * 2 + 15.99,            # 61.97   (workbook uses 15.99)
    "NiO": 58.693 + 15.999,               # 74.692
}

# cations per oxide formula unit (row 4) and oxygens per oxide (row 5)
N_CATIONS = {
    "SiO2": 1, "TiO2": 1, "Al2O3": 2, "Cr2O3": 2, "FeO": 1,
    "MnO": 1, "MgO": 1, "CaO": 1, "Na2O": 2, "NiO": 1,
}
N_OXYGENS = {
    "SiO2": 2, "TiO2": 2, "Al2O3": 3, "Cr2O3": 3, "FeO": 1,
    "MnO": 1, "MgO": 1, "CaO": 1, "Na2O": 1, "NiO": 1,
}

OXIDE_COLS = {
    "SiO2": "garnet_SiO2",
    "TiO2": "garnet_TiO2",
    "Al2O3": "garnet_Al2O3",
    "Cr2O3": "garnet_Cr2O3",
    "FeO": "garnet_FeO",
    "MnO": "garnet_MnO",
    "MgO": "garnet_MgO",
    "CaO": "garnet_CaO",
    "Na2O": "garnet_Na2O",
    "NiO": "garnet_NiO",
}


def _wt(df, oxide):
    col = OXIDE_COLS[oxide]
    if col not in df.columns:
        return np.zeros(len(df), dtype=float)
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def calculate_stoichiometric_fe3(df, verbose=True):
    """
    Add Droop and Quinn Fe3+/ΣFe columns to a copy of `df`.

    New columns
    -----------
    XFe3_droop, Fe3_droop, Fe2_droop, cation_sum_12O
    XFe3_quinn, Fe3_quinn, Fe2_quinn, valence_8cat
    """
    result = df.copy()
    n = len(result)

    # Molar proportions  wt% / MW
    mol = {ox: _wt(result, ox) / MW[ox] for ox in OXIDE_COLS}

    # Atomic proportions  mol * cations-per-oxide
    atom = {ox: mol[ox] * N_CATIONS[ox] for ox in OXIDE_COLS}

    # Oxygen units  mol * oxygens-per-oxide
    oxy = {ox: mol[ox] * N_OXYGENS[ox] for ox in OXIDE_COLS}
    oxy_sum = np.nansum(np.vstack(list(oxy.values())), axis=0)
    scale_12 = np.divide(12.0, oxy_sum, out=np.full(n, np.nan), where=oxy_sum > 0)

    # Cation units on a 12-oxygen basis
    cat12 = {ox: atom[ox] * scale_12 for ox in OXIDE_COLS}
    cation_sum = np.nansum(np.vstack(list(cat12.values())), axis=0)
    Fe_tot = cat12["FeO"]

    # Droop 1987: Fe3+ = 2*N*(1 - T/S) with N=12, T=8
    with np.errstate(divide="ignore", invalid="ignore"):
        Fe3_droop = np.where(
            (cation_sum > 8.0) & np.isfinite(cation_sum),
            2.0 * 12.0 * (1.0 - 8.0 / cation_sum),
            0.0,
        )
    Fe3_droop = np.clip(Fe3_droop, 0.0, None)
    Fe2_droop = Fe_tot - Fe3_droop
    with np.errstate(divide="ignore", invalid="ignore"):
        XFe3_droop = np.where((Fe3_droop + Fe2_droop) > 0, Fe3_droop / (Fe3_droop + Fe2_droop), np.nan)

    # 8-cation normalisation (Quinn sheet columns CC:CL).
    # Blank oxide cells must be 0, as in Excel; otherwise a missing Na2O
    # or NiO turns the whole valence sum into NaN and Fe3+ is forced to 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        cat8 = {
            ox: np.nan_to_num(cat12[ox] * 8.0 / cation_sum, nan=0.0)
            for ox in OXIDE_COLS
        }

    valence = (
        4.0 * (cat8["SiO2"] + cat8["TiO2"])
        + 3.0 * (cat8["Al2O3"] + cat8["Cr2O3"])
        + 2.0 * (cat8["FeO"] + cat8["MnO"] + cat8["MgO"] + cat8["CaO"] + cat8["NiO"])
        + 1.0 * cat8["Na2O"]
    )
    Fe3_quinn = np.where(np.isfinite(valence) & (24.0 - valence > 0), 24.0 - valence, 0.0)
    Fe_tot_8 = cat8["FeO"]
    Fe2_quinn = Fe_tot_8 - Fe3_quinn
    with np.errstate(divide="ignore", invalid="ignore"):
        XFe3_quinn = np.where((Fe3_quinn + Fe2_quinn) > 0, Fe3_quinn / (Fe3_quinn + Fe2_quinn), np.nan)

    result["cation_sum_12O"] = cation_sum
    result["Fe3_droop"] = Fe3_droop
    result["Fe2_droop"] = Fe2_droop
    result["XFe3_droop"] = XFe3_droop
    result["valence_8cat"] = valence
    result["Fe3_quinn"] = Fe3_quinn
    result["Fe2_quinn"] = Fe2_quinn
    result["XFe3_quinn"] = XFe3_quinn

    if verbose:
        n_ok = int(np.isfinite(XFe3_droop).sum())
        print(f"[stoichiometry] Droop / Quinn Fe³⁺/ΣFe computed for {n_ok} of {n} samples.")

    return result
