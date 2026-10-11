"""
stoichiometry.py
----------------
Garnet Fe3+/ΣFe from EMP oxides, translated from the Stoichiometry
sheet of SD1 - P-T-fO2.xlsm.

Charge-balance garnet Fe3+ after Droop (1987) and Quinn et al. (2016).
The two published procedures are implemented separately and give the
same result (see "Equivalence" below).

1. Droop GTR (1987) Mineral Mag 51:431-435.
   Eq. 3 (p. 432): "F = 2X(1 - T/S)", where F is "the number of Fe3+ ions
   per X oxygens", T "the ideal number of cations per formula unit" and S
   "the observed cation total per X oxygens calculated assuming all iron
   to be Fe2+"; garnet X = 12, T = 8 (Table 2).
   Implementation (p. 433):
   (ii)  "If S > T proceed. If not, ... leave all iron as Fe2+"
   (iii) F from Eq. 3
   (iv)  "Normalise the formula to T cations (i.e. multiply each number
         by T/S)"
   (v)   "Check that ... F is less than the new value of total iron
         available. If not, ... set all iron as Fe3+"
   (vi)  Fe3+ from (iii), Fe2+ "the remainder"
   Fe3+/ΣFe is therefore F / (Fe per 12 O × 8/S).  Droop's Table 3 garnet
   (S 8.268, Fe 0.818 -> Fe3+ 0.779, Fe2+ 0.013) is reproduced in
   tests/test_stoichiometry.py.

2. Quinn RJ, Valley JW, Page FZ, Fournelle JH (2016) Am Mineral
   101:1704-1707, "Method" (p. 1704):
   "(1) Calculate the formula from EPMA data normalized to eight cations.
   (2) Calculate the total charge contribution from all cations assuming
   all Fe is Fe2+. (3) (a) If the total cation charge is greater than 24
   ..., then all Fe is ferrous ... (b) If the cation charge is <24,
   ... Fe3+ = 24 - sum_i(Ci Vi)  (1)  (4) (a) If the amount of calculated
   Fe3+ is greater than SFe, then enter zero for Fe2+ and set Fe3+ to
   equal SFe. (b) ... subtract the calculated Fe3+ from the SFe to
   determine amount of Fe2+."
   Quinn et al. add no empirical or Mossbauer-based correction.

Equivalence.  For an anhydrous, vacancy-free garnet with oxygen as the
only anion, the all-Fe2+ cation charge on the 8-cation basis is
2 x 12 x 8/S, so 24 - sum(Ci Vi) = 24(1 - 8/S) = F, and both divide by
total Fe on the 8-cation basis.  Quinn et al. call their procedure
"slightly different than that of Droop (1987), but ... preferred ... due
to its simplicity": the route differs, the result does not.  XFe3_droop
and XFe3_quinn therefore agree to floating-point rounding and are
reported as a single charge-balance estimate (Table 1).

(Before 2026-10-11 the Droop branch divided F by total Fe per 12 O,
skipping step (iv), so XFe3_droop was low by a factor 8/S.)

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
    Add charge-balance Fe3+/ΣFe columns to a copy of `df`.

    New columns
    -----------
    XFe3_droop, Fe3_droop, Fe2_droop, cation_sum_12O
    XFe3_quinn, Fe3_quinn, Fe2_quinn, valence_8cat

    Fe3/Fe2 are cations per 8-cation (12-oxygen) formula.  XFe3_droop and
    XFe3_quinn are identical to rounding (see module docstring); both are
    kept so the identity can be checked.  Rows without a garnet analysis
    are NaN.
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
    # Droop step (iv): F is Fe3+ in the formula normalised to T = 8 cations,
    # so total Fe must be taken on that basis (Fe_12O × 8/S), not per 12 O.
    with np.errstate(divide="ignore", invalid="ignore"):
        Fe_tot_droop = np.where(cation_sum > 8.0, Fe_tot * 8.0 / cation_sum, Fe_tot)
    # Droop step (v): if F exceeds total Fe, all Fe is Fe3+.
    Fe3_droop = np.minimum(Fe3_droop, Fe_tot_droop)
    Fe2_droop = Fe_tot_droop - Fe3_droop
    with np.errstate(divide="ignore", invalid="ignore"):
        XFe3_droop = np.where(Fe_tot_droop > 0, Fe3_droop / Fe_tot_droop, np.nan)

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
    # Quinn step (4a): if calculated Fe3+ exceeds ΣFe, Fe3+ = ΣFe, Fe2+ = 0.
    Fe3_quinn = np.minimum(Fe3_quinn, Fe_tot_8)
    Fe2_quinn = Fe_tot_8 - Fe3_quinn
    with np.errstate(divide="ignore", invalid="ignore"):
        XFe3_quinn = np.where(Fe_tot_8 > 0, Fe3_quinn / Fe_tot_8, np.nan)

    # Rows without a garnet analysis (no oxides at all) get NaN, not the
    # artefacts of an empty formula (e.g. Fe3_quinn = 24, Fe2_quinn = -24).
    no_gt = ~(oxy_sum > 0)
    for arr in (cation_sum, Fe3_droop, Fe2_droop, XFe3_droop,
                valence, Fe3_quinn, Fe2_quinn, XFe3_quinn):
        arr[no_gt] = np.nan

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
        print(f"[stoichiometry] Charge-balance (Droop 1987 = Quinn et al. 2016) Fe³⁺/ΣFe computed for {n_ok} of {n} samples.")

    return result
