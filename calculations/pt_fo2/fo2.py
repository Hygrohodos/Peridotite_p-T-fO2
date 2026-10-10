"""
fo2.py
------
Oxygen-fugacity calculation translated from the fO2 sheet of
SD1 - P-T-fO2.xlsm (Stagno 2013 formulation).

Uses the corrected garnet Fe³⁺ / Fe²⁺ from Comparative Geothermometry
(columns BX / BY) together with CA91 pressure and TA98 temperature.

Outputs
-------
  logfO2          absolute log₁₀ fO₂
  logfO2_FMQ      FMQ buffer at the same P–T
  delta_logfO2_FMQ  logfO2 − logfO2_FMQ   (blank only if Fe³⁺ ≈ 0
                    or result is unphysical)
"""

import numpy as np
import pandas as pd


def calculate_fo2(df, verbose=True):
    """
    Compute Δlog fO₂ (FMQ) following the exact formulas on the fO2 sheet.

    Required columns
    ----------------
    Fe3_gt_corr, Fe2_gt_corr   (from fe3_correction)
    garnet_Al, garnet_Cr, garnet_Mg, garnet_Ca
    olivine_Fe, olivine_Mg
    opx_Fe, opx_Mg
    T_TA98, P_CA91
    screen_pass (optional; no longer used to blank fO₂)
    """
    result = df.copy()

    result["logfO2"] = np.nan
    result["logfO2_FMQ"] = np.nan
    result["delta_logfO2_FMQ"] = np.nan

    required = [
        "Fe3_gt_corr", "Fe2_gt_corr",
        "garnet_Al", "garnet_Cr", "garnet_Mg", "garnet_Ca",
        "olivine_Fe", "olivine_Mg",
        "opx_Fe", "opx_Mg",
        "T_TA98", "P_CA91",
    ]
    missing = [c for c in required if c not in result.columns]
    if missing:
        if verbose:
            print(f"\n[fo2] Missing columns {missing} – skipping fO₂ calculation.")
        return result

    Fe3 = result["Fe3_gt_corr"].to_numpy(dtype=float)
    Fe2 = result["Fe2_gt_corr"].to_numpy(dtype=float)
    Al  = result["garnet_Al"].to_numpy(dtype=float)
    Cr  = result["garnet_Cr"].to_numpy(dtype=float)
    Mg  = result["garnet_Mg"].to_numpy(dtype=float)
    Ca  = result["garnet_Ca"].to_numpy(dtype=float)
    Fe_ol = result["olivine_Fe"].to_numpy(dtype=float)
    Mg_ol = result["olivine_Mg"].to_numpy(dtype=float)
    Fe_opx = result["opx_Fe"].to_numpy(dtype=float)
    Mg_opx = result["opx_Mg"].to_numpy(dtype=float)
    T_C = result["T_TA98"].to_numpy(dtype=float)
    P   = result["P_CA91"].to_numpy(dtype=float)

    valid = (
        np.isfinite(Fe3) & np.isfinite(Fe2) & (Fe2 > 0)
        & np.isfinite(Al) & np.isfinite(Cr)
        & np.isfinite(Mg) & np.isfinite(Ca)
        & np.isfinite(Fe_ol) & np.isfinite(Mg_ol)
        & np.isfinite(Fe_opx) & np.isfinite(Mg_opx)
        & np.isfinite(T_C) & np.isfinite(P)
        & (T_C > 0)
    )
    if not valid.any():
        if verbose:
            print("\n[fo2] No samples with complete data for fO₂ calculation.")
        return result

    Fe3 = Fe3[valid]
    Fe2 = Fe2[valid]
    Al  = Al[valid]
    Cr  = Cr[valid]
    Mg  = Mg[valid]
    Ca  = Ca[valid]
    Fe_ol = Fe_ol[valid]
    Mg_ol = Mg_ol[valid]
    Fe_opx = Fe_opx[valid]
    Mg_opx = Mg_opx[valid]
    T_C = T_C[valid]
    P   = P[valid]
    T_K = T_C + 273.0
    R = 8.314

    eps = 1e-30

    # 1. Dodecahedral site mole fractions (E, F, G)
    denom_dod = Fe2 + Mg + Ca + eps
    XFe2 = Fe2 / denom_dod
    XMg  = Mg  / denom_dod
    XCa  = Ca  / denom_dod

    # 2. Octahedral site (H, I, J) – Excel uses Fe2 in denom for XCr/XAl
    XFe3 = Fe3 / (Fe3 + Cr + Al + eps)
    XCr  = Cr  / (Cr + Al + Fe2 + eps)
    XAl  = Al  / (Al + Cr + Fe2 + eps)

    # 3. Margules parameters (K–S)
    WCaMg = 21627 - 5.78 * T_K + 0.012 * P * 1000
    WMgCa =  9834 - 5.78 * T_K + 0.058 * P * 1000
    WCaFe =   873 - 1.69 * T_K
    WFeCa =  6773 - 1.69 * T_K + 0.03 * P * 1000
    WMgFe =  2117 + 0.07 * P * 1000
    WFeMg =  695.0
    WAlCr = 600 * 4.184
    WAlFe3 = 885 * 4.184
    WCrFe3 = 303 * 4.184
    sum_W = WCaMg + WMgCa + WCaFe + WFeCa + WMgFe + WFeMg

    # 4. Standard Gibbs terms
    dG17 = 32800 - 70 * P
    dG18 = 32800 - 70 * P
    dG19 = 50000.0
    dG20 = 50000.0

    # 5a. RTlna Ca3Al2Si3O12 (AC)
    Y = -1 * XMg * XFe3 * dG17 - XFe2 * XFe3 * dG18 - XMg * XCr * dG19 - XFe2 * XCr * dG20
    Z = (
        XFe2**2 * (WCaFe + 2 * XCa * (WFeCa - WCaFe))
        + XMg**2 * (WCaMg + 2 * XCa * (WMgCa - WCaMg))
        + 2 * XCa * XFe2 * XMg * (WFeCa + WMgCa)
        - 2 * XFe2 * XMg * (XFe2 * WMgFe + XMg * WFeMg)
        + XFe2 * XMg * (0.5 * sum_W * (1 - 2 * XCa))
    )
    AA = XFe3**2 * WAlFe3 + XCr**2 * WAlCr + XFe3 * XCr * (WAlCr + WAlFe3 - WCrFe3)
    with np.errstate(divide="ignore", invalid="ignore"):
        AB = R * T_K * np.log(np.maximum(XCa**3 * XAl**2, eps))
    AC = Y + 3 * Z + 2 * AA + AB

    # 5b. RTlna Ca3Fe2Si3O12 (AI)
    AE = (
        XMg * (1 - XFe3) * dG17
        + XFe2 * (1 - XFe3) * dG18
        - XMg * XCr * dG19
        - XFe2 * XCr * dG20
    )
    AF = Z
    AG = XAl**2 * WAlFe3 + XCr**2 * WCrFe3 + XCr * XAl * (WAlFe3 + WCrFe3 - WAlCr)
    with np.errstate(divide="ignore", invalid="ignore"):
        AH = R * T_K * np.log(np.maximum(XCa**3 * XFe3**2, eps))
    AI = AE + 3 * AF + 2 * AG + AH

    # 5c. RTlna Mg3Al2Si3O12 (AO)
    AK = (
        (1 - XMg) * XFe3 * dG17
        - XFe2 * XFe3 * dG18
        + (1 - XMg) * XCr * dG19
        - XFe2 * XCr * dG20
    )
    AL = (
        XCa**2 * (WMgCa + 2 * XMg * (WCaMg - WMgCa))
        + XFe2**2 * (WMgFe + 2 * XMg * (WFeMg - WMgFe))
        + 2 * XCa * XFe2 * XMg * (WCaMg + WFeMg)
        - 2 * XCa * XFe2 * (XCa * WFeCa + XFe2 * WCaFe)
        + XCa * XFe2 * (0.5 * sum_W) * (1 - 2 * XMg)
    )
    AM = XFe3**2 * WAlFe3 + XCr**2 * WAlCr + XFe3 * XCr * (WAlCr + WAlFe3 - WCrFe3)
    with np.errstate(divide="ignore", invalid="ignore"):
        AN = R * T_K * np.log(np.maximum(XMg**3 * XAl**2, eps))
    AO = AK + 3 * AL + AN + 2 * AM

    # 6. Olivine & opx
    XFe_ol = Fe_ol / (Fe_ol + Mg_ol + eps)
    RTlng_ol = 3700.0 * (1 - XFe_ol)**2
    with np.errstate(divide="ignore", invalid="ignore"):
        RTlna_Fa = R * T_K * np.log(np.maximum(XFe_ol, eps)) + RTlng_ol

    XFe_opx = Fe_opx / (Fe_opx + Mg_opx + eps)
    with np.errstate(divide="ignore", invalid="ignore"):
        RTlna_Fs = R * T_K * np.log(np.maximum(XFe_opx, eps))
        RTlna_En = R * T_K * np.log(np.maximum(1 - XFe_opx, eps))

    # 7. ΔG0r + RTlnK
    dG0r = (
        (-8.6623 + 21655.0 / T_K + (13.6149 - 0.01833 * T_K) * P / T_K)
        * (-np.log(10) * R * T_K)
    )
    RTlnK = 2 * AI + 2 * AO + 4 * RTlna_Fs - 2 * AC - 8 * RTlna_Fa - 6 * RTlna_En
    BA = dG0r + RTlnK

    # 8. log fO2
    with np.errstate(divide="ignore", invalid="ignore"):
        logfO2 = BA / (R * T_K * np.log(10))

    logfO2_FMQ = (
        82.75
        + 0.00484 * T_K
        - 30681.0 / T_K
        - 24.45 * np.log10(T_K)
        + 94.0 * P / T_K
        - 0.002 * P
    )
    delta = logfO2 - logfO2_FMQ

    idx = np.where(valid)[0]
    result.loc[result.index[idx], "logfO2"] = logfO2
    result.loc[result.index[idx], "logfO2_FMQ"] = logfO2_FMQ
    result.loc[result.index[idx], "delta_logfO2_FMQ"] = delta

    # Blank unphysical results (Fe3 driven to ~0 produces extreme logs)
    bad = (
        result["Fe3_gt_corr"].fillna(0).abs() < 1e-6
    ) | (result["delta_logfO2_FMQ"].abs() > 20)
    result.loc[bad, "delta_logfO2_FMQ"] = np.nan
    result.loc[bad, "logfO2"] = np.nan

    n_ok = int(np.isfinite(result["delta_logfO2_FMQ"]).sum())
    if verbose:
        print(f"\n=== Oxygen fugacity (Stagno 2013) ===")
        print(f"  ΔlogfO₂(FMQ) computed for {n_ok} samples with Fe³⁺ estimates.")
        if n_ok > 0:
            vals = result.loc[result["delta_logfO2_FMQ"].notna(), "delta_logfO2_FMQ"]
            print(f"  Range: {vals.min():.2f} to {vals.max():.2f}")
            print(f"  Mean:  {vals.mean():.2f}")

    return result
