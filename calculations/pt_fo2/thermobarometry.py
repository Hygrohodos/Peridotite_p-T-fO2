"""
thermobarometry.py
------------------
Thermometers and barometers translated from SD1 - P-T-fO2.xlsm.

Contains:
  - Taylor (1998) two-pyroxene thermometer (TA98)
  - Modified CA91 / NG85 Al-in-opx barometer
  - Supporting thermometers used for equilibrium screening
  - Iterative P–T solver that couples TA98 and CA91
"""

import numpy as np
import pandas as pd

# Gas constant used in the original workbook (cal K⁻¹ mol⁻¹)
R_cal = 1.98726


def ta98_temperature(df, P_kbar):
    """
    Taylor (1998) two-pyroxene thermometer.
    Corresponds to the Two-pyx sheet, column AV in the workbook.
    """
    P_kbar = np.asarray(P_kbar, dtype=float)

    Ti_c = df["cpx_Ti"].fillna(0).to_numpy(dtype=float)
    Al_c = df["cpx_Al"].fillna(0).to_numpy(dtype=float)
    Cr_c = df["cpx_Cr"].fillna(0).to_numpy(dtype=float)
    Fe_c = df["cpx_Fe"].fillna(0).to_numpy(dtype=float)
    Ca_c = df["cpx_Ca"].fillna(0).to_numpy(dtype=float)
    Na_c = df["cpx_Na"].fillna(0).to_numpy(dtype=float)

    Ti_o = df["opx_Ti"].fillna(0).to_numpy(dtype=float)
    Al_o = df["opx_Al"].fillna(0).to_numpy(dtype=float)
    Cr_o = df["opx_Cr"].fillna(0).to_numpy(dtype=float)
    Fe_o = df["opx_Fe"].fillna(0).to_numpy(dtype=float)
    Ca_o = df["opx_Ca"].fillna(0).to_numpy(dtype=float)
    Na_o = df["opx_Na"].fillna(0).to_numpy(dtype=float)

    Alvi_c = Al_c / 2 - Cr_c / 2 - Ti_c + Na_c / 2
    Aliv_c = Al_c / 2 + Cr_c / 2 + Ti_c - Na_c / 2
    Alvi_o = Al_o / 2 - Cr_o / 2 - Ti_o + Na_o / 2
    Aliv_o = Al_o / 2 + Cr_o / 2 + Ti_o - Na_o / 2

    Xts = Al_c + Cr_c - Na_c

    aEn_c = (1 - Ca_c - Na_c) * (1 - Alvi_c - Cr_c - Ti_c) * (1 - Aliv_c / 2) ** 2
    aEn_o = (1 - Ca_o - Na_o) * (1 - Alvi_o - Cr_o - Ti_o) * (1 - Aliv_o / 2) ** 2
    aEn_c = np.clip(aEn_c, 1e-16, None)
    aEn_o = np.clip(aEn_o, 1e-16, None)

    lnKd = np.log(aEn_c) - np.log(aEn_o)

    num = 24787.0 + 678.0 * (P_kbar / 10.0)
    den = 15.67 + 14.37 * Ti_c + 3.69 * Fe_c - 3.25 * Xts + lnKd ** 2
    T_C = num / np.clip(den, 1e-8, None) - 273.0
    return T_C


def ca91_sodic(df):
    """
    Carswell (1991, Mineral. Mag. 55, p. 24) switch from the Nickel & Green
    (1985) X(M1,Al) to the modified one: orthopyroxene Na > (Cr + Fe3+ + Ti),
    cations per 6 O, with Fe3+ = 0 (Carswell 1991, p. 25). At Na = Cr + Ti
    the two expressions give the same value.
    """
    Na = df["opx_Na"].fillna(0).to_numpy(dtype=float)
    Cr = df["opx_Cr"].fillna(0).to_numpy(dtype=float)
    Ti = df["opx_Ti"].fillna(0).to_numpy(dtype=float)
    return Na > (Cr + Ti)


def pressure_barometer_label(df):
    """'CA91' where Carswell's sodic criterion holds, otherwise 'NG85'."""
    return np.where(ca91_sodic(df), "CA91", "NG85")


def ca91_pressure(df, T_C):
    """
    Nickel & Green (1985) eq. 5 garnet–opx barometer with the Carswell
    (1991) modification of X(M1,Al) for sodic orthopyroxene ("CA91").

    NG85 applies as written for every sample. The CA91 modification is
    switched on only by Carswell's (1991, p. 24) criterion
    Na > (Cr + Fe3+ + Ti) in orthopyroxene, with Fe3+ = 0 as Carswell
    assumes for all X(M1,Al) calculations (p. 25); see ``ca91_sodic``.

        not sodic: X(M1,Al) = (Al - Cr - 2·Ti + Na)/2   NG85 = Carswell E6
        sodic:     X(M1,Al) = (Al - (Na - Cr))/2        Carswell E9

    Site fractions BB/BC/BD/BE follow NG85 (pp. 160–161) with the NG85
    X(M1,Al) (BA). For the E9 allocation the total Al on M1 (Tschermak +
    NaAlSi2O6) equals BA, so X(M1,Mg+Fe) = BB is the same in both branches.
    X(M1,Fe) = Fe/(Mg+Fe)·X(M1,Mg+Fe), as NG85 defines it (p. 161).
    Temperature in the workbook formulas is T(°C)+273 (not 273.15).
    """
    T_K = np.asarray(T_C, dtype=float) + 273.0

    Ti = df["opx_Ti"].fillna(0).to_numpy(dtype=float)
    Al = df["opx_Al"].fillna(0).to_numpy(dtype=float)
    Cr = df["opx_Cr"].fillna(0).to_numpy(dtype=float)
    Fe = df["opx_Fe"].fillna(0).to_numpy(dtype=float)
    Mn = df["opx_Mn"].fillna(0).to_numpy(dtype=float)
    Mg = df["opx_Mg"].fillna(0).to_numpy(dtype=float)
    Ca = df["opx_Ca"].fillna(0).to_numpy(dtype=float)
    Na = df["opx_Na"].fillna(0).to_numpy(dtype=float)

    Al_gt = df["garnet_Al"].fillna(0).to_numpy(dtype=float)
    Cr_gt = df["garnet_Cr"].fillna(0).to_numpy(dtype=float)
    Fe_gt = df["garnet_Fe"].fillna(0).to_numpy(dtype=float)
    Mn_gt = df["garnet_Mn"].fillna(0).to_numpy(dtype=float)
    Mg_gt = df["garnet_Mg"].fillna(0).to_numpy(dtype=float)
    Ca_gt = df["garnet_Ca"].fillna(0).to_numpy(dtype=float)

    # NG85 X(M1,Al) used for site fractions BB, BD, BE
    BA = (Al - Cr - 2.0 * Ti + Na) / 2.0
    BB = 1.0 - BA - Cr - Ti          # X(M1, Mg+Fe)
    BC = 1.0 - Ca - Na - Mn          # X(M2, Mg+Fe)
    BD = (Mg / (Mg + Fe + 1e-16)) * BB  # X(M1, Mg)
    BE = (Fe / (Mg + Fe + 1e-16)) * BB  # X(M1, Fe), NG85 p. 161 (was BC - BD, an M2 quantity)
    BF = Cr                          # X(M1, Cr)

    # X(M1,Al) in the ln and denominator terms: NG85 (BA) unless the opx is
    # sodic by Carswell's criterion, then Carswell E9.
    CJ = np.where(ca91_sodic(df), (Al - (Na - Cr)) / 2.0, BA)
    CJ = np.clip(CJ, 1e-12, None)

    sum_div = Ca_gt + Mg_gt + Fe_gt + Mn_gt + 1e-16
    BG = Ca_gt / sum_div             # X(Gt, Ca)
    BJ = Fe_gt / sum_div             # X(Gt, Fe)
    sum_AlCr = Al_gt + Cr_gt + 1e-16
    BH = Al_gt / sum_AlCr            # X(Gt, Al)
    BI = Cr_gt / sum_AlCr            # X(Gt, Cr)

    denom = 183.3 + 178.98 * CJ * (1.0 - CJ)

    ln_arg = ((1.0 - BG) ** 3 * BH ** 2) / (BB * BC ** 2 * CJ + 1e-16)
    ln_arg = np.clip(ln_arg, 1e-16, None)

    inside = (
        -R_cal * T_K * np.log(ln_arg)
        - 9000.0 * BG ** 2
        - 3400.0 * (2.0 * BI ** 2 - BD * BF)
        - BG * BI * (90853.0 - 52.1 * T_K)
        - 7590.0 * BJ * BG
        + 5157.0 * BD * BE
        + 6047.0
        - 3.23 * T_K
    )

    P = -1.0 / denom * inside
    return np.clip(P, 0.0, 80.0)


def bk90_temperature(df, P_kbar):
    """Brey & Köhler (1990) Ca-in-opx thermometer."""
    ca = np.clip(df["opx_Ca"].fillna(0).to_numpy(dtype=float), 1e-8, None)
    return (6425.0 + 26.4 * np.asarray(P_kbar)) / (-np.log(ca) + 1.843) - 273.15


def ng10_temperature(T_BK90):
    """Nimis & Grütter (2010) correction of BK90."""
    T = np.asarray(T_BK90, dtype=float)
    return -628.7 + 2.069 * T - 4.53e-4 * T ** 2


def ha84_temperature(df, P_kbar):
    """
    Harley (1984) Opx-Garnet Fe-Mg exchange thermometer.
    Corresponds to the Opx-gt sheet, column AF.
    """
    Fe_opx = df["opx_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_opx = df["opx_Mg"].fillna(0).to_numpy(dtype=float)
    Fe_gt  = df["garnet_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_gt  = df["garnet_Mg"].fillna(0).to_numpy(dtype=float)
    Mn_gt  = df["garnet_Mn"].fillna(0).to_numpy(dtype=float)
    Ca_gt  = df["garnet_Ca"].fillna(0).to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        KD = (Fe_gt / (Mg_gt + 1e-16)) / (Fe_opx / (Mg_opx + 1e-16))
    KD = np.clip(np.nan_to_num(KD, nan=1e-16, posinf=1e16, neginf=1e-16), 1e-16, None)
    Xgrs_Mn = (Ca_gt + Mn_gt) / (Ca_gt + Mg_gt + Fe_gt + Mn_gt + 1e-16)

    T_C = (3740.0 + 1400.0 * Xgrs_Mn + 22.86 * np.asarray(P_kbar)) / (
        R_cal * np.log(KD) + 1.96
    ) - 273.15
    return T_C


def kr00_temperature(df, P_kbar):
    """
    Krogh (2000) Cpx-Garnet Fe-Mg thermometer.
    Corresponds to the Cpx-gt sheet, column AW.
    """
    Fe_g = df["garnet_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_g = df["garnet_Mg"].fillna(0).to_numpy(dtype=float)
    Mn_g = df["garnet_Mn"].fillna(0).to_numpy(dtype=float)
    Ca_g = df["garnet_Ca"].fillna(0).to_numpy(dtype=float)
    Fe_c = df["cpx_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_c = df["cpx_Mg"].fillna(0).to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        KD = (Fe_g / (Mg_g + 1e-16)) / (Fe_c / (Mg_c + 1e-16))
    KD = np.clip(np.nan_to_num(KD, nan=1e-16, posinf=1e16, neginf=1e-16), 1e-16, None)

    sum_g = Mn_g + Mg_g + Ca_g + Fe_g + 1e-16
    Xca = Ca_g / sum_g
    Xmn = Mn_g / sum_g
    Xmg_num = Mg_g / (Mg_g + Fe_g + 1e-16)

    P = np.asarray(P_kbar, dtype=float)
    num = (
        1939.9
        + 3270.0 * Xca
        - 1396.0 * Xca ** 2
        + 3319.0 * Xmn
        - 3535.0 * Xmn ** 2
        + 1105.0 * Xmg_num
        - 3561.0 * Xmg_num ** 2
        + 2324.0 * Xmg_num ** 3
        + 169.4 * (P / 10.0)
    )
    den = np.log(KD) + 1.223
    T_C = num / np.clip(den, 1e-8, None) - 273.15
    return T_C


def na09_temperature(df, P_kbar):
    """
    Nakamura (2009) Cpx-Garnet thermometer, published equation.
    The (254.6 + 8.42 P) term uses the product 2·Xprp·Xalm, not the sum.
    Comparison only: not used by the screen or the Fe3+ iteration.
    """
    Fe_g = df["garnet_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_g = df["garnet_Mg"].fillna(0).to_numpy(dtype=float)
    Mn_g = df["garnet_Mn"].fillna(0).to_numpy(dtype=float)
    Ca_g = df["garnet_Ca"].fillna(0).to_numpy(dtype=float)
    Fe_c = df["cpx_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_c = df["cpx_Mg"].fillna(0).to_numpy(dtype=float)
    Al_c = df["cpx_Al"].fillna(0).to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        KD = (Fe_g / (Mg_g + 1e-16)) / (Fe_c / (Mg_c + 1e-16))
    KD = np.clip(np.nan_to_num(KD, nan=1e-16, posinf=1e16, neginf=1e-16), 1e-16, None)

    sum_g = Fe_g + Mn_g + Mg_g + Ca_g + 1e-16
    Xprp = Mg_g / sum_g
    Xalm = Fe_g / sum_g
    Xsps = Mn_g / sum_g
    Xgrs = Ca_g / sum_g

    sum_c = Mg_c + Fe_c + Al_c + 1e-16
    XMgCpx = Mg_c / sum_c
    XFeCpx = Fe_c / sum_c

    BE = 0.5 * Xgrs * (Xprp - Xalm - Xsps)
    BF = 0.5 * Xgrs * (Xprp - Xalm + Xsps)
    BG = 0.5 * (Xgrs + Xsps) * (Xprp - Xalm)

    P = np.asarray(P_kbar, dtype=float)

    num = (
        2784.0 + 14.52 * P
        + (2601.0 + 1.44 * P) * (2 * Xgrs * Xprp - BE)
        + (1183.0 + 6.98 * P) * (Xgrs ** 2 - BE)
        - 105.0 * (2 * Xgrs * Xalm + BF)
        + (814.6 + 3.61 * P) * (Xgrs ** 2 + BF)
        - (254.6 + 8.42 * P) * (2 * Xprp * Xalm - Xalm ** 2 + BG)
        - 83.6 * (Xprp ** 2 - 2 * Xprp * Xalm + BG)
        + 1388.0 * Xsps
        - 462.0 * (XMgCpx - XFeCpx)
    )
    den = (
        np.log(KD)
        + 1.431
        + 0.695 * (2 * Xgrs * Xprp + Xgrs ** 2 - 2 * BE)
        + 0.203 * (Xgrs ** 2 - 2 * Xgrs * Xalm)
        + 0.922 * Xsps
    )
    T_C = num / np.clip(den, 1e-8, None) - 273.15
    return T_C


def wz07_temperature(df, P_kbar):
    """
    Wu & Zhao (2007) Olivine-Garnet Fe-Mg thermometer.
    Corresponds to the Ol-gt sheet, column AN.
    """
    Fe_ol = df["olivine_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_ol = df["olivine_Mg"].fillna(0).to_numpy(dtype=float)
    Fe_gt = df["garnet_Fe"].fillna(0).to_numpy(dtype=float)
    Mg_gt = df["garnet_Mg"].fillna(0).to_numpy(dtype=float)
    Ca_gt = df["garnet_Ca"].fillna(0).to_numpy(dtype=float)
    Mn_gt = df["garnet_Mn"].fillna(0).to_numpy(dtype=float)

    XOlMg = Mg_ol / (Mg_ol + Fe_ol + 1e-16)
    XOlFe = Fe_ol / (Mg_ol + Fe_ol + 1e-16)
    sum_gt = Mg_gt + Fe_gt + Ca_gt + Mn_gt + 1e-16
    XGtMg = Mg_gt / sum_gt
    XGtCa = Ca_gt / sum_gt
    XGtFe = Fe_gt / sum_gt
    XGtMn = Mn_gt / sum_gt

    R = 8.3143

    AG = (  # Fea
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
    AH = (  # Feb
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
    AI = (  # Fec
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
    AJ = (  # Mga
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
    AK = (  # Mgb
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
    AL = (  # Mgc
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
        ln_term = np.log((XGtFe / (XGtMg + 1e-16)) / (XOlFe / (XOlMg + 1e-16)))
    ln_term = np.nan_to_num(ln_term, nan=0.0, posinf=0.0, neginf=0.0)

    num = (
        58136.3
        + (P * 1000.0) * (0.406 - 2 * AH + 2 * AK)
        - 2 * AI
        + 2 * AL
        - 23653.8 * (XOlFe - XOlMg)
    )
    den = 38.743 + 2 * AG - 2 * AJ + 6 * R * ln_term
    T_C = num / np.clip(den, 1e-8, None) - 273.15
    return T_C


def iterate_PT(df, max_iter=25, tol=0.1, P_start=40.0, T_start=1000.0, verbose=True):
    """
    Iteratively solve for TA98 temperature and CA91 pressure.

    Starts from an initial guess and alternates the two equations until
    both T and P change by less than `tol`.

    Returns
    -------
    T_TA98 : numpy.ndarray
    P_CA91 : numpy.ndarray
    """
    n = len(df)
    P = np.full(n, P_start, dtype=float)
    T = np.full(n, T_start, dtype=float)

    mask = (
        (df["cpx_Si"].fillna(0) > 0.5)
        & (df["opx_Si"].fillna(0) > 0.5)
        & (df["garnet_Si"].fillna(0) > 0.5)
    ).to_numpy()

    for i in range(max_iter):
        T_new = ta98_temperature(df, P)
        P_new = ca91_pressure(df, T_new)

        dT = np.nanmax(np.abs(T_new[mask] - T[mask])) if mask.any() else 0.0
        dP = np.nanmax(np.abs(P_new[mask] - P[mask])) if mask.any() else 0.0

        T = np.where(mask, T_new, np.nan)
        P = np.where(mask, P_new, np.nan)

        if dT < tol and dP < tol:
            if verbose:
                print(f"  P–T iteration converged in {i + 1} iterations "
                      f"(max ΔT={dT:.2f} °C, ΔP={dP:.2f} kbar)")
            break
    else:
        if verbose:
            print(f"  P–T iteration reached max_iter={max_iter}")

    return T, P


def calculate_all_temperatures(df, P_CA91, verbose=True):
    """
    Evaluate all supporting thermometers at the final CA91 pressure
    and attach the results to a copy of the DataFrame.
    """
    result = df.copy()
    result["T_BK90"] = bk90_temperature(result, P_CA91)
    result["T_NG10"] = ng10_temperature(result["T_BK90"])
    result["T_HA84"] = ha84_temperature(result, P_CA91)
    result["T_KR00"] = kr00_temperature(result, P_CA91)
    result["T_NA09"] = na09_temperature(result, P_CA91)
    result["T_WZ07"] = wz07_temperature(result, P_CA91)
    return result
