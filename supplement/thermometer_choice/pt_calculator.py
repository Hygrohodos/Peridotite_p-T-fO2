#!/usr/bin/env python3
"""
P–T calculator for the Hygrohodos kimberlite thermobarometry suite.

Reads a CSV of mineral chemistry (oxide wt%) and solves one thermometer
with one barometer. The equations are the ones used on the website after
the literature corrections: published Nakamura (2009), not the workbook
v1.1 substitution, and Bussweiler et al. (2017) Al-in-olivine with no Cr#
term. All probe Fe is taken as Fe2+ in the exchange thermometers.

Python 3.9+. Standard library only.

  python pt_calculator.py --list
  python pt_calculator.py chemistry.csv -t TA98 -b CA91 -o pt.csv
  python pt_calculator.py chemistry.csv -t NA09 -b BE15
  python pt_calculator.py chemistry.csv -t BU17 -b MGP
  python pt_calculator.py --self-check
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

T0 = 273.0  # site kelvin offset. BK90 pressure uses 273.15 instead.
R_CAL = 1.98726  # Nickel & Green / Carswell
R_HA = 1.9872  # Harley / Lee–Ganguly
R_J = 8.3143

OXIDES = (
    "SiO2",
    "TiO2",
    "Al2O3",
    "Cr2O3",
    "FeO",
    "MnO",
    "MgO",
    "CaO",
    "Na2O",
    "NiO",
)
OXIDE_INFO = {
    "SiO2": (60.083, 1, 2, "Si"),
    "TiO2": (79.865, 1, 2, "Ti"),
    "Al2O3": (101.961, 2, 3, "Al"),
    "Cr2O3": (151.989, 2, 3, "Cr"),
    "FeO": (71.844, 1, 1, "Fe"),
    "MnO": (70.937, 1, 1, "Mn"),
    "MgO": (40.304, 1, 1, "Mg"),
    "CaO": (56.077, 1, 1, "Ca"),
    "Na2O": (61.97, 2, 1, "Na"),
    "NiO": (74.692, 1, 1, "Ni"),
}
MINERALS = (
    ("cpx", 6),
    ("opx", 6),
    ("garnet", 12),
    ("olivine", 4),
)
PREFIX_ALIAS = {
    "cpx": "cpx",
    "clinopyroxene": "cpx",
    "opx": "opx",
    "orthopyroxene": "opx",
    "garnet": "garnet",
    "gt": "garnet",
    "grt": "garnet",
    "olivine": "olivine",
    "ol": "olivine",
}
META_KEYS = ("SampleID", "Location", "Craton", "Lithology", "Reference")
AW_AL = 26.9815385
MW_AL2O3 = 101.961
OLIVINE_MG_SPLIT = 91.5

Row = Dict[str, object]


# ---------------------------------------------------------------------------
# Numerics
# ---------------------------------------------------------------------------

def clip(x: float, lo: float, hi: float = math.inf) -> float:
    if not math.isfinite(x):
        return x
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


def round_half_even(x: float, digits: int = 4) -> float:
    """Same 4 d.p. cation rounding as the website (banker's round)."""
    if not math.isfinite(x):
        return x
    p = 10 ** digits
    n = x * p
    s = -1.0 if n < 0 else 1.0
    a = abs(n)
    f = math.floor(a + 1e-12)
    frac = a - f
    if frac > 0.5 + 1e-12:
        r = f + 1
    elif frac < 0.5 - 1e-12:
        r = f
    else:
        r = f if f % 2 == 0 else f + 1
    return (s * r) / p


def num0(row: Row, key: str) -> float:
    v = row.get(key, None)
    if isinstance(v, bool) or v is None:
        return 0.0
    if isinstance(v, (int, float)):
        return float(v) if math.isfinite(float(v)) else 0.0
    if isinstance(v, str) and v.strip() != "":
        try:
            n = float(v)
        except ValueError:
            return 0.0
        return n if math.isfinite(n) else 0.0
    return 0.0


def finite(x: float) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


# ---------------------------------------------------------------------------
# Cations
# ---------------------------------------------------------------------------

def oxides_to_cations(row: Row) -> Row:
    """Recast oxide wt% to cations per 6 (pyroxene), 12 (garnet) or 4 (olivine) oxygens."""
    out = dict(row)
    for prefix, ideal_o in MINERALS:
        # Keep a cation block already supplied with no oxides.
        if num0(out, prefix + "_Si") > 0 and num0(out, prefix + "_SiO2") <= 0:
            continue
        moles: Dict[str, float] = {}
        total_o = 0.0
        any_ox = False
        for ox in OXIDES:
            mw, n_cat, n_o, cation = OXIDE_INFO[ox]
            wt = num0(out, prefix + "_" + ox)
            if wt != 0:
                any_ox = True
            moles_ox = wt / mw
            moles[cation] = moles.get(cation, 0.0) + moles_ox * n_cat
            total_o += moles_ox * n_o
        if not any_ox:
            continue
        scale = ideal_o / total_o if total_o > 0 else 1.0
        total_cat = 0.0
        for cation, m in moles.items():
            v = round_half_even(m * scale, 4)
            out[prefix + "_" + cation] = v
            total_cat += v
        out[prefix + "_Sum_Cations"] = round_half_even(total_cat, 4)
        out[prefix + "_Norm_O"] = ideal_o
    return out


def olivine_al_ppm(row: Row) -> float:
    if num0(row, "olivine_Al_ppm") > 0:
        return num0(row, "olivine_Al_ppm")
    if num0(row, "al_ppm") > 0:
        return num0(row, "al_ppm")
    al2o3 = num0(row, "olivine_Al2O3")
    if al2o3 > 0:
        return al2o3 * 10000.0 * (2.0 * AW_AL) / MW_AL2O3
    return 0.0


def olivine_mg_number(row: Row) -> float:
    """100 Mg/(Mg+Fe) from oxide moles, else from cations."""
    mgo = num0(row, "olivine_MgO")
    feo = num0(row, "olivine_FeO")
    if mgo > 0 or feo > 0:
        mg = mgo / 40.304
        fe = feo / 71.844
    else:
        mg = num0(row, "olivine_Mg")
        fe = num0(row, "olivine_Fe")
    den = mg + fe
    if den <= 0:
        return float("nan")
    return 100.0 * mg / den


# ---------------------------------------------------------------------------
# Shared activity pieces
# ---------------------------------------------------------------------------

def px(row: Row, pre: str) -> Dict[str, float]:
    return {
        "ti": num0(row, pre + "_Ti"),
        "al": num0(row, pre + "_Al"),
        "cr": num0(row, pre + "_Cr"),
        "fe": num0(row, pre + "_Fe"),
        "mn": num0(row, pre + "_Mn"),
        "mg": num0(row, pre + "_Mg"),
        "ca": num0(row, pre + "_Ca"),
        "na": num0(row, pre + "_Na"),
    }


def a_en_wb(p: Dict[str, float]) -> float:
    m2 = p["ca"] + p["mg"] / 2 + p["fe"] / 2 + p["mn"] + p["na"]
    m1 = p["fe"] / 2 + p["al"] + p["ti"] + p["cr"] + p["mg"] / 2
    if not (m2 > 0 and m1 > 0 and p["mg"] > 0):
        return float("nan")
    return (p["mg"] / 2 / m2) * (p["mg"] / 2 / m1)


def wb_pair(row: Row) -> Optional[Dict[str, float]]:
    c = px(row, "cpx")
    o = px(row, "opx")
    ac = a_en_wb(c)
    ao = a_en_wb(o)
    if not (ac > 0 and ao > 0 and (o["fe"] + o["mg"]) > 0):
        return None
    return {"ac": ac, "ao": ao, "xFe": o["fe"] / (o["fe"] + o["mg"])}


def gt_cpx_kd(row: Row) -> Optional[Dict[str, float]]:
    fe_g = num0(row, "garnet_Fe")
    mn_g = num0(row, "garnet_Mn")
    mg_g = num0(row, "garnet_Mg")
    ca_g = num0(row, "garnet_Ca")
    fe_c = num0(row, "cpx_Fe")
    mg_c = num0(row, "cpx_Mg")
    if not (mg_g > 0 and mg_c > 0 and fe_g > 0 and fe_c > 0):
        return None
    return {
        "kd": (fe_g / mg_g) / (fe_c / mg_c),
        "feG": fe_g,
        "mnG": mn_g,
        "mgG": mg_g,
        "caG": ca_g,
    }


def opx_gt_kd(row: Row) -> Optional[Dict[str, float]]:
    fe_o = num0(row, "opx_Fe")
    mg_o = num0(row, "opx_Mg")
    fe_g = num0(row, "garnet_Fe")
    mn_g = num0(row, "garnet_Mn")
    mg_g = num0(row, "garnet_Mg")
    ca_g = num0(row, "garnet_Ca")
    total = ca_g + mg_g + fe_g + mn_g
    if not (mg_o > 0 and mg_g > 0 and fe_o > 0 and fe_g > 0 and total > 0):
        return None
    return {
        "kd": (fe_g / mg_g) / (fe_o / mg_o),
        "xCaMn": (ca_g + mn_g) / total,
        "xCa": ca_g / total,
        "xMn": mn_g / total,
    }


def opx_gt_cations(row: Row) -> Optional[Tuple[Dict[str, float], Dict[str, float]]]:
    opx = {
        "si": num0(row, "opx_Si"),
        "ti": num0(row, "opx_Ti"),
        "al": num0(row, "opx_Al"),
        "cr": num0(row, "opx_Cr"),
        "fe": num0(row, "opx_Fe"),
        "mn": num0(row, "opx_Mn"),
        "mg": num0(row, "opx_Mg"),
        "ca": num0(row, "opx_Ca"),
        "na": num0(row, "opx_Na"),
    }
    gt = {
        "al": num0(row, "garnet_Al"),
        "cr": num0(row, "garnet_Cr"),
        "fe": num0(row, "garnet_Fe"),
        "mn": num0(row, "garnet_Mn"),
        "mg": num0(row, "garnet_Mg"),
        "ca": num0(row, "garnet_Ca"),
    }
    if not (
        opx["al"] > 0
        and (opx["fe"] + opx["mg"]) > 0
        and (gt["fe"] + gt["mg"]) > 0
        and (gt["al"] + gt["cr"]) > 0
    ):
        return None
    return opx, gt


def ng_family_pressure(x_al: float, opx: Dict[str, float], gt: Dict[str, float], t_c: float) -> float:
    x_mgfe_m1 = 1 - x_al - opx["cr"] - opx["ti"]
    x_mgfe_m2 = 1 - opx["ca"] - opx["na"] - opx["mn"]
    if not (
        x_al > 0
        and x_al < 1
        and x_mgfe_m1 > 0
        and x_mgfe_m2 > 0
        and (opx["fe"] + opx["mg"]) > 0
    ):
        return float("nan")
    x_mg_m1 = (opx["mg"] / (opx["mg"] + opx["fe"])) * x_mgfe_m1
    x_fe_m1 = x_mgfe_m2 - x_mg_m1
    sum_g = gt["ca"] + gt["mg"] + gt["fe"] + gt["mn"]
    sum_al = gt["al"] + gt["cr"]
    if not (sum_g > 0 and sum_al > 0):
        return float("nan")
    x_ca = gt["ca"] / sum_g
    x_al_g = gt["al"] / sum_al
    x_cr_g = gt["cr"] / sum_al
    x_fe_g = gt["fe"] / sum_g
    if not ((1 - x_ca) > 0 and x_al_g > 0):
        return float("nan")
    tk = t_c + T0
    den = 183.3 + 178.98 * x_al * (1 - x_al)
    arg = ((1 - x_ca) ** 3 * x_al_g ** 2) / (x_mgfe_m1 * x_mgfe_m2 ** 2 * x_al)
    if not (arg > 0 and den != 0):
        return float("nan")
    inside = (
        -R_CAL * tk * math.log(arg)
        - 9e3 * x_ca ** 2
        - 3400 * (2 * x_cr_g ** 2 - x_mg_m1 * opx["cr"])
        - x_ca * x_cr_g * (90853 - 52.1 * tk)
        - 7590 * x_fe_g * x_ca
        + 5157 * x_mg_m1 * x_fe_m1
        + 6047
        - 3.23 * tk
    )
    return (-1 / den) * inside


# ---------------------------------------------------------------------------
# Thermometers. Each returns °C. P is kbar. Fe is total Fe as Fe2+.
# ---------------------------------------------------------------------------

def t_ta98(row: Row, p_kbar: float) -> float:
    """Taylor (1998) two-pyroxene. Preferred peridotite thermometer."""
    ti_c = num0(row, "cpx_Ti")
    al_c = num0(row, "cpx_Al")
    cr_c = num0(row, "cpx_Cr")
    fe_c = num0(row, "cpx_Fe")
    ca_c = num0(row, "cpx_Ca")
    na_c = num0(row, "cpx_Na")
    ti_o = num0(row, "opx_Ti")
    al_o = num0(row, "opx_Al")
    cr_o = num0(row, "opx_Cr")
    ca_o = num0(row, "opx_Ca")
    na_o = num0(row, "opx_Na")
    alvi_c = al_c / 2 - cr_c / 2 - ti_c + na_c / 2
    aliv_c = al_c / 2 + cr_c / 2 + ti_c - na_c / 2
    alvi_o = al_o / 2 - cr_o / 2 - ti_o + na_o / 2
    aliv_o = al_o / 2 + cr_o / 2 + ti_o - na_o / 2
    xts = al_c + cr_c - na_c
    a_en_c = clip((1 - ca_c - na_c) * (1 - alvi_c - cr_c - ti_c) * (1 - aliv_c / 2) ** 2, 1e-16)
    a_en_o = clip((1 - ca_o - na_o) * (1 - alvi_o - cr_o - ti_o) * (1 - aliv_o / 2) ** 2, 1e-16)
    ln_kd = math.log(a_en_c) - math.log(a_en_o)
    num = 24787 + 678 * (p_kbar / 10)
    den = 15.67 + 14.37 * ti_c + 3.69 * fe_c - 3.25 * xts + ln_kd ** 2
    return num / max(den, 1e-8) - T0


def t_bk90_ca(row: Row, p_kbar: float) -> float:
    """Brey & Köhler (1990) Ca-in-orthopyroxene."""
    ca = clip(num0(row, "opx_Ca"), 1e-8)
    return (6425 + 26.4 * p_kbar) / (-math.log(ca) + 1.843) - T0


def t_ng10(row: Row, p_kbar: float) -> float:
    """Nimis & Grütter (2010) correction of Brey & Köhler Ca-in-opx."""
    tb = t_bk90_ca(row, p_kbar)
    if not finite(tb):
        return float("nan")
    return -628.7 + 2.069 * tb - 453e-6 * tb ** 2


def t_ha84(row: Row, p_kbar: float) -> float:
    """Harley (1984) opx–garnet."""
    kd = num0(row, "garnet_Fe") / (num0(row, "garnet_Mg") + 1e-16) / (
        num0(row, "opx_Fe") / (num0(row, "opx_Mg") + 1e-16)
    )
    if not finite(kd):
        kd = 1e-16
    kd = clip(kd, 1e-16)
    fe_g = num0(row, "garnet_Fe")
    mg_g = num0(row, "garnet_Mg")
    mn_g = num0(row, "garnet_Mn")
    ca_g = num0(row, "garnet_Ca")
    x = (ca_g + mn_g) / (ca_g + mg_g + fe_g + mn_g + 1e-16)
    return (3740 + 1400 * x + 22.86 * p_kbar) / (R_HA * math.log(kd) + 1.96) - T0


def t_kr00(row: Row, p_kbar: float) -> float:
    """Krogh Ravna (2000) garnet–cpx. P term is in GPa (kbar/10)."""
    fe_g = num0(row, "garnet_Fe")
    mg_g = num0(row, "garnet_Mg")
    mn_g = num0(row, "garnet_Mn")
    ca_g = num0(row, "garnet_Ca")
    fe_c = num0(row, "cpx_Fe")
    mg_c = num0(row, "cpx_Mg")
    kd = fe_g / (mg_g + 1e-16) / (fe_c / (mg_c + 1e-16))
    if not finite(kd):
        kd = 1e-16
    kd = clip(kd, 1e-16)
    sum_g = mn_g + mg_g + ca_g + fe_g + 1e-16
    xca = ca_g / sum_g
    xmn = mn_g / sum_g
    xmg = mg_g / (mg_g + fe_g + 1e-16)
    num = (
        1939.9
        + 3270 * xca
        - 1396 * xca ** 2
        + 3319 * xmn
        - 3535 * xmn ** 2
        + 1105 * xmg
        - 3561 * xmg ** 2
        + 2324 * xmg ** 3
        + 169.4 * (p_kbar / 10)
    )
    return num / max(math.log(kd) + 1.223, 1e-8) - T0


def t_na09(row: Row, p_kbar: float) -> float:
    """Published Nakamura (2009) garnet–cpx.

    A = 0.5 Xgrs (Xprp − Xalm − Xsps)
    B = 0.5 Xgrs (Xprp − Xalm + Xsps)
    C = 0.5 (Xgrs + Xsps) (Xprp − Xalm)
    The (2 Xprp Xalm) term is a product. −273 sits outside the fraction.
    XMg and XFe in cpx are over (Al + Fe_total + Mg), with no Ca.
    """
    fe_g = num0(row, "garnet_Fe")
    mg_g = num0(row, "garnet_Mg")
    mn_g = num0(row, "garnet_Mn")
    ca_g = num0(row, "garnet_Ca")
    fe_c = num0(row, "cpx_Fe")
    mg_c = num0(row, "cpx_Mg")
    al_c = num0(row, "cpx_Al")
    kd = fe_g / (mg_g + 1e-16) / (fe_c / (mg_c + 1e-16))
    if not finite(kd):
        kd = 1e-16
    kd = clip(kd, 1e-16)
    sum_g = fe_g + mn_g + mg_g + ca_g + 1e-16
    xprp = mg_g / sum_g
    xalm = fe_g / sum_g
    xsps = mn_g / sum_g
    xgrs = ca_g / sum_g
    sum_c = mg_c + fe_c + al_c + 1e-16
    xmg = mg_c / sum_c
    xfe = fe_c / sum_c
    a = 0.5 * xgrs * (xprp - xalm - xsps)
    b = 0.5 * xgrs * (xprp - xalm + xsps)
    c = 0.5 * (xgrs + xsps) * (xprp - xalm)
    num = (
        2784
        + 14.52 * p_kbar
        + (2601 + 1.44 * p_kbar) * (2 * xgrs * xprp - a)
        + (1183 + 6.98 * p_kbar) * (xgrs ** 2 - a)
        - 105 * (2 * xgrs * xalm + b)
        + (814.6 + 3.61 * p_kbar) * (xgrs ** 2 + b)
        - (254.6 + 8.42 * p_kbar) * (2 * xprp * xalm - xalm ** 2 + c)
        - 83.6 * (xprp ** 2 - 2 * xprp * xalm + c)
        + 1388 * xsps
        - 462 * (xmg - xfe)
    )
    den = (
        math.log(kd)
        + 1.431
        + 0.695 * (2 * xgrs * xprp + xgrs ** 2 - 2 * a)
        + 0.203 * (xgrs ** 2 - 2 * xgrs * xalm)
        + 0.922 * xsps
    )
    return num / max(den, 1e-8) - T0


def _wz07_terms(xmg: float, xca: float, xfe: float, xmn: float) -> Tuple[float, float, float, float, float, float]:
    ag = (
        5.993 * xmg ** 2
        - 9.67 * xca ** 2
        - 11.006 * xfe * xmg
        + 0.674 * xfe * xca
        - 9.5725 * xmg * xca
        - 22.765 * xmg * xmn
        - 4.6665 * xca * xmn
        + 11.006 * xfe ** 2 * xmg
        - 0.674 * xfe ** 2 * xca
        - 11.986 * xmg ** 2 * xfe
        + 37.966 * xmg ** 2 * xca
        + 46.02 * xmg ** 2 * xmn
        + 19.34 * xca ** 2 * xfe
        + 25.746 * xca ** 2 * xmg
        + 46.02 * xmn ** 2 * xmg
        + 19.145 * xfe * xmg * xca
        + 45.53 * xfe * xmg * xmn
        + 9.333 * xfe * xca * xmn
        + 77.876 * xmg * xca * xmn
    )
    ah = (
        -0.034 * xmg ** 2
        + 0.135 * xca ** 2
        + 0.024 * xmn ** 2
        + 0.1 * xfe * xmg
        + 0.08 * xfe * xca
        + 0.048 * xfe * xmn
        - 0.0515 * xmg * xca
        + 0.07 * xmg * xmn
        + 0.1765 * xca * xmn
        - 0.1 * xfe ** 2 * xmg
        - 0.08 * xfe ** 2 * xca
        - 0.048 * xfe ** 2 * xmn
        + 0.068 * xmg ** 2 * xfe
        - 0.136 * xmg ** 2 * xca
        - 0.076 * xmg ** 2 * xmn
        - 0.27 * xca ** 2 * xfe
        - 0.28 * xca ** 2 * xmg
        - 0.13 * xca ** 2 * xmn
        - 0.048 * xmn ** 2 * xfe
        - 0.076 * xmn ** 2 * xmg
        - 0.13 * xmn ** 2 * xca
        + 0.103 * xfe * xmg * xca
        - 0.14 * xfe * xmg * xmn
        - 0.353 * xfe * xca * xmn
        - 0.414 * xmg * xca * xmn
    )
    ai = (
        -5672 * xmg ** 2
        + 19932 * xca ** 2
        + 1617 * xmn ** 2
        + 23244 * xfe * xmg
        - 2608 * xfe * xca
        + 3234 * xfe * xmn
        + 31326.5 * xmg * xca
        + 45841 * xmg * xmn
        + 12356 * xca * xmn
        - 23244 * xfe ** 2 * xmg
        + 2608 * xfe ** 2 * xca
        - 3234 * xfe ** 2 * xmn
        + 11344 * xmg ** 2 * xfe
        - 132228 * xmg ** 2 * xca
        - 82498 * xmg ** 2 * xmn
        - 39864 * xca ** 2 * xfe
        - 51518 * xca ** 2 * xmg
        - 2850 * xca ** 2 * xmn
        - 3234 * xmn ** 2 * xfe
        - 82498 * xmn ** 2 * xmg
        - 2850 * xmn ** 2 * xca
        - 62653 * xfe * xmg * xca
        - 91682 * xfe * xmg * xmn
        - 24712 * xfe * xca * xmn
        - 177221 * xmg * xca * xmn
    )
    aj = (
        -5.503 * xfe ** 2
        - 12.873 * xca ** 2
        - 23.01 * xmn ** 2
        + 11.986 * xfe * xmg
        - 9.5725 * xfe * xca
        - 22.765 * xfe * xmn
        - 37.966 * xmg * xca
        - 46.02 * xmg * xmn
        - 38.938 * xca * xmn
        + 11.006 * xfe ** 2 * xmg
        - 0.674 * xfe ** 2 * xca
        - 11.986 * xmg ** 2 * xfe
        + 37.966 * xmg ** 2 * xca
        + 46.02 * xmg ** 2 * xmn
        + 19.34 * xca ** 2 * xfe
        + 25.746 * xca ** 2 * xmg
        + 46.02 * xmn ** 2 * xmg
        + 19.145 * xfe * xmg * xca
        + 45.53 * xfe * xmg * xmn
        + 9.333 * xfe * xca * xmn
        + 77.876 * xmg * xca * xmn
    )
    ak = (
        0.05 * xfe ** 2
        + 0.14 * xca ** 2
        + 0.038 * xmn ** 2
        - 0.068 * xfe * xmg
        - 0.0515 * xfe * xca
        + 0.07 * xfe * xmn
        + 0.136 * xmg * xca
        + 0.076 * xmg * xmn
        + 0.207 * xca * xmn
        - 0.1 * xfe ** 2 * xmg
        - 0.08 * xfe ** 2 * xca
        - 0.048 * xfe ** 2 * xmn
        + 0.068 * xmg ** 2 * xfe
        - 0.136 * xmg ** 2 * xca
        - 0.076 * xmg ** 2 * xmn
        - 0.27 * xca ** 2 * xfe
        - 0.28 * xca ** 2 * xmg
        - 0.13 * xca ** 2 * xmn
        - 0.048 * xmn ** 2 * xfe
        - 0.076 * xmn ** 2 * xmg
        - 0.13 * xmn ** 2 * xca
        + 0.103 * xfe * xmg * xca
        - 0.14 * xfe * xmg * xmn
        - 0.353 * xfe * xca * xmn
        - 0.414 * xmg * xca * xmn
    )
    al = (
        11622 * xfe ** 2
        + 25759 * xca ** 2
        + 41249 * xmn ** 2
        - 11344 * xfe * xmg
        + 31326.5 * xfe * xca
        + 45841 * xfe * xmn
        + 132228 * xmg * xca
        + 82498 * xmg * xmn
        + 88610.5 * xca * xmn
        - 23244 * xfe ** 2 * xmg
        + 2608 * xfe ** 2 * xca
        - 3234 * xfe ** 2 * xmn
        + 11344 * xmg ** 2 * xfe
        - 132228 * xmg ** 2 * xca
        - 82498 * xmg ** 2 * xmn
        - 39864 * xca ** 2 * xfe
        - 51518 * xca ** 2 * xmg
        - 2850 * xca ** 2 * xmn
        - 3234 * xmn ** 2 * xfe
        - 82498 * xmn ** 2 * xmg
        - 2850 * xmn ** 2 * xca
        - 62653 * xfe * xmg * xca
        - 91682 * xfe * xmg * xmn
        - 24712 * xfe * xca * xmn
        - 177221 * xmg * xca * xmn
    )
    return ag, ah, ai, aj, ak, al


def t_wz07(row: Row, p_kbar: float) -> float:
    """Wu & Zhao (2007) olivine–garnet. Uncorrected (no Fe3+ iteration)."""
    mg_ol = num0(row, "olivine_Mg")
    fe_ol = num0(row, "olivine_Fe")
    mg_g = num0(row, "garnet_Mg")
    fe_g = num0(row, "garnet_Fe")
    ca_g = num0(row, "garnet_Ca")
    mn_g = num0(row, "garnet_Mn")
    x_ol_mg = mg_ol / (mg_ol + fe_ol + 1e-16)
    x_ol_fe = fe_ol / (mg_ol + fe_ol + 1e-16)
    sum_g = mg_g + fe_g + ca_g + mn_g + 1e-16
    xmg = mg_g / sum_g
    xca = ca_g / sum_g
    xfe = fe_g / sum_g
    xmn = mn_g / sum_g
    ag, ah, ai, aj, ak, al = _wz07_terms(xmg, xca, xfe, xmn)
    ln_term = math.log((xfe / (xmg + 1e-16)) / (x_ol_fe / (x_ol_mg + 1e-16)))
    if not finite(ln_term):
        ln_term = 0.0
    num = (
        58136.3
        + p_kbar * 1e3 * (0.406 - 2 * ah + 2 * ak)
        - 2 * ai
        + 2 * al
        - 23653.8 * (x_ol_fe - x_ol_mg)
    )
    den = 38.743 + 2 * ag - 2 * aj + 6 * R_J * ln_term
    return num / max(den, 1e-8) - T0


def t_bo73(row: Row, _p: float) -> float:
    """Boyd (1973) Ca-in-cpx solvus. Ca# is in mol% of Ca/(Ca+Mg). No pressure."""
    c = px(row, "cpx")
    if not (c["ca"] + c["mg"] > 0):
        return float("nan")
    ca_n = (100 * c["ca"]) / (c["ca"] + c["mg"])
    return -0.7639 * ca_n * ca_n + 33.135 * ca_n + 1097.1


def t_wb73(row: Row, _p: float) -> float:
    """Wood & Banno (1973) two-pyroxene. No pressure."""
    pair = wb_pair(row)
    if not pair:
        return float("nan")
    den = math.log(pair["ac"] / pair["ao"]) - 7.65 * pair["xFe"] + 3.88 * pair["xFe"] ** 2 - 4.6
    if not (den < 0):
        return float("nan")
    return -10202 / den - T0


def t_nw74(row: Row, _p: float) -> float:
    """Nehru & Wyllie (1974) two-pyroxene. No pressure."""
    pair = wb_pair(row)
    if not pair:
        return float("nan")
    den = math.log(pair["ac"] / pair["ao"]) - 1.72
    if not (den < 0):
        return float("nan")
    return -5006 / den - T0


def t_we77(row: Row, _p: float) -> float:
    """Wells (1977) two-pyroxene. No pressure."""
    pair = wb_pair(row)
    if not pair:
        return float("nan")
    den = 3.355 + 2.44 * pair["xFe"] - math.log(pair["ac"] / pair["ao"])
    if not (abs(den) > 1e-8):
        return float("nan")
    return 7341 / den - T0


def t_nb84(row: Row, p_kbar: float) -> float:
    """Nickel & Brey (1984) as transcribed from the workbook.

    Not a verified reading of the paper. Do not use it for interpretation.
    """
    c = px(row, "cpx")
    o = px(row, "opx")
    x_en_c = c["mg"] - c["ca"]
    x_en_o = o["mg"] - o["ca"]
    if not (x_en_c > 0 and x_en_o > 0):
        return float("nan")
    kd = x_en_c / x_en_o
    den = 0.0083143 * math.log(kd) + 0.004431 * c["ca"] ** 2 + 0.0085
    if not (abs(den) > 1e-12):
        return float("nan")
    num = (
        -7
        - p_kbar * 0.06188
        + 34 * o["ca"] ** 2
        - (21.905 - p_kbar * 0.05229) * c["ca"] ** 2
    )
    return num / den - T0


def t_bm85(row: Row, p_kbar: float) -> float:
    """Bertrand & Mercier (1985) two-pyroxene, with the Fe# iteration."""
    c = px(row, "cpx")
    o = px(row, "opx")
    if not (c["fe"] + c["mg"] > 0 and (1 - c["na"]) > 0 and (1 - o["na"]) > 0):
        return float("nan")
    fe = c["fe"] / (c["fe"] + c["mg"])
    ca_o = o["ca"] / (1 - o["na"])
    if not (ca_o < 1):
        return float("nan")
    pg = p_kbar / 10
    t = 1273.0
    for _ in range(50):
        ca_c = c["ca"] / (1 - c["na"]) + (-0.77 + 0.001 * t) * fe
        k = (1 - ca_c) / (1 - ca_o)
        if not (k > 0 and ca_c < 1):
            return float("nan")
        den = 19.31 - 8.314 * math.log(k) - 12.15 * ca_c ** 2
        if not (den > 0):
            return float("nan")
        t_new = (36273 + 399 * pg) / den
        if abs(t_new - t) < 1e-4:
            return t_new - T0
        t = t_new
    return t - T0


def t_bk90_pyx(row: Row, p_kbar: float) -> float:
    """Brey & Köhler (1990) two-pyroxene."""
    c = px(row, "cpx")
    o = px(row, "opx")
    if not (
        c["fe"] + c["mg"] > 0
        and o["fe"] + o["mg"] > 0
        and (1 - c["na"]) > 0
        and (1 - o["na"]) > 0
    ):
        return float("nan")
    xfe_c = c["fe"] / (c["fe"] + c["mg"])
    xfe_o = o["fe"] / (o["fe"] + o["mg"])
    ca_c = c["ca"] / (1 - c["na"])
    ca_o = o["ca"] / (1 - o["na"])
    kd = (1 - ca_c) / (1 - ca_o)
    if not (kd > 0):
        return float("nan")
    num = 23664 + (24.9 + 126.3 * xfe_c) * p_kbar
    den = 13.38 + math.log(kd) ** 2 + 11.59 * xfe_o
    return num / den - T0


def t_rg74(row: Row, p_kbar: float) -> float:
    """Råheim & Green (1974) garnet–cpx."""
    g = gt_cpx_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    return (3686 + 28.35 * p_kbar) / (math.log(g["kd"]) + 2.33) - T0


def t_eg79(row: Row, p_kbar: float) -> float:
    """Ellis & Green (1979) garnet–cpx."""
    g = gt_cpx_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    total = g["mgG"] + g["caG"] + g["feG"]
    if not (total > 0):
        return float("nan")
    xca = g["caG"] / total
    return (3030 + 10.86 * p_kbar + 3104 * xca) / (math.log(g["kd"]) + 1.9034) - T0


def t_po85(row: Row, p_kbar: float) -> float:
    """Powell (1985) garnet–cpx."""
    g = gt_cpx_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    total = g["mgG"] + g["caG"] + g["feG"]
    if not (total > 0):
        return float("nan")
    return (2790 + 10 * p_kbar + 3140 * g["caG"] / total) / (math.log(g["kd"]) + 1.735) - T0


def t_kr88(row: Row, p_kbar: float) -> float:
    """Krogh (1988) garnet–cpx."""
    g = gt_cpx_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    total = g["feG"] + g["mnG"] + g["mgG"] + g["caG"]
    if not (total > 0):
        return float("nan")
    xca = g["caG"] / total
    return (1879 + 10 * p_kbar + 6731 * xca - 6173 * xca * xca) / (math.log(g["kd"]) + 1.393) - T0


def t_ai94(row: Row, p_kbar: float) -> float:
    """Ai (1994) garnet–cpx. Garnet Mg# is on a 0–100 scale."""
    g = gt_cpx_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    total = g["feG"] + g["mnG"] + g["mgG"] + g["caG"]
    if not (total > 0 and (g["mgG"] + g["feG"]) > 0):
        return float("nan")
    xca = g["caG"] / total
    mg100 = (100 * g["mgG"]) / (g["mgG"] + g["feG"])
    num = -1629 * xca * xca + 3648.55 * xca - 6.59 * mg100 + 1987.98 + 17.66 * p_kbar
    return num / (math.log(g["kd"]) + 1.076) - T0


def t_lg88(row: Row, p_kbar: float) -> float:
    """Lee & Ganguly (1988) opx–garnet."""
    g = opx_gt_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    num = 1971 + 11.91 * p_kbar + (1 / R_HA) * (3e3 * g["xCa"] + 3e3 * g["xMn"])
    return num / (math.log(g["kd"]) + 0.96) - T0


def t_ng09(row: Row, p_kbar: float) -> float:
    """Nimis & Grütter (2009) opx–garnet."""
    g = opx_gt_kd(row)
    if not g or not (g["kd"] > 0):
        return float("nan")
    return (1215 + 17.4 * p_kbar + 1495 * (g["xCa"] + g["xMn"])) / (math.log(g["kd"]) + 0.732) - T0


def t_ss81(row: Row, _p: float) -> float:
    """Sachtleben & Seck (1981) as transcribed from the workbook, on opx Ca cations per 6 oxygens.

    The expression is almost flat across mantle orthopyroxene. It is not a
    verified reading of the paper and is not used for interpretation.
    """
    ca = num0(row, "opx_Ca")
    if not (ca > 0):
        return float("nan")
    return 812 + 163.1 * ca - 10 * ca * ca


def _ow79_parts(row: Row) -> Optional[Dict[str, float]]:
    mg_o = num0(row, "olivine_Mg")
    fe_o = num0(row, "olivine_Fe")
    fe_g = num0(row, "garnet_Fe")
    mg_g = num0(row, "garnet_Mg")
    ca_g = num0(row, "garnet_Ca")
    sum_g = fe_g + mg_g + ca_g
    if not ((mg_o + fe_o) > 0 and sum_g > 0 and fe_o > 0 and mg_g > 0):
        return None
    xmg = mg_o / (mg_o + fe_o)
    xfe = fe_o / (fe_o + mg_o)
    xg_mg = mg_g / sum_g
    xg_fe = fe_g / sum_g
    xg_ca = ca_g / sum_g
    kd = (xmg * xg_fe) / (xfe * xg_mg)
    if not (kd > 0):
        return None
    return {"xMg": xmg, "xFe": xfe, "xgMg": xg_mg, "xgFe": xg_fe, "xgCa": xg_ca, "kd": kd}


def _ow79_at(parts: Dict[str, float], p_kbar: float, w: float) -> float:
    num = (
        902
        + w
        + (parts["xMg"] - parts["xFe"]) * (498 + 1.51 * (p_kbar - 30))
        - 98 * (parts["xgMg"] - parts["xgFe"])
        + 1347 * parts["xgCa"]
    )
    return num / (math.log(parts["kd"]) + 0.357) - T0


def t_ow79(row: Row, p_kbar: float) -> float:
    """O'Neill & Wood (1979) olivine–garnet, constant volume term."""
    parts = _ow79_parts(row)
    if not parts:
        return float("nan")
    return _ow79_at(parts, p_kbar, 25)


def t_ow79v(row: Row, p_kbar: float) -> float:
    """O'Neill & Wood (1979) with the ΔV correction of O'Neill (1980)."""
    parts = _ow79_parts(row)
    if not parts:
        return float("nan")
    t0 = _ow79_at(parts, p_kbar, 25)
    if not finite(t0):
        return float("nan")
    dt = t0 + T0 - 1073
    p = p_kbar
    w = (
        -462.5 * ((1.0191 + dt) * 287e-7) * (p - 263e-6 * p * p - 29.76)
        - 262.4 * (1.0292 + dt) * 45e-6 * (p - 39e-5 * p * p - 29.65)
        + 454 * (1.02 + dt) * 284e-7 * (p - 236e-6 * p * p - 29.79)
        + 278.3 * (1.0234 + dt) * 23e-6 * (p - 45e-5 * p * p - 29.6)
    )
    return _ow79_at(parts, p_kbar, w)


def t_nt00(row: Row, p_kbar: float) -> float:
    """Nimis & Taylor (2000) single-cpx enstatite activity."""
    c = px(row, "cpx")
    a_en = (1 - c["ca"] - c["na"]) * (1 - 0.5 * (c["al"] + c["cr"] + c["na"]))
    if not (a_en > 0):
        return float("nan")
    num = 23166 + 39.28 * p_kbar
    den = 13.25 + 15.35 * c["ti"] + 4.5 * c["fe"] - 1.55 * (c["al"] + c["cr"] - c["na"]) + math.log(a_en) ** 2
    if not (abs(den) > 1e-8):
        return float("nan")
    return num / den - T0


def t_as16(row: Row, p_kbar: float) -> float:
    """Ashchepkov empirical correction of Nimis & Taylor (2000). Uses cpx oxide wt%."""
    t = t_nt00(row, p_kbar)
    fe = row.get("cpx_FeO", None)
    try:
        fe_v = float(fe) if fe is not None and str(fe).strip() != "" else float("nan")
    except (TypeError, ValueError):
        fe_v = float("nan")
    if not (finite(t) and finite(fe_v)):
        return float("nan")
    al = num0(row, "cpx_Al2O3")
    cr = num0(row, "cpx_Cr2O3")
    na = num0(row, "cpx_Na2O")
    tk = t + T0
    return -1e-6 * tk * tk + 0.95 * tk + 33.75 * (al + cr) - 55 - 17.9 * (fe_v - 2.5) - 10.9 * na - T0


def t_bu17(row: Row, p_kbar: float) -> float:
    """Bussweiler et al. (2017) Al-in-olivine. Cr# is not used.

    T (°C) = (11245 + 46 P) / (13.68 − ln Al_ppm) − 273
    Al_ppm may be a column, or is taken from olivine Al2O3 wt%.
    """
    al = olivine_al_ppm(row)
    if not (al > 0 and finite(p_kbar)):
        return float("nan")
    den = 13.68 - math.log(al)
    if abs(den) < 1e-8:
        return float("nan")
    t = (11245 + 46 * p_kbar) / den - T0
    return t if finite(t) else float("nan")


# ---------------------------------------------------------------------------
# Barometers. Each returns kbar at a Celsius temperature.
# ---------------------------------------------------------------------------

def p_ca91(row: Row, t_c: float) -> float:
    """Carswell (1991) Al-in-opx, the Nickel & Green (1985) recast. Clipped to 0–80 kbar."""
    if not finite(t_c):
        return float("nan")
    tk = t_c + T0
    al = num0(row, "opx_Al")
    cr = num0(row, "opx_Cr")
    ti = num0(row, "opx_Ti")
    fe = num0(row, "opx_Fe")
    mn = num0(row, "opx_Mn")
    mg = num0(row, "opx_Mg")
    ca = num0(row, "opx_Ca")
    na = num0(row, "opx_Na")
    al_gt = num0(row, "garnet_Al")
    cr_gt = num0(row, "garnet_Cr")
    fe_gt = num0(row, "garnet_Fe")
    mn_gt = num0(row, "garnet_Mn")
    mg_gt = num0(row, "garnet_Mg")
    ca_gt = num0(row, "garnet_Ca")
    ba = (al - cr - 2 * ti + na) / 2
    bb = 1 - ba - cr - ti
    bc = 1 - ca - na - mn
    bd = (mg / (mg + fe + 1e-16)) * bb
    be = bc - bd
    bf = cr
    raw_cj = (al - (na - cr)) / 2 if na > cr + ti else ba
    cj = clip(raw_cj, 1e-12)
    sum_div = ca_gt + mg_gt + fe_gt + mn_gt + 1e-16
    bg = ca_gt / sum_div
    bj = fe_gt / sum_div
    sum_alcr = al_gt + cr_gt + 1e-16
    bh = al_gt / sum_alcr
    bi = cr_gt / sum_alcr
    denom = 183.3 + 178.98 * cj * (1 - cj)
    ln_arg = clip(((1 - bg) ** 3 * bh ** 2) / (bb * bc ** 2 * cj + 1e-16), 1e-16)
    inside = (
        -R_CAL * tk * math.log(ln_arg)
        - 9e3 * bg ** 2
        - 3400 * (2 * bi ** 2 - bd * bf)
        - bg * bi * (90853 - 52.1 * tk)
        - 7590 * bj * bg
        + 5157 * bd * be
        + 6047
        - 3.23 * tk
    )
    return clip((-1 / denom) * inside, 0, 80)


def p_ng85(row: Row, t_c: float) -> float:
    """Nickel & Green (1985) Al-in-opx. XAl = (Al − Cr − 2 Ti + Na) / 2."""
    g = opx_gt_cations(row)
    if not g or not finite(t_c):
        return float("nan")
    opx, gt = g
    return ng_family_pressure((opx["al"] - opx["cr"] - 2 * opx["ti"] + opx["na"]) / 2, opx, gt, t_c)


def p_ta98al(row: Row, t_c: float) -> float:
    """Taylor (1998) Al-in-opx. XAl = (Al − Cr + Ti + Na) / 2."""
    g = opx_gt_cations(row)
    if not g or not finite(t_c):
        return float("nan")
    opx, gt = g
    return ng_family_pressure((opx["al"] - opx["cr"] + opx["ti"] + opx["na"]) / 2, opx, gt, t_c)


def p_bk90(row: Row, t_c: float) -> float:
    """Brey & Köhler (1990) Al-in-opx. Published assembly. Kelvin is T + 273.15."""
    g = opx_gt_cations(row)
    if not g or not finite(t_c):
        return float("nan")
    opx, gt = g
    x_jd = opx["na"] - opx["cr"] - 2 * opx["ti"]
    x_al_m1 = (opx["al"] + x_jd) / 2
    x_al_ts = x_al_m1 if x_jd < 0 else (opx["al"] - x_jd) / 2
    x_mf_m2 = 1 - opx["ca"] - opx["na"] - opx["mn"]
    x_mf_m1 = 1 - x_al_m1 - opx["cr"] - opx["ti"]
    sum_g = gt["ca"] + gt["mg"] + gt["fe"] + gt["mn"]
    sum_al = gt["al"] + gt["cr"]
    if not (
        x_al_ts > 0
        and x_mf_m1 > 0
        and x_mf_m2 > 0
        and (opx["fe"] + opx["mg"]) > 0
        and sum_g > 0
        and sum_al > 0
    ):
        return float("nan")
    x_ca = gt["ca"] / sum_g
    x_al_g = gt["al"] / sum_al
    x_cr_g = gt["cr"] / sum_al
    x_fe_g = gt["fe"] / sum_g
    if not ((1 - x_ca) > 0 and x_al_g > 0):
        return float("nan")
    k = ((1 - x_ca) ** 3 * x_al_g ** 2) / (x_mf_m1 * x_mf_m2 ** 2 * x_al_ts)
    if not (k > 0):
        return float("nan")
    mg_no = opx["mg"] / (opx["mg"] + opx["fe"])
    x_mg_m1 = x_mf_m1 * mg_no
    x_fe_m1 = x_mf_m1 * (1 - mg_no)
    tk = t_c + 273.15
    c1 = (
        -8.3143 * tk * math.log(k)
        - 5510
        + 88.91 * tk
        - 19 * tk ** 1.2
        + 3 * x_ca ** 2 * 82458
        + x_mg_m1 * x_fe_m1 * (80942 - 46.7 * tk)
        - 3 * x_fe_g * x_ca * 17793
        - x_ca * x_cr_g * (1.164e6 - 420.4 * tk)
        - x_fe_g * x_cr_g * (-1.25e6 + 565 * tk)
    )
    c2 = (
        -0.832
        - 8.78e-5 * (tk - 298)
        + 3 * x_ca ** 2 * 3.305
        - x_ca * x_cr_g * 13.45
        + x_fe_g * x_cr_g * 10.5
    )
    c3 = 1.66e-3
    disc = c2 ** 2 + (4 * c3 * c1) / 1000
    if disc < 0:
        return float("nan")
    return (-c2 - math.sqrt(disc)) / (2 * c3)


def p_hg82(row: Row, t_c: float) -> float:
    """Harley & Green (1982) Al-in-opx."""
    g = opx_gt_cations(row)
    if not g or not finite(t_c):
        return float("nan")
    opx, gt = g
    x_al = opx["al"] / 2
    x_fe_g = gt["fe"] / (gt["fe"] + gt["mg"])
    x_ca = gt["ca"] / (gt["fe"] + gt["mg"] + gt["ca"])
    x_fe = opx["fe"] / (opx["fe"] + opx["mg"])
    vol = -(183.3 + 178.98 * x_al * (1 - x_al))
    if not (vol != 0 and x_al * (1 - x_al) > 0 and (1 - x_ca) > 0):
        return float("nan")
    k = (x_al * (1 - x_al)) / (1 - x_ca)
    if not (k > 0):
        return float("nan")
    tk = t_c + T0
    return (1 / vol) * (
        (R_HA * math.log(k) - 2.93) * tk
        + 5650
        + 5157 * (1 - x_al) * (1 - 2 * x_al) * x_fe
        - 6300 * (x_ca * x_fe_g + x_ca ** 2)
    )


def p_ha84(row: Row, t_c: float) -> float:
    """Harley (1984) Al-in-opx barometer."""
    g = opx_gt_cations(row)
    if not g or not finite(t_c):
        return float("nan")
    opx, gt = g
    x_al = opx["al"] / 2
    x_fe = opx["fe"] / (opx["fe"] + opx["mg"])
    x_mg = opx["mg"] / (opx["fe"] + opx["mg"])
    x_fe_g = gt["fe"] / (gt["fe"] + gt["mg"])
    x_mg_g = gt["mg"] / (gt["fe"] + gt["mg"])
    sum3 = gt["fe"] + gt["mg"] + gt["ca"]
    if not (sum3 > 0 and x_al * (1 - x_al) > 0 and x_mg > 0):
        return float("nan")
    x_ca = gt["ca"] / sum3
    k2 = x_mg_g ** 3 / (x_mg ** 3 * x_al * (1 - x_al))
    if not (k2 > 0):
        return float("nan")
    tk = t_c + T0
    return (1 / 206.74) * (
        R_HA * tk * math.log(k2)
        + tk * (0.15 + 0.001507 * (tk - 970))
        - 2467
        - (2458 * (1e3 / tk) - 1261) * (2 * x_fe ** 2)
        - (3525 * (1e3 / tk) - 1667) * (1 - 2 * x_mg * (1 - x_al)) * (1 - x_mg) * (1 - x_al)
        - (4.75 * tk - 6680) * (2 * (1 - x_mg) * (1 - x_al))
        + 920 * ((1 - 2 * x_al) * (1 - x_al) * (1 - x_mg))
        + (5436 - 2.45 * tk) * (1 - x_mg) * (x_mg * (1 - x_al) + x_al)
        + 5700 * (x_ca * (x_ca + x_fe_g))
    )


def p_wo74(row: Row, t_c: float) -> float:
    """Wood (1974) Al-in-opx, as implemented (R in cal, T in °C in the R T term)."""
    g = opx_gt_cations(row)
    if not g or not finite(t_c):
        return float("nan")
    opx, gt = g
    x_al = opx["al"] / 2
    x_al_100 = x_al * 100
    dv = 5e-4 * x_al_100 ** 2 - 0.061 * x_al_100 - 7.7701
    x_fe_g = gt["fe"] / (gt["fe"] + gt["mg"])
    x_mg_g = gt["mg"] / (gt["fe"] + gt["mg"])
    ao = (x_fe_g + x_mg_g) ** 3
    x_fe = opx["fe"] / (opx["fe"] + opx["mg"])
    if not (dv != 0 and x_al * (1 - x_al) > 0 and ao > 0):
        return float("nan")
    return (
        ((R_HA * t_c) / dv) * math.log((x_al * (1 - x_al)) / ao)
        + (7012 - 3.89 * t_c) / dv
        - 10.45 * x_fe * (1 - 2 * x_al)
    )


def _pos(x: float) -> float:
    return x if x > 0 else 0.0


def _droop_cpx_fe3(si, ti, al, cr, fe, mn, mg, ca, na) -> float:
    charge = (
        4 * (_pos(si) + _pos(ti))
        + 3 * (_pos(al) + _pos(cr))
        + 2 * (_pos(fe) + _pos(mn) + _pos(mg) + _pos(ca))
        + _pos(na)
    )
    fe3 = 12 - charge
    if fe3 < 0:
        return 0.0
    return min(fe3, _pos(fe))


def _cpx_sites(si, ti, al, fe, mg, ca, na, mn=0.0, fe3=0.0) -> Dict[str, float]:
    si, ti, al = _pos(si), _pos(ti), _pos(al)
    fe, mg, ca, na, mn = _pos(fe), _pos(mg), _pos(ca), _pos(na), _pos(mn)
    fe3 = _pos(min(fe3, fe))
    fe2 = _pos(fe - fe3)
    al_iv = _pos(2 - si)
    if al_iv > al:
        al_iv = al
    al_m1 = _pos(al - al_iv)
    m1_fixed = al_m1 + ti + fe3
    fe_mg_m1 = _pos(1 - m1_fixed)
    fe_mg_tot = fe2 + mg
    if fe_mg_m1 > fe_mg_tot:
        fe_mg_m1 = fe_mg_tot
    fe_mg_m2 = fe_mg_tot - fe_mg_m1
    r_fe = fe2 / fe_mg_tot if fe_mg_tot > 0 else 0.0
    r_mg = 1 - r_fe
    return {
        "XAl_M1": al_m1,
        "XFe_M1": fe_mg_m1 * r_fe,
        "XMg_M1": fe_mg_m1 * r_mg,
        "XCa_M2": ca,
        "XNa_M2": na,
        "XMg_M2": fe_mg_m2 * r_mg,
        "XFe_M2": fe_mg_m2 * r_fe,
        "XAl_IV": min(max(al_iv / 2, 0), 1),
        "XSi_IV": min(max(si / 2, 0), 1),
    }


def _garnet_rtlna(xfe: float, xmg: float, xca: float, tk: float) -> Tuple[float, float]:
    r = 8.314
    w_camg = 21627 - 5.78 * tk
    w_mgca = 9834 - 5.78 * tk
    w_cafe = 873 - 1.69 * tk
    w_feca = 6773 - 1.69 * tk
    w_mgfe = 2117
    w_femg = 695
    w_tern = ((w_mgca - w_camg) + (w_cafe - w_feca) + (w_femg - w_mgfe)) / 2
    xca, xmg, xfe = _pos(xca), _pos(xmg), _pos(xfe)
    rtln_gr = (
        0.5 * w_camg * xmg * (1 - xca + xmg + 2 * xca * (xca - xmg - 1))
        + 0.5 * w_mgca * xmg * (1 - xca - xmg - 2 * xca * (xca - xmg - 1))
        + 0.5 * w_cafe * xfe * (1 - xca + xfe + 2 * xca * (xca - xfe - 1))
        + 0.5 * w_feca * xfe * (1 - xca - xfe - 2 * xca * (xca - xfe - 1))
        + w_mgfe * xmg * xfe * (xmg - xfe - 0.5)
        + w_femg * xfe * xmg * (xfe - xmg - 0.5)
        + w_tern * xmg * xfe * (1 - 2 * xca)
    )
    rtln_py = (
        0.5 * w_mgca * xca * (1 - xmg + xca + 2 * xmg * (xmg - xca - 1))
        + 0.5 * w_camg * xca * (1 - xmg - xca - 2 * xmg * (xmg - xca - 1))
        + 0.5 * w_mgfe * xfe * (1 - xmg + xfe + 2 * xmg * (xmg - xfe - 1))
        + 0.5 * w_femg * xfe * (1 - xmg - xfe - 2 * xmg * (xmg - xfe - 1))
        + w_cafe * xca * xfe * (xca - xfe - 0.5)
        + w_feca * xfe * xca * (xfe - xca - 0.5)
        + w_tern * xca * xfe * (1 - 2 * xmg)
    )
    rtln_gr *= 3
    rtln_py *= 3

    def rtlnx3(x: float) -> float:
        return -1e12 if x <= 0 else 3 * r * tk * math.log(x)

    return rtln_gr + rtlnx3(xca), rtln_py + rtlnx3(xmg)


def _cpx_rtlna(sites: Dict[str, float], tk: float) -> Tuple[float, float]:
    r = 8.314
    dg_i, dg_ii, dg_iii = -21000.0, -4680.0, -5010.0
    w_mgal, w_feal, w_mgfe = 4849, 3573, 0
    w_cana, w_mgna, w_fena = 11567, -35873, 59319
    w_camg, w_cafe, w_alsi = 24920, -49275, 7038
    xal = sites["XAl_M1"]
    xfe = sites["XFe_M1"]
    xmg = sites["XMg_M1"]
    xca2 = sites["XCa_M2"]
    xna2 = sites["XNa_M2"]
    xmg2 = sites["XMg_M2"]
    xfe2 = sites["XFe_M2"]
    xal_iv = sites["XAl_IV"]
    xsi_iv = sites["XSi_IV"]
    x_di = max(xca2 * xmg * xsi_iv ** 2, 1e-30)
    x_cats = max(4 * xca2 * xal * xal_iv * xsi_iv, 1e-30)
    rt = r * tk
    rtln_m1_di = w_mgal * xal ** 2 + w_mgfe * xfe ** 2 + xfe * xal * (w_mgfe + w_mgal - w_feal)
    rtln_m1_cats = w_mgal * xmg ** 2 + w_feal * xfe ** 2 + xmg * xfe * (w_mgal + w_feal - w_mgfe)
    rtln_m2 = (
        w_camg * xmg2 ** 2
        + w_cana * xna2 ** 2
        + w_cafe * xfe2 ** 2
        + xmg2 * xfe2 * (w_camg + w_cafe - w_mgfe)
        + xmg2 * xna2 * (w_camg + w_cana - w_mgna)
        + xna2 * xfe2 * (w_cana + w_cafe - w_fena)
    )
    recip = xfe * (xfe2 * dg_ii + xmg2 * dg_iii)
    rtln_di = rt * math.log(x_di) + rtln_m2 + rtln_m1_di + w_alsi * xal_iv ** 2 - xal * (xfe2 * dg_i) - recip
    rtln_cats = (
        rt * math.log(x_cats) + rtln_m2 + rtln_m1_cats + w_alsi * xsi_iv ** 2 + (1 - xal) * (xfe2 * dg_i) - recip
    )
    return rtln_di, rtln_cats


def p_be15(row: Row, t_c: float) -> float:
    """Beyer et al. (2015) garnet–cpx. Published activities, not the workbook CQ cell.

    Garnet fractions are Fe–Mg–Ca with all Fe as Fe2+. Clinopyroxene Fe3+ for
    the activity model is Droop charge balance. The thermometer is separate
    and is not given that Fe3+.
    """
    if not finite(t_c):
        return float("nan")
    fe_g = num0(row, "garnet_Fe")
    mg_g = num0(row, "garnet_Mg")
    ca_g = num0(row, "garnet_Ca")
    si_c = num0(row, "cpx_Si")
    if not (si_c > 0.5 and (fe_g + mg_g + ca_g) > 0 and (num0(row, "cpx_Fe") + num0(row, "cpx_Mg")) > 0):
        return float("nan")
    s = fe_g + mg_g + ca_g
    xfe, xmg, xca = fe_g / s, mg_g / s, ca_g / s
    fe3 = _droop_cpx_fe3(
        si_c,
        num0(row, "cpx_Ti"),
        num0(row, "cpx_Al"),
        num0(row, "cpx_Cr"),
        num0(row, "cpx_Fe"),
        num0(row, "cpx_Mn"),
        num0(row, "cpx_Mg"),
        num0(row, "cpx_Ca"),
        num0(row, "cpx_Na"),
    )
    sites = _cpx_sites(
        si_c,
        num0(row, "cpx_Ti"),
        num0(row, "cpx_Al"),
        num0(row, "cpx_Fe"),
        num0(row, "cpx_Mg"),
        num0(row, "cpx_Ca"),
        num0(row, "cpx_Na"),
        num0(row, "cpx_Mn"),
        fe3,
    )
    tk = t_c + T0
    rtln_gr, rtln_py = _garnet_rtlna(xfe, xmg, xca, tk)
    rtln_di, rtln_cats = _cpx_rtlna(sites, tk)
    rtlnk = rtln_cats + rtln_di - (2 / 3) * rtln_gr - (1 / 3) * rtln_py
    dh, ds, dv = 16500.0, 23.5, 0.719
    return -(dh - tk * ds + rtlnk) / dv / 1000


def p_mgproj(row: Row, _t_c: float) -> float:
    """Olivine Mg# pressure projection. Not a mineral barometer.

    Split at Mg# 91.5:
      P_SGP  = −16.78 Mg# + 1575.6
      P_MSGP = −14.49 Mg# + 1356.9
    with P in kbar. Independent of temperature.
    """
    mg = olivine_mg_number(row)
    if not finite(mg):
        return float("nan")
    return -16.78 * mg + 1575.6 if mg >= OLIVINE_MG_SPLIT else -14.49 * mg + 1356.9


# ---------------------------------------------------------------------------
# Registries
# ---------------------------------------------------------------------------

ThermoFn = Callable[[Row, float], float]
BaroFn = Callable[[Row, float], float]


class Calibration:
    def __init__(self, code: str, cite: str, needs: str, fn: Callable, note: str = ""):
        self.code = code
        self.cite = cite
        self.needs = needs
        self.fn = fn
        self.note = note


THERMOMETERS: Dict[str, Calibration] = {}
BAROMETERS: Dict[str, Calibration] = {}


def _add_t(code, cite, needs, fn, note=""):
    THERMOMETERS[code] = Calibration(code, cite, needs, fn, note)


def _add_b(code, cite, needs, fn, note=""):
    BAROMETERS[code] = Calibration(code, cite, needs, fn, note)


_add_t("TA98", "Taylor (1998) two-pyroxene", "cpx + opx", t_ta98, "Preferred peridotite thermometer.")
_add_t("NG10", "Nimis & Grütter (2010) Ca-in-opx", "opx", t_ng10, "Polynomial in Brey & Köhler Ca-in-opx.")
_add_t("BK90", "Brey & Köhler (1990) Ca-in-opx", "opx", t_bk90_ca)
_add_t("BK90pyx", "Brey & Köhler (1990) two-pyroxene", "cpx + opx", t_bk90_pyx)
_add_t("BM85", "Bertrand & Mercier (1985) two-pyroxene", "cpx + opx", t_bm85)
_add_t("HA84", "Harley (1984) opx–garnet", "opx + garnet", t_ha84)
_add_t("LG88", "Lee & Ganguly (1988) opx–garnet", "opx + garnet", t_lg88)
_add_t("NG09", "Nimis & Grütter (2009) opx–garnet", "opx + garnet", t_ng09)
_add_t("NA09", "Nakamura (2009) garnet–cpx", "cpx + garnet", t_na09, "Published equation, including Mn.")
_add_t("KR00", "Krogh Ravna (2000) garnet–cpx", "cpx + garnet", t_kr00)
_add_t("KR88", "Krogh (1988) garnet–cpx", "cpx + garnet", t_kr88)
_add_t("AI94", "Ai (1994) garnet–cpx", "cpx + garnet", t_ai94, "Garnet Mg# is 0–100.")
_add_t("EG79", "Ellis & Green (1979) garnet–cpx", "cpx + garnet", t_eg79)
_add_t("PO85", "Powell (1985) garnet–cpx", "cpx + garnet", t_po85)
_add_t("RG74", "Råheim & Green (1974) garnet–cpx", "cpx + garnet", t_rg74)
_add_t("WZ07", "Wu & Zhao (2007) olivine–garnet", "olivine + garnet", t_wz07, "Raw thermometer. No Fe3+ iteration.")
_add_t("OW79", "O'Neill & Wood (1979) olivine–garnet", "olivine + garnet", t_ow79)
_add_t("OW79v", "O'Neill & Wood (1979) with ΔV", "olivine + garnet", t_ow79v)
_add_t("NT00", "Nimis & Taylor (2000) single cpx", "cpx", t_nt00)
_add_t("AS16", "Ashchepkov correction of NT00", "cpx oxides", t_as16, "Needs cpx oxide wt%, not only cations.")
_add_t("BO73", "Boyd (1973) Ca-in-cpx", "cpx", t_bo73, "No pressure term. Ca# in mol%.")
_add_t("WB73", "Wood & Banno (1973) two-pyroxene", "cpx + opx", t_wb73, "No pressure term.")
_add_t("WE77", "Wells (1977) two-pyroxene", "cpx + opx", t_we77, "No pressure term.")
_add_t("NW74", "Nehru & Wyllie (1974) two-pyroxene", "cpx + opx", t_nw74, "No pressure term.")
_add_t("BU17", "Bussweiler et al. (2017) Al-in-olivine", "olivine Al", t_bu17, "No Cr# term. Needs Al ppm or olivine Al2O3.")
_add_t("SS81", "Sachtleben & Seck (1981), workbook form", "opx", t_ss81, "NOT VERIFIED. Do not interpret.")
_add_t("NB84", "Nickel & Brey (1984), workbook form", "cpx + opx", t_nb84, "NOT VERIFIED. Do not interpret.")

_add_b("CA91", "Carswell (1991) Al-in-opx", "opx + garnet", p_ca91, "Preferred peridotite barometer. P clipped to 0–80 kbar.")
_add_b("NG85", "Nickel & Green (1985) Al-in-opx", "opx + garnet", p_ng85)
_add_b("TA98Al", "Taylor (1998) Al-in-opx", "opx + garnet", p_ta98al)
_add_b("BK90", "Brey & Köhler (1990) Al-in-opx", "opx + garnet", p_bk90, "Uses 273.15. Published assembly.")
_add_b("HG82", "Harley & Green (1982) Al-in-opx", "opx + garnet", p_hg82)
_add_b("HA84", "Harley (1984) Al-in-opx", "opx + garnet", p_ha84)
_add_b("WO74", "Wood (1974) Al-in-opx", "opx + garnet", p_wo74, "Comparison only. The R T term uses °C, matching the website; do not interpret a result of hundreds of kbar.")
_add_b("BE15", "Beyer et al. (2015) garnet–cpx", "cpx + garnet", p_be15, "Preferred eclogite barometer, with NA09.")
_add_b("MGP", "Olivine Mg# projection", "olivine", p_mgproj, "Not a mineral barometer. Split at Mg# 91.5.")


# ---------------------------------------------------------------------------
# Iteration
# ---------------------------------------------------------------------------

def iterate_pair(
    row: Row,
    thermo: str,
    baro: str,
    p_start: Optional[float] = None,
    t_start: float = 1000.0,
    max_iter: Optional[int] = None,
    tol: Optional[float] = None,
) -> Dict[str, object]:
    """Solve T(P) and P(T) together.

    Peridotite-style pairs use the website defaults (start 40 kbar / 1000 °C,
    25 iterations, 0.1 tolerance). Beyer (BE15) uses the eclogite defaults
    (start 50 kbar, 1000 iterations, 0.001 tolerance).
    """
    thermo_u = thermo.upper()
    baro_u = baro.upper()
    # Accept a few friendly aliases.
    aliases_t = {"NAKAMURA": "NA09", "TAYLOR": "TA98", "BUSSWEILER": "BU17", "ALOL": "BU17"}
    aliases_b = {"CARSWELL": "CA91", "BEYER": "BE15", "MG": "MGP"}
    thermo_u = aliases_t.get(thermo_u, thermo_u)
    baro_u = aliases_b.get(baro_u, baro_u)
    if thermo_u not in THERMOMETERS:
        raise KeyError("Unknown thermometer %s. Use --list." % thermo)
    if baro_u not in BAROMETERS:
        raise KeyError("Unknown barometer %s. Use --list." % baro)
    if baro_u == "BE15":
        p0 = 50.0 if p_start is None else p_start
        n_max = 1000 if max_iter is None else max_iter
        tol_v = 0.001 if tol is None else tol
    else:
        p0 = 40.0 if p_start is None else p_start
        n_max = 25 if max_iter is None else max_iter
        tol_v = 0.1 if tol is None else tol
    t_fn = THERMOMETERS[thermo_u].fn
    p_fn = BAROMETERS[baro_u].fn
    p = p0
    t = t_start
    converged = False
    it = 0
    for it in range(1, n_max + 1):
        t_new = t_fn(row, p)
        p_new = p_fn(row, t_new)
        if not (finite(t_new) and finite(p_new)):
            t, p = t_new, p_new
            break
        if abs(t_new - t) < tol_v and abs(p_new - p) < tol_v:
            t, p = t_new, p_new
            converged = True
            break
        t, p = t_new, p_new
    return {
        "thermometer": thermo_u,
        "barometer": baro_u,
        "T_C": t,
        "P_kbar": p,
        "converged": converged,
        "iterations": it,
    }


def evaluate_all(row: Row, t_c: float, p_kbar: float) -> Dict[str, float]:
    """Every thermometer at one pressure, and every barometer at one temperature."""
    out: Dict[str, float] = {}
    for code, cal in THERMOMETERS.items():
        try:
            out["T_" + code] = cal.fn(row, p_kbar)
        except (ValueError, ZeroDivisionError, OverflowError):
            out["T_" + code] = float("nan")
    for code, cal in BAROMETERS.items():
        try:
            out["P_" + code] = cal.fn(row, t_c)
        except (ValueError, ZeroDivisionError, OverflowError):
            out["P_" + code] = float("nan")
    return out


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def _canon_header(h: str) -> str:
    s = h.strip().lstrip("\ufeff")
    key = s.lower().replace(" ", "").replace("-", "")
    meta = {
        "sampleid": "SampleID",
        "sample": "SampleID",
        "id": "SampleID",
        "location": "Location",
        "locality": "Location",
        "craton": "Craton",
        "lithology": "Lithology",
        "rock": "Lithology",
        "reference": "Reference",
        "ref": "Reference",
        "olivinealppm": "olivine_Al_ppm",
        "olalppm": "olivine_Al_ppm",
        "alppm": "olivine_Al_ppm",
        "al_ppm": "olivine_Al_ppm",
    }
    if key in meta:
        return meta[key]
    if "_" not in s:
        return s
    pre, rest = s.split("_", 1)
    pre_l = pre.lower()
    if pre_l not in PREFIX_ALIAS:
        return s
    prefix = PREFIX_ALIAS[pre_l]
    rest_key = rest.strip()
    for ox in OXIDES:
        if rest_key.lower() == ox.lower():
            return prefix + "_" + ox
    cations = {"si": "Si", "ti": "Ti", "al": "Al", "cr": "Cr", "fe": "Fe", "mn": "Mn", "mg": "Mg", "ca": "Ca", "na": "Na", "ni": "Ni"}
    if rest_key.lower() in cations:
        return prefix + "_" + cations[rest_key.lower()]
    if rest_key.lower() in ("alppm", "al_ppm"):
        return "olivine_Al_ppm"
    return prefix + "_" + rest_key


def load_csv(path: str) -> List[Row]:
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise SystemExit("CSV has no header row.")
        mapping = {h: _canon_header(h) for h in reader.fieldnames}
        rows: List[Row] = []
        for raw in reader:
            row: Row = {}
            for src, dst in mapping.items():
                val = raw.get(src)
                if val is None:
                    continue
                text = val.strip()
                if text == "" or text.lower() in ("na", "nan", "null", "none"):
                    continue
                try:
                    row[dst] = float(text)
                except ValueError:
                    row[dst] = text
            if any(str(v).strip() for v in row.values()):
                rows.append(oxides_to_cations(row))
    return rows


def _cell(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if not math.isfinite(v):
            return ""
        return repr(v) if abs(v) >= 1e6 else ("%s" % v)
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def write_csv(path: Optional[str], rows: Sequence[Dict[str, object]], fields: Sequence[str]) -> None:
    fh = sys.stdout if path in (None, "-") else open(path, "w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(fh, fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _cell(row.get(k)) for k in fields})
    finally:
        if fh is not sys.stdout:
            fh.close()


def solve_table(
    rows: Sequence[Row],
    thermo: str,
    baro: str,
    also_all: bool = False,
    p_start: Optional[float] = None,
    t_start: float = 1000.0,
    max_iter: Optional[int] = None,
    tol: Optional[float] = None,
) -> Tuple[List[Dict[str, object]], List[str]]:
    out_rows: List[Dict[str, object]] = []
    extra: List[str] = []
    if also_all:
        extra = ["T_" + c for c in THERMOMETERS] + ["P_" + c for c in BAROMETERS]
    for row in rows:
        solved = iterate_pair(row, thermo, baro, p_start, t_start, max_iter, tol)
        rec: Dict[str, object] = {}
        for k in META_KEYS:
            if k in row:
                rec[k] = row[k]
        rec.update(solved)
        if also_all and finite(float(solved["T_C"])) and finite(float(solved["P_kbar"])):
            rec.update(evaluate_all(row, float(solved["T_C"]), float(solved["P_kbar"])))
        elif also_all:
            rec.update({k: float("nan") for k in extra})
        out_rows.append(rec)
    fields = [k for k in META_KEYS if any(k in r for r in rows)]
    fields += ["thermometer", "barometer", "T_C", "P_kbar", "converged", "iterations"]
    fields += extra
    return out_rows, fields


def print_catalog() -> None:
    print("Thermometers  (pass the code to -t)")
    print("%-8s  %-42s  %s" % ("code", "calibration", "needs"))
    for code, cal in THERMOMETERS.items():
        mark = "  **" if cal.note.startswith("NOT") else ""
        print("%-8s  %-42s  %s%s" % (code, cal.cite, cal.needs, mark))
        if cal.note:
            print("          " + cal.note)
    print()
    print("Barometers  (pass the code to -b)")
    print("%-8s  %-42s  %s" % ("code", "calibration", "needs"))
    for code, cal in BAROMETERS.items():
        print("%-8s  %-42s  %s" % (code, cal.cite, cal.needs))
        if cal.note:
            print("          " + cal.note)
    print()
    print("Preferred pairs:  -t TA98 -b CA91     peridotite")
    print("                  -t NA09 -b BE15     eclogite")
    print("                  -t BU17 -b MGP      olivine xenocryst")
    print("** SS81 and NB84 are workbook transcriptions, not verified papers.")


def prompt_code(title: str, options: Dict[str, Calibration]) -> str:
    print(title)
    codes = list(options)
    for i, code in enumerate(codes, 1):
        print("  %2d  %-8s  %s" % (i, code, options[code].cite))
    while True:
        raw = input("Number or code: ").strip()
        if raw.upper() in options:
            return raw.upper()
        if raw.isdigit() and 1 <= int(raw) <= len(codes):
            return codes[int(raw) - 1]
        print("Not in the list.")


# ---------------------------------------------------------------------------
# Self-check: Jericho 22-4 against the website gold values
# ---------------------------------------------------------------------------

JERICHO_22_4: Row = {
    "SampleID": "22-4",
    "Location": "Jericho",
    "Lithology": "Sp-gt peridotite",
    "Reference": "Kopylova et al. 1999",
    "cpx_SiO2": 54.4,
    "cpx_Al2O3": 1.31,
    "cpx_Cr2O3": 1.13,
    "cpx_FeO": 1.97,
    "cpx_MnO": 0.05,
    "cpx_MgO": 16.9,
    "cpx_CaO": 22.2,
    "cpx_Na2O": 1.19,
    "opx_SiO2": 57.6,
    "opx_Al2O3": 0.49,
    "opx_Cr2O3": 0.17,
    "opx_FeO": 5.01,
    "opx_MnO": 0.11,
    "opx_MgO": 36.1,
    "opx_CaO": 0.21,
    "opx_Na2O": 0.03,
    "opx_NiO": 0.08,
    "garnet_SiO2": 41.0,
    "garnet_Al2O3": 20.0,
    "garnet_Cr2O3": 5.18,
    "garnet_FeO": 8.86,
    "garnet_MnO": 0.49,
    "garnet_MgO": 18.0,
    "garnet_CaO": 6.36,
    "olivine_SiO2": 40.9,
    "olivine_FeO": 8.16,
    "olivine_MnO": 0.11,
    "olivine_MgO": 50.7,
    "olivine_NiO": 0.39,
}


def self_check() -> int:
    row = oxides_to_cations(dict(JERICHO_22_4))
    solved = iterate_pair(row, "TA98", "CA91")
    p = float(solved["P_kbar"])
    t = float(solved["T_C"])
    cation_tol = {
        "cpx_Si": (1.9869, 5e-4),
        "opx_Al": (0.0198, 5e-4),
        "garnet_Fe": (0.539, 5e-4),
        "olivine_Mg": (1.8366, 5e-4),
    }
    fail = 0
    for name, (gold, tol) in cation_tol.items():
        got = float(row[name])
        if abs(got - gold) > tol:
            print("FAIL", name, got, "expected", gold)
            fail += 1
        else:
            print("ok  ", name, got)
    temp_gold = {
        "T_TA98": (t, 794.1390116862908),
        "P_CA91": (p, 37.157514567054065),
        "T_BK90": (t_bk90_ca(row, p), 830.7960811926337),
        "T_NG10": (t_ng10(row, p), 777.546467765717),
        "T_HA84": (t_ha84(row, p), 806.8593371807335),
        "T_KR00": (t_kr00(row, p), 815.7227988053044),
        "T_NA09": (t_na09(row, p), 831.5849965450154),
        "T_WZ07": (t_wz07(row, p), 866.0329903377851),
    }
    for name, (got, gold) in temp_gold.items():
        if abs(got - gold) > 0.05:
            print("FAIL", name, got, "expected", gold, "d", got - gold)
            fail += 1
        else:
            print("ok  ", name, round(got, 3))
    if fail:
        print("self-check failed:", fail)
        return 1
    print("self-check ok  Jericho 22-4  T=%.1f °C  P=%.2f kbar" % (t, p))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Solve a thermometer–barometer pair from a mineral-chemistry CSV."
    )
    parser.add_argument("csv", nargs="?", help="Mineral chemistry CSV (oxide wt%).")
    parser.add_argument("-t", "--thermometer", help="Thermometer code. See --list.")
    parser.add_argument("-b", "--barometer", help="Barometer code. See --list.")
    parser.add_argument("-o", "--output", help="Output CSV. Default: print to the terminal.")
    parser.add_argument("--list", action="store_true", help="Print the calibrations and exit.")
    parser.add_argument("--interactive", action="store_true", help="Ask which pair to use.")
    parser.add_argument("--also-all", action="store_true", help="Also report every calibration at the solved P and T.")
    parser.add_argument("--pressure", type=float, default=None, help="Starting pressure in kbar (default 40, or 50 with BE15).")
    parser.add_argument("--temperature", type=float, default=1000.0, help="Starting temperature in °C (default 1000).")
    parser.add_argument("--max-iter", type=int, default=None)
    parser.add_argument("--tol", type=float, default=None, help="Stop when |ΔT| and |ΔP| are both below this.")
    parser.add_argument("--self-check", action="store_true", help="Check Jericho 22-4 against the website.")
    args = parser.parse_args(argv)

    if args.list:
        print_catalog()
        return 0
    if args.self_check:
        return self_check()
    if not args.csv:
        parser.error("give a CSV, or use --list / --self-check")

    thermo = args.thermometer
    baro = args.barometer
    if args.interactive or (not thermo and not baro and sys.stdin.isatty()):
        print_catalog()
        print()
        thermo = prompt_code("Choose a thermometer", THERMOMETERS)
        baro = prompt_code("Choose a barometer", BAROMETERS)
    if not thermo or not baro:
        parser.error("choose both -t THERMOMETER and -b BAROMETER (or --interactive)")

    rows = load_csv(args.csv)
    if not rows:
        raise SystemExit("No data rows in %s" % args.csv)
    solved, fields = solve_table(
        rows,
        thermo,
        baro,
        also_all=args.also_all,
        p_start=args.pressure,
        t_start=args.temperature,
        max_iter=args.max_iter,
        tol=args.tol,
    )
    write_csv(args.output, solved, fields)
    if args.output and args.output != "-":
        n_ok = sum(1 for r in solved if r.get("converged"))
        print("Wrote %d rows (%d converged) to %s" % (len(solved), n_ok, args.output), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
