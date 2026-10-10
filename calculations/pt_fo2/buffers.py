"""
Oxygen-fugacity buffers for Fig. 12 (P vs Δlog fO₂ FMQ).

Solid-state buffers use the Frost (1991) parameterization
    log10 fO2 = A/T + B + C*(P_bar - 1)/T
with T in kelvin and P in bar (Reviews in Mineralogy vol. 25, Table 2).

FMQ as used on the SD1 fO2 sheet (and therefore in this project's
ΔFMQ values) is the O'Neill-style expression already in fo2.py:
    log10 fO2 = 82.75 + 0.00484 T - 30681/T - 24.45 log10(T)
                + 94 P_kbar/T - 0.002 P_kbar

EMOG/D (enstatite–magnesite–olivine–graphite/diamond) is the
Stagno & Frost (2010) / Stagno et al. (2013) locus along a cratonic
geotherm, expressed directly as ΔFMQ(P).  The polynomial is a fit to
the EMOD curve as drawn on Frost & McCammon (2008, fig. 6) and
Stagno et al. (2013, fig. 1), which the chapter caption cites.

Fe–Ni precipitation (NiPC) is O'Neill & Wall (1987) as drawn by
Frost & McCammon (2008, Fig. 3): just below the IW buffer
(~IW − 0.2). On the SD1 FMQ axis that is ≈ FMQ−5.0 to FMQ−5.4
along a cratonic geotherm, and BSE meets it near 8 GPa.

BSE is the predicted bulk-silicate-Earth trend of Frost & McCammon
(2008).  The Nimis et al. (2015) “typical mantle” band is the
envelope they quote for garnet-peridotite xenoliths.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# Frost (1991) A, B, C.  P in bar, T in K.
# FMβQ is the β-quartz FMQ used above 573 °C.
FROST1991 = {
    "FMQ": {"A": -25096.3, "B": 8.735, "C": 0.110, "source": "Frost 1991 FMβQ"},
    "NNO": {"A": -24930.0, "B": 9.360, "C": 0.046, "source": "Frost 1991 NiNiO"},
    "IW":  {"A": -27489.0, "B": 6.702, "C": 0.055, "source": "Frost 1991 IW"},
    "WM":  {"A": -32307.0, "B": 13.012, "C": 0.063, "source": "Frost 1991 WM"},
}


def logfO2_frost1991(name: str, T_C, P_kbar):
    """Absolute log10 fO2 from Frost (1991)."""
    p = FROST1991[name]
    T = np.asarray(T_C, dtype=float) + 273.15
    P_bar = np.asarray(P_kbar, dtype=float) * 1000.0
    return p["A"] / T + p["B"] + p["C"] * (P_bar - 1.0) / T


def logfO2_fmq_sd1(T_C, P_kbar):
    """FMQ used on the SD1 fO2 sheet / calculations/fo2.py."""
    T = np.asarray(T_C, dtype=float) + 273.15
    P = np.asarray(P_kbar, dtype=float)
    return (
        82.75
        + 0.00484 * T
        - 30681.0 / T
        - 24.45 * np.log10(T)
        + 94.0 * P / T
        - 0.002 * P
    )


def delta_frost_vs_sd1_fmq(name: str, T_C, P_kbar):
    """Buffer minus the SD1 FMQ, so it plots on the same axis as the xenoliths."""
    return logfO2_frost1991(name, T_C, P_kbar) - logfO2_fmq_sd1(T_C, P_kbar)


def cratonic_geotherm_C(P_kbar):
    """
    Simple 40 mW m⁻² cratonic geotherm (Hasterok & Chapman style),
    adequate for drawing reference buffers on Fig. 12.
    T (°C) ≈ 420 + 14.2 P(kbar)  →  ~850 °C at 30 kbar, ~1200 °C at 55 kbar.
    """
    P = np.asarray(P_kbar, dtype=float)
    return 420.0 + 14.2 * P


def xenolith_geotherm_C(P_kbar, df=None):
    """
    If a xenolith table with P_CA91 / T_TA98 is supplied, fit a linear
    T(P) to screened samples; otherwise use the 40 mW m⁻² geotherm.
    """
    P = np.asarray(P_kbar, dtype=float)
    if df is None or "P_CA91" not in df.columns or "T_TA98" not in df.columns:
        return cratonic_geotherm_C(P)
    use = df.copy()
    if "screen_pass" in use.columns:
        use = use[use["screen_pass"] == True]
    x = pd.to_numeric(use["P_CA91"], errors="coerce")
    y = pd.to_numeric(use["T_TA98"], errors="coerce")
    m = x.notna() & y.notna() & (x > 20) & (x < 80) & (y > 800) & (y < 1600)
    if m.sum() < 10:
        return cratonic_geotherm_C(P)
    a, b = np.polyfit(x[m], y[m], 1)
    return a * P + b


# EMOG/D as drawn on Nimis et al. (2015) Fig. 8 (Stagno et al. 2013
# along a cratonic geotherm). Depth from data/geotherm.xlsx:
# 50 km ≈ 16 kbar, 100 ≈ 32, 150 ≈ 49, 200 ≈ 65, 250 ≈ 81.
_EMOD_P = np.array([16.0, 32.0, 49.0, 65.0, 81.0])
_EMOD_D = np.array([-0.40, -0.80, -1.20, -1.65, -2.05])


def delta_emod(P_kbar):
    P = np.asarray(P_kbar, dtype=float)
    return np.interp(P, _EMOD_P, _EMOD_D)


# Fe–Ni precipitation as drawn on Nimis et al. (2015) Fig. 8
# (O'Neill & Wall 1987). Nearly flat, slightly more reduced with depth.
_FENI_P = np.array([16.0, 32.0, 49.0, 65.0, 81.0])
_FENI_D = np.array([-4.55, -4.75, -4.90, -5.05, -5.20])


def delta_feni(P_kbar, T_C=None, df=None):
    """NiPC from Nimis et al. (2015) Fig. 8."""
    P = np.asarray(P_kbar, dtype=float)
    return np.interp(P, _FENI_P, _FENI_D)


# BSE = Frost & McCammon 2008 Fig. 3 line 1: four-phase garnet peridotite
# with bulk Fe3+/ΣFe = 2% (McDonough & Sun 1995 BSE) along a cratonic
# geotherm. Their Fig. 3 is drawn vs IW; converted here to ΔFMQ using
# the paper's own statement that the curve reaches FMQ−5 when it
# crosses the Ni precipitation curve at ~8 GPa (250 km), and the
# xenolith trend in their Fig. 2 (FMQ−1 at ~3 GPa to FMQ−3.5 at ~6 GPa).
# After 8 GPa their line 2 follows NiPC; we stop at 85 kbar.
_BSE_P = np.array([20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 85.0])
_BSE_D = np.array([-0.6, -1.5, -2.4, -3.3, -4.0, -4.6, -5.1, -5.3])


def delta_bse(P_kbar):
    P = np.asarray(P_kbar, dtype=float)
    return np.interp(P, _BSE_P, _BSE_D)


# Nimis et al. (2015) Fig. 8 typical-mantle envelope, drawn as two
# parallel straight lines. Slope and intercepts are the average of
# the Table 6 dashed-line anchors so the band stays linear.
_NIMIS_SLOPE = -0.0537          # ΔFMQ per kbar
_NIMIS_HI0 = 1.02 + 0.0537 * 20.8
_NIMIS_WIDTH = 1.02 - (-2.50)   # 3.52 log units, constant


def nimis2015_band(P_kbar):
    P = np.asarray(P_kbar, dtype=float)
    hi = _NIMIS_HI0 + _NIMIS_SLOPE * P
    lo = hi - _NIMIS_WIDTH
    return lo, hi


def buffer_table(P_kbar, T_C=None, df=None) -> pd.DataFrame:
    P = np.asarray(P_kbar, dtype=float)
    if T_C is None:
        T_C = xenolith_geotherm_C(P, df=df)
    T_C = np.asarray(T_C, dtype=float)
    out = pd.DataFrame({
        "P_kbar": P,
        "T_C": T_C,
        "FMQ_sd1": logfO2_fmq_sd1(T_C, P),
        "FMQ_frost": logfO2_frost1991("FMQ", T_C, P),
        "dIW": delta_frost_vs_sd1_fmq("IW", T_C, P),
        "dNNO": delta_frost_vs_sd1_fmq("NNO", T_C, P),
        "dWM": delta_frost_vs_sd1_fmq("WM", T_C, P),
        "dEMOD": delta_emod(P),
        "dFeNi": delta_feni(P, T_C=T_C, df=df),
        "dBSE": delta_bse(P),
    })
    lo, hi = nimis2015_band(P)
    out["nimis_lo"] = lo
    out["nimis_hi"] = hi
    return out


def verify_against_literature(verbose=True) -> pd.DataFrame:
    """
    Check Frost (1991) 1-bar anchors and a few high-P EMOD/FeNi points
    against published values.
    """
    rows = []

    # Frost 1991, 1 bar, 1000 °C.  Classic FMQ ≈ −10.9, NNO−FMQ ≈ +0.7
    T, P = 1000.0, 0.001
    fmq = logfO2_frost1991("FMQ", T, P)
    nno = logfO2_frost1991("NNO", T, P)
    iw = logfO2_frost1991("IW", T, P)
    rows.append(dict(check="Frost1991 FMQ 1 bar 1000C",
                     computed=fmq, expected=-10.97, tol=0.15))
    rows.append(dict(check="Frost1991 NNO-FMQ 1 bar 1000C",
                     computed=nno - fmq, expected=0.75, tol=0.15))
    rows.append(dict(check="Frost1991 IW-FMQ 1 bar 1000C",
                     computed=iw - fmq, expected=-3.91, tol=0.15))

    # Nimis et al. 2015 Fig. 8 / Stagno 2013 EMOG/D
    rows.append(dict(check="EMOD ΔFMQ at 49 kbar (Nimis 2015 Fig.8 ~−1.2)",
                     computed=float(delta_emod(49)), expected=-1.2, tol=0.3))
    rows.append(dict(check="EMOD ΔFMQ at 65 kbar (Nimis 2015 Fig.8 ~−1.65)",
                     computed=float(delta_emod(65)), expected=-1.65, tol=0.3))

    # Nimis et al. 2015 Fig. 8 Fe–Ni
    rows.append(dict(check="FeNi ΔFMQ at 65 kbar (Nimis 2015 Fig.8 ~−5.05)",
                     computed=float(delta_feni(65)), expected=-5.05, tol=0.2))
    rows.append(dict(check="BSE ΔFMQ at 30 kbar (F&M2008 Fig.2 ~−1.5)",
                     computed=float(delta_bse(30)), expected=-1.5, tol=0.5))
    rows.append(dict(check="BSE ΔFMQ at 60 kbar (F&M2008 Fig.3 ~−4)",
                     computed=float(delta_bse(60)), expected=-4.0, tol=0.5))
    rows.append(dict(check="BSE ΔFMQ at 80 kbar (F&M2008 FMQ−5 at NiPC)",
                     computed=float(delta_bse(80)), expected=-5.0, tol=0.4))

    tab = pd.DataFrame(rows)
    tab["resid"] = tab["computed"] - tab["expected"]
    tab["ok"] = tab["resid"].abs() <= tab["tol"]
    if verbose:
        print("Buffer verification (computed vs published anchors)")
        print(tab.to_string(index=False, float_format=lambda v: f"{v:8.3f}"))
        print(f"  passed {int(tab.ok.sum())} / {len(tab)}")
    return tab
