"""Analytical-quality screen of Nimis and Grütter (2010, p. 413) in calculations/pt_fo2/screening.py."""
import sys
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from calculations.pt_fo2.cations import convert_oxides_to_cations  # noqa: E402
from calculations.pt_fo2.screening import analytical_quality, oxide_total  # noqa: E402

# A good garnet lherzolite assemblage (wt%); each test changes one analysis.
GOOD = {
    "opx_SiO2": 57.6, "opx_TiO2": 0.05, "opx_Al2O3": 0.80, "opx_Cr2O3": 0.30, "opx_FeO": 4.60, "opx_MnO": 0.11,
    "opx_MgO": 35.3, "opx_CaO": 0.60, "opx_Na2O": 0.15, "opx_NiO": 0.10,
    "cpx_SiO2": 54.9, "cpx_TiO2": 0.10, "cpx_Al2O3": 2.20, "cpx_Cr2O3": 1.50, "cpx_FeO": 2.20, "cpx_MnO": 0.08,
    "cpx_MgO": 16.9, "cpx_CaO": 19.8, "cpx_Na2O": 1.90, "cpx_NiO": 0.05,
    "garnet_SiO2": 41.9, "garnet_TiO2": 0.10, "garnet_Al2O3": 20.0, "garnet_Cr2O3": 4.50, "garnet_FeO": 7.00,
    "garnet_MnO": 0.35, "garnet_MgO": 20.5, "garnet_CaO": 5.20, "garnet_Na2O": 0.04, "garnet_NiO": 0.01,
}


def _quality(**changes):
    row = dict(GOOD, **changes)
    return bool(analytical_quality(convert_oxides_to_cations(pd.DataFrame([row]), verbose=False)).iloc[0])


def test_good_analysis_passes():
    assert _quality()


def test_oxide_total_limits_use_reported_total_where_present():
    assert not _quality(opx_Total=98.4) and not _quality(cpx_Total=101.6) and not _quality(garnet_Total=101.6)
    assert _quality(opx_Total=98.5) and _quality(garnet_Total=101.5)
    assert not _quality(cpx_SiO2=52.0)                       # no reported total: oxide sum 96.7 wt%
    df = pd.DataFrame([dict(GOOD, garnet_Total=500.05)])     # not an oxide total: the oxide sum is used
    assert oxide_total(df, "garnet").iloc[0] == pytest.approx(sum(v for k, v in GOOD.items() if k.startswith("garnet_")))


def test_pyroxene_cation_sum_on_6_oxygens():
    assert not _quality(opx_SiO2=59.6, opx_MgO=33.3)         # Si-rich, Mg-poor opx: cation sum < 3.990
    assert not _quality(cpx_SiO2=56.9, cpx_MgO=14.9)


def test_garnet_cr_or_mg_number():
    assert _quality(garnet_Cr2O3=0.5, garnet_Al2O3=24.0)                                       # low Cr, Mg# 0.84
    assert not _quality(garnet_Cr2O3=0.5, garnet_Al2O3=24.0, garnet_MgO=8.0, garnet_FeO=19.5)  # low Cr, Mg# < 0.6
    assert _quality(garnet_Cr2O3=4.5, garnet_MgO=8.0, garnet_FeO=19.5, garnet_Al2O3=20.0)      # Cr2O3 > 1 wt%


def test_quality_screen_runs_before_four_phase():
    """Quality is an independent first screen; it can pass when olivine is missing."""
    from calculations.pt_fo2.screening import apply_equilibrium_screening
    row = dict(GOOD, olivine_SiO2=0.0, T_TA98=1200, T_NG10=1200, T_HA84=1200, T_WZ07=1200)
    df, stats = apply_equilibrium_screening(convert_oxides_to_cations(pd.DataFrame([row]), verbose=False), verbose=False)
    assert bool(df["screen_quality"].iloc[0])
    assert not bool(df["screen_4phases"].iloc[0])
    assert not bool(df["screen_pass"].iloc[0])
    assert stats["n_quality"] == 1 and stats["n_4phases"] == 0 and stats["n_quality_4phases"] == 0
