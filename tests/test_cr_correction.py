"""Eq. 4 (v2 method): TA98 + garnet-Cr2O3 correction of TA98 - WZ07 with a flat low-T level (calculations/pt_fo2/cr_correction.py)."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from calculations import run_full_pipeline  # noqa: E402
from calculations.pt_fo2 import cr_correction as crc  # noqa: E402
from calculations.pt_fo2.fe3_correction import _wz07corr_fe2only  # noqa: E402
from calculations.pt_fo2.params import PARAMS  # noqa: E402
import fit_cr_correction as fcc  # noqa: E402


@pytest.fixture(scope="module")
def table():
    assert PARAMS.fe3.wz07_correction == "cr2o3"
    t, _ = run_full_pipeline(fcc.INPUT, verbose=False)
    return t


def test_params_are_the_module_coefficients():
    fp = PARAMS.fe3
    assert (fp.cr_a, fp.cr_b, fp.cr_c, fp.cr_ref, fp.cr_pivot_c) == (crc.CR_A, crc.CR_B, crc.CR_C, crc.CR_REF, crc.CR_PIVOT)
    assert (fp.cr_min_wt, fp.cr_max_wt, fp.cr_tmin_c, fp.cr_tmax_c) == (1.3, 9.8, 1124.0, 1394.0)
    assert (fp.cr_d_low_c, fp.cr_low_t_c) == (crc.CR_D_LOW, crc.T_LOW) == (-113.24224691422943, 1124.0)
    assert fp.wz07_offset_c == 0 and fp.wz07_offset_slope == 0


def test_level_is_the_wls_intercept_folded_into_eq4():
    """The adopted level is CR_A_FIT (not -102) and dT_match is +10 C on the Eq. 4 scale."""
    fp = PARAMS.fe3
    assert crc.CR_A == crc.CR_A_FIT == fp.cr_a == pytest.approx(-92.7734, abs=1e-4)
    assert fp.dt_match_c == crc.DT_MATCH == 10.0 and fp.emp_fit_dt_c == 10.0
    assert crc.predicted_discrepancy(crc.CR_PIVOT, 2.0) == pytest.approx(crc.CR_A_FIT)
    # the pipeline subtracts the whole correction, level included, from WZ07
    from calculations.pt_fo2.fe3_correction import wz07_offset
    from calculations.pt_fo2.level_anchor import curve
    for T, cr in ((1000.0, 5.0), (1123.9, 2.0), (1124.0, 2.0), (1207.0, 8.0), (1450.0, 1.0)):
        assert -wz07_offset(fp, T, cr) == pytest.approx(crc.predicted_discrepancy(T, cr))
        assert curve(T, cr2o3=cr) == pytest.approx(crc.predicted_discrepancy(T, cr))


def test_formula():
    T, cr = 1300.0, 2.0
    t = T - crc.CR_PIVOT
    assert crc.predicted_discrepancy(T, cr) == pytest.approx(crc.CR_A + crc.CR_B * t + crc.CR_C * (cr - crc.CR_REF) * t)
    # rounded form quoted in the paper: -92.77 - 0.285 t + 0.067 (Cr2O3 - 5.36) t, t = TA98 - 1207.02
    assert crc.predicted_discrepancy(T, cr) == pytest.approx(-92.77 - 0.285 * t + 0.067 * (cr - 5.36) * t, abs=0.5)
    assert crc.corrected_wz07(1000.0, T, cr) == pytest.approx(1000.0 + crc.predicted_discrepancy(T, cr))


def test_low_t_branch_is_flat():
    """Below 1124 C: T_WZ07corr = T_WZ07 - 113.24, independent of TA98 and Cr2O3."""
    f = crc.predicted_discrepancy
    for T in (900.0, 1000.0, 1123.99):
        for cr in (0.5, 2.0, 5.36, 9.0, np.nan, None):
            assert f(T, cr) == pytest.approx(-113.24224691422943, abs=1e-9)
    assert crc.corrected_wz07(1000.0, 1050.0, 3.0) == pytest.approx(1000.0 - 113.24, abs=0.01)
    # at and above 1124 C the Eq. 4 slope applies (a step at the boundary, no interpolation)
    t = 1124.0 - crc.CR_PIVOT
    assert f(1124.0, 2.0) == pytest.approx(crc.CR_A + crc.CR_B * t + crc.CR_C * (2.0 - crc.CR_REF) * t)
    arr = f(np.array([1000.0, 1124.0, 1300.0]), np.array([2.0, 2.0, 2.0]))
    assert arr[0] == pytest.approx(crc.CR_D_LOW) and arr[1] == pytest.approx(f(1124.0, 2.0)) and arr[2] == pytest.approx(f(1300.0, 2.0))
    # cr_d_low_c = None restores the fe606c1 behaviour (slope held at 1124 C)
    from dataclasses import replace
    old = replace(PARAMS.fe3, cr_d_low_c=None)
    assert f(1000.0, 5.0, old) == pytest.approx(f(1124.0, 5.0, old))


def test_clamping():
    f = crc.predicted_discrepancy
    assert f(1500.0, 5.0) == pytest.approx(f(1394.0, 5.0))
    assert f(1300.0, 0.5) == pytest.approx(f(1300.0, 1.3))
    assert f(1300.0, 12.0) == pytest.approx(f(1300.0, 9.8))
    assert f(1300.0, np.nan) == pytest.approx(f(1300.0, crc.CR_REF))
    assert f(1300.0, None) == pytest.approx(f(1300.0, crc.CR_REF))
    fl = crc.flags(pd.DataFrame({"SampleID": list("abcde"), "T_TA98": [1000, 1200, 1500, 1200, 1200],
                                 "garnet_Cr2O3": [5, 5, 5, 0.5, np.nan]}))
    assert fl.T_low.tolist() == [True, False, False, False, False]
    assert fl.T_high.tolist() == [False, False, True, False, False]
    assert fl.Cr_low.tolist() == [False, False, False, True, False]
    assert fl.no_Cr.tolist() == [False, False, False, False, True]
    assert fl.flag.tolist() == [True, False, True, True, True]


def test_fit_reproduces_coefficients(table):
    r = fcc.fit(table)
    assert r["n_natural"] == 37 and r["n_experiments"] == 8
    assert r["a"]["value"] == pytest.approx(crc.CR_A_FIT, abs=1e-6)
    assert r["b"]["value"] == pytest.approx(crc.CR_B, abs=1e-9)
    assert r["c"]["value"] == pytest.approx(crc.CR_C, abs=1e-9)
    assert r["cr_ref_wt"] == pytest.approx(crc.CR_REF, abs=1e-12)
    assert r["c"]["p"] < 0.001
    assert r["d_low"]["n"] == 8
    assert r["d_low"]["value"] == pytest.approx(crc.CR_D_LOW, abs=1e-6)
    # slope check against the experiments: predicted -0.5566 vs measured weighted -0.5547
    assert r["slope_check"]["predicted_slope"] == pytest.approx(-0.5566, abs=5e-4)
    assert r["slope_check"]["predicted_slope"] == pytest.approx(r["slope_check"]["measured_weighted_slope_13_runs"], abs=0.005)


def test_iteration_stops_at_the_cr_correction(table):
    """Every screened sample above the floor stops at the first 0.001 step that meets T_TA98 - T_WZ07corr <= +10 C,
    T_WZ07corr = WZ07 + C(T, Cr)."""
    s = table[table.screen_pass == True]
    x = pd.to_numeric(s.XFe3_gt_iter, errors="coerce").to_numpy(float)
    num = lambda c: pd.to_numeric(s[c], errors="coerce").to_numpy(float)
    T, cr = num("T_TA98"), num("garnet_Cr2O3")
    C = crc.predicted_discrepancy(T, cr) + 10.0
    wz = lambda xf: _wz07corr_fe2only(num("olivine_Mg"), num("olivine_Fe"), num("garnet_Mg"), num("garnet_Fe") * (1 - xf),
                                      num("garnet_Ca"), num("garnet_Mn"), num("P_CA91"), offset_c=0.0)
    assert np.all(T - wz(x) <= C + 1e-9)
    above = x > PARAMS.fe3.xfe3_start + 1e-9
    assert above.sum() > 200
    assert np.all((T - wz(x - PARAMS.fe3.dxfe3))[above] > C[above])
    assert np.allclose(num("T_WZ07corr"), wz(x) + C - 10.0)
    low = T < 1124.0
    assert low.sum() == 231 and ((x > PARAMS.fe3.xfe3_start + 1e-9) & low).sum() > 100


def test_screened_flag_counts(table):
    s = table[table.screen_pass == True]
    f = crc.flags(s)
    assert len(s) == 424
    assert (f.T_low.sum(), f.T_high.sum(), f.Cr_low.sum(), f.Cr_high.sum(), f.no_Cr.sum()) == (231, 6, 2, 2, 0)
