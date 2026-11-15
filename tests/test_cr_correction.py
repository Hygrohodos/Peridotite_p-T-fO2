"""Revision v7_cr: TA98 + garnet-Cr2O3 correction of TA98 - WZ07 (calculations/pt_fo2/cr_correction.py)."""
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
    assert fp.dt_match_c == fp.cr_a and fp.wz07_offset_c == 0 and fp.wz07_offset_slope == 0


def test_formula():
    assert crc.predicted_discrepancy(crc.CR_PIVOT, 2.0) == pytest.approx(crc.CR_A)
    T, cr = 1300.0, 2.0
    t = T - crc.CR_PIVOT
    assert crc.predicted_discrepancy(T, cr) == pytest.approx(crc.CR_A + crc.CR_B * t + crc.CR_C * (cr - crc.CR_REF) * t)
    # rounded form quoted in the paper: -102 - 0.285 t + 0.067 (Cr2O3 - 5.36) t
    assert crc.predicted_discrepancy(T, cr) == pytest.approx(-102 - 0.285 * t + 0.067 * (cr - 5.36) * t, abs=0.5)
    assert crc.corrected_wz07(1000.0, T, cr) == pytest.approx(1000.0 + crc.predicted_discrepancy(T, cr))


def test_clamping():
    f = crc.predicted_discrepancy
    assert f(1000.0, 5.0) == pytest.approx(f(1124.0, 5.0))
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
    # slope check against the experiments: predicted -0.5566 vs measured weighted -0.5547
    assert r["slope_check"]["predicted_slope"] == pytest.approx(-0.5566, abs=5e-4)
    assert r["slope_check"]["predicted_slope"] == pytest.approx(r["slope_check"]["measured_weighted_slope_13_runs"], abs=0.005)


def test_iteration_stops_at_the_cr_correction(table):
    """Every screened sample above the floor stops at the first 0.001 step that meets T_TA98 - WZ07 <= C(T, Cr)."""
    s = table[table.screen_pass == True]
    x = pd.to_numeric(s.XFe3_gt_iter, errors="coerce").to_numpy(float)
    num = lambda c: pd.to_numeric(s[c], errors="coerce").to_numpy(float)
    T, cr = num("T_TA98"), num("garnet_Cr2O3")
    C = crc.predicted_discrepancy(T, cr)
    wz = lambda xf: _wz07corr_fe2only(num("olivine_Mg"), num("olivine_Fe"), num("garnet_Mg"), num("garnet_Fe") * (1 - xf),
                                      num("garnet_Ca"), num("garnet_Mn"), num("P_CA91"), offset_c=0.0)
    assert np.all(T - wz(x) <= C + 1e-9)
    above = x > PARAMS.fe3.xfe3_start + 1e-9
    assert above.sum() > 200
    assert np.all((T - wz(x - PARAMS.fe3.dxfe3))[above] > C[above])


def test_screened_flag_counts(table):
    s = table[table.screen_pass == True]
    f = crc.flags(s)
    assert len(s) == 424
    assert (f.T_low.sum(), f.T_high.sum(), f.Cr_low.sum(), f.Cr_high.sum(), f.no_Cr.sum()) == (231, 6, 2, 2, 0)
