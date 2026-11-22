"""Revision v3 (option C): flat WZ07 match level anchored to the median natural Slave calibration discrepancy.
Kept as a regression test of the v3 anchor rule: run with the v3 parameters (no Cr2O3 correction). With the
analytical-quality screen the anchor set is the 23 screened Slave samples (median -99.48 C; v3 had 30 and -96.04 C);
the default Eq. 4 correction is tested in tests/test_cr_correction.py."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
import make_figures as mf  # noqa: E402
from calculations.pt_fo2.params import PARAMS  # noqa: E402
from calculations.pt_fo2 import level_anchor as la  # noqa: E402
from calculations.pt_fo2 import calibration as calib  # noqa: E402


V3_LEVEL = -99.48   # v3 rule (Slave median of D) on the quality-screened anchor set


@pytest.fixture(autouse=True, scope="module")
def v3_params():
    fp = PARAMS.fe3
    saved = (fp.wz07_correction, fp.dt_match_c)
    fp.wz07_correction, fp.dt_match_c = "experimental", V3_LEVEL
    yield
    fp.wz07_correction, fp.dt_match_c = saved


@pytest.fixture(scope="module")
def df():
    d = pd.read_csv(mf.FULL, encoding="utf-8-sig", low_memory=False)
    return d.rename(columns={d.columns[0]: "SampleID"})


def test_params_level_is_the_anchor(df):
    a = la.anchor_level(df)
    assert a["n"] == 23
    assert PARAMS.fe3.dt_match_c == pytest.approx(a["median_D_C"], abs=0.05)
    assert a["median_r_C"] == pytest.approx(a["median_D_C"])        # no curve: r = D
    assert a["ci95_bootstrap_D_C"][0] < a["median_D_C"] < a["ci95_bootstrap_D_C"][1]


def test_no_slope_no_offset():
    fp = PARAMS.fe3
    assert fp.wz07_offset_c == 0.0 and fp.wz07_offset_slope == 0.0
    assert np.allclose(la.curve([900.0, 1200.0, 1500.0]), 0.0)
    assert not fp.apply_empirical_correction


def test_anchor_uses_all_calibration_and_no_validation_samples(df):
    m = la.anchor_set(df)
    assert not (m & calib.validation_mask(df)).any()
    assert (m == (calib.calibration_mask(df) & calib.measured_fe3(df).notna())).all()
    _, _, blue, red, _, _ = calib.error_groups(df)
    assert (m & (blue | red)).sum() == len(la.anchor_level(df)["error_groups"]) > 0   # reported, not excluded


def test_match_level_is_slave_median_at_pivot(df):
    level, pivot = la.match_level(df)
    a = la.anchor_level(df)
    assert pivot == pytest.approx(a["median_T_TA98_C"])
    assert level == pytest.approx(a["median_r_C"] + float(la.curve(pivot)))


def test_anchor_is_median_of_calibration_discrepancy(df):
    """By construction about half the calibration samples stop above / below their measured Fe3+."""
    m = la.anchor_set(df)
    D = la.discrepancy_measured(df)[m]
    assert ((D <= PARAMS.fe3.dt_match_c).mean()) == pytest.approx(0.5, abs=0.04)
