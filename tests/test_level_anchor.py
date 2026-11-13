"""Revision v3 (option C): flat WZ07 match level anchored to the median natural Slave calibration discrepancy."""
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


@pytest.fixture(scope="module")
def df():
    d = pd.read_csv(mf.FULL, encoding="utf-8-sig", low_memory=False)
    return d.rename(columns={d.columns[0]: "SampleID"})


def test_params_level_is_the_anchor(df):
    a = la.anchor_level(df)
    assert a["n"] == 30
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
    assert df.loc[m, "SampleID"].astype(str).isin(la.BLUE_IDS + la.RED_IDS).sum() == 5   # not excluded


def test_error_groups_unchanged_at_adopted_level(df):
    _, _, blue, red, rem, _ = calib.error_groups(df)
    assert sorted(df.SampleID[blue].astype(str)) == sorted(la.BLUE_IDS)
    assert sorted(df.SampleID[red].astype(str)) == sorted(la.RED_IDS)


def test_anchor_is_median_of_calibration_discrepancy(df):
    """By construction about half the calibration samples stop above / below their measured Fe3+."""
    m = la.anchor_set(df)
    D = la.discrepancy_measured(df)[m]
    assert ((D <= PARAMS.fe3.dt_match_c).mean()) == pytest.approx(0.5, abs=0.04)
